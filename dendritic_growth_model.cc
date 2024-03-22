/*
  © Copyright 2008 Randal A. Koene <randalk@netmorph.org>
  
  With design assistance from J. van Pelt & A. van Ooyen, and support
  from the Netherlands Organization for Scientific Research (NWO)
  Program Computational Life Sciences grant CLS2003 (635.100.005) and
  from the EC Marie Curie Research and Training Network (RTN)
  NEURoVERS-it 019247.

  This file is part of NETMORPH.

  NETMORPH is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  NETMORPH is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.
*/
// dendritic_growth_model.cc
// Randal A. Koene, 20041118, 20051109

/* Documentation: The Process of Branching/Turning

The following steps are involved in the branching and turning process of
dendritic or axonal arbor growth:
1) Probabilistic calculations in a grow() function of a model derived from
   dendritic_growth_model (which may address growth in dendrites or axons)
   determine when (a) a branch occurs, (b) a turn occurs, or (c) an existing
   turn is converted into a new branch.
Note A: The SELECT_TYPED_BRANCH macro enables branching with two different
        sets of angle constraints.
Note B: The fibre_segment::branches_at_continuation_nodes() function with
        a Delayed_Branching_Model derived parameter object contains the
        code that determines the possible conversion of a turn into a
        branch (option (c) above).
        In case (c), the recursive fibre_segment function calls
        Delayed_Branching_Model::branch_selection().
2) A van_Pelt_dendritic_growth_model derived branch(), dirmodel->direction(),
   turn(), turn_free() function may be called for options (a) and (b) above.
   When turn to branch conversion is selected, in a model derived from
   Delayed_Branching_Model::branch_selection() then the
   Delayed_Branching_Model::continuation_node_to_branch() function is called.
3) The branch(), dirmodel->direction(), turn() or turn_free() functions
   calculate spatial information for the new fibre segments and may use the
   veer() function in dendritic_growth_model.cc.
   These then call terminal_segment::turn_without_update_of_segment_vector()
   or terminal_segment::branch() functions.
   The Delayed_Branching_Model::continuation_node_to_branch() function
   calculates spatial information for the added fibre segment.
4) The terminal_segment::turn() or
   terminal_segment::turn_without_update_of_segment_vector() functions call
   fibre_segment::turn(), then link new terminal segment objects into the
   list of terminal segments and remove the old one.
   The terminal_segment::branch() funcion calls fibre_segment::branch(),
   then creates two new terminal segments and adds them to the list of
   terminal segments, after which the old terminal segment is usually removed.
   The Delayed_Branching_Model::continuation_node_to_branch() also creates
   and links in another terminal_segment object for the new branch, then
   calls fibre_segment::continuation_node_to_branch().
5) The fibre_segment::turn() function creates a new fibre_segment object
   for branch1 and sets up initial location information (branch1->P0=P1).
   The fibre_segment::branch() function does the same, but for both branch1
   and branch2. fibre_segment::continuation_node_to_branch() makes the branch
   fibre segment.

Note: Many of the functions below expect that sampled_output points to a
      valid object.

NOTE: THE FIXED STEP CONVENTION HERE IS TO COMPUTE FROM t TO t+dt. THIS IS
      NOT THE SAME CONVENTION AS IN THE GENERATION-FRAMEWORK PAPER, WHERE
      COMPUTATIONS GO FROM t-dt TO t.
 */

#include <math.h>
#include "diagnostic.hh"
#include "network.hh"
#include "Sampled_Output.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"
#include "axon_direction_model.hh"
#include "BigString.hh"
#include "branching_models.hh"

#define DAYSECONDS 86400.0
#define DAYSECONDSSQR (86400.0*86400.0)
#define DAYSECONDSCUB (86400.0*86400.0*86400.0)

#define INCLUDE_BRANCH_ANGLE_SAMPLING
#define INCLUDE_TURN_ANGLE_SAMPLING

#define USE_REMAINDER_FUNCTION_FOR_ANGLE_DATA_STANDARDIZATION

#ifdef USE_REMAINDER_FUNCTION_FOR_ANGLE_DATA_STANDARDIZATION
#define STANDARDIZE_ANGLE_DATA(angle) \
  angle = remainder(angle,2.0*M_PI);\
  if (angle<-M_PI) angle += 2.0*M_PI;\
  else if (angle>M_PI) angle -= 2.0*M_PI;\
  angle = fabs(angle)
#else
#define STANDARDIZE_ANGLE_DATA(angle) \
  while (angle<-M_PI) angle += 2.0*M_PI;\
  while (angle>M_PI) angle -= 2.0*M_PI;\
  angle = fabs(angle)
#endif


// Some locally shared (cache) variables
#ifdef VECTOR3D
double cosphiSx, cosphi, sinphiSz, sinphi, costheta, sintheta;
#endif

int branch_at_continuation_node_counter = 0;

bool growing_dendrites;

// *** The following cache variables should be allocated as unique variables
// wherever turnanglemin and turnanglemax are allocated.
//double turn_a, turn_max_y, turn_SQb, turn_negb, turn_twoa;
#ifdef VECTOR2D
//double turn_twomax_y;
#endif

#ifdef TEST_FOR_NAN
unsigned long int DBcounter = 0;
#endif

#ifdef VECTOR3D
void prepare_parent_coordinate_system(spatial & parentacoords) {
  // Uses the spherical coordinates representation of a "parent" vector to establish
  // conversion variables that allow a "daughter" direction interpreted within the
  // parent vector axis coordinate system to be converted to the universal cartesian
  // coordinate system.
  sinphi = sin(parentacoords.Z());
  cosphi = cos(parentacoords.Z());
  sintheta = sin(parentacoords.Y());
  costheta = cos(parentacoords.Y());
}

void daughter_to_network_coordinate_system(spatial & daughter, spatial & daughterxyz) {
  // This converts daughter spatial data from spherical coordinates in the
  // parent coordinate system to spherical coordinates in the network
  // coordinate system and initializes the daughter length.
  // Parameters:
  //   daughter = direction vector in the parent vector axis system
  //   daughterxyz = result, direction vector in the universal cartesian coordinate system, in spherical coordinates
  // The resulting direction vector has length zero.
  // Used by:
  //   axon_direction_model.cc:direction_model_base::turn() in 3D
  //   dendritic_growth_model.cc:Delayed_Branching_Model::continuation_node_to_branch() in 3D
  //   dendritic_growth_model.cc:veer() in 3D
  //   dendritic_growth_model.cc:van_Pelt_dendritic_growth_model::branch() in 3D
  //   dendritic_growth_model.cc:turn() in 3D
  daughter.convert_from_spherical();
  cosphiSx = cosphi*daughter.X();
  sinphiSz = sinphi*daughter.Z();
  //   Sx(in x,y,z) = cos(phi)cos(theta)Sx-sin(theta)Sy+sin(phi)cos(theta)Sz
  daughterxyz.set_X((costheta*cosphiSx) - (sintheta*daughter.Y()) + (costheta*sinphiSz));
  //   Sy(in x,y,z) = cos(phi)sin(theta)Sx+cos(theta)Sy+sin(phi)sin(theta)Sz
  daughterxyz.set_Y((sintheta*cosphiSx) + (costheta*daughter.Y()) + (sintheta*sinphiSz));
  //   Sz(in x,y,z) = -sin(phi)Sx+0Sy+cos(phi)Sz
  daughterxyz.set_Z((cosphi*daughter.Z()) - (sinphi*daughter.X()));
  daughterxyz.convert_to_spherical(); // ***IMMEDIATELY BEFORE THIS CARTESIAN MODIFICATION OF THE DIRECTION IN UNIVERSAL COORDINATES IS POSSIBLE
  daughterxyz.set_X(0.0); // initial length zero
}

double veer(terminal_segment * ts, double veerangle, spatial & S1xyz, double rollangle) {
  // The abstract "veer" function is supplied so that it can be used
  // both for branch functions and those turn functions that do not involve
  // the Samsonovich hypothesis for dendritic turns.
  // If the ROLL angle given is less than zero, then a roll angle is drawn
  // from a uniform random distribution.
  // The ROLL angle is returned.
  // This process is described in the figure
  // ~/src/nnmodels/nibr/3D-branch-turn-angles.fig.
  // Spherical angles around the w axis, determined by constrained randomness
  if (rollangle<0.0) rollangle = 2.0*M_PI*X_turn.get_rand_real1();
  spatial S1(1.0,rollangle,veerangle);
  // Double rotation transformation T S, for phi and theta of the parent
  prepare_parent_coordinate_system(ts->AngularCoords());
  daughter_to_network_coordinate_system(S1,S1xyz);
  // Do not apply environment pressure and initial length here. That is done
  // upon return from this function, as the specific call warrants, just as
  // upon return from the daughter_to_network_coordinate_system() function.
  return rollangle;
}
#endif

void Delayed_Branching_Model::branch_selection(fibre_segment * fs) {
  // For more information, see the comments in
  // fibre_segment::branches_at_continuation_nodes().
  // [***INCOMPLETE] The base class version of this function is merely
  // and example and is not yet based on a rational model.
  if (X_branch.get_rand_real1()<=Btpercn) continuation_node_to_branch(fs,branchanglemax,branchanglemin);
}

void Delayed_Branching_Model::continuation_node_to_branch(fibre_segment * fs, double branchanglemin, double branchanglemax) {
  fibre_segment * existingbranch = fs->Branch1();
  if (!existingbranch) existingbranch = fs->Branch2();
  if (!existingbranch) {
    warning("Warning: fibre_segment was terminal_segment in Delayed_Branching_Model::continuation_node_to_branch()\n");
    return;
  }
  // [***INCOMPLETE] The code here may have to test whether the existing
  // branch would be considered a main or a daughter branch upon conversion
  // of the turn to a branch. To do this, a parent_phi() request function
  // may be needed.
  // Currently, the new branch is treated as the daughter branch.
  double angletoparent = X_branch.get_rand_range_real1(branchanglemin,branchanglemax); // in existing parent-daughter plane
  spatial acoords;
  fibre_segment * newbranch;
#ifdef VECTOR3D
#ifdef TEST_FOR_NAN
  cout << "DB(" << DBcounter << "), fs used to get acoords: P0.X()=" << fs->P0.X() << ", P1.X()=" << fs->P1.X() << '\n'; cout.flush();
  spatial a = fs->P1;
  a -= fs->P0;
  cout << "acoords vector in Euler: (" << a.X() << ',' << a.Y() << ',' << a.Z() << ")\n"; cout.flush();
#endif
  fs->get_angles_and_length(acoords);
#ifdef TEST_FOR_NAN
  cout << "DELAYED BRANCH acoords: [get:] length=" << acoords.X() << ", theta=" << acoords.Y() << ", phi=" << acoords.Z(); cout.flush();
#endif
  prepare_parent_coordinate_system(acoords);
#ifdef TEST_FOR_NAN
  cout << " [prep:] length=" << acoords.X() << ", theta=" << acoords.Y() << ", phi=" << acoords.Z() << '\n'; cout.flush();
#endif
  double newparenttheta = existingbranch->parent_theta()+M_PI;
  spatial S2(1.0,newparenttheta,angletoparent);
  daughter_to_network_coordinate_system(S2,acoords);
  // [***NOTE] The ROLL angle (parent_theta()) must be stored
  // wherever a new fibre segment is created.
  // Create a new terminal segment object
  extra_fibre_data efd(newparenttheta);
  if ((newbranch=fs->continuation_node_to_branch(efd))==NULL) {
#endif
#ifdef VECTOR2D
  existingbranch->get_angles_and_length(acoords);
  double existingangle = acoords.Y();
  while (existingangle<0.0) existingangle += (2.0*M_PI);
  fs->get_angles_and_length(acoords);
  double parentangle = acoords.Y();
  while (parentangle<0.0) parentangle += (2.0*M_PI);
  if (existingangle>=parentangle) angletoparent = -angletoparent;
  acoords.set_all(0.0,parentangle+angletoparent);
  // Create a new terminal segment object
  if ((newbranch=fs->continuation_node_to_branch())==NULL) {
#endif
    warning("Warning: Unable to create branch fibre_segment in Delayed_Branching_Model::continuation_node_to_branch()\n");
    return;
  }
#ifdef TEST_FOR_NAN
  String nanstr(newbranch->P0.X(),"%f");
  bool newisNAN = (nanstr==String("nan"));
  if (newisNAN) { cout << "New branch fiber at Delayed Branch is NAN!\n"; cout.flush(); }
#ifdef VECTOR3D
  cout << " [result:] length=" << acoords.X() << ", theta=" << acoords.Y() << ", phi=" << acoords.Z() << '\n'; cout.flush();
#endif
#endif
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
  // Here we create a new terminal segment and insure it has
  // a direction model.
  // Determine centrifugal order:
  int centord = 1;
  for (fibre_segment * fsp = fs->parent; (fsp); fsp = fsp->parent) if ((fsp->branch1) && (fsp->branch2)) centord++;
  terminal_segment * ts = new terminal_segment(*arbor,*newbranch,acoords,centord);
  register PLLRoot<terminal_segment> * tsroot = arbor->TerminalSegments();
  register direction_model_base * dirmodel = tsroot->head()->DirectionModel();
  // ALL HANDLERS FOR MODEL PROPAGATION SHOULD BE CALLED HERE (Look for similar notes such as this elsewhere, e.g. in terminal_segment::branch().)
  if (dirmodel) ts->set_direction_model(dirmodel->clone());
  ts->set_elongation_model(tsroot->head()->ElongationModel()->clone(ts));
  ts->set_ERI_model(tsroot->head()->ElongationRateInitializationModel()->clone());
  ts->set_TSBM_model(tsroot->head()->TSBModel()->clone());
  ts->set_branch_angle_model(tsroot->head()->BranchAngleModel()->clone());
  tsroot->link_after(ts);
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,0,sizeof(terminal_segment));
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_BRANCH_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_branch_angles,dendrites_branch_angles};
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    // [***NOTE] I need to determine which is which if I am to collect data
    // about slightly and strongly deviating daughter segments separately.
    STANDARDIZE_ANGLE_DATA(angletoparent);
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(angletoparent);
  }
#endif
#endif
}

van_Pelt_dendritic_growth_model::van_Pelt_dendritic_growth_model(bool _dendrites, TSTM tm, double tf, double tr, bool bt): dendrites(_dendrites), competition_with_all_trees(gcc_within_tree), branchesatturns(bt) {
  set_fixed_time_step_size(van_Pelt_SUGGESTED_dt);
}

van_Pelt_dendritic_growth_model::van_Pelt_dendritic_growth_model(bool _dendrites, Command_Line_Parameters & clp): dendrites(_dendrites), competition_with_all_trees(gcc_within_tree), branchesatturns(false) {
  set_fixed_time_step_size(van_Pelt_SUGGESTED_dt);
  parse_CLP(clp);
}

const char neuritetypeletter[2] = {'A','D'};
const char neuritetypeword[2][9] = {"axon","dendrite"};

void van_Pelt_dendritic_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  //set_fixed_time_step_size(_dt);
  if (fibreswithturns) {
    if ((n=clp.Specifies_Parameter(String(neuritetypeletter[dendrites])+"branchesatturns"))>=0) branchesatturns = (downcase(clp.ParValue(n))==String("true"));
  }
}

String van_Pelt_dendritic_growth_model::report_parameters() {
  String res("van_Pelt or Polynomial O1/O3 ");
  res += neuritetypeword[dendrites];
  res += " growth model:\n";
  return res;
}

void van_Pelt_dendritic_growth_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
}

void van_Pelt_dendritic_growth_model::initialize(network * net, Connection_Statistics_Root * cstats, double & t) {
  if (!net) error("nibr Error: No network in van_Pelt_dendritic_growth_model::initialize()\n");
  //if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_dendritic_growth_model::initialize()\n");
  if (outattr_show_progress) progress("Initializing first "+String(neuritetypeword[dendrites])+" segments...\n");
  // initialize the fibre structures
  // [***NOTE] The actual length of the root terminal segments can be modified during the prepost constructor, based
  // on region and natural set specific L0_min and L0_max specifications. (See TL#200808071040.1.)
  if (dendrites) PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->initialize_input_structure(1.0,1.0);
  else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->initialize_output_structure(1.0,1.0);
  // grow the fibre structures
  firstiteration = initiallengthatactualfirstbranch; // only branch at the first iteration if "branch at initlength" is true
}

// Some functions that simplify dendrite specific elongation requests made in
// multiple growth models below:
void van_Pelt_dendritic_growth_model::request_elongation(fibre_structure * first_fs) {
  // first_fs is the head element of the list of fibre_structure roots of a neuron
  PLL_LOOP_FORWARD(fibre_structure,first_fs,1) {
#ifdef ELONGATION_EVENT_PREDICTION
    e->devactivity->queue(_arbor_development_activity_elongation);
#else
    e->batch_elongate();
#endif
  }
}

// Growth "batch" functions for dendrites:

// Note that at time t = 0.0 there was no dt step preceding, so the same
// branching calculation does not apply. Unless we know at what time in the
// cell's life the initial segments were completed, we also don't know how
// to apply the calculation using times t1 and t2.
// How to regard t at that time determines which t1 and t2 values should be
// applied. If t=0, i.e. t2=t1, or very small then there is no or almost no
// expectation of a branch - so that the initial segment is in effect
// lengthened. The second implemented function always assumes that at least
// dt time has passed, so that it makes no sense when used in combination
// with t1=0, unless we are to consider negative time in which the initial
// segment was grown.

void van_Pelt_dendritic_growth_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // growing in fixed time steps
  // <A NAME="van_Pelt_dendritic_growth_plus_turns_model::grow"> </A>
  // Note that this is actually a growth "batch" function, as it directly
  // calls branch model functions and requests elongation updates for all
  // terminal segments.
  // Making turns at continuation points is included as a simulation option.
  growing_dendrites = dendrites;
  fibre_structure * first_fs;
  if (!net) error("nibr Error: No network in van_Pelt_dendritic_growth_model::grow()\n");
  //if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_dendritic_growth_model::grow()\n");
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    if (dendrites) first_fs = e->InputStructure()->head();
    else first_fs = e->OutputStructure()->head();
    request_elongation(first_fs); // See TL#200603151008.1.
    // *** request_branching(e); // This is where the arbor-specific branching models will be called (see branching_models.hh). (Update 20080314:) Unless called by the TSBM and then cached.
    PLL_LOOP_FORWARD_NESTED(fibre_structure,first_fs,1,ps) {
      terminal_segment * ts_next; bool DBM_parse_structure = fibreswithturns; // only parse once per fiber structure per simulation step if fibreswithturns==true
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	double p_i = 1.0; // in case of "branch at initlength"
	if (!firstiteration) {
	  p_i = ts->TSBModel()->branching_probability(*ts);
	  if ((branchesatturns) && (DBM_parse_structure)) { // Note: See TL#200808060250.1 this method of calling Delayed Branching Models.
	    // [***NOTE] Further consideration of delayed branching models should include making
	    // such branching part of the overall branching statistics that are subject to arbor
	    // wide constraints.
	    Delayed_Branching_Model dbm(e,ps,p_i); // [***NOTE] This could flexibly choose a DBM model and allocate or call it here.
	    ps->branches_at_continuation_nodes(&dbm);
	    DBM_parse_structure=false;
	  }
	}
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  fibre_segment * fs = ts->TerminalSegment();
	  double r = X_branch.get_rand_real1();
	  if (r<=p_i) ts->branch(NULL); // possible bifurcation point (branching)
	  else
	    if (fibreswithturns)
	      if (ts->TSTModel()->turn_event(r,*ts,p_i)) ts->turn(e); // possible coninuation point (turning)
	      else fs = NULL;
	    else fs = NULL;
	  if (fs) { // The fiber segment is no longer a terminal segment
	    if (net->Seek_Candidate_Synapses()) {
	      if (dendrites) {
		if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning("Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_model::grow() (at t="+String(t,"%.1f)\n"));
	      } else {
		if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning("Warning: Unable to add segment to axonal spatial subsets in van_Pelt_dendritic_growth_model::grow() (at t="+String(t,"%.1f)\n"));
	      }
	    }
	    if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);    
	  }
	}
      }
    }
  }
  firstiteration = false;
}
// [***NOTE] There is a small issue here. At present the branching probability updates per arbor on every
//           iteration, regardless whether a particular segment has been unable to branch, due to the
//           minimum node separation restriction. After the minimum separation is reached, calculation of
//           the probability continues, as if the time interval since the last such valid calculation was
//           smaller than the time interval since creation of the last node. It is presently unclear if
//           this is the explicit intent. (20080314)

bool van_Pelt_dendritic_growth_model::IsComplete(double & t) {
  t += _dt;
  return (t>max_growth_time);
}

void van_Pelt_dendritic_growth_model::postop(network * net, Connection_Statistics_Root * cstats, double & t) {
  // Note that the fibre_segment coordinates are updated here for each
  // terminal segment.
  if (dendrites) {
    if (net->Seek_Candidate_Synapses()) {
      PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
	PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
	  PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
	    ts->update_segment_vector();
	    if (!net->spatial_segment_subset()->add_dendritic_segment(ts->TerminalSegment())) warning("Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_model::postop() (at t="+String(t,"%.1f)\n"));
	  }
	}
	if (outattr_show_progress) progress('_');
      }
    } else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
	     PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps)
	     PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) ts->update_segment_vector();
  } else {
    if (net->Seek_Candidate_Synapses()) {
      PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
	PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
	  PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
	    ts->update_segment_vector();
	    if (!net->spatial_segment_subset()->add_axonal_segment(ts->TerminalSegment())) warning("Warning: Unable to add segment to axonal spatial subsets in van_Pelt_dendritic_growth_model::postop() (at t="+String(t,"%.1f)\n"));
	  }
	}
	if (outattr_show_progress) progress('_');
      }
    } else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
	     PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps)
	     PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) ts->update_segment_vector();
  }
}

#ifdef TESTING_RANDOM_DISTRIBUTION
unsigned long int randomdist[11] = {0,0,0,0,0,0,0,0,0,0,0};
#endif

/*//#define TEST_AXON_POSTOP
void van_Peltlike_axonal_growth_model::postop(network * net, Connection_Statistics_Root * cstats, double & t) {
  // Note that the fibre_segment coordinates are updated here for each
  // terminal segment.
  if (net->Seek_Candidate_Synapses()) {
#ifdef TEST_SSS_DEPTH
	  cout << "Min. memory consumption per SSS: " << sizeof(Spatial_Segment_Subset)+(sizeof(fibre_segment_ptr)*spatialsegmentsubsetsizelimit) << " bytes\n"; cout.flush();
#endif
#ifdef TEST_AXON_POSTOP
    int tap_n = 0;
#endif
    PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
#ifdef TEST_AXON_POSTOP
    int tap_a = 0;
#endif
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
#ifdef TEST_AXON_POSTOP
    int tap_t = 0;
#endif
	PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
	  ts->update_segment_vector();
	  if (!net->spatial_segment_subset()->add_axonal_segment(ts->TerminalSegment())) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Peltlike_axonal_growth_model::postop() (at t=%.1f)\n"));
#ifdef TEST_SSS_DEPTH
	  cout << "SSS#" << sss_number << '\n'; cout.flush();
#endif
#ifdef TEST_AXON_POSTOP
	  cout << "N/A/T = " << tap_n << '/' << tap_a << '/' << tap_t << '\n'; cout.flush();
	  tap_t++;
#endif
	}
#ifdef TEST_AXON_POSTOP
	tap_a++;
#endif
      }
#ifdef TEST_AXON_POSTOP
      tap_n++;
#endif
      if (outattr_show_progress) progress('_');
    }
#ifdef TEST_SSS_DEPTH
    cout << "Final SSS#" << sss_number << " min. memory consumption per SSS: " << sizeof(Spatial_Segment_Subset)+(sizeof(fibre_segment_ptr)*spatialsegmentsubsetsizelimit) << " bytes\n"; cout.flush();
#endif
  } else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
	   PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps)
	     PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) ts->update_segment_vector();
}
*/
