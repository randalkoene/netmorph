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
// fibre_structure.cc
// Randal A. Koene, 20041119

//#include "nibr.hh"
#include "Fig_Object.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"
#include "axon_direction_model.hh"
//#include "branching_models.hh"
#include "Network_Generated_Statistics.hh"
#include "Sampled_Output.hh"

general_fibre_structure_parameters_interface general_fibre_structure_parameters;

region_structure_initialization_ptr * structure_initialization_region_subset_schemas = NULL;

nodegenesis_data * NodeGenesis_Data = NULL;

const char fibre_structure_name[NUM_fs][9] = {
  "axon",
  "dendrite"
};
fibre_structure_class_id parsing_fs_type = axon_fs; // Modify this wherever it is important to keep track (e.g. neuron::net_Txt())

terminal_segment * most_recent_branch1_ts = NULL;
terminal_segment * most_recent_branch2_ts = NULL;

void general_fibre_structure_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("fibrediameter"))>=0) fibrediameter = (downcase(clp.ParValue(n))==String("true"));
}

String general_fibre_structure_parameters_interface::report_parameters() {
  String res("Fiber structure parameters:\n  fibres have ");
  if (fibrediameter) res += "specific ";
  else res += "uniform ";
  res += "diameters.\n";
  return res;
}

#ifdef VECTOR3D
extra_fibre_data root_efd(0.0);

bool fibre_segment::branch(extra_fibre_data & efd1, extra_fibre_data & efd2) {
#endif
#ifdef VECTOR2D
bool fibre_segment::branch() {
#endif
  // Note: This leaves P1 of the branches undefined! If this is undesired
  // the the new fibre_segments should be constructed with the call
  // new fibre_segment(P1,P1).
  // The current behavior is designed to work efficiently when branches
  // are created by terminal_segment::branch().
  if (branch1) return false;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
#ifdef VECTOR3D
  branch1 = new fibre_segment(*n,this,efd1);
  branch2 = new fibre_segment(*n,this,efd2);
#endif
#ifdef VECTOR2D
  branch1 = new fibre_segment(*n,this);
  branch2 = new fibre_segment(*n,this);
#endif
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,2,sizeof(*branch1));
  branch1->P0 = P1;
  branch2->P0 = P1;
  return true;
}

#ifdef VECTOR3D
bool fibre_segment::turn(extra_fibre_data & efd) {
#endif
#ifdef VECTOR2D
bool fibre_segment::turn() {
#endif
  if (branch1) return false;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
#ifdef VECTOR3D
  branch1 = new fibre_segment(*n,this,efd);
#endif
#ifdef VECTOR2D
  branch1 = new fibre_segment(*n,this);
#endif
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,1,sizeof(*branch1));
  branch1->P0 = P1;
  return true;
}

#ifdef VECTOR3D
fibre_segment * fibre_segment::continuation_node_to_branch(extra_fibre_data & efd) {
#endif
#ifdef VECTOR2D
fibre_segment * fibre_segment::continuation_node_to_branch() {
#endif
  if (branch1 && branch2) return NULL;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
#ifdef VECTOR3D
  fibre_segment * newbranch = new fibre_segment(*n,this,efd);
#endif
#ifdef VECTOR2D
  fibre_segment * newbranch = new fibre_segment(*n,this);
#endif
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,1,sizeof(*newbranch));
  newbranch->P0 = P1;
#ifdef TEST_FOR_NAN
  String nanstr(newbranch->P0.X(),"%f");
  bool newisNAN = (nanstr==String("nan"));
  if (newisNAN) { 
    cout << "New branch DBM P0 is set to NAN!\n";
    String oldnanstr(P1.X(),"%f");
    if (oldnanstr==String("nan")) {
      cout << "The P1 from which it inherits is also NAN!\n";
      String oldPOnanstr(P0.X(),"%f");
      if (oldPOnanstr==String("nan")) cout << "The P0 of the existing branch was also NAN!\n";
      if (branch1) {
	String branchnanstr(branch1->P0.X(),"%f");
	if (branchnanstr==String("nan")) cout << "Its pre-existing continuation on branch1 had NAN P0!\n";
      }
      if (branch2) {
	String branchnanstr(branch2->P0.X(),"%f");
	if (branchnanstr==String("nan")) cout << "Its pre-existing continuation on branch2 had NAN P0!\n";
      }
    }
    cout.flush();
  }  
#endif
  if (!branch1) branch1 = newbranch;
  else branch2 = newbranch;
  return newbranch;
}

void fibre_segment::branches_at_continuation_nodes(Delayed_Branching_Model * dbm) {
  // This function can be used both for presynaptic and postsynaptic structures
  // and it searchest through a tree for continuation nodes, possibly causing
  // a branch to appear at such a node.
  // Method:
  // 1) A Delayed_Branching_Model is created before the call to this function.
  //    The corresponding object (dbm) contains all the necessary data and
  //    functions with which to make the selections based on information about
  //    (a) a specific continuation node and (b) delayed branching events that
  //    have taken place during this call. The information in dbm can include
  //    a pointer to the neuron that is involved.
  // 2) A call to this recursive function progresses through all the branches
  //    (that are not terminal segments).
  // 3) For each continuation node encountered (i.e. a fibre segment that has
  //    exactly one branch), the selection procedure is enacted by passing
  //    that fibre segment to dbm->branch_selection().
  // [*** NOTE] A method similar to this, or even the same one could be
  // applied to terminal segments when cleaning up the code or rewriting the
  // program.
#ifdef TEST_FOR_NAN
  DBcounter++;
#endif
  if ((branch1 || branch2) && (!(branch1 && branch2))) dbm->branch_selection(this);
  if (branch1) branch1->branches_at_continuation_nodes(dbm);
  if (branch2) branch2->branches_at_continuation_nodes(dbm);
}

double fibre_segment::cartesian_soma_center_distance2() {
  // Returns the square of the cartesian distance between the center of
  // the neuron soma and the end point of the fiber segment.
  return distance2(P1,n->Pos());
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
double fibre_segment::cartesian_to_fiber_ratio(fibre_collected_data & fcd) {
  // Calculates the ratio of cartesian distance and fiber distance, given
  // the fiber segment before a bifurcation and this fiber segment.
  // This takes into account the possibility that there was no previous
  // bifurcation, in which case the cartesian distance is calculated from
  // the center of the soma, from which the cell radius is then substracted.
  // Note that if fcd.lensincebifurcation==0.0 then the cartesian distance
  // must also be 0.0. To avoid a "nan" result in such cases, return 1.0.
  if (fcd.lensincebifurcation<=0.0) return 1.0;
  if (fcd.prevbifurcation) return P1.distance(fcd.prevbifurcation->P1)/fcd.lensincebifurcation;
  double d = P1.distance(n->Pos());
  d -= n->Radius();
  return d/fcd.lensincebifurcation;
}

double fibre_segment::soma_to_term_cartesian_to_fiber_ratio(double len) {
  if (len<=0.0) return 1.0;
  double d = P1.distance(n->Pos());
  d -= n->Radius();
  return d/len;
}

void fibre_segment::fiber_location(double dist, spatial & loc) {
  // [***INCOMPLETE] Implement this.
}

void fibre_segment::Collect_Data(fibre_collected_data & fcd) {
  double len = Length();
  fcd.arborlength += len;
  fcd.lensincesoma += len;
  fcd.lensincebifurcation += len;
  if ((!branch1) && (!branch2)) { // terminal segment
    if (fcd.ispresynaptic) {
      fcd.nsb->netstats[axons_turns_between_bifurcations].Collect_Data((double) fcd.turnssincebifurcation);
      fcd.nsb->netstats[axons_term_len_since_bifurcation].Collect_Data(fcd.lensincebifurcation);
      if (fcd.nsb->netstats[axons_cartratio_bifurcationtoterm].Collect()) fcd.nsb->netstats[axons_cartratio_bifurcationtoterm].Collect_Data(cartesian_to_fiber_ratio(fcd));
      if (fcd.nsb->netstats[axons_cartratio_somatoterm].Collect()) fcd.nsb->netstats[axons_cartratio_somatoterm].Collect_Data(soma_to_term_cartesian_to_fiber_ratio(fcd.lensincesoma));
      if (fcd.nsb->netstats[axons_term_len_since_soma].Collect()) fcd.nsb->netstats[axons_term_len_since_soma].Collect_Data(fcd.lensincesoma);
    } else {
      fcd.nsb->netstats[dendrites_turns_between_bifurcations].Collect_Data((double) fcd.turnssincebifurcation);
      fcd.nsb->netstats[dendrites_term_len_since_bifurcation].Collect_Data(fcd.lensincebifurcation);
      if (fcd.nsb->netstats[dendrites_cartratio_bifurcationtoterm].Collect()) fcd.nsb->netstats[dendrites_cartratio_bifurcationtoterm].Collect_Data(cartesian_to_fiber_ratio(fcd));
      if (fcd.nsb->netstats[dendrites_cartratio_somatoterm].Collect()) fcd.nsb->netstats[dendrites_cartratio_somatoterm].Collect_Data(soma_to_term_cartesian_to_fiber_ratio(fcd.lensincesoma));
      if (fcd.nsb->netstats[dendrites_term_len_since_soma].Collect()) fcd.nsb->netstats[dendrites_term_len_since_soma].Collect_Data(fcd.lensincesoma);
    }
    return; // also see terminal_segment::Collect_Data();
  }
  if ((branch1) && (branch2)) { // bifurcates at end
    if (fcd.ispresynaptic) {
      fcd.nsb->netstats[axons_len_between_bifurcations].Collect_Data(fcd.lensincebifurcation);
      fcd.nsb->netstats[axons_turns_between_bifurcations].Collect_Data((double) fcd.turnssincebifurcation);
      if (fcd.nsb->netstats[axons_cartratio_between_bifurcations].Collect()) fcd.nsb->netstats[axons_cartratio_between_bifurcations].Collect_Data(cartesian_to_fiber_ratio(fcd));
      if (fcd.nsb->netstats[axons_seven_micron_branch_angles].Collect()) {
	spatial p,b1, b2;
	// [***INCOMPLETE] Implement the fiber_location() functions.
	fiber_location(-7.0,p); // a location 7 micrometers before bifurcation
	branch1->fiber_location(7.0,b1); // a location 7 micrometers along
	branch2->fiber_location(7.0,b2); // a location 7 micrometers along
	// Compute the angles (see http://rak.minduploading.org:8080/caspan/Members/randalk/results/validating-network-structure).
	// [***INCOMPLETE] Add the angle calculation and storage code here.
	// *** subtract P1 from p, b1, and b2
	// *** compute the dot product p1*b1 and p1*b2
	// *** divide each by all lengths (e.g. |p1||b1|), normalizing
	// *** obtain arccos of the results
	// [***INCOMPLETE] Deal with location distances less than 7 microns
	// as indicated in TL#200605111045.
      }
    } else {
      fcd.nsb->netstats[dendrites_len_between_bifurcations].Collect_Data(fcd.lensincebifurcation);
      fcd.nsb->netstats[dendrites_turns_between_bifurcations].Collect_Data((double) fcd.turnssincebifurcation);
      if (fcd.nsb->netstats[dendrites_cartratio_between_bifurcations].Collect()) fcd.nsb->netstats[dendrites_cartratio_between_bifurcations].Collect_Data(cartesian_to_fiber_ratio(fcd));
      if (fcd.nsb->netstats[dendrites_seven_micron_branch_angles].Collect()) {
	spatial p,b1, b2;
	// [***INCOMPLETE] Implement the fiber_location() functions.
	fiber_location(-7.0,p); // a location 7 micrometers before bifurcation
	branch1->fiber_location(7.0,b1); // a location 7 micrometers along
	branch2->fiber_location(7.0,b2); // a location 7 micrometers along
	// Compute the angles (see http://rak.minduploading.org:8080/caspan/Members/randalk/results/validating-network-structure).
	// [***INCOMPLETE] Add the angle calculation and storage code here.
	// [***INCOMPLETE] Deal with location distances less than 7 microns
	// as indicated in TL#200605111045.
      }
    }
    fcd.branch_at_end(*this);
  } else fcd.turnssincebifurcation++; // turn at end
  if (fcd.ispresynaptic) fcd.nsb->netstats[axons_intersegs_len].Collect_Data(len);
  else fcd.nsb->netstats[dendrites_intersegs_len].Collect_Data(len);
    if (branch1) {
    double cached_lensincesoma = fcd.lensincesoma;
    branch1->Collect_Data(fcd);
    if (branch2) { // reset for next branch
      fcd.lensincesoma = cached_lensincesoma;
      fcd.branch_at_end(*this);
    }
  }
  if (branch2) branch2->Collect_Data(fcd);
}
#endif

int fibre_segment::count_segments() {
  int res = 1;
  if (branch1) res += branch1->count_segments();
  if (branch2) res += branch2->count_segments();
  return res;
}

int fibre_segment::count_terminal_segments() {
  int res = 0;
  if (branch1) res += branch1->count_terminal_segments();
  if (branch2) res += branch2->count_terminal_segments();
  if (!res) return 1; // if there are no branches then this is a terminal segment
  return res;
}

void fibre_segment::count_segment_types(unsigned long & continuation, unsigned long & bifurcation, unsigned long terminal) {
  // Note that this function purposely does not initialize the return variables.
  if (branch1) {
    if (branch2) bifurcation++;
    else continuation++;
  } else {
    if (branch2) continuation++;
    else terminal++;
  }
  if (branch1) branch1->count_segment_types(continuation,bifurcation,terminal);
  if (branch2) branch2->count_segment_types(continuation,bifurcation,terminal);
}

double fibre_segment::sum_of_intermediate_segment_lengths(network_statistics_data * internsd) {
  if (!has_daughters()) return 0.0; // this is a terminal segment
  double res = Length();
  if (internsd) {
    internsd->mean += res;
    internsd->n += 1.0;
    if (res<internsd->min) internsd->min = res;
    if (res>internsd->max) internsd->max = res;
  }
  if (branch1) res += branch1->sum_of_intermediate_segment_lengths(internsd);
  if (branch2) res += branch2->sum_of_intermediate_segment_lengths(internsd);
  return res;
}

void fibre_segment::move_add(spatial & addvec) {
  P0 += addvec;
  P1 += addvec;
  if (branch1) branch1->move_add(addvec);
  if (branch2) branch2->move_add(addvec);
}

void fibre_segment::fanin_rot() {
#define FANIN_ANGLE_SLICE (M_PI/50.0)
  P0.convert_to_spherical();
  P0.set_Y(0.0);
  P0.convert_from_spherical();
  P1.convert_to_spherical();
  P1.set_Y(0.0);
#ifdef VECTOR3D
  double phi = P1.Z();
#endif
#ifdef VECTOR2D
  double phi = 0.0;
#endif
  P1.convert_from_spherical();
  if ((net_fanin_histogram.collect) && (!branch1) && (!branch2)) {
    unsigned int i = (int) (phi/FANIN_ANGLE_SLICE);
    if (i>=50) i=49;
    net_fanin_histogram.len[i] += (P1.len() - N()->Radius());
    net_fanin_histogram.n[i]++;
  }
  if (branch1) branch1->fanin_rot();
  if (branch2) branch2->fanin_rot();
}

//#define TEST_FIG_SIZES
Fig_Object * fibre_segment::net_Fig() {
  /* Documentation:
     Use of the "figattr_fibres_nobox" parameter can optionally avoid the
     computational overhead of fig_in_zoom() tests, and by defining a box
     around one neuron, but not around fibre structure a similar result can
     be achieved as with the connections structure output (see nibr.cc).
   */
  //cout << '['; cout.flush();
  double x1, y1, x2, y2;
  if (iszerolength()) warning_fibre_segment_net_Fig_zero_length++;
  Fig_Line * fibre_segment_fig = NULL;
  if (figattr_fibres_nobox || (fig_in_zoom(P0) && fig_in_zoom(P1))) {
#ifdef TEST_FOR_NAN
    if (isnan(P0.X())) { cout << "DID DETECT NAN!\n"; cout.flush(); }
#endif
    P0.plane_mapped(x1,y1);
    P1.plane_mapped(x2,y2);
    DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
    if ((general_fibre_structure_parameters.fibrediameter) && (diameter>0.0)) { // include a representation of fiber thickness
      long x[4], y[4]; double r = 0.5*diameter;
      if (abs(x2-x1)>abs(y2-y1)) { // [***INCOMPLETE] This does not produce perfectly correct proportional representations of fiber diameters.
	x[0] = FIGSCALED(x1);
	x[1] = x[0];
	x[2] = FIGSCALED(x2);
	x[3] = x[2];
	y[0] = FIGSCALED((y1-r));
	y[1] = FIGSCALED((y1+r));
	y[2] = FIGSCALED((y2+r));
	y[3] = FIGSCALED((y2-r));
      } else {
	x[0] = FIGSCALED((x1-r));
	x[1] = FIGSCALED((x1+r));
	x[2] = FIGSCALED((x2+r));
	x[3] = FIGSCALED((x2-r));
	y[0] = FIGSCALED(y1);
	y[1] = y[0];
	y[2] = FIGSCALED(y2);
	y[3] = y[2];
      }
      fibre_segment_fig = new Fig_Polygon(0,1,_figvar_connection_color,_figvar_connection_color,_figvar_connection_depth,20,0.0,0,0,x,y,4);
    } else {
#endif
      long X1(FIGSCALED(x1)), Y1(FIGSCALED(y1)), X2(FIGSCALED(x2)), Y2(FIGSCALED(y2));
#define TESTPOSTMA
#ifdef TESTPOSTMA
#ifdef VECTOR3D
      if ((X1==-20) && (Y1==-710) && (X2==-2422) && (Y2==-2441)) cout << "FOUND ERROR LINE GENERATION:\nP0=(" << P0.X() << ',' << P0.Y() << ',' << P0.Z() << "), P1=(" << P1.X() << ',' << P1.Y() << ',' << P1.Z() << ")\n";
#endif
#endif
      fibre_segment_fig = new Fig_Line(0,1,_figvar_connection_color,7,_figvar_connection_depth,-1,0.0,0,0,X1,Y1,X2,Y2);
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
    }
#endif
    DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(*fibre_segment_fig));
    if (_figgroup) _figgroup->Add_Fig_Object(fibre_segment_fig);
  }
  //refresh_branch_xy(x2,y2);
  if (branch1) branch1->net_Fig();
  if (branch2) branch2->net_Fig();
  return fibre_segment_fig;
}

Txt_Object * fibre_segment::bifurcationnode_net_Txt(long parentnodeindex) {
  // This fiber is followed by two daughter fibers after a bifurcation node.
  // Output data about the bifurcation node.
  // Format: <node index>, <fiber segment label>, <axon/dendrite>, <x>, <y>, <z>, <neuron label>, <parent node index>, <simulated creation time>
  Txt_nodeindex++;
  long nodeindexcache = Txt_nodeindex; // To remember this after depth-first parse along one branch
  (*Txt_bifurcationnodelist) += String(Txt_nodeindex); // Txt_nodeindex is the index of this bifurcation node
  (*Txt_bifurcationnodelist) += ',';
  (*Txt_bifurcationnodelist) += String((long) this); // This is a memory location label that can be crossreferenced from synapse text output
  (*Txt_bifurcationnodelist) += ',';
  (*Txt_bifurcationnodelist) += fibre_structure_name[parsing_fs_type];
  double x=0.0, y=0.0, z=0.0;
  P1.get_all(x,y,z); // the bifurcation node is at the end-point of this fiber (equal to branch1->P0 and branch2->P0)
  (*Txt_bifurcationnodelist) += String(x,",%f");
  (*Txt_bifurcationnodelist) += String(y,",%f");
  (*Txt_bifurcationnodelist) += String(z,",%f,");
  (*Txt_bifurcationnodelist) += String((long) n);
  (*Txt_bifurcationnodelist) += ',';
  (*Txt_bifurcationnodelist) += String(parentnodeindex);
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  (*Txt_bifurcationnodelist) += String(diameter,",%f");
#endif
  if (NodeGenesis_Data) (*Txt_bifurcationnodelist) += String(NodeGenesis_Data->find_t_genesis(branch1),",%f"); // same as that of branch2
  (*Txt_bifurcationnodelist) += '\n';
  branch1->net_Txt(nodeindexcache);
  branch2->net_Txt(nodeindexcache);
  return NULL;
}

Txt_Object * fibre_segment::continuationnode_net_Txt(long parentnodeindex) {
  // This fiber is followed by a daughter fiber after a continuation node.
  // Output data about the continuation node.
  // Format: <node index>, <fiber segment label>, <axon/dendrite>, <x>, <y>, <z>, <neuron label>, <parent node index>, <simulated creation time>
  Txt_nodeindex++;
  fibre_segment * fs = branch1; if (!fs) fs = branch2; // whichever continues the segment of fiber
  (*Txt_continuationnodelist) += String(Txt_nodeindex); // Txt_nodeindex is the index of this continuation node
  (*Txt_continuationnodelist) += ',';
  (*Txt_continuationnodelist) += String((long) this); // This is a memory location label that can be crossreferenced from synapse text output
  (*Txt_continuationnodelist) += ',';
  (*Txt_continuationnodelist) += fibre_structure_name[parsing_fs_type];
  double x=0.0, y=0.0, z=0.0;
  P1.get_all(x,y,z); // the continuation node is at the end-point of this fiber (equal to fs->P0)
  (*Txt_continuationnodelist) += String(x,",%f");
  (*Txt_continuationnodelist) += String(y,",%f");
  (*Txt_continuationnodelist) += String(z,",%f,");
  (*Txt_continuationnodelist) += String((long) n);
  (*Txt_continuationnodelist) += ',';
  (*Txt_continuationnodelist) += String(parentnodeindex);
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  (*Txt_continuationnodelist) += String(diameter,",%f");
#endif
  if (NodeGenesis_Data) (*Txt_continuationnodelist) += String(NodeGenesis_Data->find_t_genesis(fs),",%f");
  (*Txt_continuationnodelist) += '\n';
  fs->net_Txt(Txt_nodeindex);
  return NULL;
}

Txt_Object * fibre_segment::net_Txt(long parentnodeindex) {
  // It is easier for a fiber segment to determine type for a daughter
  // segment than for itself, so that is the way this function works.
  // Format: <node index>, <fiber segment lable>, <axon/dendrite>, <x>, <y>, <z>, <neuron label>, <parent node index>, <simulated creation time>
  if ((branch1) && (branch2)) return bifurcationnode_net_Txt(parentnodeindex); // the fiber segment ends in a bifurcation node
  if ((branch1) || (branch2)) return continuationnode_net_Txt(parentnodeindex); // the fiber segment ends in a continuation node
  // This is a terminal fiber segment.
  // Output data about the growth cone.
  Txt_nodeindex++;
  (*Txt_terminalgrowthconelist) += String(Txt_nodeindex); // Txt_nodeindex is the index of this growth cone
  (*Txt_terminalgrowthconelist) += ',';
  (*Txt_terminalgrowthconelist) += String((long) this); // This is a memory location label that can be crossreferenced from synapse text output
  (*Txt_terminalgrowthconelist) += ',';
  (*Txt_terminalgrowthconelist) += fibre_structure_name[parsing_fs_type];
  double x=0.0, y=0.0, z=0.0;
  P1.get_all(x,y,z); // end-point of this fiber segment
  (*Txt_terminalgrowthconelist) += String(x,",%f");
  (*Txt_terminalgrowthconelist) += String(y,",%f");
  (*Txt_terminalgrowthconelist) += String(z,",%f,");
  (*Txt_terminalgrowthconelist) += String((long) n);
  (*Txt_terminalgrowthconelist) += ',';
  (*Txt_terminalgrowthconelist) += String(parentnodeindex);
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  (*Txt_terminalgrowthconelist) += String(diameter,",%f");
#endif
  if (NodeGenesis_Data) (*Txt_terminalgrowthconelist) += String(eq->T(),",%f"); // growth cone position at this time
  (*Txt_terminalgrowthconelist) += '\n';
  return NULL;
}

VRML_Object * fibre_segment::net_VRML() {
  spatial center, rotaxis;
  double theta, len;
  if (figattr_fibres_nobox || (fig_in_zoom(P0) && fig_in_zoom(P1))) {
    X3D_cylinder_center(P0,P1,center);
    if (X3D_SFrotation_to_vector(P0,P1,rotaxis,theta,len)) { // only show pieces of non-zero length
      (*VRML_fiberlist) += "<Transform center='0 0 0' rotation='";
      double x, y, z;
#ifdef VECTOR2D
      z = 0.0;
#endif
      rotaxis.get_all(x,y,z);
      (*VRML_fiberlist) += String(x,"%.2f");
      (*VRML_fiberlist) += String(y," %.2f");
      (*VRML_fiberlist) += String(z," %.2f");
      (*VRML_fiberlist) += String(theta," %.2f");
      (*VRML_fiberlist) += "' translation='";
      center.get_all(x,y,z);
      (*VRML_fiberlist) += String(x,"%.2f");
      (*VRML_fiberlist) += String(y," %.2f");
      (*VRML_fiberlist) += String(z," %.2f");
      (*VRML_fiberlist) += "'>\n<ProtoInstance name='NeuritePiece'>\n<fieldValue name='neuriteheight' value='";
      (*VRML_fiberlist) += String(len,"%.2f");
      (*VRML_fiberlist) += "'/>\n";
      if (parsing_fs_type==axon_fs) (*VRML_fiberlist) += "<fieldValue name='neuritecolor' value='.8 .8 .0'/>\n";
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
      (*VRML_fiberlist) += "<fieldValue name='neuriteradius' value='";
      (*VRML_fiberlist) += String(diameter,"%.2f");
      (*VRML_fiberlist) += "'/>\n";
#endif
      (*VRML_fiberlist) += "</ProtoInstance>\n</Transform>\n";
      //    (*VRML_fiberlist) += "<Transform translation=''>\n";
      // *** I can put spheres at the ends of the cylinders, but I can integrate this and all the transforms in the PROTO
      //    (*VRML_fiberlist) += "</Transform>\n";
    }
  }
  if (branch1) branch1->net_VRML();
  if (branch2) branch2->net_VRML();
  return NULL;
}

terminal_segment_ptr * fibre_structure::array_of_terminal_segments() {
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,terminalarraylength,-sizeof(terminal_segment_ptr));
  delete[] terminalsegmentsarray;
  terminalarraylength = terminalsegments.length();
  if (terminalarraylength==0) return NULL;
  terminalsegmentsarray = new terminal_segment_ptr[terminalarraylength];
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,terminalarraylength,sizeof(terminal_segment_ptr));
  int i = 0;
  PLL_LOOP_FORWARD(terminal_segment,terminalsegments.head(),1) terminalsegmentsarray[i++] = e;
  return terminalsegmentsarray;
}

double fibre_structure::total_arbor_length(double * totterminalseglength, network_statistics_data * internsd) {
  // Returns the sum of the lengths of all fibre segments in the arbor
  // defined by this fibre_structure.
  // [*** NOTE] This can be done more efficiently by adding to an arbor
  // length variable each time a segment is allocated and managing changes
  // to the variable whenever segments are modified (e.g. end points change
  // or the segment is pruned). This function is a reliable way to refresh
  // such a variable.
  // If totterminalseglength is provided then the sum of the lengths of
  // all terminal segments is stored in that variable.
  double totterminallen = 0.0;
  PLL_LOOP_FORWARD(terminal_segment,terminalsegments.head(),1) totterminallen += e->Length();
  if (totterminalseglength) *totterminalseglength = totterminallen;
  // Now the lengths of all brances must be obtained.
  double res = sum_of_intermediate_segment_lengths(internsd);
  cache = res + totterminallen; // value is cached for future references
  return cache;
}

void fibre_structure::center_of_gravity(spatial & cg) {
  int numtrees = 0;
  PLL_LOOP_FORWARD(fibre_structure,this,1) {
    cg += e->P0;
    numtrees++;
  }
  cg /= (double) numtrees;
}

#ifndef ELONGATION_EVENT_PREDICTION
void fibre_structure::batch_elongate() {
  // This function does not schedule events on a queue, instead calling event
  // handlers directly for all terminal segments. (See TL#200603071856.1.)
#ifdef TESTING_ELONGATION_TOTAL
  usedarborelongationfraction = 0.0;
#endif
#ifdef TEST_FOR_NAN
  terminal_segment * e_next = NULL; // Use this in case e self-destructs during the elongate() call (e.g. by branching in BESTLNN_pyramidal_AD_terminal_segment_elongation_model::elongate().
  for (terminal_segment * e = TerminalSegments()->head(); (e); e = e_next) {
    e_next = e->Next(); // cached for safety
    bool NANbefore = (e->TerminalSegment()->P0.X()==NAN);
    e->ElongationModel()->elongate(e);
    bool NANafter = (e->TerminalSegment()->P0.X()==NAN);
    if ((NANafter) && (!NANbefore)) { cout << "batch_elongate caused NAN\n"; cout.flush(); }
  }
#else
  terminal_segment * e_next = NULL; // Use this in case e self-destructs during the elongate() call (e.g. by branching in BESTLNN_pyramidal_AD_terminal_segment_elongation_model::elongate().
  for (terminal_segment * e = TerminalSegments()->head(); (e); e = e_next) {
    e_next = e->Next(); // cached for safety
    e->ElongationModel()->elongate(e);
  }
#endif
#ifdef TESTING_ELONGATION_TOTAL
  if (usedarborelongationfraction<0.9999) cout << "TESTING_ELONGATION_TOTAL: Used arbor elongation fraction = " << usedarborelongationfraction << '\n';
#endif
}
#endif

double fibre_structure::sum_of_perturbed_expected_elongations() {
  if (eq->T()<=sum_l_i_cache_t) return sum_l_i_cache;
  sum_l_i_cache = 0.0;
  PLL_LOOP_FORWARD(terminal_segment,TerminalSegments()->head(),1) sum_l_i_cache += e->ElongationModel()->perturbed_expected_elongation();
  sum_l_i_cache_t = eq->T();
  return sum_l_i_cache;     
}

int fibre_structure::count_terminal_segments_in_PLL() {
  int res = 0;
  PLL_LOOP_FORWARD(fibre_structure,head(),1) res += e->count_terminal_segments();
  return res;
}

void fibre_structure::move_add(spatial & addvec) {
  this->fibre_segment::move_add(addvec);
}

void fibre_structure::fanin_rot() {
  this->fibre_segment::fanin_rot();
}

Fig_Object * fibre_structure::net_Fig() {
  //cout << '('; cout.flush();
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
  Fig_Group * fibre_structure_fig = new Fig_Group(); // *** can add info in comment here;
  DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(Fig_Group));
  _figgroup = fibre_structure_fig;
  //  cout << "X(" << this << ")=" << x << ", Y=" << y << '\n';
  _figvar_connection_color = colnum; // set display color per fiber structure
  this->fibre_segment::net_Fig();
  _figgroup = NULL;
  //cout << ')'; cout.flush();
  return fibre_structure_fig;
}

Txt_Object * fibre_structure::net_Txt() {
  // Output data bout the root point of the fiber structure.
  // Format: <node index>, <fiber segment label>, <axon/dendrite>, <x>, <y>, <z>, <neuron label>, <simulated creation time>
  Txt_nodeindex++;
  (*Txt_fiberrootlist) += String(Txt_nodeindex); // Txt_nodeindex is the index of this root node
  (*Txt_fiberrootlist) += ',';
  (*Txt_fiberrootlist) += String((long) this); // This is a memory location label that can be crossreferenced from synapse text output
  (*Txt_fiberrootlist) += ',';
  (*Txt_fiberrootlist) += fibre_structure_name[parsing_fs_type];
  double x=0.0, y=0.0, z=0.0;
  P0.get_all(x,y,z); // begin-point of this fiber segment is the root point of the fiber structure
  (*Txt_fiberrootlist) += String(x,",%f");
  (*Txt_fiberrootlist) += String(y,",%f");
  (*Txt_fiberrootlist) += String(z,",%f,");
  (*Txt_fiberrootlist) += String((long) n);
  if (NodeGenesis_Data) (*Txt_fiberrootlist) += String(NodeGenesis_Data->find_t_genesis(this),",%f\n");
  fibre_segment::net_Txt(Txt_nodeindex);
  return NULL;
}

void nodegenesis_data::add(fibre_segment * s, double t) {
  if (dataidx>=datalen) {
    unsigned int newdatalen = datalen + 10240;
    fibre_segment_ptr * nptr = new fibre_segment_ptr[newdatalen];
    double * d = new double[newdatalen];
    for (unsigned int i = 0; i<dataidx; i++) {
      nptr[i] = nodelist[i];
      d[i] = data[i];
    }
    datalen = newdatalen;
    delete[] nodelist;
    delete[] data;
    nodelist = nptr;
    data = d;
  }
  nodelist[dataidx] = s;
  data[dataidx] = t;
  dataidx++;
}

double nodegenesis_data::find_t_genesis(fibre_segment * s) {
  for (unsigned int i = 0; i<dataidx; i++) if (nodelist[i] == s) return data[i];
  return -1.0;
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void terminal_segment::Collect_Data(fibre_collected_data & fcd) {
  // Be sure to Collect_Data() from the fibre_segment tree before a call to
  // this function in order to assure valid cache values!
  // This function does some preparatory work for the
  // fibre_segment::Collect_Data() function, and some statistical data about
  // terminal segments is also collected in fibre_segment::Collect_Data()
  // to avoid the need for additional cached data.
  update_segment_vector();
  double len = Length();
  if (fcd.ispresynaptic) {
    fcd.nsb->netstats[axons_termsegs_len].Collect_Data(len);
  } else {
    fcd.nsb->netstats[dendrites_termsegs_len].Collect_Data(len);
  }
}
#endif

// [***NOTE] This could become part of (a) the initialization of terminal segment specific elongaiton models if straightforward,
//           or (b) the main function of a protocol-compliant Branch_Distribution model, either of which may be combined with
//           a phenomenologically justified manner of selecting initial elongation rates after a bifurcation (as proposed by
//           Jaap van Pelt). In an implementation as a stand-alone model, it would also be appropriate to include the ability
//           to choose a different PDF with which to select X1 and X2.
double Distribute_Elongation(terminal_segment & ts1, terminal_segment & ts2, double distribute) {
  // This determines the distribution of remaininig elongation to daughter branches
  // after a bifurcation, when the bifurcation falls within a preceding fixed
  // simulation time step.
  // [***INCOMPLETE] According to the paper, there should be a recalculation here, determining from remainder
  //                 at exactly which stimulated time the branch event took place, and what the actual remainder
  //                 should be, considering the new number of terminal segments.
  //                 If the time of the bifurcation is determined (either by computing the time from the location,
  //                 or by directly picking a time within the step and computing a location from that), then
  //                 the elongation models of ts1 and ts2 can be used to compute the actual total amount of
  //                 elongation since that time.
  // This function returns the proportion prop = X1/(X1+X2). That proportion suffices for computations, such as
  // the distribution of elongation and the balance of forces for branch angle determination, since X2/X1 = (1-prop)/prop,
  // with prop = 0 when X1 = 0. The only reason why we need to draw to random numbers here is to determine which daughter
  // branch receives the greater proportion. Even when nothing is distributed, a proportion value between 0 and 1 is returned.
  // [See TL#200709111456 and preceding entries.]
  // [UPDATE 20080522] This function now *adds* the remainder instead of setting it. This allows the distribution to be called in
  // any relative order to other initialization tasks during branching (see TL#200805190817.1).
  double X1 = X_branch.get_rand_real1(); // [***INCOMPLETE] The paper mentions use of a branch ratio PDF, PDF_br, here.
  if (distribute<=0.0) return X1;
  double X2 = X_branch.get_rand_real1();
  double X1plusX2 = X1+X2;
  if (X1plusX2==0.0) X1plusX2 = 1.0;
  double prop = X1/X1plusX2;
  double x1 = distribute*prop;
  ts1.AngularCoords().add_X(x1); // [UPDATE 20080522] Changed set_X() to add_X(). Be careful about x coordinate initialization!
  ts2.AngularCoords().add_X(distribute-x1); // [UPDATE 20080522] Changed set_X() to add_X(). Be careful about x coordinate initialization!
  return prop;
}

 bool terminal_segment::branch(PLLRoot<terminal_segment> * prevts) {
  // If prevts!=NULL then the terminal_segment object is maintained for
  // possible further branches and is linked to prevts for clean-up
  // outside the fibre_structure object.
  // The fibre_segment coordinates for this terminal_segment are updated
  // to correspond with the current angle and length of the terminal_segment
  // before branching processes are done.
  // The fibre_segment branches receive initial segment coordinates according
  // to the end point of the parent fibre_segment (in the branch() function).
  // Note that new terminal segments  are placed at the head of the list of
  // terminal segments so that they are not processed as previously existing
  // terminal segments when created during a looping growth function.
  // The initial length of a new terminal segment is determined by the
  // handle_branching() function of an Elongation Model.
  // The direction of a new terminal segment is determined by a Branch Angle
  // model.
  // The length and angle are then used to specify new terminal segments.
  // The direction of a new terminal segments may need to be modified according
  // to environment pressures. Here, this is dealt with in the constructor calls.
  // Presently, branch angles are determined relative to the parent fiber, and
  // absolute angular vectors are then derived in this function for the daughter
  // terminal segments.
  // [***NOTE] If Branch Angle Models will also include models that
  // apply effects of the environment to branching angles, which may be expressed
  // more readily in absolute coordinates, then the derivation of absolute
  // angular vectors may also be moved into specific Branch Angle Models that
  // compute parent-relative angles. Similarly, it might then be necessary to
  // obtain the ROLL angle (Stheta) from a Branch Angle Model.
#define Pitch() Theta()
 
  // 1. Determine actual bifurcation node location in fixed step approximation
  double remainder = 0.0;
  if (branchsomewhereinsegment) remainder = randomize_last_segment_elongation(Arbor()->BranchingModel()->Min_Node_Interval());
#ifdef VECTOR3D
  // 2. Compute a ROLL angle around the parent fiber segment (this terminal segment)
  double Stheta = 2.0*M_PI*X_turn.get_rand_real1();
#endif
  // 3. Prepare daughter fibre_segment objects
  update_segment_vector();
#ifdef VECTOR3D
  extra_fibre_data efd1(Stheta), efd2(Stheta+M_PI);
  if (!terminalsegment->branch(efd1,efd2)) return false;
#endif
#ifdef VECTOR2D
  if (!terminalsegment->branch()) return false;
#endif
  // 4. Add daughter terminal_segment objects to the set of terminal segments
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
  terminal_segment * t1 = new terminal_segment(*arbor,*(terminalsegment->Branch1()),centrifugalorder+1); //,acoords1);
  terminal_segment * t2 = new terminal_segment(*arbor,*(terminalsegment->Branch2()),centrifugalorder+1); //,acoords2);
  Root()->link_after(t1);
  Root()->link_after(t2);
  most_recent_branch1_ts = t1;
  most_recent_branch2_ts = t2;
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,2,sizeof(terminal_segment));
  // ALL HANDLERS FOR MODEL PROPAGATION SHOULD BE CALLED HERE (Look for similar notes such as this elsewhere, e.g. in Delayed_Branching_Model::continuation_node_to_branch().)
  // 5. Handle the bifurcation in the Elongation Model
  elmodel->handle_branching(*t1,*t2);
  // 6. Handle the bifurcation in the Elongation Rate Initialization, TSB and TST Models
  erimodel->handle_branching(*t1,*t2);
  tsbmodel->handle_branching(*t1,*t2);
  tstmodel->handle_branching(*t1,*t2);
  // 7. Handle the bifurcation in the Branch Angle Model
  bamodel->handle_branching(*t1,*t2);
  // 8. Determine the distribution of initial lengths, returned in rho (x) of t1->angularcoords and t2->angularcoords
  Distribute_Elongation(*t1,*t2,remainder); // [***NOTE] Environmental pressures are not tested here, since these are still relative angular coordinates.
  // 9. Determine the distribution of initial elongation rates, set in t1/t2->ElongationModel()->base_parameters.l_i_cache
  erimodel->initialize(t1,t2);
  // 8. Determine branching PITCH angles, returned in theta (y) of t1->angularcoords and t2->angularcoords
  bamodel->branch_angle(t1,t2);
  // 9. Handle the bifurcation in the Direction Model
  if (dirmodel) dirmodel->handle_branching(*t1,*t2); // pass on the direction model
  // 10. Compute the angular vectors of the daughter terminal segments
#ifdef VECTOR3D
  // This process is described in the figure ~/src/nnmodels/nibr/3D-branch-turn-angles.fig.
  spatial S1absoluteangular, S2absoluteangular;
  veer(this,t1->Pitch(),S1absoluteangular,Stheta); // for branch 1
  spatial S2(1.0,efd2,t2->Pitch());
  daughter_to_network_coordinate_system(S2,S2absoluteangular); // for branch 2
#endif
#ifdef VECTOR2D
  // [***NOTE] It must not be predictable, which of the two daughters will have the greater angle.
  //           Angles must both be positive or both be negative.
  spatial S1absoluteangular(0.0,Theta()+t1->Pitch());
  spatial S2absoluteangular(0.0,Theta()-t2->Pitch());
#endif
  // Collect angle sample data during branch generation
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_BRANCH_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_branch_angles,dendrites_branch_angles};
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    // [***NOTE] I need to determine which is which if I am to collect data
    // about slightly and strongly deviating daughter segments separately.
    STANDARDIZE_ANGLE_DATA(t1->Pitch());
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(t1->Pitch());
    STANDARDIZE_ANGLE_DATA(t2->Pitch());
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(t2->Pitch());
  }
#endif
#endif
  // [***INCOMPLETE] In places such as this, I may have to set
  //                 elmodel and dirmodel to NULL before remove(), just in case
  //                 a destructor is defined that would clean those up if they
  //                 were not propagated to a new terminal segment.
  // 11. Set the angle and length data in the terminal_segment branch objects, using the appropriate call functions
  S1absoluteangular.set_X(t1->Length());
  S2absoluteangular.set_X(t2->Length());
  t1->set_angularcoords(S1absoluteangular); // Here, environmental pressures are also applied.
  t2->set_angularcoords(S2absoluteangular); // Here, environmental pressures are also applied.
  t1->update_segment_vector();
  t2->update_segment_vector();
  // 12. Remove this terminal_segment object
  if (prevts) {
    isolate();
    prevts->link_before(this);
  } else {
    remove();
    DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,-1,sizeof(terminal_segment));
  }
  return true;
}

bool terminal_segment::turn(neuron * n) {
  // 1. Determine actual bifurcation node location in fixed step approximation
  double remainder = 0.0;
  if (branchsomewhereinsegment) remainder = randomize_last_segment_elongation(Arbor()->BranchingModel()->Min_Node_Interval());
  DirectionModel()->direction(this,n,remainder);
  return true;
}

#ifdef VECTOR3D
bool terminal_segment::turn_without_update_of_segment_vector(spatial & acoords, extra_fibre_data & efd) {
#endif
#ifdef VECTOR2D
bool terminal_segment::turn_without_update_of_segment_vector(spatial & acoords) {
#endif
  // create a turn in the structure
  // The direction of a new terminal segments may need to be modified according
  // to environment pressures. Here, this is dealt with in the constructor call.
#ifdef VECTOR3D
  if (!terminalsegment->turn(efd)) return false;
#endif
#ifdef VECTOR2D
  if (!terminalsegment->turn()) return false;
#endif
  // update the set of terminal segments
  // Note that they are placed at the head of the list of terminal
  // segments so that they are not processed as previously existing
  // terminal segments when created during a looping growth function.
  DIAGNOSTIC_BEFORE_ALLOCATION(new_fibre_segment);
  terminal_segment * ts = new terminal_segment(*arbor,*(terminalsegment->Branch1()),acoords,centrifugalorder);
  Root()->link_after(ts);
  // ALL HANDLERS FOR MODEL PROPAGATION SHOULD BE CALLED HERE (Look for similar notes such as this elsewhere, e.g. in Delayed_Branching_Model::continuation_node_to_branch().)
  if (dirmodel) dirmodel->handle_turning(*ts);
  elmodel->handle_turning(*ts);
  erimodel->handle_turning(*ts);
  tsbmodel->handle_turning(*ts);
  tstmodel->handle_turning(*ts);
  bamodel->handle_turning(*ts);
  remove();
  // [***NOTE] A new history element for the direction model is not accounted
  // for in the diagnostics below.
  DIAGNOSTIC_AFTER_ALLOCATION(new_fibre_segment,0,sizeof(terminal_segment));
  return true;
}

#ifdef ENABLE_FIXED_STEP_SIMULATION
double terminal_segment::randomize_last_segment_elongation(double mininterval) {
  // This function compensates for the fixed step nature of the simulation,
  // making turns/branches possible at points between steps.
  // Note that this function (a) can avoid zero-length fibre segments
  // [Update 20070509/20070624: Now Obsolete: "and
  // (b) favors fibre segment lengths near the maximum, both in accordance
  // with model presuppositions described in TL#200511221231.1.", removed
  // subtraction by product with sqrt(random_double())]
  // Updated as described in TL#200603151008.1.
  // Update 200704102307 (see related task log entries): This does not actually prevent nodes from being
  // created when the fiber length happens to be zero.
  // Update 200706240729: In combination with the new strict ts->Length() detection, this does
  // enforce the minimum node distance. Also, the selection of a point within the most recent
  // elongation is now done with a uniform random selection, since there is no reason to
  // give any part of the elongation more or less significance. This way, the combination of
  // this function with the implemented probabilistic branch or turn selection is a proper
  // approximation of the phenomenological model.
  // Note that the fixedstepelongation variable is set only when adding length, not when
  // setting terminal segment length. Consequently, initialization of terminal segments with
  // initial lengths (using the L0_min and L0_max parameters) does not affect fixedstepelongation
  // and that initial length cannot be reduced (even on the first iteration).
  // In addition to updating the terminal segment, this function returns the length that has
  // been subtracted.
  double termseglength = angularcoords.X();
  if (termseglength<=mininterval) return 0.0;
  double x = fixedstepelongation;
  if ((termseglength-x)<mininterval) x = termseglength - mininterval;
  x *= X_elongate.get_rand_real1(); // [***NOTE] Possibly replace this with a X_branchorturn in order not to affect any of the other random generators.
  angularcoords.set_X(termseglength-x);
#ifdef DEBUGGING_ELONGATION
  if (Length()<0.0) debugging.randomize_last_segment_elongation_negative++;
  if (Length()==0.0) debugging.randomize_last_segment_elongation_zero_length++;
  if ((Length()<=0.0) && (fixedstepelongation<0.0)) debugging.randomize_last_segment_elongation_fixedstepelongation_negative++;
  if ((Length()<=0.0) && (fixedstepelongation==0.0)) debugging.randomize_last_segment_elongation_fixedstepelongation_zero++;
#endif
  return x;
}
#endif // ENABLE_FIXED_STEP_SIMULATION

void structure_initialization::init(natural_schema_parent_set n, double _L0_min, double _L0_max) {
  L0_min = _L0_min;
  L0_max = _L0_max;
  colnum = definedcolors;
  definedcolors++;
  // Some specific defaults that depend on the type of fiber structure
  if (colortable) switch (n) {
  case all_axons_sps: case all_multipolar_axons_sps: case all_bipolar_axons_sps: case all_pyramidal_axons_sps: colortable->set_color(colnum,colortable->coldef(CT_axon_excitatory)); break;;
  case all_interneuron_axons_sps: colortable->set_color(colnum,colortable->coldef(CT_axon_inhibitory)); break;;
  default: colortable->set_color(colnum,colortable->coldef(CT_dendrites));
  }
}

String structure_initialization::report_parameters() {
  String res("L0_min=");
  res += String(L0_min,"%.3f")+" L0_max="+String(L0_max,"%.3f");
  res += " color(#" + String(colnum) + ")=0x";
  if (colortable) res += RGBstr(colortable->get_color(colnum));
  return res;
}

void structure_initialization::initialize(fibre_structure * fs_root) {
  if (!fs_root) {
    warning("Warning: Fiber structure root segment pointer (fs_root) is NULL in structure_initialization::initialize().\n");
    return;
  }
  if ((L0_min<0.0) || (L0_max<0.0) || (L0_max<L0_min)) error("Error: Invalid (or negative) initial segment length range [L0_min,L0_max] for a fiber structure root segment of the neuron with numerical ID "+String((long) fs_root->N()->numerical_ID())+'\n');
  fs_root->TerminalSegments()->head()->set_length(X_misc.get_rand_range_real1(L0_min,L0_max));
  if (colnum<1) warning("Warning: The color with which to display neurite fiber of the neuron with numerical ID "+String((long) fs_root->N()->numerical_ID())+" was not set.\n");
  fs_root->colnum = colnum;
}

structure_initialization * structure_initialization_selection(String & label, Command_Line_Parameters & clp, structure_initialization & superior_set, structure_initialization & strucinit) {
  // Note: The result pointer is not necessarily a pointer to strucinit,
  // since it is used here primarily to indicate whether a valid schema
  // initialization was accomplished.
  // Note: This is protocol compliant in terms of set specification, but
  // unlike specifications for models, structure initialization data is
  // allocated directly into existing array elements, instead setting
  // pointers to clones.
  int n;
  // Determine if specific structure initialization information is provided for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"L0"))>=0) {
    int L0veclen;
    double * L0vec = str2vec(clp.ParValue(n),&L0veclen);
    if (L0veclen<1) warning("Warning: Missing or non-numerical value for "+label+"L0.\n");
    else {
      strucinit.L0_min = L0vec[0];
      if (L0veclen>1) strucinit.L0_max = L0vec[1];
      else strucinit.L0_max = L0vec[0];
    }
    delete L0vec;
  }
  if (strucinit.L0_min<0.0) {// Try to obtain missing parameter values from the designated superior set
    strucinit.L0_min = superior_set.L0_min;
    strucinit.L0_max = superior_set.L0_max;
  }
  if ((n=clp.Specifies_Parameter(label+"color"))>=0) {
    if (strucinit.colnum<1) {
      strucinit.colnum = definedcolors;
      definedcolors++;
    }
    if (colortable) colortable->set_color(strucinit.colnum,strtoul(clp.ParValue(n),NULL,0));
  }
  if (strucinit.colnum<1) { // Try to obtain missing parameter values from the designated superior set
    strucinit.colnum = superior_set.colnum;
  }
  if ((strucinit.L0_min<0.0) || (strucinit.L0_max<0.0) || (strucinit.colnum<1)) return NULL;
  return &strucinit;  
}
