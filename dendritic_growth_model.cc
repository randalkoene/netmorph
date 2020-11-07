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

#define DAYSECONDS 86400.0
#define DAYSECONDSSQR (86400.0*86400.0)
#define DAYSECONDSCUB (86400.0*86400.0*86400.0)

#define TESTING_SIMPLE_RANDOMIZATION

#define INCLUDE_BRANCH_ANGLE_SAMPLING
#define INCLUDE_TURN_ANGLE_SAMPLING

#define USE_REMAINDER_FUNCTION_FOR_ANGLE_DATA_STANDARDIZATION

#ifdef LEGACY_BRANCHING_PROTOCOL
#define SELECT_TYPED_BRANCH(two_branch_types) \
{ \
  if (!two_branch_types) branch(ts,branchanglemin,branchanglemax,branchangleothermax); /* main-daughter branching */ \
  else { /* if this is the only terminal segment */ \
    if ((ps->TerminalSegments()->head()==ps->TerminalSegments()->tail()) && (X_branch.get_rand_real1() <= branchIIprobability)) branch(ts,branchIIanglemin,branchIIanglemax,branchIIangleothermax); /* fork branch */ \
    else branch(ts,branchanglemin,branchanglemax,branchangleothermax); /* main-daughter branching */ \
  } \
}
#endif

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
double turn_a, turn_max_y, turn_SQb, turn_negb, turn_twoa;
#ifdef VECTOR2D
double turn_twomax_y;
#endif

#ifdef TEST_FOR_NAN
unsigned long int DBcounter = 0;
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
    cout << "Warning: fibre_segment was terminal_segment in Delayed_Branching_Model::continuation_node_to_branch()\n";
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
    cout << "Warning: Unable to create branch fibre_segment in Delayed_Branching_Model::continuation_node_to_branch()\n";
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
  terminal_segment * ts = new terminal_segment(*arbor,*newbranch,acoords);
  register PLLRoot<terminal_segment> * tsroot = arbor->TerminalSegments();
  register direction_model_base * dirmodel = tsroot->head()->DirectionModel();
  // ALL HANDLERS FOR MODEL PROPAGATION SHOULD BE CALLED HERE (Look for similar notes such as this elsewhere, e.g. in terminal_segment::branch().)
  if (dirmodel) ts->set_direction_model(dirmodel->clone());
  ts->set_elongation_model(tsroot->head()->ElongationModel()->clone());
  ts->set_ERI_model(tsroot->head()->ElongationRateInitializationModel()->clone());
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

void van_Pelt_dendritic_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("Dgrowth_E"))>=0) _E = atof(clp.ParValue(n));
#ifdef STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION
  if ((n=clp.Specifies_Parameter("Dgrowth_tau"))>=0) { _tau = atof(clp.ParValue(n)); _c = _B_inf / _tau; }
#else
  if ((n=clp.Specifies_Parameter("Dgrowth_tau"))>=0) _tau = atof(clp.ParValue(n));
#endif
  //if ((n=clp.Specifies_Parameter("Dgrowth_n_inf"))>=0) _n_inf = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_B_inf"))>=0) { _B_inf = atof(clp.ParValue(n)); _c = _B_inf / _tau; }
  if ((n=clp.Specifies_Parameter("Dgrowth_c"))>=0) { _c = atof(clp.ParValue(n)); _B_inf = _c*_tau; }
  //if ((n=clp.Specifies_Parameter("Dgrowth_nu0"))>=0) _nu0 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Dgrowth_F"))>=0) _F = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_L0_min"))>=0) _L0_min = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_L0_max"))>=0) _L0_max = atof(clp.ParValue(n));
#ifdef LEGACY_BRANCHING_PROTOCOL
  if ((n=clp.Specifies_Parameter("Dmin_node_interval"))>=0) min_node_interval = atof(clp.ParValue(n));
#endif
  if ((n=clp.Specifies_Parameter("Dcompetition_with_all_trees"))>=0) {
    if (downcase(clp.ParValue(n))==String("true")) competition_with_all_trees = gcc_within_all_AorD;
    else competition_with_all_trees = gcc_within_tree;
  }
  if ((n=clp.Specifies_Parameter("Dcompetition_with_whole_neuron"))>=0) if (downcase(clp.ParValue(n))==String("true")) competition_with_all_trees = gcc_whole_neuron;
  set_fixed_time_step_size(_dt);
}

String van_Pelt_dendritic_growth_model::local_report_parameters() {
  String res(_E,"%.3f _tau = ");
  res += String(_tau,"%.3f _B_inf = ") + String(_B_inf,"%.3f _c = ") + String(_c,"%.3f _L0_min = ") + String(_L0_min,"%.3f _L0_max = ") + String(_L0_max,"%.3f min_node_interval = ") + String(_c_dt,"%f)");
  // + String(_nu0,"%f.3 _F = ") + String(_F,"%f.3 ")
  res += " E competes ";
  switch (competition_with_all_trees) {
  case gcc_whole_neuron: res += "with whole neuron"; break;;
  case gcc_within_all_AorD: res += "with all axons or with all dendrites"; break;;
  default: res += "within the axon/dendrite arbor"; // gcc_within_tree
  };
  return res;
}

String van_Pelt_dendritic_growth_model::report_parameters() {
  return "van_Pelt_dendritic_growth_model: _E = " + local_report_parameters() + '\n';
}

van_Pelt_dendritic_growth_model::van_Pelt_dendritic_growth_model(double E, double tau, double B_inf): _E(E), _tau(tau), _B_inf(B_inf), _c(B_inf/tau), _L0_min(9.0), _L0_max(11.0),
#ifdef LEGACY_BRANCHING_PROTOCOL
min_node_interval(1.0),
#endif
competition_with_all_trees(gcc_within_tree) {
  set_fixed_time_step_size(van_Pelt_SUGGESTED_dt);
  e_1 = exp(-1.0);
}

van_Pelt_dendritic_growth_model::van_Pelt_dendritic_growth_model(Command_Line_Parameters & clp): _E(0.051), _tau(319680.0), _B_inf(1.124), _L0_min(9.0), _L0_max(11.0),
#ifdef LEGACY_BRANCHING_PROTOCOL
min_node_interval(1.0),
#endif
competition_with_all_trees(gcc_within_tree) {
  _c = _B_inf / _tau;
  e_1 = exp(-1.0);
  _dt = van_Pelt_SUGGESTED_dt;
  parse_CLP(clp);
}

void van_Pelt_dendritic_growth_plus_turns_model::local_parse_CLP(Command_Line_Parameters & clp) {
  int n;
  // apical dendrite branching model parameters
  if ((n=clp.Specifies_Parameter("Dgrowth_E_apicaltrunc"))>=0) _E_apicaltrunc = atof(clp.ParValue(n));
#ifdef STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION
  if ((n=clp.Specifies_Parameter("Dgrowth_tau_apicaltrunc"))>=0) { _tau_apicaltrunc = atof(clp.ParValue(n)); _c_apicaltrunc = _B_inf_apicaltrunc / _tau_apicaltrunc; }
#else
  if ((n=clp.Specifies_Parameter("Dgrowth_tau_apicaltrunc"))>=0) _tau_apicaltrunc = atof(clp.ParValue(n));
#endif
  if ((n=clp.Specifies_Parameter("Dgrowth_B_inf_apicaltrunc"))>=0) { _B_inf_apicaltrunc = atof(clp.ParValue(n)); _c_apicaltrunc = _B_inf_apicaltrunc / _tau_apicaltrunc; }
  if ((n=clp.Specifies_Parameter("Dgrowth_c_apicaltrunc"))>=0) { _c_apicaltrunc = atof(clp.ParValue(n)); _B_inf_apicaltrunc = _c_apicaltrunc*_tau_apicaltrunc; }
  if ((n=clp.Specifies_Parameter("Dgrowth_E_apicaltuft"))>=0) _E_apicaltuft = atof(clp.ParValue(n));
#ifdef STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION
  if ((n=clp.Specifies_Parameter("Dgrowth_tau_apicaltuft"))>=0) { _tau_apicaltuft = atof(clp.ParValue(n)); _c_apicaltuft = _B_inf_apicaltuft / _tau_apicaltuft; }
#else
  if ((n=clp.Specifies_Parameter("Dgrowth_tau_apicaltuft"))>=0) _tau_apicaltuft = atof(clp.ParValue(n));
#endif
  if ((n=clp.Specifies_Parameter("Dgrowth_B_inf_apicaltuft"))>=0) { _B_inf_apicaltuft = atof(clp.ParValue(n)); _c_apicaltuft = _B_inf_apicaltuft / _tau_apicaltuft; }
  if ((n=clp.Specifies_Parameter("Dgrowth_c_apicaltuft"))>=0) { _c_apicaltuft = atof(clp.ParValue(n)); _B_inf_apicaltuft = _c_apicaltuft*_tau_apicaltuft; }
  if ((n=clp.Specifies_Parameter("apical_tufting_trigger_time"))>=0) apical_tufting_trigger_time = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("apical_tufting_trigger_distance"))>=0) {
    apical_tufting_trigger_distance = atof(clp.ParValue(n));
    apical_tufting_trigger_distance *= apical_tufting_trigger_distance; // squared for use with distance2()
  }
  //_c_dt_apicaltrunc = 1.0 - exp(-_dt/_tau_apicaltrunc); // this is akin to the set_fixed_time_step_size() call above
  //_c_dt_apicaltuft = 1.0 - exp(-_dt/_tau_apicaltuft); // this is akin to the set_fixed_time_step_size() call above
  // Updated to match the t-dt to t interval order of the paper
  _c_dt_apicaltrunc = exp(_dt/_tau_apicaltrunc)-1.0; // this is akin to the set_fixed_time_step_size() call above
  _c_dt_apicaltuft = exp(_dt/_tau_apicaltuft)-1.0; // this is akin to the set_fixed_time_step_size() call above
  // direction change model parameters
  if ((n=clp.Specifies_Parameter("Dgrowth.turn_model"))>=0) { // +++++
    if (downcase(clp.ParValue(n))==String("branch_coupled")) turnmodel = tm_branch_coupled;
    else if (downcase(clp.ParValue(n))==String("linear_rate")) turnmodel = tm_linear_rate;
    else if (downcase(clp.ParValue(n))==String("each_dt")) turnmodel = tm_each_dt;
    else cout << "Warning: Unknown turn model '" << clp.ParValue(n) << "'\n";
  }
  if ((n=clp.Specifies_Parameter("Dgrowth_tf"))>=0) _turnfactor = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dturn_rate"))>=0) turn_rate = atof(clp.ParValue(n)); // +++++
  if ((n=clp.Specifies_Parameter("Dturn_separation"))>=0) turn_rate = 1.0/atof(clp.ParValue(n)); // see explanation in DirectionChange_Threshold()
  if ((n=clp.Specifies_Parameter("Dbranchesatturns"))>=0) branchesatturns = (downcase(clp.ParValue(n))==String("true"));
  double max_x = turnanglemax - turnanglemin;
  if (max_x<=0.0) error("Error: turnanglemax must be greater than turnanglemin\n");
  double b = 1.0;
  turn_a = -b/max_x;
  turn_max_y = 0.5*b*max_x;
  turn_SQb = b*b;
  turn_negb = -b;
  turn_twoa = 2.0*turn_a;
#ifdef VECTOR2D
  turn_twomax_y = 2.0*turn_max_y;
#endif
}

void van_Pelt_dendritic_growth_plus_turns_model::parse_CLP(Command_Line_Parameters & clp) {
  van_Pelt_dendritic_growth_model::parse_CLP(clp);
  local_parse_CLP(clp);
}

String van_Pelt_dendritic_growth_plus_turns_model::local_report_parameters() {
  String res(van_Pelt_dendritic_growth_model::local_report_parameters() + String(_turnfactor," _turnfactor = %.3f")
	     + String(_E_apicaltrunc," _E_apicaltrunc = %f")
	     + String(_tau_apicaltrunc," _tau_apicaltrunc = %f")
	     + String(_B_inf_apicaltrunc," _B_inf_apicaltrunc = %f")
	     + String(_c_apicaltrunc," _c_apicaltrunc = %f")
	     + String(_E_apicaltuft," _E_apicaltuft = %f")
	     + String(_tau_apicaltuft," _tau_apicaltuft = %f")
	     + String(_B_inf_apicaltuft," _B_inf_apicaltuft = %f")
	     + String(_c_apicaltuft," _c_apicaltuft = %f")
	     + " (_c_dt_apicaltrunc=" + String(_c_dt_apicaltrunc,"%f,")
	     + " _c_dt_apicaltuft=" + String(_c_dt_apicaltuft,"%f)"));
  if (apical_tufting_trigger_time<DBL_MAX) res += " apical_tufting_trigger_time = " + String(apical_tufting_trigger_time,"%.3f");
  else res += " apical_tufting_trigger_time = NONE" ;
  if (apical_tufting_trigger_distance>0.0) res += " apical_tufting_trigger_distance = " + String(sqrt(apical_tufting_trigger_distance),"%.3f");
  else res += " apical_tufting_trigger_distance = NONE" ;
  switch (turnmodel) { // +++++
  case tm_branch_coupled: res += " turnmodel = branch_coupled"; break;;
  case tm_linear_rate: res += " turnmodel = linear_rate"; break;;
  default: res += " turnmodel = each_dt";
  }
  res += String(turn_rate," turn_rate = %f.3"); // +++++
  return res;
}

String van_Pelt_dendritic_growth_plus_turns_model::report_parameters() {

  return "van_Pelt_dendritic_growth_plus_turns_model: _E = " + local_report_parameters() + '\n';
}

van_Pelt_dendritic_growth_plus_turns_model::van_Pelt_dendritic_growth_plus_turns_model(double E, double tau, double B_inf, double E_apicaltrunc, double tau_apicaltrunc, double B_inf_apicaltrunc, double E_apicaltuft, double tau_apicaltuft, double B_inf_apicaltuft, double tufttriggertime, double tufttriggerdistance, turn_model_type tm, double tf, double tr, bool bt): van_Pelt_dendritic_growth_model(E,tau,B_inf), _E_apicaltrunc(E_apicaltrunc), _tau_apicaltrunc(tau_apicaltrunc), _B_inf_apicaltrunc(B_inf_apicaltrunc), _E_apicaltuft(E_apicaltuft), _tau_apicaltuft(tau_apicaltuft), _B_inf_apicaltuft(B_inf_apicaltuft), apical_tufting_trigger_time(tufttriggertime), apical_tufting_trigger_distance(tufttriggerdistance), turnmodel(tm), _turnfactor(tf), turn_rate(tr), branchesatturns(bt) {
  _c_dt_apicaltrunc = exp(_dt/_tau_apicaltrunc)-1.0; // this is akin to the set_fixed_time_step_size() call above
  _c_dt_apicaltuft = exp(_dt/_tau_apicaltuft)-1.0; // this is akin to the set_fixed_time_step_size() call above
}

van_Pelt_dendritic_growth_plus_turns_model::van_Pelt_dendritic_growth_plus_turns_model(Command_Line_Parameters & clp): van_Pelt_dendritic_growth_model(clp), _E_apicaltrunc(0.051), _tau_apicaltrunc(400000.0), _B_inf_apicaltrunc(2.0), _E_apicaltuft(0.051), _tau_apicaltuft(400000.0), _B_inf_apicaltuft(5.0), apical_tufting_trigger_time(DBL_MAX), apical_tufting_trigger_distance(-1.0), turnmodel(tm_linear_rate), _turnfactor(9.25), turn_rate(0.1), branchesatturns(false) {
  _c_apicaltrunc = _B_inf_apicaltrunc / _tau_apicaltrunc;
  _c_apicaltuft = _B_inf_apicaltuft / _tau_apicaltuft;
  local_parse_CLP(clp);
}

void Polynomial_O1_dendritic_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  van_Pelt_dendritic_growth_plus_turns_model::local_parse_CLP(clp);
  if ((n=clp.Specifies_Parameter("Dgrowth_E"))>=0) _E = atof(clp.ParValue(n));
  //double n_inf, if ((n=clp.Specifies_Parameter("Dgrowth_n_inf"))>=0) _n_inf = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_c"))>=0) _c = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Dgrowth_nu0"))>=0) _nu0 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Dgrowth_F"))>=0) _F = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_L0_min"))>=0) _L0_min = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_L0_max"))>=0) _L0_max = atof(clp.ParValue(n));
  set_fixed_time_step_size(_dt);
}

String Polynomial_O1_dendritic_growth_model::report_parameters() {
  String res("Polynomial_O1_dendritic_growth_model: _E = ");
  res += String(_E,"%f.3 _c = ") + String(_c,"%f.3 _L0_min = ") + String(_L0_min,"%f.3 _L0_max = ") + String(_L0_max,"%f.3 _turnfactor = ") + String(_turnfactor,"%f.3\n");
  // + String(_nu0,"%f.3 _F = ") + String(_F,"%f.3 _L0_min = ")
  return res;
}

Polynomial_O1_dendritic_growth_model::Polynomial_O1_dendritic_growth_model(double E, double c, double tf): van_Pelt_dendritic_growth_plus_turns_model(E,1.0,1.0,E,1.0,1.0,E,1.0,1.0,DBL_MAX,-1.0,tm_each_dt,tf,0.05) {
  _c = c;
  _L0_min = 0.0;
  _L0_max = 0.0;
  set_fixed_time_step_size(van_Pelt_SUGGESTED_dt);
  e_1 = exp(-1.0);
}

Polynomial_O1_dendritic_growth_model::Polynomial_O1_dendritic_growth_model(Command_Line_Parameters & clp) {
  // default values are obtained by polynomial approximation of mirrored data
  _E = 1.0; _c = 0.156204/DAYSECONDS; _L0_min = 0.0; _L0_max = 0.0; _turnfactor = 9.25;
  _dt = van_Pelt_SUGGESTED_dt;
  e_1 = exp(-1.0);
  parse_CLP(clp);
}

void Polynomial_O3_dendritic_growth_model::local_parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("Dgrowth_c1"))>=0) _c1 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Dgrowth_c2"))>=0) _c2 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Dgrowth_nu1"))>=0) _nu1 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Dgrowth_nu2"))>=0) _nu2 = atof(clp.ParValue(n));
}

void Polynomial_O3_dendritic_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  Polynomial_O1_dendritic_growth_model::parse_CLP(clp);
  local_parse_CLP(clp);
  set_fixed_time_step_size(_dt);
}

String Polynomial_O3_dendritic_growth_model::report_parameters() {
  String res("Polynomial_O3_dendritic_growth_model: _E = ");
  res += String(_E,"%f.3 _c = ") + String(_c,"%f.3 _c1 = ") + String(_c1,"%f.3 _c2 = ") + String(_c2,"%f.3 _L0_min = ") + String(_L0_min,"%f.3 _L0_max = ") + String(_L0_max,"%f.3 _turnfactor = ") + String(_turnfactor,"%f.3\n");
  // + String(_nu0,"%f.3 _nu1 = ") + String(_nu1,"%f.3 _nu2 = ") + String(_nu2,"%f.3 _F = ") + String(_F,"%f.3 _L0_min = ") 
  return res;
}

Polynomial_O3_dendritic_growth_model::Polynomial_O3_dendritic_growth_model(double E, double c, double c1, double c2, double tf): Polynomial_O1_dendritic_growth_model(E,c,tf), _c1(c1), _c2(c2) {
  _L0_min = 0.0; // set by L(0) cubic approximation of dendrite length for mirrored data
  _L0_max = 0.0;
}

Polynomial_O3_dendritic_growth_model::Polynomial_O3_dendritic_growth_model(Command_Line_Parameters & clp): Polynomial_O1_dendritic_growth_model(clp), _c1(0.0/DAYSECONDSSQR), _c2(0.000039/DAYSECONDSCUB) {
  _c = 0.145511/DAYSECONDS;
  local_parse_CLP(clp);
}

double van_Pelt_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  //double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
  return _c*_tau*(exp(-t1/_tau)-exp(-t2/_tau)); //*nt_E;
}

double Polynomial_O1_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  //double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
  return _c*(t2-t1);//*nt_E; // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c =  0.145156 / (24.0*60.0*60.0))
}

double Polynomial_O3_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  double t1sq = t1*t1, t2sq = t2*t2;
  return (_c*(t2-t1)) + (_c1*(t2sq-t1sq)) + (_c2*((t2sq*t2)-(t1sq*t1))); // set the c parameters by dividing daily fits to Ger Ramaker's data by (24.0*60.0*60.0)
}

void van_Pelt_dendritic_growth_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
  //_c_dt = 1.0 - exp(-dt/_tau);
  // Updated to match the t-dt to t interval as in the paper
  _c_dt = exp(dt/_tau)-1.0;
  //_dt_nu0 = dt*_nu0;
}

void van_Pelt_dendritic_growth_plus_turns_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
  //_c_dt = 1.0 - exp(-dt/_tau);
  //_c_dt_apicaltrunc = 1.0 - exp(-dt/_tau_apicaltrunc); // this is akin to the set_fixed_time_step_size() call above
  //_c_dt_apicaltuft = 1.0 - exp(-dt/_tau_apicaltuft); // this is akin to the set_fixed_time_step_size() call above
  // Updated to match the t-dt to t interval as in the paper
  _c_dt = exp(dt/_tau)-1.0;
  _c_dt_apicaltrunc = exp(dt/_tau_apicaltrunc)-1.0; // this is akin to the set_fixed_time_step_size() call above
  _c_dt_apicaltuft = exp(dt/_tau_apicaltuft)-1.0; // this is akin to the set_fixed_time_step_size() call above
  //_dt_nu0 = dt*_nu0;
}

void Polynomial_O1_dendritic_growth_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
  _c_dt = _c*dt;
  //_dt_nu0 = dt*_nu0;
}

double van_Pelt_dendritic_growth_model::expected_number_of_branches(neuron * n, double t) {
  // Returns the integral of p_s(t) over the time interval _dt using an
  // n(t)^(-E) that counts all terminal segments of all dendritic arbors of
  // neuron n if competition_with_all_trees==gcc_within_all_AorD. The competition does not
  // take into account the value of n(t) after the update interval, i.e.
  // assumes that it is a constant.
  // Otherwise, returns the integral of the baseline component D(t) over the
  // time interval _dt, leaving it up to the calling function to multiply
  // with an appropriate n(t) dependence.
  double res = _B_inf*_c_dt*exp(-t/_tau);
  if (competition_with_all_trees==gcc_within_tree) return res;
  if (competition_with_all_trees==gcc_within_all_AorD) {
    double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
    return res*nt_E;
  }
  double nt_E = exp(-_E*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
  return res*nt_E;
}

double van_Pelt_dendritic_growth_plus_turns_model::expected_number_of_branches_APICAL(neuron * n, double t) {
  // A variant of the above that uses apical dendrite specific parameters and
  // switches them in accordance with a development time trigger.
  if (t>=apical_tufting_trigger_time) {
    double res = _c_apicaltuft*_tau_apicaltuft*_c_dt_apicaltuft*exp(-t/_tau_apicaltuft);
    if (competition_with_all_trees==gcc_within_tree) return res;
    if (competition_with_all_trees==gcc_within_all_AorD) {
      double nt_E = exp(-_E_apicaltuft*log((double) n->total_input_terminal_segments()));
      return res*nt_E;
    }
    double nt_E = exp(-_E_apicaltuft*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
    return res*nt_E;
  } else {
    double res = _c_apicaltrunc*_tau_apicaltrunc*_c_dt_apicaltrunc*exp(-t/_tau_apicaltrunc);
    if (competition_with_all_trees==gcc_within_tree) return res;
    if (competition_with_all_trees==gcc_within_all_AorD) {
      double nt_E = exp(-_E_apicaltrunc*log((double) n->total_input_terminal_segments()));
      return res*nt_E;
    }
    double nt_E = exp(-_E_apicaltrunc*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
    return res*nt_E;
  }
}

double van_Pelt_dendritic_growth_plus_turns_model::expected_number_of_branches_APICAL_TRUNC(neuron * n, double t) {
  // A variant of the above that uses apical dendrite TRUNC specific parameters.
  double res = _c_apicaltrunc*_tau_apicaltrunc*_c_dt_apicaltrunc*exp(-t/_tau_apicaltrunc);
  if (competition_with_all_trees==gcc_within_tree) return res;
  if (competition_with_all_trees==gcc_within_all_AorD) {
    double nt_E = exp(-_E_apicaltrunc*log((double) n->total_input_terminal_segments()));
    return res*nt_E;
  }
  double nt_E = exp(-_E_apicaltrunc*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
  return res*nt_E;
}

double van_Pelt_dendritic_growth_plus_turns_model::expected_number_of_branches_APICAL_TUFT(neuron * n, double t) {
  // A variant of the above that uses apical dendrite TUFT specific parameters.
  double res = _c_apicaltuft*_tau_apicaltuft*_c_dt_apicaltuft*exp(-t/_tau_apicaltuft);
  if (competition_with_all_trees==gcc_within_tree) return res;
  if (competition_with_all_trees==gcc_within_all_AorD) {
    double nt_E = exp(-_E_apicaltuft*log((double) n->total_input_terminal_segments()));
    return res*nt_E;
  }
  double nt_E = exp(-_E_apicaltuft*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
  return res*nt_E;
}

double Polynomial_O1_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1) {
  // [*** NOTE] The O3 polynomial fit to bifurcation in Ger Ramakers' data
  // does not use the t1 parameter.
  //double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
  return _c_dt;//*nt_E; // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c =  0.145156 / (24.0*60.0*60.0))
}

double Polynomial_O3_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  return expected_number_of_branches(n,t1,t1+_dt);
}

double van_Pelt_dendritic_growth_plus_turns_model::DirectionChange_Threshold(terminal_segment * ts) { // +++++
  // Returns the probability that a turn happened in the most recently elongated piece of
  // a neurite. This depends on the turning rate in expected number of turns per micron.
  // If the elongated length is not cached, then the function interprets the turning rate
  // as the expected number of turns per second.
  // the pointer ts must be non-NULL.
  // When we are able to calculate the probability of a turn as a function of elongated
  // length, i.e. when ENABLE_FIXED_STEP_SIMULATION is defined, then the following
  // equations hold: The probability of a turn is equal to the turn rate expressed as
  // the average number of turns per unit length multiplied by the elongated length.
  // Also, the probability of a turn is equal to the elongated length divided by the
  // average elongated length per turn (the average distance between turns).
  // The first way is specified by A/Dturn_rate, the second way by A/Dturn_separation.
  // If turn_rate*len = len / turn_separation, then turn_rate = 1 / turn_separation.
#ifdef ENABLE_FIXED_STEP_SIMULATION
  double P_turn = turn_rate*ts->elongated_length();
#else
  double P_turn = turn_rate*_dt;
#endif
  if (P_turn>=1.0) cout << "Warning: One or more turning points are expected (dendrite P_turn=" << String(P_turn,"%.2f") << "), perhaps development time steps should be smaller\n";
  return P_turn;
} // +++++

void van_Pelt_dendritic_growth_model::initialize(network * net, Connection_Statistics_Root * cstats, double & t) {
  if (!net) error("nibr Error: No network in van_Pelt_dendritic_growth_model::initialize()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_dendritic_growth_model::initialize()\n");
  if (outattr_show_progress) { cout << "Initializing first dendrite segments...\n"; cout.flush(); }
  // initialize the fibre structures
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->initialize_input_structure(_L0_min,_L0_max);
  // grow the fibre structures
  firstiteration = true;
}

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

#ifdef LEGACY_BRANCHING_PROTOCOL

void van_Pelt_dendritic_growth_model::branch(terminal_segment * ts, double anglemin, double anglemax, double othermax) {
  /* Documentation:
     The statistical data collected below (sampled_output) is about the
     veer angles (phi) from the parent direction, as measured in the plane
     that contains the parent and daughter segments.
   */
#ifdef VECTOR3D
  // This process is described in the figure
  // ~/src/nnmodels/nibr/3D-branch-turn-angles.fig.
  // Branch 1
  spatial S1xyz;
  double randomangle = X_branch.get_rand_range_real1(anglemin,anglemax);
  double rollangle = veer(ts,randomangle,S1xyz);
  // Branch 2
  // Opposite side of w in terms of theta, phi set by constrained randomness
  extra_fibre_data efd1(rollangle), efd2(rollangle+M_PI);
  double otherrandomangle = othermax*X_branch.get_rand_real1();
  spatial S2(1.0,efd2,otherrandomangle), S2xyz;
  daughter_to_network_coordinate_system(S2,S2xyz);
  ts->branch(NULL,S1xyz,S2xyz,efd1,efd2);
#endif
#ifdef VECTOR2D
  double anglespread = anglemax - anglemin;
  double randomangle = anglespread*((2.0*X_branch.get_rand_real1()) - 1.0);
  double otherrandomangle = othermax*X_branch.get_rand_real1();
  if (randomangle<0.0) randomangle -= anglemin;
  else {
    randomangle += anglemin;
    otherrandomangle = -otherrandomangle;
  }
  spatial S1xyz(0.0,ts->AngularCoords().Y()+otherrandomangle); // initial length zero
  spatial S2xyz(0.0,ts->AngularCoords().Y()+randomangle); // initial length zero
  ts->branch(NULL,S1xyz,S2xyz);
#endif
  // Collect angle sample data during branch generation
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_BRANCH_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_branch_angles,dendrites_branch_angles};
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    // [***NOTE] I need to determine which is which if I am to collect data
    // about slightly and strongly deviating daughter segments separately.
    STANDARDIZE_ANGLE_DATA(otherrandomangle);
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(otherrandomangle);
    STANDARDIZE_ANGLE_DATA(randomangle);
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(randomangle);
  }
#endif
#endif
}

#endif // LEGACY_BRANCHING_PROTOCOL

// Some functions that simplify dendrite specific elongation requests made in
// multiple growth models below:
void van_Pelt_dendritic_growth_model::request_elongation(neuron * e) {
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
#ifdef ELONGATION_EVENT_PREDICTION
      ps->devactivity->queue(_arbor_development_activity_elongation);
#else
      ps->batch_elongate();
#endif
    }
}

// Growth "batch" functions for dendrites:

void van_Pelt_dendritic_growth_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // growing in fixed time steps
  // Note that this is actually a growth "batch" function, as it directly
  // calls branch model functions and requests elongation updates for all
  // terminal segments.
  growing_dendrites = true;
  if (!net) error("nibr Error: No network in van_Pelt_dendritic_growth_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_dendritic_growth_model::grow()\n");
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
  double Bt = 0.0, intDt;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    //processing_neuron = e;
    request_elongation(e); // See TL#200603151008.1.
    fibre_structure * apical_dendrite = NULL;
    if (e->TypeID()==PYRAMIDAL) apical_dendrite = e->InputStructure()->tail(); // an apical dendrite to be treated differently
    // compute probability of terminal segment branching in this time step
    intDt = expected_number_of_branches(e,t);
    if (competition_with_all_trees!=gcc_within_tree) Bt = intDt; // See expected_number_of_branches(e,t) (the e,t parameter one, not the e,t1,t2 parameter one).
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Bt = 1.0; // insure branching
      else if (competition_with_all_trees==gcc_within_tree) Bt = intDt*exp(-_E*log((double) ps->count_terminal_segments()));
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	//if ((!strict_node_separation) || (ts->Length()>0.0)) {
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  if (X_branch.get_rand_real1()<=Bt) {
	    fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
	    double remainder = 0.0;
#ifdef ENABLE_FIXED_STEP_SIMULATION
	    // Since the equation for Bt gives us the probability or
	    // number of branches that should happen up to this time in the
	    // growth simulation, I assume that branching must happen at
	    // the latest at that time, but could happen at any time before
	    // up to the initiation time of the most recent time step.
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	    if (branchsomewhereinsegment) {
	      if (firstiteration) {
		if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
	      } else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	    }
#else
	    if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	    SELECT_TYPED_BRANCH(dendrite_two_branch_types);
#else
	    ts->branch(NULL);
#endif
	    if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_model::grow() (at t=%.1f)\n"));
	    if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	  }
	}
      }
    }
    // elongate terminal segments
    // This elongates terminal segments for the next iteration of branching
    // evaluations.
    // [***MOVED] request_elongation(e); See TL#200603151008.1.
  }
  firstiteration = false;
}

void turn_free(terminal_segment * ts) {
  // This function applies when the Samsonovich hypothesis for dendritic turns
  // is not used.
  // This is one of the functions that may be used when the chosen direction
  // model is "legacy_dm".
  double dremainder = 0.0;
  if (branchsomewhereinsegment) dremainder = ts->randomize_last_segment_elongation(ts->Arbor()->BranchingModel()->Min_Node_Interval());
  ts->update_segment_vector(); // *** perhaps this can be done in ts->turn().
#ifdef VECTOR3D
  spatial S1xyz;
  // See TL#200511261312.2.
  double y = turn_max_y*X_turn.get_rand_real1();
  double randomangle = (turn_negb + sqrt(turn_SQb + (turn_twoa*y)))/turn_a;
  randomangle += turnanglemin;
  //double randomangle = random_double(turnanglemin,turnanglemax);
  double rollangle = veer(ts,randomangle,S1xyz);
  extra_fibre_data efd(rollangle);
  S1xyz.set_X(dremainder);
  ts->turn_without_update_of_segment_vector(S1xyz,efd);
#endif
#ifdef VECTOR2D
  // See TL#200511261312.2.
  double r = turn_twomax_y*X_turn.get_rand_real1();
  double y = r;
  if (r>turn_max_y) y -= turn_max_y;
  double randomangle = (turn_negb + sqrt(turn_SQb + (turn_twoa*y)))/turn_a;
  if (r>turn_max_y) randomangle += turnanglemin;
  else randomangle = -turnanglemin - randomangle;
  //double anglespread = turnanglemax - turnanglemin;
  //double randomangle = anglespread*((2.0*random_double()) - 1.0);
  //if (randomangle<0.0) randomangle -= turnanglemin;
  //else randomangle += turnanglemin;
  spatial S1xyz(dremainder,ts->AngularCoords().Y()+randomangle);
  ts->turn_without_update_of_segment_vector(S1xyz);
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_TURN_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_turn_angles,dendrites_turn_angles};
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    STANDARDIZE_ANGLE_DATA(randomangle);
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(randomangle);
  }
#endif
#endif
}

//void van_Pelt_dendritic_growth_plus_turns_model::turn(terminal_segment * ts, neuron * e) {
void turn(terminal_segment * ts, neuron * e) {
  // The normal from the neuron to P1 determines the preferred angles
  // This is one of the functions that may be used when the chosen direction
  // model is "legacy_dm".
  double dremainder = 0.0;
  if (branchsomewhereinsegment) dremainder = ts->randomize_last_segment_elongation(ts->Arbor()->BranchingModel()->Min_Node_Interval());
  ts->update_segment_vector(); // done here, so don't repeat at ts->turn...
  spatial normal(ts->TerminalSegment()->P1); normal -= e->Pos();
  normal.convert_to_spherical();
#ifdef VECTOR3D
  // *** An alternative here is to add a small random deviation directly to
  // the angles of the normal and to use the resulting angles, even though
  // it is not quite the same as a constrained random deviation with any
  // orientation around the normal as computed below.
  double Stheta = 2.0*M_PI*X_turn.get_rand_real1();
  spatial S1(1.0,Stheta,turnmaxdeviationfromnormal*X_turn.get_rand_real1()), S1xyz;
  prepare_parent_coordinate_system(normal);
  daughter_to_network_coordinate_system(S1,S1xyz);
  // [***INCOMPLETE] It may be that this Stheta is not directly usable for
  // later conversion of turns to branches, as it is derived for a deviation
  // from the normal vector from the cell body.
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_TURN_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_turn_angles,dendrites_turn_angles};
  double angletoparent = 0.0;
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) { // calculate the angle with the parent
    // 1) parent and daughter relative to the origin
    spatial p(ts->AngularCoords()), d(S1xyz);
    // 2) unit vectors
    p.set_X(1.0); d.set_X(1.0);
    p.convert_from_spherical(); d.convert_from_spherical();
    // 3) dot product to obtain cosine of angle between vectors
    angletoparent = p * d;
    // 4) inverse cosine to obtain angle
    // [***NOTE] Beware of possible values outside the -1.0 to 1.0 range!
    angletoparent = acos(angletoparent);
  }
#endif
#endif
  extra_fibre_data efd(Stheta);
  S1xyz.set_X(dremainder);
  ts->turn_without_update_of_segment_vector(S1xyz,efd);
#endif
#ifdef VECTOR2D
  double absoluteangle = normal.Y()+(turnmaxdeviationfromnormal*((2.0*X_turn.get_rand_real1()) - 1.0));
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_TURN_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_turn_angles,dendrites_turn_angles};
  double angletoparent = absoluteangle - ts->AngularCoords().Y();
#endif
#endif
  normal.set_Y(absoluteangle);
  normal.set_X(dremainder);
  ts->turn_without_update_of_segment_vector(normal);
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_TURN_ANGLE_SAMPLING
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    STANDARDIZE_ANGLE_DATA(angletoparent);
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(angletoparent);
  }
#endif
#endif
}

void van_Pelt_dendritic_growth_plus_turns_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // <A NAME="van_Pelt_dendritic_growth_plus_turns_model::grow"> </A>
  // growing in fixed time steps
  // this is based on the growth function of the van_Pelt_dendritic_growth_model,
  // but includes turns at continuation points in dendrites.
  // Note that the reduce_rand_calls option is not yet included, although
  // one may be added if a clear benefit is apparent.
  growing_dendrites = true;
  if (!net) error("nibr Error: No network in van_Pelt_dendritic_growth_plus_turns_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_dendritic_growth_plus_turns_model::grow()\n");
  double Bt = 0.0, intDt, Bt_plus_turns = 0.0, intDt_plus_turns, nt_E;
  double BASALBt = 0.0, BASALintDt = 0.0, BASALBt_plus_turns = 0.0, BASALintDt_plus_turns;
  double APICALBt = 0.0, APICALintDt = 0.0, APICALBt_plus_turns = 0.0, APICALintDt_plus_turns = 0.0, APICAL_E = 0.0;
#ifdef ENABLE_APICAL_TUFTING_TRIGGER_DISTANCE
  double ATTD_TUFT_intDt = 0.0, ATTD_TUFT_intDt_plus_turns = 0.0, ATTD_TRUNC_intDt = 0.0, ATTD_TRUNC_intDt_plus_turns = 0.0, ATTD_TUFT_Bt = 0.0, ATTD_TUFT_Bt_plus_turns = 0.0, ATTD_TRUNC_Bt = 0.0, ATTD_TRUNC_Bt_plus_turns = 0.0;
#endif
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    //processing_neuron = e;
    request_elongation(e); // See TL#200603151008.1.
    // compute probability of terminal segment branching in this time step
    fibre_structure * apical_dendrite = NULL;
    if (e->TypeID()==PYRAMIDAL) apical_dendrite = e->InputStructure()->tail(); // APICAL DENDRITE to be treated differently
    if (apical_dendrite) {
      APICALintDt = expected_number_of_branches_APICAL(e,t);
      APICALintDt_plus_turns = APICALintDt + (_turnfactor*APICALintDt);
      if (competition_with_all_trees!=gcc_within_tree) {
	APICALBt = APICALintDt;
	APICALBt_plus_turns = APICALintDt_plus_turns;
	if (APICALintDt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per APICAL dendrite terminal segment (APICALBt_plus_turns=" << String(APICALintDt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
      }
      APICAL_E = apical_trunc_or_tuft_E(t);
#ifdef ENABLE_APICAL_TUFTING_TRIGGER_DISTANCE
      if (apical_tufting_trigger_distance>0.0) {
	ATTD_TUFT_intDt = expected_number_of_branches_APICAL_TUFT(e,t);
	ATTD_TUFT_intDt_plus_turns = ATTD_TUFT_intDt + (_turnfactor*ATTD_TUFT_intDt);
	ATTD_TRUNC_intDt = expected_number_of_branches_APICAL_TRUNC(e,t);
	ATTD_TRUNC_intDt_plus_turns = ATTD_TRUNC_intDt + (_turnfactor*ATTD_TRUNC_intDt);
	if (competition_with_all_trees!=gcc_within_tree) {
	  ATTD_TUFT_Bt = ATTD_TUFT_intDt;
	  ATTD_TUFT_Bt_plus_turns = ATTD_TUFT_intDt_plus_turns;
	  if (ATTD_TUFT_intDt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per APICAL dendrite terminal segment (ATTD_TUFT_intDt_plus_turns=" << String(ATTD_TUFT_intDt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
	  ATTD_TRUNC_Bt = ATTD_TRUNC_intDt;
	  ATTD_TRUNC_Bt_plus_turns = ATTD_TRUNC_intDt_plus_turns;
	  if (ATTD_TRUNC_intDt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per APICAL dendrite terminal segment (ATTD_TRUNC_intDt_plus_turns=" << String(ATTD_TRUNC_intDt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
	}
      }
#endif
    }
    BASALintDt = expected_number_of_branches(e,t);
    BASALintDt_plus_turns = BASALintDt + (_turnfactor*BASALintDt);
    if (competition_with_all_trees!=gcc_within_tree) {
      BASALBt = BASALintDt;
      BASALBt_plus_turns = BASALintDt_plus_turns;
      if (BASALintDt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per dendrite terminal segment (BASALBt_plus_turns=" << String(BASALintDt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
    }
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
      // [***INCOMPLETE] Take the desired statistic (i.e. the calculated
      // branching probability) into account for possible branches at
      // continuation nodes (and apply this method to other grow() functions
      // above as well.
      // A specific Delayed Branching Model can be chosen here.
      //Delayed_Branching_Model * dbm = new Delayed_Branching_Model(e,ps,Btperts);
      if (apical_dendrite==ps) {
	Bt = APICALBt; intDt = APICALintDt; Bt_plus_turns = APICALBt_plus_turns; intDt_plus_turns = APICALintDt_plus_turns;
      } else {
	Bt = BASALBt; intDt = BASALintDt; Bt_plus_turns = BASALBt_plus_turns; intDt_plus_turns = BASALintDt_plus_turns;
      }
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Bt = 1.0; // insure branching
      else {
	if (competition_with_all_trees==gcc_within_tree) {
	  if (apical_dendrite==ps) {
	    register double logpscount = log((double) ps->count_terminal_segments());
	    nt_E = exp(-APICAL_E*logpscount);
#ifdef ENABLE_APICAL_TUFTING_TRIGGER_DISTANCE
	    if (apical_tufting_trigger_distance>0.0) {
	      double attdtuftnt_E = exp(-_E_apicaltuft*logpscount);
	      double attdtruncnt_E = exp(-_E_apicaltrunc*logpscount);
	      ATTD_TUFT_Bt = ATTD_TUFT_intDt*attdtuftnt_E;
	      ATTD_TUFT_Bt_plus_turns = ATTD_TUFT_intDt_plus_turns*attdtuftnt_E;
	      ATTD_TRUNC_Bt = ATTD_TRUNC_intDt*attdtruncnt_E;
	      ATTD_TRUNC_Bt_plus_turns = ATTD_TRUNC_intDt_plus_turns*attdtruncnt_E;
	    }
#endif
	  }
	  else nt_E = exp(-_E*log((double) ps->count_terminal_segments()));
	  Bt = intDt*nt_E;
	  Bt_plus_turns = intDt_plus_turns*nt_E;
	  if (Bt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per dendrite terminal segment (Bt_plus_turns=" << String(Bt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
	}
      }
      Delayed_Branching_Model dbm(e,ps,Bt);
      if (branchesatturns) ps->branches_at_continuation_nodes(&dbm);
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
#ifdef ENABLE_APICAL_TUFTING_TRIGGER_DISTANCE
	if ((apical_dendrite==ps) && (apical_tufting_trigger_distance>0.0)) { // Use distance criterion to determine apical trunc or tuft
	  if (!ts->IsTuft()) ts->Check_Tuft_Distance_Threshold(apical_tufting_trigger_distance);
	  if (ts->IsTuft()) { // use tuft parameter values
	    Bt = ATTD_TUFT_Bt;
	    Bt_plus_turns = ATTD_TUFT_Bt_plus_turns;
	  } else { // use trunc parameter values
	    Bt = ATTD_TRUNC_Bt;
	    Bt_plus_turns = ATTD_TRUNC_Bt_plus_turns;
	  }
	}
#endif
	ts_next = ts->Next(); // preserved in a separate variable
	//if ((!strict_node_separation) || (ts->Length()>0.0)) {
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  double r = X_branch.get_rand_real1();
	  if (turnmodel==tm_branch_coupled) { // +++++
	    if (r<=Bt_plus_turns) { // a branch or a continuation point (turn)
	      // identify segment and finalize its geometric properties
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	      if (r<=Bt) SELECT_TYPED_BRANCH(dendrite_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	      if (r<=Bt) ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	      else { // continuation point (turn)
		if (ts->DirectionModel()) ts->turn(e); // The test is for legacy purposes only.
		else if (dendrite_Samsonovich_hypothesis) turn(ts,e);
		else turn_free(ts);
	      }
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    }
	  } else { // Apply an uncoupled turn model // +++++
	    if (r<=Bt) {
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	      SELECT_TYPED_BRANCH(dendrite_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	      ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    } else if ((turnmodel==tm_each_dt) || (r<=DirectionChange_Threshold(ts))) {
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
#endif // LEGACY_BRANCHING_PROTOCOL
	      if (ts->DirectionModel()) ts->turn(e); // The test is for legacy purposes only.
	      else if (dendrite_Samsonovich_hypothesis) turn(ts,e);
	      else turn_free(ts);
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    }
	  } // +++++
	}
      }
    }
  }
  firstiteration = false;
}

void Polynomial_O1_dendritic_growth_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // <A NAME="Polynomial_O1_dendritic_growth_model::grow"> </A>
  // Growing in fixed time steps.
  // This is based on a linear polynomial growth function for which default
  // parameters are based on polynomial curve fitting to data by Ramakers.
  // Turns at continuation points are included.
  growing_dendrites = true;
  if (!net) error("nibr Error: No network in Polynomial_O1_dendritic_growth_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in Polynomial_O1_dendritic_growth_model::grow()\n");
  double Bt, Btperts, Bt_plus_turns, nt_E, lognt;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) { //processing_neuron = e;
    fibre_structure * apical_dendrite = NULL;
    if (e->TypeID()==PYRAMIDAL) apical_dendrite = e->InputStructure()->tail(); // an apical dendrite to be treated differently
    request_elongation(e); // See TL#200603151008.1.
    // compute probability of terminal segment branching in this time step
    Bt = expected_number_of_branches(e,t); // without arbor specific division by n
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
      lognt = log((double) ps->count_terminal_segments());
      nt_E = exp(-_E*lognt);
      Btperts = Bt*nt_E; // probability per terminal segment of this arbor
      Bt_plus_turns = Btperts + (_turnfactor*Btperts);
      if (Bt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per dendrite terminal segment (Bt_plus_turns=" << String(Bt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Btperts = 1.0; // insure branching
      // [***INCOMPLETE] Take the desired statistic (i.e. the calculated
      // branching probability) into account for possible branches at
      // continuation nodes (and apply this method to other grow() functions
      // above as well.
      // A specific Delayed Branching Model can be chosen here.
      //Delayed_Branching_Model * dbm = new Delayed_Branching_Model(e,ps,Btperts);
      Delayed_Branching_Model dbm(e,ps,Btperts);
      if (branchesatturns) ps->branches_at_continuation_nodes(&dbm);
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  double r = X_branch.get_rand_real1();
	  if (r<=Bt_plus_turns) { // a branch or a continuation point (turn)
	    // identify segment and finalize its geometric properties
	    fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	    double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	    if (branchsomewhereinsegment) {
	      if (firstiteration) {
		if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
	      } else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	    }
#else
	    if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	    if (r<=Btperts) SELECT_TYPED_BRANCH(dendrite_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	    if (r<=Btperts) ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
            else { // continuation point (turn)
	      if (ts->DirectionModel()) ts->turn(e);
	      else if (dendrite_Samsonovich_hypothesis) turn(ts,e);
	      else turn_free(ts);
	    }
	    if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_dendritic_segment(fs)) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	    if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	  }
	}
      }
    }
  }
  firstiteration = false;
}

bool van_Pelt_dendritic_growth_model::IsComplete(double & t) {
  t += _dt;
  return (t>max_growth_time);
}

void van_Pelt_dendritic_growth_model::postop(network * net, Connection_Statistics_Root * cstats, double & t) {
  // Note that the fibre_segment coordinates are updated here for each
  // terminal segment.
  if (net->Seek_Candidate_Synapses()) {
    PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
      //processing_neuron = e;
      PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps) {
	PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
	  ts->update_segment_vector();
	  if (!net->spatial_segment_subset()->add_dendritic_segment(ts->TerminalSegment())) warning(String(t,"Warning: Unable to add segment to dendritic spatial subsets in van_Pelt_dendritic_growth_model::postop() (at t=%.1f)\n"));
	}
      }
      if (outattr_show_progress) { cout << '_'; cout.flush(); }
    }
  } else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
	   PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,ps)
	     PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) ts->update_segment_vector();
}

void van_Peltlike_axonal_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("Agrowth_E"))>=0) _E = atof(clp.ParValue(n));
#ifdef STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION
  if ((n=clp.Specifies_Parameter("Agrowth_tau"))>=0) { _tau = atof(clp.ParValue(n)); _c = _B_inf / _tau; }
#else
  if ((n=clp.Specifies_Parameter("Agrowth_tau"))>=0) _tau = atof(clp.ParValue(n));
#endif
  //if ((n=clp.Specifies_Parameter("Agrowth_n_inf"))>=0) _n_inf = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_B_inf"))>=0) { _B_inf = atof(clp.ParValue(n)); _c = _B_inf / _tau; }
  if ((n=clp.Specifies_Parameter("Agrowth_c"))>=0) { _c = atof(clp.ParValue(n)); _B_inf = _c*_tau; }
  if ((n=clp.Specifies_Parameter("Agrowth_L0_min"))>=0) _L0_min = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_L0_max"))>=0) _L0_max = atof(clp.ParValue(n));
#ifdef LEGACY_BRANCHING_PROTOCOL
  if ((n=clp.Specifies_Parameter("Amin_node_interval"))>=0) min_node_interval = atof(clp.ParValue(n));
#endif
  if ((n=clp.Specifies_Parameter("Acompetition_with_all_trees"))>=0) {
    if (downcase(clp.ParValue(n))==String("true")) competition_with_all_trees = gcc_within_all_AorD;
    else competition_with_all_trees = gcc_within_tree;
  }
  if ((n=clp.Specifies_Parameter("Acompetition_with_whole_neuron"))>=0) if (downcase(clp.ParValue(n))==String("true")) competition_with_all_trees = gcc_whole_neuron;
  set_fixed_time_step_size(_dt);
}

String van_Peltlike_axonal_growth_model::local_report_parameters() {
  return van_Pelt_dendritic_growth_model::local_report_parameters();
}

String van_Peltlike_axonal_growth_model::report_parameters() {
  return "van_Peltlike_axonal_growth_model: _E = " + local_report_parameters() + '\n';
}

van_Peltlike_axonal_growth_model::van_Peltlike_axonal_growth_model(Command_Line_Parameters & clp) {
  _E = 0.051; _tau = 432000.0; _B_inf = 5.31; _L0_min = 18.0; _L0_max = 22.0;
  _c = _B_inf / _tau;
  e_1 = exp(-1.0);
  _dt = van_Pelt_SUGGESTED_dt;
  parse_CLP(clp);
}

void van_Peltlike_axonal_growth_plus_turns_model::local_parse_CLP(Command_Line_Parameters & clp) {
  int n;
  // direction change model parameters
  if ((n=clp.Specifies_Parameter("Agrowth.turn_model"))>=0) {
    if (downcase(clp.ParValue(n))==String("branch_coupled")) turnmodel = tm_branch_coupled;
    else if (downcase(clp.ParValue(n))==String("linear_rate")) turnmodel = tm_linear_rate;
    else if (downcase(clp.ParValue(n))==String("each_dt")) turnmodel = tm_each_dt;
    else cout << "Warning: Unknown turn model '" << clp.ParValue(n) << "'\n";
  }
  if ((n=clp.Specifies_Parameter("Agrowth_tf"))>=0) _turnfactor = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Aturn_rate"))>=0) turn_rate = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Aturn_separation"))>=0) turn_rate = 1.0/atof(clp.ParValue(n)); // see explanation in DirectionChange_Threshold()
  if ((n=clp.Specifies_Parameter("Abranchesatturns"))>=0) branchesatturns = (downcase(clp.ParValue(n))==String("true"));
  double max_x = turnanglemax - turnanglemin;
  if (max_x<=0) error("Error: turnanglemax must be greater than turnanglemin\n");
  double b = 1.0;
  turn_a = -b/max_x;
  turn_max_y = 0.5*b*max_x;
  turn_SQb = b*b;
  turn_negb = -b;
  turn_twoa = 2.0*turn_a;
#ifdef VECTOR2D
  turn_twomax_y = 2.0*turn_max_y;
#endif
}

void van_Peltlike_axonal_growth_plus_turns_model::parse_CLP(Command_Line_Parameters & clp) {
  van_Peltlike_axonal_growth_model::parse_CLP(clp);
  local_parse_CLP(clp);
}

String van_Peltlike_axonal_growth_plus_turns_model::local_report_parameters() {
  String res(van_Peltlike_axonal_growth_model::local_report_parameters() + " _turnfactor = " + String(_turnfactor,"%f.3"));
  switch (turnmodel) {
  case tm_branch_coupled: res += " turnmodel = branch_coupled"; break;;
  case tm_linear_rate: res += " turnmodel = linear_rate"; break;;
  default: res += " turnmodel = each_dt";
  }
  res += String(turn_rate," turn_rate = %f.3");
  return res;
}

String van_Peltlike_axonal_growth_plus_turns_model::report_parameters() {

  return "van_Peltlike_axonal_growth_plus_turns_model: _E = " + local_report_parameters() + '\n';
}

van_Peltlike_axonal_growth_plus_turns_model::van_Peltlike_axonal_growth_plus_turns_model(Command_Line_Parameters & clp): van_Peltlike_axonal_growth_model(clp), turnmodel(tm_linear_rate), _turnfactor(9.25), turn_rate(0.1) {
  local_parse_CLP(clp);
}

void Polynomial_O1_axonal_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  van_Peltlike_axonal_growth_plus_turns_model::local_parse_CLP(clp);
  if ((n=clp.Specifies_Parameter("Agrowth_E"))>=0) _E = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Agrowth_n_inf"))>=0) _n_inf = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_c"))>=0) _c = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Agrowth_nu0"))>=0) _nu0 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Agrowth_F"))>=0) _F = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_L0_min"))>=0) _L0_min = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_L0_max"))>=0) _L0_max = atof(clp.ParValue(n));
  set_fixed_time_step_size(_dt);
}

String Polynomial_O1_axonal_growth_model::report_parameters() {
  String res("Polynomial_O1_axonal_growth_model: _E = ");
  res += String(_E,"%f.3 _c = ") + String(_c,"%f.3 _L0_min = ")  + String(_L0_min,"%f.3 _L0_max = ") + String(_L0_max,"%f.3 _turnfactor = ") + String(_turnfactor,"%f.3\n");
  //+ String(_nu0,"%f.3 _F = ") + String(_F,"%f.3 _L0_min = ")
  return res;
}

Polynomial_O1_axonal_growth_model::Polynomial_O1_axonal_growth_model(double E, double c, double tf): van_Peltlike_axonal_growth_plus_turns_model(E,1.0,1.0,tm_linear_rate,tf,0.1) {
  _c = c;
  _L0_min = 0.0;
  _L0_max = 0.0;
  set_fixed_time_step_size(van_Pelt_SUGGESTED_dt);
  e_1 = exp(-1.0);
}

Polynomial_O1_axonal_growth_model::Polynomial_O1_axonal_growth_model(Command_Line_Parameters & clp) {
  // default values are obtained by polynomial approximation of mirrored data
  _E = 1.0; _c = 5.018568/DAYSECONDS; _L0_min = 0.0; _L0_max = 0.0; _turnfactor = 9.25;
  _dt = van_Pelt_SUGGESTED_dt;
  e_1 = exp(-1.0);
  parse_CLP(clp);
}

void Polynomial_O3_axonal_growth_model::local_parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("Agrowth_c1"))>=0) _c1 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("Agrowth_c2"))>=0) _c2 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Agrowth_nu1"))>=0) _nu1 = atof(clp.ParValue(n));
  //if ((n=clp.Specifies_Parameter("Agrowth_nu2"))>=0) _nu2 = atof(clp.ParValue(n));
}

void Polynomial_O3_axonal_growth_model::parse_CLP(Command_Line_Parameters & clp) {
  Polynomial_O1_axonal_growth_model::parse_CLP(clp);
  local_parse_CLP(clp);
  set_fixed_time_step_size(_dt);
}

String Polynomial_O3_axonal_growth_model::report_parameters() {
  String res("Polynomial_O3_axonal_growth_model: _E = ");
  res += String(_E,"%f.3 _c = ") + String(_c,"%f.3 _c1 = ") + String(_c1,"%f.3 _c2 = ") + String(_c2,"%f.3 _L0_min = ") + String(_L0_min,"%f.3 _L0_max = ") + String(_L0_max,"%f.3 _turnfactor = ") + String(_turnfactor,"%f.3\n");
  // + String(_nu0,"%f.3 _nu1 = ") + String(_nu1,"%f.3 _nu2 = ") + String(_nu2,"%f.3 _F = ") + String(_F,"%f.3 _L0_min = ") 
  return res;
}

Polynomial_O3_axonal_growth_model::Polynomial_O3_axonal_growth_model(double E, double c, double c1, double c2, double tf): Polynomial_O1_axonal_growth_model(E,c,tf), _c1(c1), _c2(c2) {
  _L0_min = 0.0; // L(0) cubic approximation with mirrored data
  _L0_max = 0.0;
}

Polynomial_O3_axonal_growth_model::Polynomial_O3_axonal_growth_model(Command_Line_Parameters & clp): Polynomial_O1_axonal_growth_model(clp), _c1(0.0/DAYSECONDSSQR), _c2(0.001076/DAYSECONDSCUB) {
  _c = 4.724645/DAYSECONDS;
  local_parse_CLP(clp);
}

double van_Peltlike_axonal_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  double nt_E = exp(-_E*log((double) n->total_output_terminal_segments()));
  return _c*_tau*(exp(-t1/_tau)-exp(-t2/_tau))*nt_E;
}

double Polynomial_O1_axonal_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  return _c*(t2-t1); // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c =  4.844813 / (24.0*60.0*60.0))
}

double Polynomial_O3_axonal_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  double t1sq = t1*t1, t2sq = t2*t2;
  return (_c*(t2-t1)) + (_c1*(t2sq-t1sq)) + (_c2*((t2sq*t2)-(t1sq*t1)));
}

double van_Peltlike_axonal_growth_plus_turns_model::DirectionChange_Threshold(terminal_segment * ts) {
  // Returns the probability that a turn happened in the most recently elongated piece of
  // a neurite. This depends on the turning rate in expected number of turns per micron.
  // If the elongated length is not cached, then the function interprets the turning rate
  // as the expected number of turns per second.
  // the pointer ts must be non-NULL.
  // When we are able to calculate the probability of a turn as a function of elongated
  // length, i.e. when ENABLE_FIXED_STEP_SIMULATION is defined, then the following
  // equations hold: The probability of a turn is equal to the turn rate expressed as
  // the average number of turns per unit length multiplied by the elongated length.
  // Also, the probability of a turn is equal to the elongated length divided by the
  // average elongated length per turn (the average distance between turns).
  // The first way is specified by A/Dturn_rate, the second way by A/Dturn_separation.
  // If turn_rate*len = len / turn_separation, then turn_rate = 1 / turn_separation.
#ifdef ENABLE_FIXED_STEP_SIMULATION
  double P_turn = turn_rate*ts->elongated_length();
#else
  double P_turn = turn_rate*_dt;
#endif
  if (P_turn>=1.0) cout << "Warning: One or more turning points are expected (axon P_turn=" << String(P_turn,"%.2f") << "), perhaps development time steps should be smaller\n";
  return P_turn;
}

void Polynomial_O1_axonal_growth_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
  _c_dt = _c*dt;
  //_dt_nu0 = dt*_nu0;
}

double van_Peltlike_axonal_growth_model::expected_number_of_branches(neuron * n, double t) {
  // Returns the integral of p_s(t) over the time interval _dt using an
  // n(t)^(-E) that counts all terminal segments of all axonal arbors of
  // neuron n if competition_with_all_trees==true. The competition does not
  // take into account the value of n(t) after the update interval, i.e.
  // assumes that it is a constant.
  // Otherwise, returns the integral of the baseline component D(t) over the
  // time interval _dt, leaving it up to the calling function to multiply
  // with an appropriate n(t) dependence.
  double res = _B_inf*_c_dt*exp(-t/_tau);
  if (competition_with_all_trees==gcc_within_tree) return res;
  if (competition_with_all_trees==gcc_within_all_AorD) {
    double nt_E = exp(-_E*log((double) n->total_output_terminal_segments()));
    return res*nt_E;
  }
  double nt_E = exp(-_E*log((double) (n->total_input_terminal_segments()+n->total_output_terminal_segments()))); // gcc_whole_neuron
  return res*nt_E;
}

double Polynomial_O3_axonal_growth_model::expected_number_of_branches(neuron * n, double t1) {
  return expected_number_of_branches(n,t1,t1+_dt);
}

double Polynomial_O1_axonal_growth_model::expected_number_of_branches(neuron * n, double t1) {
  // [*** NOTE] The O3 polynomial fit to bifurcation in Ger Ramakers' data
  // does not use the t1 parameter.
  return _c_dt; // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c = 4.844813  / (24.0*60.0*60.0))
}

void van_Peltlike_axonal_growth_model::initialize(network * net, Connection_Statistics_Root * cstats, double & t) {
  if (!net) error("nibr Error: No network in van_Peltlike_axonal_growth_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Peltlike_axonal_growth_model::grow()\n");
  if (outattr_show_progress) { cout << "Initializing first axonal segments...\n"; cout.flush(); }
  // initialize the fibre structures
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->initialize_output_structure(_L0_min,_L0_max);
  // grow the fibre structures
  firstiteration = true;
}

void van_Peltlike_axonal_growth_model::request_elongation(neuron * e) {
  PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
#ifdef ELONGATION_EVENT_PREDICTION
      ps->devactivity->queue(_arbor_development_activity_elongation);
#else
      ps->batch_elongate();
#endif
  }
}

#ifdef TESTING_RANDOM_DISTRIBUTION
unsigned long int randomdist[11] = {0,0,0,0,0,0,0,0,0,0,0};
#endif
void van_Peltlike_axonal_growth_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // growing in fixed time steps
  growing_dendrites = false;
  double Bt = 0.0, intDt;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    //processing_neuron = e;
    request_elongation(e); // See TL#200603151008.1.
    // compute probability of terminal segment branching in this time step
    intDt = expected_number_of_branches(e,t);
    if (competition_with_all_trees!=gcc_within_tree) Bt = intDt;
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Bt = 1.0; // insure branching
      else if (competition_with_all_trees==gcc_within_tree) Bt = intDt*exp(-_E*log((double) ps->count_terminal_segments()));
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  if (X_branch.get_rand_real1()<=Bt) {
	    fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	    double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	    if (branchsomewhereinsegment) {
	      if (firstiteration) {
		if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
	      } else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	    }
#else
	    if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	    SELECT_TYPED_BRANCH(axon_two_branch_types);
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
            ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	    if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Peltlike_axonal_growth_model::grow() (at t=%.1f)\n"));
	    if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	  }
	}
      }
    }
    // elongate terminal segments
  }
  firstiteration = false;
}

void van_Peltlike_axonal_growth_plus_turns_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // <A NAME="van_Pelt_axonal_growth_plus_turns_model::grow"> </A>
  // growing in fixed time steps
  // this is based on the growth function of the van_Pelt_axonal_growth_model,
  // but includes turns at continuation points in axons.
  // Note that the reduce_rand_calls option is not yet included, although
  // one may be added if a clear benefit is apparent.
  growing_dendrites = false;
  if (!net) error("nibr Error: No network in van_Pelt_axonal_growth_plus_turns_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in van_Pelt_axonal_growth_plus_turns_model::grow()\n");
  double Bt = 0.0, intDt, Bt_plus_turns = 0.0, intDt_plus_turns, nt_E;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) { //processing_neuron = e;
    request_elongation(e); // See TL#200603151008.1.
    // compute probability of terminal segment branching in this time step
    intDt = expected_number_of_branches(e,t);
    intDt_plus_turns = intDt + (_turnfactor*intDt);
    if (competition_with_all_trees!=gcc_within_tree) {
      Bt = intDt;
      Bt_plus_turns = intDt_plus_turns;
      if (intDt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per axon terminal segment (Bt_plus_turns=" << String(intDt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
    }
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
      // [***INCOMPLETE] Take the desired statistic (i.e. the calculated
      // branching probability) into account for possible branches at
      // continuation nodes (and apply this method to other grow() functions
      // above as well.
      // A specific Delayed Branching Model can be chosen here.
      //Delayed_Branching_Model * dbm = new Delayed_Branching_Model(e,ps,Btperts);
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Bt = 1.0; // insure branching
      else {
	if (competition_with_all_trees==gcc_within_tree) {
	  nt_E = exp(-_E*log((double) ps->count_terminal_segments()));
	  Bt = intDt*nt_E;
	  Bt_plus_turns = intDt_plus_turns*nt_E;
	  if (Bt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per axonal terminal segment (Bt_plus_turns=" << String(Bt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
	}
      }
      Delayed_Branching_Model dbm(e,ps,Bt);
      if (branchesatturns) ps->branches_at_continuation_nodes(&dbm);
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  double r = X_branch.get_rand_real1();
	  if (turnmodel==tm_branch_coupled) {
	    if (r<=Bt_plus_turns) { // a branch or a continuation point (turn)
	      // identify segment and finalize its geometric properties
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	      if (r<=Bt) SELECT_TYPED_BRANCH(axon_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	      if (r<=Bt) ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	      else { // continuation point (turn)
		if (ts->DirectionModel()) ts->turn(e);
		else if (axon_Samsonovich_hypothesis) turn(ts,e);
		else turn_free(ts);
	      }
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Peltlike_axonal_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    }
	  } else { // Apply an uncoupled turn model
	    if (r<=Bt) {
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	      SELECT_TYPED_BRANCH(axon_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	      ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Peltlike_axonal_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    } else if ((turnmodel==tm_each_dt) || (r<=DirectionChange_Threshold(ts))) {
	      fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	      double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	      if (branchsomewhereinsegment) {
		if (firstiteration) {
		  if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
		} else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	      }
#else
	      if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
#endif // LEGACY_BRANCHING_PROTOCOL
	      if (ts->DirectionModel()) ts->turn(e);
	      else if (axon_Samsonovich_hypothesis) turn(ts,e);
	      else turn_free(ts);
	      if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Peltlike_axonal_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	      if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	    }
	  }
	}
      }
    }
  }
  firstiteration = false;
}

void Polynomial_O1_axonal_growth_model::grow(network * net, Connection_Statistics_Root * cstats, double & t) {
  // <A NAME="Polynomial_O1_axonal_growth_model::grow"> </A>
  // Growing in fixed time steps.
  // This is based on a linear polynomial growth function for which default
  // parameters are based on polynomial curve fitting to data by Ramakers.
  // Turns at continuation points are included.
  growing_dendrites = false;
  if (!net) error("nibr Error: No network in Polynomial_O1_axonal_growth_model::grow()\n");
  if (!cstats) error("nibr Error: No Connection_Statistics_Root in Polynomial_O1_axonal_growth_model::grow()\n");
  double Bt, Btperts, Bt_plus_turns, nt_E, lognt;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) { //processing_neuron = e;
    request_elongation(e); // See TL#200603151008.1.
    // compute probability of terminal segment branching in this time step
    Bt = expected_number_of_branches(e,t); // without arbor specific division by n
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps) {
      lognt = log((double) ps->count_terminal_segments());
      nt_E = exp(-_E*lognt);
      Btperts = Bt*nt_E; // probability per terminal segment of this arbor
      Bt_plus_turns = Btperts + (_turnfactor*Btperts);
      if (Bt_plus_turns>=1.0) cout << "Warning: One or more bifurcations or turning points are expected per axon terminal segment (Bt_plus_turns=" << String(Bt_plus_turns,"%.2f") << "), perhaps development time steps should be smaller\n";
      if ((firstiteration) && (initiallengthatactualfirstbranch)) Btperts = 1.0; // insure branching
      // [***INCOMPLETE] Take the desired statistic (i.e. the calculated
      // branching probability) into account for possible branches at
      // continuation nodes (and apply this method to other grow() functions
      // above as well.
      // A specific Delayed Branching Model can be chosen here.
      //Delayed_Branching_Model * dbm = new Delayed_Branching_Model(e,ps,Btperts);
      Delayed_Branching_Model dbm(e,ps,Btperts);
      if (branchesatturns) ps->branches_at_continuation_nodes(&dbm);
      terminal_segment * ts_next;
      for (terminal_segment * ts = ps->TerminalSegments()->head(); (ts); ts = ts_next) {
	ts_next = ts->Next(); // preserved in a separate variable
	if (ts->Length()>=ps->BranchingModel()->Min_Node_Interval()) { // [***NOTE] This is now "strict" [20070624]
	  double r = X_branch.get_rand_real1();
	  if (r<=Bt_plus_turns) { // a branch or a continuation point (turn)
	    // identify segment and finalize its geometric properties
	    fibre_segment * fs = ts->TerminalSegment();
#ifdef LEGACY_BRANCHING_PROTOCOL
#ifdef ENABLE_FIXED_STEP_SIMULATION
	    double remainder = 0.0;
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
	    if (branchsomewhereinsegment) {
	      if (firstiteration) {
		if (!initiallengthatactualfirstbranch) remainder = ts->randomize_last_segment_elongation(0.0); // can branch at L0 or beyond
	      } else remainder = ts->randomize_last_segment_elongation(min_node_interval);
	    }
#else
	    if ((branchsomewhereinsegment) && (!firstiteration)) remainder = ts->randomize_last_segment_elongation(min_node_interval); // doesn't reduce the size of the initial segment length
#endif
#endif
	    if (r<=Btperts) SELECT_TYPED_BRANCH(axon_two_branch_types) // bifurcation point (branch)
#else // NOT(LEGACY_BRANCHING_PROTOCOL)
	    if (r<=Btperts) ts->branch(NULL);
#endif // LEGACY_BRANCHING_PROTOCOL
	    else { // continuation point (turn)
	      if (ts->DirectionModel()) ts->turn(e);
	      else if (axon_Samsonovich_hypothesis) turn(ts,e);
	      else turn_free(ts);
	    }
	    if (net->Seek_Candidate_Synapses()) if (!net->spatial_segment_subset()->add_axonal_segment(fs)) warning(String(t,"Warning: Unable to add segment to axonal spatial subsets in van_Pelt_axonal_growth_plus_turns_model::grow() (at t=%.1f)\n"));
	    if (sampled_output->max_res()) sampled_output->Max_Resolution_Sample(t,fs);
	  }
	}
      }
    }
  }
  firstiteration = false;
}

//#define TEST_AXON_POSTOP
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
      //processing_neuron = e;
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
      if (outattr_show_progress) { cout << '_'; cout.flush(); }
    }
#ifdef TEST_SSS_DEPTH
    cout << "Final SSS#" << sss_number << " min. memory consumption per SSS: " << sizeof(Spatial_Segment_Subset)+(sizeof(fibre_segment_ptr)*spatialsegmentsubsetsizelimit) << " bytes\n"; cout.flush();
#endif
  } else PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
	   PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps)
	     PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) ts->update_segment_vector();
}
