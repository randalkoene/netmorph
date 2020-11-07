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
// dendritic_growth_model.hh
// Randal A. Koene, 20041118, 20060311
//
// Classes that contain fibre growth models and functions that apply
// those models.
// With the specification of separate fibre_elongation_model.hh classes,
// this set of classes now deals with the specification of branching and
// turning models. The grow() functions do include the calls to specific
// elongation models.

#ifndef __DENDRITIC_GROWTH_MODEL_HH
#define __DENDRITIC_GROWTH_MODEL_HH

#include "neuron.hh"
//#include "network.hh" // include in .cc if needed to avoid problems
#include "Connection_Statistics.hh"
#include "Command_Line_Parameters.hh"

#ifdef TEST_FOR_NAN
extern unsigned long int DBcounter;
#endif
//#define TESTING_RANDOM_DISTRIBUTION
#ifdef TESTING_RANDOM_DISTRIBUTION
extern unsigned long int randomdist[];
#endif

class network; // forward declaration resolved in network.hh

//class Connection_Statistics_Root; // forward declaration resolved in a nibr.hh-like header
//class neuron; // forward declaration resolved in a neuron.hh header
//class network; // forward declaration resolved in a nibr.hh-like header

class Delayed_Branching_Model {
// [***INCOMPLETE] The following factor should be provided elsewhere!
#define at_continuation_node_factor 0.1
protected:
  neuron * n;
  fibre_structure * arbor;
  double Btpercn; // branching probability per continuation node
public:
  Delayed_Branching_Model(neuron * N, fibre_structure * A, double Bprob): n(N), arbor(A), Btpercn(Bprob*at_continuation_node_factor) {}
  virtual ~Delayed_Branching_Model() {}
  virtual void branch_selection(fibre_segment * fs);
  virtual void continuation_node_to_branch(fibre_segment * fs, double branchanglemin, double branchanglemax);
};

class dendritic_growth_model: public CLP_Modifiable {
// this is intended as a template
public:
  dendritic_growth_model() {}
  virtual ~dendritic_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2) = 0;
  virtual void set_fixed_time_step_size(double dt) = 0;
  virtual double get_fixed_time_step_size() = 0;
  virtual double suggested_fixed_time_step_size() = 0;
  virtual double expected_number_of_branches(neuron * n, double t1) = 0;
  virtual void initialize(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  virtual bool IsComplete(double & t) = 0;
  virtual void postop(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  // [***INCOMPLETE] a synapse development function should select candidate
  // synapses that become actual synapses
  virtual void parse_CLP(Command_Line_Parameters & clp) = 0;
  virtual String report_parameters() = 0;
};

class van_Pelt_dendritic_growth_model: public dendritic_growth_model {
// see [<A HREF="../../../doc/html/planning/bibliography.html#PELT:GROWTHFUNCTIONS">PELT:GROWTHFUNCTIONS</A>]
#define van_Pelt_SUGGESTED_dt 8640.0
protected:
  // parameters without default values
  // [***IMPORTANT] When turning this into properly modular models, please keep in mind that the
  // apical dendrite trunc and tuft specific data and functions implemented for the "plus_turns" model
  // below should also be available without turns!
  double _E; // factor for resource competition between terminating segments
  double _tau; // branching rate time constant
  //double _n_inf; // asymptotic value of the number of terminating segments
  double _B_inf; // asymptotic expected number of branching events at a single terminal segment
  double _c; // scaling factor of the branching rate function
  // parameters with default values
  double _L0_min; // minimum sum of the lengths of initial segments
  double _L0_max; // maximum sum of the lengths of initial segments
#ifdef LEGACY_BRANCHING_PROTOCOL
  double min_node_interval; // desired minimum fiber length between nodes
#endif
  growth_cone_competition_type competition_with_all_trees;
  // cached values
  double _dt, _c_dt;
  double e_1; // exp(-1.0)
  bool firstiteration;
#ifdef LEGACY_BRANCHING_PROTOCOL
  void branch(terminal_segment * ts, double anglemin, double anglemax, double othermax);
#endif
  String local_report_parameters();
  van_Pelt_dendritic_growth_model() {} // not intended for direct declaration
  virtual void request_elongation(neuron * e);
public:
  van_Pelt_dendritic_growth_model(double E, double tau, double B_inf);
  van_Pelt_dendritic_growth_model(Command_Line_Parameters & clp);
  virtual ~van_Pelt_dendritic_growth_model() {}
  void set_initial_length_spread(double lmin, double lmax) { _L0_min=lmin; _L0_max=lmax; }
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual void set_fixed_time_step_size(double dt);
  virtual double get_fixed_time_step_size() { return _dt; };
  virtual double suggested_fixed_time_step_size() { return van_Pelt_SUGGESTED_dt; }
  virtual double expected_number_of_branches(neuron * n, double t);
  virtual void initialize(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual bool IsComplete(double & t);
  virtual void postop(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

enum turn_model_type { tm_each_dt, tm_linear_rate, tm_branch_coupled, tm_NUM }; // +++++

class van_Pelt_dendritic_growth_plus_turns_model: public van_Pelt_dendritic_growth_model {
// this includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>]
// some typical _turnfactors:
//   CA3 cell basal dendrites 9.54
//   CA3 cell apical dendrites 8.96
//   CA1 cell basal dendrites 8.88
//   CA1 cell apical dendrites 7.17
protected:
  // parameters without default values
  double _E_apicaltrunc; // factor for resource competition between terminating segments
  double _tau_apicaltrunc; // branching rate time constant
  double _B_inf_apicaltrunc; // asymptotic expected number of branching events at a single terminal segment
  double _c_apicaltrunc; // scaling factor of the branching rate function
  double _c_dt_apicaltrunc;
  double _E_apicaltuft; // factor for resource competition between terminating segments
  double _tau_apicaltuft; // branching rate time constant
  double _B_inf_apicaltuft; // asymptotic expected number of branching events at a single terminal segment
  double _c_apicaltuft; // scaling factor of the branching rate function
  double _c_dt_apicaltuft;
  double apical_tufting_trigger_time; // optional threshold for time-triggered tufting of apical dendrites
  double apical_tufting_trigger_distance; // optional threshold for distance-triggered tufting of apical dendrite terminal segments (if > 0.0)
  turn_model_type turnmodel; // selected turn model to apply // +++++
  double _turnfactor; // how many times more likely does a turn occur than a branch
  double turn_rate; // turn rate expressed in expected number of turns per micron or per second // +++++
  bool branchesatturns;
  //  void turn(terminal_segment * ts, neuron * e);
  void local_parse_CLP(Command_Line_Parameters & clp);
  String local_report_parameters();
  van_Pelt_dendritic_growth_plus_turns_model() {} // not intended for direct declaration
public:
  van_Pelt_dendritic_growth_plus_turns_model(double E, double tau, double B_inf, double E_apicaltrunc, double tau_apicaltrunc, double B_inf_apicaltrunc, double E_apicaltuft, double tau_apicaltuft, double B_inf_apicaltuft, double tufttriggertime, double tufttriggerdistance, turn_model_type tm, double tf, double tr, bool bt = false);
  van_Pelt_dendritic_growth_plus_turns_model(Command_Line_Parameters & clp);
  virtual ~van_Pelt_dendritic_growth_plus_turns_model() {}
  double apical_trunc_or_tuft_E(double t) {
    // Chooses a parameter for apical dendrites dependent on a development time trigger.
    if (t>=apical_tufting_trigger_time) return _E_apicaltuft;
    return _E_apicaltrunc;
  }
  virtual void set_fixed_time_step_size(double dt);
  virtual double expected_number_of_branches_APICAL(neuron * n, double t);
  virtual double expected_number_of_branches_APICAL_TRUNC(neuron * n, double t);
  virtual double expected_number_of_branches_APICAL_TUFT(neuron * n, double t);
  virtual double DirectionChange_Threshold(terminal_segment * ts); // +++++
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class Polynomial_O1_dendritic_growth_model: public van_Pelt_dendritic_growth_plus_turns_model {
// This growth function is derived from a first order polynomial
// approximation to (smoothing of) dendritic growth measurements.
// This includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>].
public:
  Polynomial_O1_dendritic_growth_model(double E, double c, double tf);
  Polynomial_O1_dendritic_growth_model(Command_Line_Parameters & clp);
  virtual ~Polynomial_O1_dendritic_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual void set_fixed_time_step_size(double dt);
  virtual double expected_number_of_branches(neuron * n, double t1);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class Polynomial_O3_dendritic_growth_model: public Polynomial_O1_dendritic_growth_model {
// This growth function is derived from a first order polynomial
// approximation to (smoothing of) dendritic growth measurements.
// This includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>].
protected:
  double _c1, _c2; // square and cubic coefficients of branching function
  void local_parse_CLP(Command_Line_Parameters & clp);
public:
  Polynomial_O3_dendritic_growth_model(double E, double c, double c1, double c2, double tf);
  Polynomial_O3_dendritic_growth_model(Command_Line_Parameters & clp);
  virtual ~Polynomial_O3_dendritic_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual double expected_number_of_branches(neuron * n, double t1);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class van_Peltlike_axonal_growth_model: public van_Pelt_dendritic_growth_model {
protected:
  String local_report_parameters();
  van_Peltlike_axonal_growth_model() {} // not intended for direct declaration
  virtual void request_elongation(neuron * e);
public:
  van_Peltlike_axonal_growth_model(double E, double tau, double B_inf): van_Pelt_dendritic_growth_model(E,tau,B_inf) { _L0_min = 18.0; _L0_max = 22.0; }
  van_Peltlike_axonal_growth_model(Command_Line_Parameters & clp);
  virtual ~van_Peltlike_axonal_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual double expected_number_of_branches(neuron * n, double t);
  virtual void initialize(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void postop(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class van_Peltlike_axonal_growth_plus_turns_model: public van_Peltlike_axonal_growth_model {
// Here the preferred direction rectifying effect of continuation points
// (turns) is a presupposition for axon segments.
// (Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>].)
protected:
  // parameters without default values
  turn_model_type turnmodel; // selected turn model to apply
  double _turnfactor; // how many times more likely does a turn occur than a branch
  double turn_rate; // turn rate expressed in expected number of turns per micron or per second // +++++
  bool branchesatturns;
  //  void turn(terminal_segment * ts, neuron * e);
  void local_parse_CLP(Command_Line_Parameters & clp);
  String local_report_parameters();
  van_Peltlike_axonal_growth_plus_turns_model() {} // not intended for direct declaration
public:
  van_Peltlike_axonal_growth_plus_turns_model(double E, double tau, double B_inf, turn_model_type tm, double tf, double tr, bool bt = false): van_Peltlike_axonal_growth_model(E,tau,B_inf), turnmodel(tm), _turnfactor(tf), turn_rate(tr), branchesatturns(bt) {}
  van_Peltlike_axonal_growth_plus_turns_model(Command_Line_Parameters & clp);
  virtual ~van_Peltlike_axonal_growth_plus_turns_model() {}
  double DirectionChange_Threshold(terminal_segment * ts);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class Polynomial_O1_axonal_growth_model: public van_Peltlike_axonal_growth_plus_turns_model {
// This growth function is derived from a first order polynomial
// approximation to (smoothing of) axonal growth measurements.
// This includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>].
//protected:
//  virtual void request_elongation(neuron * e, double t);
public:
  Polynomial_O1_axonal_growth_model(double E, double c, double tf);
  Polynomial_O1_axonal_growth_model(Command_Line_Parameters & clp);
  virtual ~Polynomial_O1_axonal_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual void set_fixed_time_step_size(double dt);
  virtual double expected_number_of_branches(neuron * n, double t1);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class Polynomial_O3_axonal_growth_model: public Polynomial_O1_axonal_growth_model {
// This growth function is derived from a first order polynomial
// approximation to (smoothing of) axonal growth measurements.
// This includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>].
protected:
  double _c1, _c2; // square and cubic coefficients of branching function
  //double _nu1, _nu2; // square and cubic coefficients of elongation function
  void local_parse_CLP(Command_Line_Parameters & clp);
public:
  Polynomial_O3_axonal_growth_model(double E, double c, double c1, double c2, double tf);
  Polynomial_O3_axonal_growth_model(Command_Line_Parameters & clp);
  virtual ~Polynomial_O3_axonal_growth_model() {}
  virtual double expected_number_of_branches(neuron * n, double t1, double t2);
  virtual double expected_number_of_branches(neuron * n, double t1);
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

// functions

#ifdef VECTOR3D
//extern double Stheta; // export this, just as veer() is exported

void prepare_parent_coordinate_system(spatial & parentacoords);
void daughter_to_network_coordinate_system(spatial & daughter, spatial & daughterxyz);
double veer(terminal_segment * ts, double veerangle, spatial & S1xyz, double rollangle = -1.0);
#endif

#endif
