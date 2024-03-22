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
#include "turning_models.hh"

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
  //virtual double expected_number_of_branches(neuron * n, double t1, double t2) = 0;
  virtual void set_fixed_time_step_size(double dt) = 0;
  virtual double get_fixed_time_step_size() = 0;
  virtual double suggested_fixed_time_step_size() = 0;
  //virtual double expected_number_of_branches(neuron * n, double t1) = 0;
  virtual void initialize(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  virtual bool IsComplete(double & t) = 0;
  virtual void postop(network * net, Connection_Statistics_Root * cstats, double & t) = 0;
  // [***INCOMPLETE] a synapse development function should select candidate
  // synapses that become actual synapses
  virtual void parse_CLP(Command_Line_Parameters & clp) = 0;
  virtual String report_parameters() = 0;
};

// I am joining together the with and without turns classes and using a global
// turns flag to conditionally execute turns related code. There are more
// elegant, protocol compliant ways to do this, for example by simply creating
// a turning model that always returns a turn probabilty of 0 to represent the
// no turns case. I am not doing that right now, because I am using the no-turns
// case as a way to conserve computational resources! (See TL#200808061810.1.)

class van_Pelt_dendritic_growth_model: public dendritic_growth_model {
// see [<A HREF="../../../doc/html/planning/bibliography.html#PELT:GROWTHFUNCTIONS">PELT:GROWTHFUNCTIONS</A>]
// this includes the preferred direction rectifying effect of continuation
// points (turns) in dendrite segments, as described in the paper by
// Samsonovich and Ascoli, see [<A HREF="../../../doc/html/planning/bibliography.html#SAMS:DENDRDIR">SAMS:DENDRDIR</A>]
// some typical _turnfactors:
//   CA3 cell basal dendrites 9.54
//   CA3 cell apical dendrites 8.96
//   CA1 cell basal dendrites 8.88
//   CA1 cell apical dendrites 7.17
#define van_Pelt_SUGGESTED_dt 8640.0
protected:
  bool dendrites; // dendrites or axons flag
  //double _L0_min; // minimum sum of the lengths of initial segments
  //double _L0_max; // maximum sum of the lengths of initial segments
  growth_cone_competition_type competition_with_all_trees;
  TSTM turnmodel; // selected turn model to apply
  double _turnfactor; // how many times more likely does a turn occur than a branch
  double turn_rate; // turn rate expressed in expected number of turns per micron or per second
  bool branchesatturns;
  // cached values
  double _dt;
  bool firstiteration;
  virtual void request_elongation(fibre_structure * first_fs);
  van_Pelt_dendritic_growth_model() {} // not intended for direct declaration
public:
  van_Pelt_dendritic_growth_model(bool _dendrites, TSTM tm, double tf, double tr, bool bt = false);
  van_Pelt_dendritic_growth_model(bool _dendrites, Command_Line_Parameters & clp);
  virtual ~van_Pelt_dendritic_growth_model() {}
  virtual void set_fixed_time_step_size(double dt);
  virtual double get_fixed_time_step_size() { return _dt; };
  virtual double suggested_fixed_time_step_size() { return van_Pelt_SUGGESTED_dt; }
  virtual void initialize(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual void grow(network * net, Connection_Statistics_Root * cstats, double & t);
  virtual bool IsComplete(double & t);
  virtual void postop(network * net, Connection_Statistics_Root * cstats, double & t);
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
