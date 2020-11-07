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
// branching_models.hh
// Randal A. Koene, 20070509

/* Documentation:
   This implements base and derived classes for models of bifurcation selection
   at the arbor level and for branch angle determination at the level of
   growth cones or terminal segments.

   These classes stive to adhere to the NETMORPH interface protocol for
   model and parameter selection according to the smallest matching set.

   (Involved DIL entries: 20070612153815.1.)
 */

#ifndef __BRANCHING_MODELS_HH
#define __BRANCHING_MODELS_HH

#include "templates.hh"
//#include "fibre_structure.hh"
#include "spatial.hh"
#include "Command_Line_Parameters.hh"
#include "global.hh"

#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
class terminal_segment;
#endif
#ifndef __NEURON_HH
class neuron;
#endif

class branching_model_base {
  // Each axon/dendrite arbor is designated a branching_model object, which
  // is responsible for determining the momentary branching rate, and
  // therefore the probability that a bifurcation occurs at a growth cone.
  // In an event-driven implementation of the simulation, the branching_model
  // objects can also select growth cones at which branching occurs and
  // send those events to the queue.
  // Selection of a branching_model is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all branching_model objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // branching_model objects.
protected:
  branching_model_base * contributing;
  double contributingweight;
  struct branching_model_base_parameters {
    double min_node_interval;
    branching_model_base_parameters(double _min_node_interval): min_node_interval(_min_node_interval) {}
    } base_parameters;
  virtual ~branching_model_base() { if (contributing) contributing->delete_shared(); }
public:
  branching_model_base(branching_model_base * bmcontrib, double & bmweight, branching_model_base & schema);
  branching_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual branching_model_base * clone() = 0;
  virtual void delete_shared() { delete this; }
  //  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2) = 0; // Not needed at Arbor level.
  //  virtual void handle_turning(terminal_segment & ts) = 0; // Not needed at Arbor level.
  virtual double predict_branching(double weight = 1.0) = 0;
  virtual double predict(double weight, double & probability);
  virtual double branching_probability();
  virtual String report_parameters_specific() = 0;
  String report_parameters();
  double Min_Node_Interval() { return base_parameters.min_node_interval; }
};

typedef branching_model_base * branching_model_base_ptr;

class van_Pelt_branching_model: public branching_model_base {
  // This branching model applies an exponential branching function and
  // competition between growth cones to compute the probability that
  // a bifurcation event occurred in a preceding simulation time interval.
  // This model is described in the gfneuroinformatics paper, and originally
  // in [PELT:GROWTHFUNCTIONS].
protected:
  /* 
  [***INCOMPLETE] This class is currently used only to convey local information
                  about minimum node intervals.
  */
public:
  van_Pelt_branching_model(branching_model_base * bmbcontrib, double & bmbweight, van_Pelt_branching_model & schema);
  van_Pelt_branching_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual branching_model_base * clone();
  //  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2); // [***NOTE] Could use this to calculate a new N-dependent temp.variable.
  //  virtual void handle_turning(terminal_segment & ts);
  virtual double predict_branching(double weight = 1.0);
  virtual String report_parameters_specific();
};

class branch_angle_model_base {
  // Each terminal segment is designated a branch_angle_model object, which
  // is responsible for determining the branch angles when a bifurcation occurs.
  // Selection of a branch_angle_model is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all branch_angle_model objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // branch_angle_model objects.
protected:
  branch_angle_model_base * contributing;
  double contributingweight;
  // The PDF is used for perturbation and can be specified per arbor
  // [***NOTE] Efficiency may be improved if different terminal segments can have pointers to
  //           shared PDF objects. Use profiling to see if this is worth implementing.
  //           Or, see the REDUCED_MEMORY_PTSEM option in terminal_segment_elongation_model_base.
  probability_distribution_function * pdf; // this also needs cloning and deletion
  void bifurcate(terminal_segment * ts1, double branchangle1, terminal_segment * ts2, double branchangle2);
  /*  struct branch_angle_model_base_parameters {
    sometype something;
    branch_angle_model_base_parameters(sometime _something): something(_something) {}
    } base_parameters; */
  virtual ~branch_angle_model_base() { if (contributing) contributing->delete_shared(); if (pdf) delete pdf; }
public:
  branch_angle_model_base(branch_angle_model_base * bamcontrib, double & bamweight, branch_angle_model_base & schema);
  branch_angle_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual branch_angle_model_base * clone() = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2) = 0;
  virtual void handle_turning(terminal_segment & ts) = 0;
  virtual double predict(double weight, double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2) = 0;
  virtual void branch_angle(terminal_segment * ts1, terminal_segment * ts2) = 0; // fibre_segment * parent, 
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef branch_angle_model_base * branch_angle_model_base_ptr;

class balanced_forces_branch_angle_model: public branch_angle_model_base {
  // This branch angle model uses the predictions of elongation models to
  // obtain initial elongation rates of daughter branches. Those elongation
  // rates, aside from being used to determine the ratios of remaining
  // elongated neurite to receive as initial lengths, are also used to
  // represent relative forces along the two branches. Those forces are
  // used to determine the angles that each of the two branches makes
  // with the projection of the parent direction.
protected:
  void predict_branch_angles(double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2, double weight = 1.0);
public:
  balanced_forces_branch_angle_model(branch_angle_model_base * bambcontrib, double & bambweight, balanced_forces_branch_angle_model & schema);
  balanced_forces_branch_angle_model(String & thislabel, String & label, Command_Line_Parameters & clp): branch_angle_model_base(thislabel,label,clp) {}
  virtual branch_angle_model_base * clone();
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual double predict(double weight, double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2);
  virtual void branch_angle(terminal_segment * ts1, terminal_segment * ts2); // fibre_segment * parent, 
  virtual String report_parameters_specific();
};

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
extern branching_model default_branching_model; // identifier, changed when a different universal model is set
extern String universal_branching_model_root; // Used only to report if chaining at the universal set level.
extern branching_model_base_ptr branching_model_subset_schemas[];
extern branch_angle_model default_branch_angle_model; // identifier, changed when a different universal model is set
extern String universal_branch_angle_model_root; // Used only to report if chaining at the universal set level.
extern branch_angle_model_base_ptr branch_angle_model_subset_schemas[];

// functions
branching_model_base * branching_model_selection(String & label, Command_Line_Parameters & clp, branching_model_base * superior_set, branching_model_base_ptr & bmbptr);
branch_angle_model_base * branch_angle_model_selection(String & label, Command_Line_Parameters & clp, branch_angle_model_base * superior_set, branch_angle_model_base_ptr & bambptr);

#endif
// __BRANCHING_MODELS_HH

