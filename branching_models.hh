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

   These classes strive to adhere to the NETMORPH interface protocol for
   model and parameter selection according to the smallest matching set.

   (Involved DIL entries: 20070612153815.1.)
 */

#ifndef __BRANCHING_MODELS_HH
#define __BRANCHING_MODELS_HH

#include <include/templates.hh>
//#include "fibre_structure.hh"
#include "spatial.hh"
#include "Command_Line_Parameters.hh"
#include "global.hh"

#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
class terminal_segment;
class fibre_structure;
#endif
#ifndef __NEURON_HH
class neuron;
#endif

extern double Btestarray[];

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
  // Note that the "_arbor" parameter in the clone() function allows specification
  // of a different arbor, for instance when cloning from a schema to set up a
  // branching model for an arbor.
protected:
  branching_model_base * contributing;
  double contributingweight;
  struct branching_model_base_parameters {
    fibre_structure * arbor;
    double min_node_interval;
    branching_model_base_parameters(fibre_structure * _arbor, double _min_node_interval): arbor(_arbor), min_node_interval(_min_node_interval) {}
    } base_parameters;
public:
  virtual ~branching_model_base() { if (contributing) contributing->delete_shared(); }// NOTE: On 2024-03-19, moved this from protected to public so that fibre_structure.hh can destruct this.
  branching_model_base(branching_model_base * bmcontrib, double & bmweight, branching_model_base & schema);
  branching_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual branching_model_base * clone(fibre_structure * _arbor = NULL) = 0;
  virtual void delete_shared() { delete this; } // This is never actually shared, see the note in van_Pelt_branching_model.
  //  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2) = 0; // Not needed at Arbor level.
  //  virtual void handle_turning(terminal_segment & ts) = 0; // Not needed at Arbor level.
  virtual double predict_branching(double weight = 1.0) = 0;
  virtual double predict(double weight, double & probability);
  virtual double branching_probability();
  virtual String report_parameters_specific() = 0;
  String report_parameters();
  double Min_Node_Interval() { return base_parameters.min_node_interval; }
  fibre_structure * Arbor() { return base_parameters.arbor; }
  void set_arbor(fibre_structure * _arbor) { base_parameters.arbor = _arbor; }
public:
  // cache variables, mostly intended for use by TSBMs
  // Note that if several chained TSBMs make use of cache variables at the
  // arbor level, then special care needs to be taken to avoid interference.
  // The easiest solution with little overhead (overhead only at the arbor
  // level) seems to be to hard-code each type of TSBM that needs such
  // cache to access a different cache variable.
  double cache1; // established for the van_Pelt_TSBM
};

typedef branching_model_base * branching_model_base_ptr;

/*! \brief Arbor branching model that implements the arbor-wide constraint function of the phenomenological model published by van Pelt and Uylings (2003), <em>Brain and Mind</em>.

  This branching model applies an exponential branching function and
  competition between growth cones to compute the probability that
  a bifurcation event occurred in a preceding simulation time interval.
  This model is described in the NETMORPH paper by Koene et al. and
  originally in van Pelt and Uylings (2003), <em>Brain and Mind</em>.

  In the default case, the branching probability \f$p_{i,j}\f$ of terminal
  segment j in the i-th time interval is calculated as
  \f[
  p_{i,j} = {n_i}^{-E} B_{\infty} \exp{(-t_i/\tau)} (\exp{(\Delta t/\tau)}-1)2^{-S\gamma_j}C_{n_i},
  \f]
  where \f$n_i\f$ is the number of terminal segments in the arbor
  (dendrite or axon tree), \f$E\f$ is a competition parameter,
  \f$\gamma_j\f$ is the centrifugal order of the terminal segment, \f$S\f$
  is a parameter that modulates the order-dependency of the branching
  probability, \f$C_n\f$ is a normalization constant. Assuming an
  exponential function \f$D(t) = c e^{-t/\tau}\f$ then
  \f$B(t) = \int_0^t{D(t)dt}\f$.
  
 */
class van_Pelt_branching_model: public branching_model_base {
  // Note: Due to the caching of the calculated probability, an instance of
  // this model cannot be shared by multiple arbor objects. Making use of
  // this fact, I allow the branching model to remember the arbor that it
  // belongs to.
protected:
  // cache variables
  double creation_t; // The time of neurogenesis.
  double last_t; // The last t for which model calculations were done, initialized from eq->T().
  double probability; // cached probability value
protected:
  struct van_Pelt_branching_model_parameters {
    double B_inf;
    double tau;
    double E;
    growth_cone_competition_type competition_with_all_trees; // select scope of N(t)
    van_Pelt_branching_model_parameters(double _B_inf, double _tau, double _E, growth_cone_competition_type _competition_with_all_trees): B_inf(_B_inf), tau(_tau), E(_E), competition_with_all_trees(_competition_with_all_trees) {}
  } parameters;
public:
  van_Pelt_branching_model(branching_model_base * bmbcontrib, double & bmbweight, van_Pelt_branching_model & schema);
  van_Pelt_branching_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual branching_model_base * clone(fibre_structure * _arbor = NULL);
  //  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2); // [***NOTE] Could use this to calculate a new N-dependent temp.variable.
  //  virtual void handle_turning(terminal_segment & ts);
  virtual double predict_branching(double weight = 1.0);
  virtual String report_parameters_specific();
};

class TSBM_base {
  // Each terminal segment is designated a TSBM object, which
  // is responsible for determining the branching probability for
  // the specific terminal segment, which may involve local data,
  // such as cetrifugal order or environmental influences.
  // Selection of a TSBM is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all branching_model objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // branching_model objects.
protected:
  // cache variables
  terminal_segment * cache_ts;
protected:
  TSBM_base * contributing;
  double contributingweight;
  //struct TSBM_base_parameters {
  //  TSBM_base_parameters() {}
  //} base_parameters;
  virtual ~TSBM_base() { if (contributing) contributing->delete_shared(); }
public:
  TSBM_base(TSBM_base * tsbmcontrib, double & tsbmweight, TSBM_base & schema);
  TSBM_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSBM_base * clone(terminal_segment * _ts = NULL) = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual double predict_branching() = 0; //double weight = 1.0) = 0;
  virtual double predict(double weight, double & probability);
  virtual double branching_probability(terminal_segment & ts);
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef TSBM_base * TSBM_base_ptr;

class van_Pelt_TSBM: public TSBM_base {
  // This TSBM modifies the arbor branching model probability in
  // accordance with centrifugal order and the S parameter.
protected:
  // cache variables
  double binSgamma;
protected:
  struct van_Pelt_TSBM_parameters {
    double S;
    van_Pelt_TSBM_parameters(double _Sinit): S(_Sinit) {}
  } parameters;
  void init_cache(terminal_segment * _ts);
public:
  van_Pelt_TSBM(TSBM_base * tsbmbcontrib, double & tsbmbweight, van_Pelt_TSBM & schema);
  van_Pelt_TSBM(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSBM_base * clone(terminal_segment * _ts = NULL);
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  //virtual void handle_turning(terminal_segment & ts);
  virtual double predict_branching(); //double weight = 1.0);
  virtual String report_parameters_specific();
};

class van_Pelt_specBM_TSBM: public van_Pelt_TSBM {
  // This TSBM is designed for use in oblique and tuft branches of pyramidal apical dendrites.
  // Instead of referencing the arbor branching model for overall branching resources and
  // the consequent overall probability of branching, this TSBM references a separate
  // branching model that is set when a first oblique or tuft branch is created at
  // bifurcation from the apical trunk. Subsequent branches share the same
  // branching model. See TL#200807210612.5.
protected:
  struct van_Pelt_specBM_TSBM_parameters {
    branching_model_base * bm;
    van_Pelt_specBM_TSBM_parameters(branching_model_base * _bm): bm(_bm) {}
  } specparameters;
public:
  van_Pelt_specBM_TSBM(TSBM_base * tsbmbcontrib, double & tsbmbweight, van_Pelt_specBM_TSBM & schema);
  van_Pelt_specBM_TSBM(String & thislabel, String & label, Command_Line_Parameters & clp, terminal_segment * ts = NULL);
  virtual TSBM_base * clone(terminal_segment * _ts = NULL);
  //virtual double predict_branching(double weight = 1.0);
  virtual double branching_probability(terminal_segment & ts);
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
  probability_distribution_function * parallelogramanglepdf; // this also needs cloning and deletion
  void predict_branch_angles(double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2, double weight = 1.0);
  virtual ~balanced_forces_branch_angle_model() { if (parallelogramanglepdf) delete parallelogramanglepdf; }
public:
  balanced_forces_branch_angle_model(branch_angle_model_base * bambcontrib, double & bambweight, balanced_forces_branch_angle_model & schema);
  balanced_forces_branch_angle_model(String & thislabel, String & label, Command_Line_Parameters & clp);
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
extern TSBM default_TSBM; // identifier, changed when a different universal model is set
extern String universal_TSBM_root; // Used only to report if chaining at the universal set level.
extern branch_angle_model default_branch_angle_model; // identifier, changed when a different universal model is set
extern String universal_branch_angle_model_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
typedef branching_model_base_ptr * region_branching_model_base_ptr;
extern region_branching_model_base_ptr * branching_model_region_subset_schemas;
// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
typedef TSBM_base_ptr * region_TSBM_base_ptr;
extern region_TSBM_base_ptr * TSBM_region_subset_schemas;
// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
typedef branch_angle_model_base_ptr * region_branch_angle_model_base_ptr;
extern region_branch_angle_model_base_ptr * branch_angle_model_region_subset_schemas;

// functions
branching_model_base * branching_model_selection(String & label, Command_Line_Parameters & clp, branching_model_base * superior_set, branching_model_base_ptr & bmbptr);
TSBM_base * TSBM_selection(String & label, Command_Line_Parameters & clp, TSBM_base * superior_set, TSBM_base_ptr & tsbmbptr, terminal_segment * ts = NULL);
branch_angle_model_base * branch_angle_model_selection(String & label, Command_Line_Parameters & clp, branch_angle_model_base * superior_set, branch_angle_model_base_ptr & bambptr);

#endif
// __BRANCHING_MODELS_HH

