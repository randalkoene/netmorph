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
// fibre_elongation_model.hh
// Randal A. Koene, 20051221

/* Documentation:
   This implements classes for models of neural fiber elongation
   during network development. These models seek to preserve the
   desired elongation statistics for axon and dendrite arbor.

   As terminal segments elongate, they are the obvious place in
   which to link to a fiber elongation model.

   The models here strive to abide by the guidelines for the
   network generation framework. For example, as a consequence of
   those guidelines, the models offer two types of constructors:
   (1) A constructor that determines if a chained model should
   be sought and linked, and that searches for parameter
   specifications at a given resolution. (2) A constructor that
   copies parameters from a schema and clones chained model
   objects. The first constructor is used when schemas are
   created. The second constructor is used when model objects
   are cloned according to the most subset specific schema
   (or directly from other model objects at the same resolution,
   as in the case of model propagation during branching at
   terminal segments).

   Some more documenation of the parameter specification method
   is available in neuron.cc.

   When creating a new model BASE, keep in mind the following:
   A. Find specifications, as at the top of neuron.cc.
   B. Allocate to arbor and terminal segments in prepost_structure.cc.
   C. Have a model selection function.
   D. Call handle_branching from terminal_segment::branch().
   E. Call handle_turning from terminal_segment::turn() (actually in turn_without_update_of_segment_vector()).
 */

#ifndef __FIBRE_ELONGATION_MODEL_HH
#include "event.hh"
#include "Command_Line_Parameters.hh"
#include "global.hh"
class elongation_rate_initialization_model_base;
#define __FIBRE_ELONGATION_MODEL_HH

#ifndef __FIBRE_STRUCTURE_HH
class terminal_segment;
class fibre_structure;
#endif

#ifdef TESTING_ELONGATION_TOTAL
extern double usedarborelongationfraction;
#endif

  /* To do:
     - Separate the tasks of (a) obtaining the amount of arbor elongation
       and (b) obtaining the proportions allocated to specific terminal
       segments by creating actual separate model classes. Each task could
       in fact be carried out by a chain of models.
       See also notes starting with TL#200602171257.1. [PARTLY DONE]
     - Link the models of task (a) to the pre/postsynaptic structures,
       see where best to link (i) model selection and (ii) calls to model
       functions for the models of task (b).
     - I can allow calls to (a) to occur once per cycle for each arbor,
       while calls to (b) can appear as needed (see current growth
       functions). Note that the assignment of proportions may also occur
       just once if that simplifies the distribution of elongation, i.e.
       (a) can call a proportion determining function of (b) in all
       terminal segments.
     - Create at least one functional derived class for the statistical
       approach.
     - In that class, obtain the mean elongation of an arbor for a specific
       time step. Also obtain the standard deviation and use a good
       probability density function for the random selection of the actual
       arbor elongation.
     - Divvie up the amount of fiber available among terminal segments,
       using plausible models for speeding up and slowing down of elongation.
       This is done in terms of a proportion of the available fiber, and
       that proportion can be determined by a combination of the chained
       elongation models.
     - Check to see if the explit handling in separate classes of the
       separate model tasks has any performance disadvantages.
     - Check to see that all guidelines of the framework were met.
   */

class arbor_elongation_model_base {
  // Each pre/postsynaptic arbor is designated an elongation model, since
  // the resources provided for elongation are hypothesized to be limited
  // and distributed per arbor. Some derived classes of this class
  // may allow sharing of the same object for many
  // arbors, in which case specific information is obtained when elongation
  // functions are called. Otherwise, each assigned elongation model object
  // is unique. If an elongation model object is shared by multiple arbors,
  // then such a model must define a delete_shared() function that counts
  // how many arbors are still using the object, so that it is only deleted
  // when no longer in use.
  // Elongation models to use can be chosen at general and increasingly
  // specific levels: network, neuron, pre/postsynaptic structure.
  // This base object enables chaining of multiple elongation models.
protected:
  arbor_elongation_model_base * contributing;
  double contributingweight;
  double prev_t; // equivalent of model t_0 (see TL#200603140355.2)
  probability_distribution_function * pdf; // this also needs cloning and deletion
  struct arbor_elongation_model_base_parameters {
    double dL; // elongation resources for this update
    arbor_elongation_model_base_parameters(double _dL): dL(_dL) {}
  } base_parameters;
  virtual ~arbor_elongation_model_base() { if (contributing) contributing->delete_shared(); if (pdf) delete pdf; }
public:
  arbor_elongation_model_base(arbor_elongation_model_base * aemcontrib, double & aemweight, arbor_elongation_model_base & schema);
  arbor_elongation_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual void delete_shared() { delete this; }
  void reset(double t = 0.0);
  virtual arbor_elongation_model_base * clone() = 0;
  int total_terminal_segments(fibre_structure * arbor, growth_cone_competition_type scope);
  virtual double predict_elongation(fibre_structure * arbor) = 0;
  virtual double predict(double weight, fibre_structure * arbor);
  virtual double elongation(fibre_structure * arbor);
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef arbor_elongation_model_base * arbor_elongation_model_base_ptr;

class van_Pelt_arbor_elongation_model: public arbor_elongation_model_base {
  // This is an elongation model that deals specifically with establishing
  // the resources applied to elongation in an arbor, regardless of
  // the specific distribution to terminal segments. It uses the model
  // functions by van Pelt and Uylings to determine the elongation over a
  // time interval of development.
protected:
  struct van_Pelt_arbor_elongation_model_parameters {
    double _nu0; // initial mean elongation rate
    double _F; // parameter for the elongation rate dependency on the distribution of resources to multiple terminal segments (growth cones)
    growth_cone_competition_type _competition_with_all_trees; // select scope of N(t)
    //van_Pelt_arbor_elongation_model_parameters() {}
    van_Pelt_arbor_elongation_model_parameters(double nu0, double F, growth_cone_competition_type competition_with_all_trees): _nu0(nu0), _F(F), _competition_with_all_trees(competition_with_all_trees) {}
  } parameters;
public:
  van_Pelt_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, van_Pelt_arbor_elongation_model & schema);
  van_Pelt_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual arbor_elongation_model_base * clone();
  virtual double predict_elongation(fibre_structure * arbor);
  //void copy_parameters(van_Pelt_arbor_elongation_model_parameters & p) { p = parameters; }
  //double get_nu0() { return _nu0; }
  //double get_F() { return _F; }
  virtual String report_parameters_specific();
};

class Polynomial_O1_arbor_elongation_model: public arbor_elongation_model_base {
  // This arbor elongation model attempts to match elongation statistics as
  // in the van Pelt model, but is compatible with the polynomial approach to
  // branching. There is one significant difference in the elongation
  // calculation: The van Pelt model adjusts elongation according to the number
  // of terminal segments in all arbors of one kind (dendrite, axon) combined.
  // The Polynomial models adjust elongation according to the number of
  // terminal segments in the same arbor.
protected:
  // Note that I choose to allocate local parameters instead of inheriting
  // van_Pelt_arbor_elongation_model, so that I can specify different default
  // parameter values.
  struct Polynomial_O1_elongation_model_parameters {
    double _nu0; // initial mean elongation rate
    double _F; // parameter for the elongation rate dependency on the distribution of resources to multiple terminal segments (growth cones)
    Polynomial_O1_elongation_model_parameters(double nu0, double F): _nu0(nu0), _F(F) {}
  } parameters;
public:
  Polynomial_O1_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, Polynomial_O1_arbor_elongation_model & schema);
  Polynomial_O1_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual arbor_elongation_model_base * clone();
  virtual double predict_elongation(fibre_structure * arbor);
  virtual String report_parameters_specific();
};

class Polynomial_O3_arbor_elongation_model: public arbor_elongation_model_base {
  // Here the assumed polynomial is of order 3.
protected:
  // Note that I choose to allocate local parameters instead of inheriting
  // Polynomial_O1_arbor_elongation_model, so that I can specify different
  // default parameter values.
  struct Polynomial_O3_elongation_model_parameters {
    double _nu0; // initial mean elongation rate
    double _F; // parameter for the elongation rate dependency on the distribution of resources to multiple terminal segments (growth cones)
    double _nu1, _nu2; // square and cubic coefficients of elongation function
    Polynomial_O3_elongation_model_parameters(double nu0, double nu1, double nu2, double F): _nu0(nu0), _F(F), _nu1(nu1), _nu2(nu2) {}
  } parameters;
public:
  Polynomial_O3_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, Polynomial_O3_arbor_elongation_model & schema);
  Polynomial_O3_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual arbor_elongation_model_base * clone();
  virtual double predict_elongation(fibre_structure * arbor);
  virtual String report_parameters_specific();
};

class terminal_segment_elongation_model_base {
  // Models specific for terminal segments may be shared, in which case it
  // is sensible to designate them per pre/postsynaptic arbor.
  // Some derived classes of this calss
  // may allow sharing of the same object for many
  // arbors, in which case specific information is obtained when elongation
  // functions are called. Otherwise, each assigned elongation model object
  // is unique. If an elongation model object is shared by multiple arbors,
  // then such a model must define a delete_shared() function that counts
  // how many arbors are still using the object, so that it is only deleted
  // when no longer in use.
  // Elongation models to use can be chosen at general and increasingly
  // specific levels: network, neuron, pre/postsynaptic structure.
  // This base object enables chaining of multiple elongation models.
  friend class elongation_rate_initialization_model_base;
protected:
  terminal_segment_elongation_model_base * contributing;
  double contributingweight;
  double prev_t; // time of last elongation update
  double l_i_cache_t; // time of last perturbed expected elongation update
#ifndef REDUCED_MEMORY_PTSEM
  probability_distribution_function * pdf; // this also needs cloning and deletion
#endif
  struct terminal_segment_elongation_model_base_parameters {
    // The main purpose of l_i_cache: To cache the result of perturbed_expected_elongation(),
    //   so that multiple calculations per simulation time are avoided. This is done in
    //   the l_i_cache variable of the model at the HEAD of a chain of models.
    // Additional uses of l_i_cache: The cache value is READ in perturbed_expected_elongation()
    //   only if no function calls are needed, i.e. if the cached value is returned. Otherwise,
    //   the cached values is only WRITTEN to in perturbed_expected_elongation(). That means,
    //   it can be safely used as a local cache in functions called during calculations, such
    //   as predict_elongate(). Used as a local cache, l_i_cache will reflect the overall
    //   COMBINED result to functions in a model at the HEAD of a chain (since any local
    //   caching will be overwritten in perturbed_expected_elongation()), but will reflect a
    //   purely local value to functions in a model elsewhere in a chain. This difference can
    //   ALTER the implicit model response. E.g. the BESTL_TSEM will attempt to maintain the
    //   same combined elongation result if at the head of a chain, but will attempt to
    //   contribute the same expected elongation quota if elsewhere in a chain.
    //   [See TL#200709150740.1.]
    double l_i_cache; // perturbed expected elongation
    terminal_segment_elongation_model_base_parameters(double l_i): l_i_cache(l_i) {}
  } base_parameters;
  virtual ~terminal_segment_elongation_model_base() {
    if (contributing) contributing->delete_shared();
#ifndef REDUCED_MEMORY_PTSEM
    if (pdf) delete pdf;
#endif
  }
public:
  terminal_segment_elongation_model_base(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, terminal_segment_elongation_model_base & schema);
  terminal_segment_elongation_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual void delete_shared() { delete this; }
  void reset(double t = 0.0);
  virtual terminal_segment_elongation_model_base * clone() = 0;
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual double perturbed_expected_elongation();
  virtual double predict_elongate() = 0;
  virtual double predict(double weight);
  virtual void  elongate(terminal_segment * ts);
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef terminal_segment_elongation_model_base * terminal_segment_elongation_model_base_ptr;

class simple_terminal_segment_elongation_model: public terminal_segment_elongation_model_base {
  // This elongation model has no memory and evenly divides resources,
  // prior to random perturbation (see generation-framework paper).
  // This is the equivalent of the legacy elongation model.
public:
  simple_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, simple_terminal_segment_elongation_model & schema);
  simple_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

class inertia_terminal_segment_elongation_model: public terminal_segment_elongation_model_base {
  // This is an elongation model that deals specifically with the
  // allocation of a proportion of available elongation resources to a
  // particular terminal segment. This model implements a hypothesis of
  // physical inertia in the distribution of growth over terminal segments,
  // which may speed up or slow down for each terminal segment, but do so
  // in a reasonably gradual manner. This hypothesis is related to the
  // hypothesis at the base of the tension direction model.
  // This may be sensibly combined with a model that can react more
  // abruptly to features that are encountered in the environment.
public:
  inertia_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, inertia_terminal_segment_elongation_model & schema);
  inertia_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

class second_order_terminal_segment_elongation_model: public terminal_segment_elongation_model_base {
  // This second order model maintains a parameter of elongation rate (or 
  // speed) and a parameter of acceleration. Together, these define the
  // relative state of terminal segment elongation (relative, since the
  // resulting elongation is scaled in proportion to the desired elongation
  // of other terminal segments that share the same resources).
  // A coefficient of inertia defines how strongly perturbations of the
  // acceleration affect the system.
protected:
  struct second_order_terminal_segment_elongation_model_parameters {
    double l_i_acceleration;
    double inertia;
    second_order_terminal_segment_elongation_model_parameters(double liacc, double _inertia): l_i_acceleration(liacc), inertia(_inertia) {}
  } parameters;
public:
  second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, second_order_terminal_segment_elongation_model & schema);
  second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

class constrained_second_order_terminal_segment_elongation_model: public second_order_terminal_segment_elongation_model {
  // This model is identical to second_order_terminal_segment_elongation_model
  // in all respects, except that the relative elongation speeds are
  // constrained between the values 0.0 and 1.0, thereby moderating the
  // differences between elongation rates at competing terminal segments.
public:
  constrained_second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, constrained_second_order_terminal_segment_elongation_model & schema);
  constrained_second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
};

class initialized_CSO_terminal_segment_elongation_model: public constrained_second_order_terminal_segment_elongation_model {
  // This model is identical to the constrained_second_order_terminal_
  // segment_elongation_model in all respects, except that the initial
  // acceleration of the relative elongation speed of each new fiber can
  // be specified. A specifion of 0.0 can result in delayed elongation
  // (appearing as delayed branching), while negative specifications are
  // interpreted as "clone the parent or schema current acceleration".
  // If the initial acceleration parameter is not specified then the default
  // behavior is to clone the parent or schema current acceleration.
  // Note the special case when initial acceleration is specified as a
  // negative value and an object constructor receives no parent or schema
  // information (i.e. it is itself a schema), in which case the initial
  // value is determined by the base class constructor (see comments in the
  // fibre_elongation_model.cc source code for constructors of this class).
protected:
  struct initialized_CSO_terminal_segment_elongation_model_parameters {
    double initial_l_i_acceleration;
    initialized_CSO_terminal_segment_elongation_model_parameters(double initialliacceleration): initial_l_i_acceleration(initialliacceleration) {}
  } icsoparameters;
public:
  initialized_CSO_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, initialized_CSO_terminal_segment_elongation_model & schema);
  initialized_CSO_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual String report_parameters_specific();
};

class decaying_second_order_terminal_segment_elongation_model: public second_order_terminal_segment_elongation_model {
  // This model is identical to second_order_terminal_segment_elongation_model
  // in all respects, except that the acceleration of elongation speed decays
  // toward zero unless it is sufficiently perturbed from zero.
protected:
  struct decaying_second_order_terminal_segment_elongation_model_parameters {
    double t_decay;
    decaying_second_order_terminal_segment_elongation_model_parameters(double tdecay): t_decay(tdecay) {}
  } dsoparameters;
public:
  decaying_second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, decaying_second_order_terminal_segment_elongation_model & schema);
  decaying_second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

class delayed_branch_terminal_segment_elongation_model: public terminal_segment_elongation_model_base {
  // This terminal segment elongation model can temporarily control
  // elongation of a new terminal segment after a bifurcation. This enables
  // delayed branching, while simply generating all branches according to
  // global stochastic functions.
  // Note: When used, the delay is automatically applied to "daughter"
  // branches at bifurcations, since those receive a new cloned copy of
  // the parent models, so that the delay is initialized in the
  // constructor.
  // See DIL#20050927095951.1.
protected:
  struct delayed_branch_terminal_segment_elongation_model_parameters {
    // The delay_time_constant parameter specifies the time constant of an
    // exponential function that governs the probability of specific
    // delay duration in this non-mechanistic delayed branching model.
    // [***INCOMPLETE] Adjust this to the way it is described in the
    // paper!
    double delay_time_constant;
    delayed_branch_terminal_segment_elongation_model_parameters(double d_t_c): delay_time_constant(d_t_c) {}
  } parameters;
    // The elongation_commencement variable is the computed time after
    // which elongation will no longer be delayed.
  double elongation_commencement;
public:
  delayed_branch_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, delayed_branch_terminal_segment_elongation_model & schema);
  delayed_branch_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

#ifdef INCLUDE_BRANCHWISE_COMPETING_ICSO_TSEM
class branchwise_competing_ICSO_terminal_segment_elongation_model: public initialized_CSO_terminal_segment_elongation_model {
  // In this model, terminal segments compete for resources on a branch-by-branch bases. This model is
  // planned in DIL#20060323042401.1 and described in ~/src/nnmodels/nibr/elongation-competing-per-bifurcation.pdf.
};
#endif

class BESTL_terminal_segment_elongation_model: public terminal_segment_elongation_model_base {
  // This elongation model assumes that elongation rates between bifurcation points are
  // fixed at the value in the parameter cache. A new elongation rate is drawn only
  // after bifurcation. This class overloads the handle_branching() function.
  // In the simplest form, no perturbation is applied at all at updates. This is
  // very similar to the BESTL model used in [PELT:VARIABILITY].
public:
  BESTL_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, BESTL_terminal_segment_elongation_model & schema);
  BESTL_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual terminal_segment_elongation_model_base * clone();
  virtual double predict_elongate();
  virtual String report_parameters_specific();
};

class elongation_rate_initialization_model_base {
  // Specific models that determine the initial elongation rate of
  // terminal segments after a bifurcation. The cached perturbed elongation
  // quota terminal_segment::base_parameters.l_i_cache is set.
  // Note that this is not necessary during growth prior to bifurcation,
  // since the quota has no effect then. All elongation resources determined
  // at the arbor level are used for elongation at the single growth cone.
  //
  // Models specific for terminal segments may be shared, in which case it
  // is sensible to designate them per pre/postsynaptic arbor.
  // Some derived classes of this calss
  // may allow sharing of the same object for many
  // arbors, in which case specific information is obtained when elongation
  // functions are called. Otherwise, each assigned elongation model object
  // is unique. If an elongation model object is shared by multiple arbors,
  // then such a model must define a delete_shared() function that counts
  // how many arbors are still using the object, so that it is only deleted
  // when no longer in use.
  // Elongation models to use can be chosen at general and increasingly
  // specific levels: network, neuron, pre/postsynaptic structure.
  // This base object enables chaining of multiple elongation models.
protected:
  elongation_rate_initialization_model_base * contributing;
  double contributingweight;
  //probability_distribution_function * pdf; // this also needs cloning and deletion
  /* struct elongation_rate_initialization_model_base_parameters {
    sometype something; // perturbed expected elongation
    elongation_rate_initialization_model_base_parameters(sometype _something): something(_something) {}
    } base_parameters; */
  double Get_Mean_and_Max_Quotas(terminal_segment * t, double & maxquota);
  virtual ~elongation_rate_initialization_model_base() {
    if (contributing) contributing->delete_shared();
    //if (pdf) delete pdf;
  }
public:
  elongation_rate_initialization_model_base(elongation_rate_initialization_model_base * ericontrib, double & eriweight, elongation_rate_initialization_model_base & schema);
  elongation_rate_initialization_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual void delete_shared() { delete this; }
  virtual elongation_rate_initialization_model_base * clone() = 0;
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0) = 0;
  virtual double predict(double weight, double & predictedquota1, double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2);
  virtual void initialize(terminal_segment * ts1, terminal_segment * ts2);
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef elongation_rate_initialization_model_base * elongation_rate_initialization_model_base_ptr;

class length_distribution_eri_model: public elongation_rate_initialization_model_base {
  // This model uses the initial lengths (parts of remainders) of the two terminal segments,
  // when available, as a hint to their relative initial elongation quotas.
protected:
  probability_distribution_function * pdf; // this also needs cloning and deletion
  virtual ~length_distribution_eri_model() {
    if (pdf) delete pdf;
  }
public:
  length_distribution_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, length_distribution_eri_model & schema);
  length_distribution_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual elongation_rate_initialization_model_base * clone();
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0);
  virtual String report_parameters_specific();
};

class pure_stochastic_eri_model: public elongation_rate_initialization_model_base {
  // This model does not take into account initial lengths.
protected:
  probability_distribution_function * pdf; // this also needs cloning and deletion
  virtual ~pure_stochastic_eri_model() {
    if (pdf) delete pdf;
  }
public:
  pure_stochastic_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, pure_stochastic_eri_model & schema);
  pure_stochastic_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual elongation_rate_initialization_model_base * clone();
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0);
  virtual String report_parameters_specific();
};

class zero_eri_model: public elongation_rate_initialization_model_base {
  // This model assumes initialization with a zero quota. Elongation depends on
  // modification or perturbation of that initial elongation rate quota.
public:
  zero_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, zero_eri_model & schema);
  zero_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual elongation_rate_initialization_model_base * clone();
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0);
  virtual String report_parameters_specific();
};

class unitary_eri_model: public elongation_rate_initialization_model_base {
  // This model assumes initialization with a 1.0 quota.
  struct unitary_eri_model_parameters {
    double initialquota; // perturbed expected elongation
    unitary_eri_model_parameters(double _initialquota): initialquota(_initialquota) {}
    } unitary_parameters;
public:
  unitary_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, unitary_eri_model & schema);
  unitary_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual elongation_rate_initialization_model_base * clone();
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0);
  virtual String report_parameters_specific();
};

class continue_defaults_eri_model: public elongation_rate_initialization_model_base {
  // This model continues the branch that received the longest remainder of preceding
  // elongation with the quota of the parent fiber, while the other branch is initialized
  // according to default values (often zero, 1.0 in the case of BESTL TSEM).
public:
  continue_defaults_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, continue_defaults_eri_model & schema);
  continue_defaults_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual elongation_rate_initialization_model_base * clone();
  virtual void predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight = 1.0);
  virtual String report_parameters_specific();
};

class terminal_segment_elongation_event: public Event {
protected:
  terminal_segment * ts;
public:
  terminal_segment_elongation_event(terminal_segment * _ts): ts(_ts) {}
  terminal_segment_elongation_event(double _t, terminal_segment * _ts): Event(_t), ts(_ts) {}
  virtual void event();
};

// global variables
extern arbor_elongation_model general_arbor_elongation_model;
extern String general_arbor_elongation_model_root;
extern terminal_segment_elongation_model general_terminal_segment_elongation_model;
extern String general_terminal_segment_elongation_model_root;
// [Now using the regular protocol natural set arrays.] extern arbor_elongation_model_base * general_arbor_elongation_model_schema;
// [Now using the regular protocol natural set arrays.] extern terminal_segment_elongation_model_base * general_terminal_segment_elongation_model_schema;
extern arbor_elongation_model_base_ptr general_arbor_elongation_model_base_subset_schemas[];
extern terminal_segment_elongation_model_base_ptr general_terminal_segment_elongation_model_base_subset_schemas[];
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
extern elongation_rate_initialization_model default_elongation_rate_initialization_model; // identifier, changed when a different universal model is set
extern String universal_elongation_rate_initialization_model_root; // Used only to report if chaining at the universal set level.
extern elongation_rate_initialization_model_base_ptr elongation_rate_initialization_model_subset_schemas[];

// Indices for a few natural subsets, see TL#200603120618.2
/* Now using the regular protocol natural set arrays.*/
#define ALL_AXONS_AEM                      0
#define ALL_DENDRITES_AEM                  1
#define ALL_PYRAMIDAL_AXONS_AEM            2
#define ALL_PYRAMIDAL_DENDRITES_AEM        3
#define ALL_INTERNEURON_AXONS_AEM          4
#define ALL_INTERNEURON_DENDRITES_AEM      5
#define ALL_APICAL_PYRAMIDAL_DENDRITES_AEM 6
#define NUM_NATURAL_SUBSETS_AEM            7
#define ALL_AXONS_TSEM                      0
#define ALL_DENDRITES_TSEM                  1
#define ALL_PYRAMIDAL_AXONS_TSEM            2
#define ALL_PYRAMIDAL_DENDRITES_TSEM        3
#define ALL_INTERNEURON_AXONS_TSEM          4
#define ALL_INTERNEURON_DENDRITES_TSEM      5
#define ALL_APICAL_PYRAMIDAL_DENDRITES_TSEM 6
#define NUM_NATURAL_SUBSETS_TSEM            7

// functions

arbor_elongation_model_base * arbor_elongation_model_selection(String & label, Command_Line_Parameters & clp, arbor_elongation_model_base * superior_set, arbor_elongation_model_base_ptr & aembptr);
terminal_segment_elongation_model_base * terminal_segment_elongation_model_selection(String & label, Command_Line_Parameters & clp, terminal_segment_elongation_model_base * superior_set, terminal_segment_elongation_model_base_ptr & tsembptr);
elongation_rate_initialization_model_base * elongation_rate_initialization_model_selection(String & label, Command_Line_Parameters & clp, elongation_rate_initialization_model_base * superior_set, elongation_rate_initialization_model_base_ptr & eribptr);

#endif
