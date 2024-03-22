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
// turning_models.hh
// Randal A. Koene, 20070509

/* Documentation:
   This implements base and derived classes for models of bifurcation selection
   at the arbor level and for branch angle determination at the level of
   growth cones or terminal segments.

   These classes stive to adhere to the NETMORPH interface protocol for
   model and parameter selection according to the smallest matching set.
 */

#ifndef __TURNING_MODELS_HH
#define __TURNING_MODELS_HH

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

class turning_model_base {
  // [***INCOMPLETE] Currently, NETMORPH uses only the terminal segment
  // turning model.
  // Each axon/dendrite arbor is designated a turning_model object, which
  // is responsible for determining the momentary turning rate, and
  // therefore the probability that a change of direction occurs at a growth cone.
  // In an event-driven implementation of the simulation, the turning_model
  // objects can also select growth cones at which turning occurs and
  // send those events to the queue.
  // Selection of a turning_model is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all turning_model objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // turning_model objects.
protected:
  turning_model_base * contributing;
  double contributingweight;
  /*  struct turning_model_base_parameters {
    sometype something;
    turning_model_base_parameters(sometime _something): something(_something) {}
    } base_parameters; */
  virtual ~turning_model_base() { if (contributing) contributing->delete_shared(); }
public:
  turning_model_base(turning_model_base * tmcontrib, double & tmweight, turning_model_base & schema);
  turning_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual turning_model_base * clone() = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_turning(terminal_segment & ts1, terminal_segment & ts2) = 0;
  virtual void handle_turning(terminal_segment & ts) = 0;
  //  virtual spatial predict(double weight, terminal_segment * ts, neuron * e) = 0;
  //  virtual void turning(terminal_segment * ts, neuron * e) = 0;
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef turning_model_base * turning_model_base_ptr;

class TSTM_base {
  // Each terminal segment is designated a TSTM object, the
  // terminal segment turning model, which determines if a
  // turn event takes place at a specific terminal segment.
  // This can involve local data. For example in the case of
  // mechanistic models the presence of another object in the
  // near vicinity ahead of the growth cone could elicit a
  // change of direction. The subsequent direction is determined
  // by a direction model.
  // Selection of a TSTM is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all TSTM objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // TSTM objects.
protected:
  // cache variables
  terminal_segment * cache_ts;
  double cache_p_i;
protected:
  TSTM_base * contributing;
  double contributingweight;
  //struct TSTM_base_parameters {
  //  TSTM_base_parameters() {}
  //} base_parameters;
  virtual ~TSTM_base() { if (contributing) contributing->delete_shared(); }
public:
  TSTM_base(TSTM_base * tstmcontrib, double & tstmweight, TSTM_base & schema);
  TSTM_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSTM_base * clone() = 0; //terminal_segment * _ts = NULL) = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual double predict_turn() = 0; //double weight = 1.0) = 0;
  virtual double predict(double weight, double & probability);
  virtual bool turn_event(double r_selected, terminal_segment & ts, double p_i);
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef TSTM_base * TSTM_base_ptr;

class none_TSTM: public TSTM_base {
  // This never signals a turn event. This model can be used
  // to selectively disable turns on select fibre structures.
public:
  none_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, none_TSTM & schema);
  none_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSTM_base * clone(); //terminal_segment * _ts = NULL);
  //virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  //virtual void handle_turning(terminal_segment & ts);
  virtual double predict_turn(); // double weight = 1.0);
  virtual String report_parameters_specific();
};

class each_dt_TSTM: public TSTM_base {
  // This TSTM signals a turn event whenever the simulation
  // time exceeds the last simulation time at which a turn
  // event was determined. In effect, one turn is elicited
  // per simulation time step.
protected:
  struct each_dt_TSTM_parameters {
    double t_prev;
    each_dt_TSTM_parameters(double _t = 0.0): t_prev(_t) {}
  } parameters;
public:
  each_dt_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, each_dt_TSTM & schema);
  each_dt_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSTM_base * clone(); //terminal_segment * _ts = NULL);
  //virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  //virtual void handle_turning(terminal_segment & ts);
  virtual double predict_turn(); // double weight = 1.0);
  virtual String report_parameters_specific();
};

class linear_rate_TSTM: public TSTM_base {
  // This TSTM probabilistically signals a turn event at a
  // linear rate.
protected:
  // cache variables
protected:
  struct linear_rate_TSTM_parameters {
    double turn_rate;
    linear_rate_TSTM_parameters(double _turnrate): turn_rate(_turnrate) {}
  } parameters;
public:
  linear_rate_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, linear_rate_TSTM & schema);
  linear_rate_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSTM_base * clone(); //terminal_segment * _ts = NULL);
  //virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  //virtual void handle_turning(terminal_segment & ts);
  virtual double predict_turn(); //double weight = 1.0);
  virtual String report_parameters_specific();
};

class branch_coupled_TSTM: public TSTM_base {
  // This TSTM probabilistically signals a turn event with
  // a likelihood that is specified factor of the branching
  // probability.
protected:
  struct branch_coupled_TSTM_parameters {
    double turnfactor;
    branch_coupled_TSTM_parameters(double _turnfactor): turnfactor(_turnfactor) {}
  } parameters;
public:
  branch_coupled_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, branch_coupled_TSTM & schema);
  branch_coupled_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual TSTM_base * clone(); //terminal_segment * _ts = NULL);
  //virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  //virtual void handle_turning(terminal_segment & ts);
  virtual double predict_turn(); //double weight = 1.0);
  virtual String report_parameters_specific();
};

// The turn angle is determined by a Direction Models in the axon_direction_model.hh/cc files.

// functions
turning_model_base * turning_model_selection(String & label, Command_Line_Parameters & clp, turning_model_base_ptr & tmbptr);

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
extern turning_model default_turning_model; // identifier, changed when a different universal model is set
extern String universal_turning_model_root; // Used only to report if chaining at the universal set level.
extern TSTM default_TSTM; // identifier, changed when a different universal model is set
extern String universal_TSTM_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
typedef turning_model_base * turning_model_base_ptr;
typedef TSTM_base_ptr * region_TSTM_base_ptr;
extern region_TSTM_base_ptr * TSTM_region_subset_schemas;

// functions
TSTM_base * TSTM_selection(String & label, Command_Line_Parameters & clp, TSTM_base * superior_set, TSTM_base_ptr & tstmbptr);

#endif
// __TURNING_MODELS_HH

