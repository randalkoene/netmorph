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
// axon_direction_model.hh
// Randal A. Koene, 20051021

/* Documentation:
   This implements base and derived classes for models of axon direction during
   network development that are based on curvature statistics (rather than
   causal hypotheses about environmental interactions).
   The base model provides a reference to the history of an axon's growth.

   An obvious place in which to link to an axon direction model is in the
   terminal segment objects of axons.
 */

#ifndef __AXON_DIRECTION_MODEL_HH
#define __AXON_DIRECTION_MODEL_HH

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
class direction_history_selection;

// When speed is prioretized over space then histories contain actual
// data with which a direction can be calculated. Otherwise, histories
// are largely restricted to pointers to the fibre segment directly
// after the last branch and the number of segments to skip (i.e. not
// to include during calculations), with which the necessary data can
// be computed.
// [***NOTE] This is now defined in the Makefile.
//#define ADM_PRIORITIZE_SPEED_OVER_SPACE

enum direction_history_selection_type { no_dhs, random_truncation_dhs };

class general_direction_model_parameters_interface: public CLP_Modifiable {
  // This is used to set and provide general (default) parameter settings
  // used in functions of direction models and related models such as
  // direction history selection models.
public:
  general_direction_model_parameters_interface(): direction_history_selection_model(no_dhs), random_truncation_dh_probability(0.1), random_truncation_dh_minfraction(0.5), random_truncation_dh_maxfraction(1.0) {}
  virtual ~general_direction_model_parameters_interface() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  direction_history_selection_type direction_history_selection_model;
  double random_truncation_dh_probability;
  double random_truncation_dh_minfraction;
  double random_truncation_dh_maxfraction;
};

extern general_direction_model_parameters_interface general_direction_model_parameters;

class axon_direction_history {
  friend class random_truncation_dh_selection;
  // This is a utility class that can autonomously handle data in a preferred
  // manner. Derivatives may make use of arrays, linked lists, or may even
  // calculate required data on-demand.
  // The history tracks those parts of a fibre that may affect the direction
  // of growth (in the statistical/phenomenological model). The intention is
  // that elements dropped from the history are not recovered. The actual
  // elements that are taken into account and their weighting depends on the
  // model used to compute the effect of the history on direction of growth,
  // which may involve the entire history, a certain number of turns, a
  // certain distance from the growth cone, etc.
  // Whether random or non-random events cause elements to be dropped from
  // the history, or whether they only affect the portion or weighting of
  // elements taken into account determines if a shortening of the history
  // that affects the direction of growth can be reversed by subsequent
  // events.
  // The history does not contain information about the fibre segment
  // associated with the terminal segment object, since the history is
  // only needed when the direction must be determined for a new
  // terminal segment, so that the associated fibre segment data can not
  // provide contributing information.
  // The update_effective_history() function provides an opportunity to
  // update the selection of a portion of the segment history that is
  // considered effective, if a direction history selection model has been
  // specified in dhs.
protected:
  direction_history_selection * dhs; // an optional model of history modification
  double effective_length; // in terms of fibre length
public:
  axon_direction_history();
  virtual ~axon_direction_history() {}
  virtual void add(fibre_segment * fs, double len) = 0;
  virtual void remove_tail() = 0;
  virtual void clear() { effective_length = 0.0; }
  virtual void restart_at_nearest() = 0;
  virtual void next_further() = 0;
  virtual double segment_length() = 0;
  virtual double segment_radius_squared() = 0;
  virtual spatial & segment_unit_direction() = 0;
  virtual void update_effective_history() = 0;
};

class axon_direction_element: public PLLHandle<axon_direction_element> {
public:
  axon_direction_element(double _l): l(_l) {}
  axon_direction_element(fibre_segment * _fs);
  spatial u;
  double l;
  //double r; // This needs updating each time another element is added.
};

class pll_axon_direction_history: public axon_direction_history {
protected:
  PLLRoot<axon_direction_element> adr;
  axon_direction_element * segidx;
  axon_direction_element nullelement;
  double remaining;
public:
  pll_axon_direction_history(): segidx(&nullelement), nullelement(-1.0) {}
  virtual ~pll_axon_direction_history() {}
  virtual void add(fibre_segment * fs, double len) {
    adr.link_before(new axon_direction_element(fs));
    effective_length += len;
  }
  virtual void remove_tail() { if (adr.tail()) adr.tail()->remove(); }
  virtual void clear();
  virtual void restart_at_nearest() {
    segidx = adr.tail();
    if (!segidx) segidx = &nullelement;
    remaining = effective_length;
  }
  virtual void next_further() {
    remaining -= segidx->l;
    if (remaining<=0.0) segidx = &nullelement;
    else {
      segidx = segidx->Prev();
      if (!segidx) segidx = &nullelement;
    }
  }
  virtual double segment_length() {
    if (segidx->l < remaining) return segidx->l;
    else return remaining;
  }
  // [***INCOMPLETE] Segment radius is temporarily a fixed value.
  virtual double segment_radius_squared() { return 1.0; }
  virtual spatial & segment_unit_direction() { return segidx->u; }
  virtual void update_effective_history();
};

class parse_axon_direction_history: public axon_direction_history {
protected:
  fibre_segment * firstafterbranch;
  int skipafterbranch;
  int segmentsafterbranch;
  int segidx;
  fibre_segment * segidxfs;
  spatial u;
  fibre_segment * seek(int idx);
  double remaining;
public:
  parse_axon_direction_history(fibre_segment & fab): firstafterbranch(&fab), skipafterbranch(0), segmentsafterbranch(0), segidx(-1), segidxfs(NULL) {}
  ~parse_axon_direction_history() {}
  virtual void add(fibre_segment * fs, double len) { segmentsafterbranch++; effective_length += len; }
  virtual void remove_tail() { if (segmentsafterbranch>skipafterbranch) skipafterbranch++; }
  virtual void clear(fibre_segment & fab);
  virtual void restart_at_nearest() {
    if (segmentsafterbranch>skipafterbranch) {
      segidx = segmentsafterbranch-1;
      segidxfs = seek(segidx);
    } else segidx = -1;
    remaining = effective_length; 
  }
  virtual void next_further();
  virtual double segment_length();
  // [***INCOMPLETE] Segment radius is temporarily a fixed value.
  virtual double segment_radius_squared() { return 1.0; }
  virtual spatial & segment_unit_direction();
  virtual void update_effective_history();
};

class direction_history_selection {
  // This is the general class of models that vary the effective segment
  // history used to determine the predicted direction.
public:
  direction_history_selection() {}
  virtual ~direction_history_selection() {}
  virtual void history_reset() = 0;
  virtual void effective_history(axon_direction_history * adh) = 0;
};

class random_truncation_dh_selection: public direction_history_selection {
  // This simple segment history selection model uses only a uniform random
  // function to decide that the segment history should be truncated and
  // another random value to determine at which fibre length the segment
  // history should be truncated.
public:
  random_truncation_dh_selection() {}
  virtual void history_reset() {}
  virtual void effective_history(axon_direction_history * adh);
};

class direction_model_base {
  // Unless the alternative (legacy) direction models are selected, in which
  // case a NULL pointer is assigned, each terminal_segment object is
  // designated an axon_direction_model object. That way, there is no need
  // to specify the specific history to use each time an event must be
  // handled.
  // The selection of the type of axon_direction_model to use (if any) can
  // be done at the level of all neurons (globally), per neuron or per
  // pre/postsynaptic structure.
  // [***INCOMPLETE] If I want to insure that all direction models define
  // the correct set of constructors and go about procedures correctly
  // then I should create a template for direction models that all inherit
  // this direction_model_base.
protected:
  // The following is shared among direction models in order to produce
  // random veering about a predicted direction. If there should also be
  // direction models with no such randomness then create a
  // direction_model_base_base without these elements.
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  linear_pdf admpdf; // [***INCOMPLETE] The pdf is not yet an option.
#endif
  direction_model_base * contributing;
  double contributingweight;
  void turn(terminal_segment * ts, spatial & predicted, double initlen);
  struct direction_model_base_parameters {
    double veeranglemin, veeranglemax; // perturbed expected elongation
    direction_model_base_parameters(double _veeranglemin, double _veeranglemax): veeranglemin(_veeranglemin), veeranglemax(_veeranglemax) {}
  } base_parameters;
  virtual ~direction_model_base() { if (contributing) contributing->delete_shared(); }
public:
  direction_model_base(direction_model_base * dmcontrib, double & dmweight, direction_model_base & schema);
  direction_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual direction_model_base * clone() = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2) = 0;
  virtual void handle_turning(terminal_segment & ts) = 0;
  virtual spatial predict(double weight, terminal_segment * ts, neuron * e) = 0;
  virtual void direction(terminal_segment * ts, neuron * e, double initlen = 0.0) = 0;
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef direction_model_base * direction_model_base_ptr;

class radial_direction_model: public direction_model_base {
  // This is the simplest model of growth cone direction, which
  // makes only the optional assumption (based on a hypothesis
  // by Samsonovich) that dendrite growth cones tend to return
  // to a growth direction that is radial to the soma.
protected:
  struct radial_dm_parameters {
    bool Samsonovich_hypothesis;
    radial_dm_parameters(bool samshyp): Samsonovich_hypothesis(samshyp) {}
  } radial_parameters;
  void predict_direction(spatial & predicted, neuron * n, terminal_segment * ts);
public:
  radial_direction_model(direction_model_base * dmbcontrib, double & dmbweight, radial_direction_model & schema);
  radial_direction_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  ~radial_direction_model() {}
  virtual direction_model_base * clone();
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual spatial predict(double weight, terminal_segment * ts, neuron * e);
  virtual void direction(terminal_segment * ts, neuron * e, double initlen = 0.0);
  virtual String report_parameters_specific();
};

class tension_direction_model: public direction_model_base {
  // This default model is based on the following hypotheses:
  // 1. The effect that each segment in the segment history has on the
  //    predicted direction of growth can be summed through a vector summation
  //    in which individual components have an amplitude that is represented
  //    by the length of their vectors.
  // 2. The amplitude of a specific segment's effect on the predicted direction
  //    can be calculated by a model that is analogous to gravitation, where
  //    the distance is measured from the "center of mass" of a segment to the
  //    tip of the growth cone along the path of the intervening fibre
  //    segments.
  // 3. A separately chosen model updates the effective segment history and
  //    the optiona ability to retrieve truncated segments.
  // (See further comments in axon_direction_model::direction().)
protected:
  axon_direction_history * adh;
  struct tension_dm_parameters {
    double history_power;
    tension_dm_parameters(double hp): history_power(hp) {}
  } tension_parameters;
  void predict_direction(spatial & predicted);
public:
  //tension_direction_model(axon_direction_history & _adh): adh(&_adh) {}
  tension_direction_model(axon_direction_history & _adh, direction_model_base * dmbcontrib, double & dmbweight,  tension_direction_model & schema);
  tension_direction_model(axon_direction_history & _adh, String & thislabel, String & label, Command_Line_Parameters & clp);
  ~tension_direction_model() { delete adh; }
  virtual direction_model_base * clone();
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual spatial predict(double weight, terminal_segment * ts, neuron * e);
  virtual void direction(terminal_segment * ts, neuron * e, double initlen = 0.0);
  virtual String report_parameters_specific();
};

class cell_attraction_direction_model: public direction_model_base {
  // This is a very simple direction model that uses attraction associations
  // beteween specific cells to predict the direction of growth.
  // This model can be used in conjunction with the axon_direction_model.
protected:
  void predict_direction(spatial & predicted, neuron * n, spatial & growthcone);
public:
  //cell_attraction_direction_model() {}
  cell_attraction_direction_model(direction_model_base * dmbcontrib, double & dmbweight, cell_attraction_direction_model & schema);
  cell_attraction_direction_model(String & thislabel, String & label, Command_Line_Parameters & clp): direction_model_base(thislabel,label,clp) {}
  virtual direction_model_base * clone();
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual spatial predict(double weight, terminal_segment * ts, neuron * e);
  virtual void direction(terminal_segment * ts, neuron * e, double initlen = 0.0);
  virtual String report_parameters_specific();
};

class vector_direction_model: public direction_model_base {
  // This growth cone direction model expects the direction
  // of growth to be along a specified vector.
  // This model can be useful, for example in a chain of
  // direction models that guide apical dendrite trunk growth
  // for demonstration purposes.
protected:
  struct vector_dm_parameters {
    spatial direction;
    vector_dm_parameters() {}
    vector_dm_parameters(spatial d): direction(d) {}
  } vector_parameters;
  void predict_direction(spatial & predicted, neuron * n, terminal_segment * ts);
public:
  vector_direction_model(direction_model_base * dmbcontrib, double & dmbweight, vector_direction_model & schema);
  vector_direction_model(String & thislabel, String & label, Command_Line_Parameters & clp);
  ~vector_direction_model() {}
  virtual direction_model_base * clone();
  virtual void handle_branching(terminal_segment & ts1, terminal_segment & ts2);
  virtual void handle_turning(terminal_segment & ts);
  virtual spatial predict(double weight, terminal_segment * ts, neuron * e);
  virtual void direction(terminal_segment * ts, neuron * e, double initlen = 0.0);
  virtual String report_parameters_specific();
};

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
extern direction_model default_direction_model; // identifier, changed when a different universal model is set
extern String universal_direction_model_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
typedef direction_model_base_ptr * region_direction_model_base_ptr;
extern region_direction_model_base_ptr * direction_model_region_subset_schemas;

// functions
direction_model_base * direction_model_selection(String & label, Command_Line_Parameters & clp, direction_model_base * superior_set, direction_model_base_ptr & dmbptr);

#endif
