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
// axon_direction_model.cc
// Randal A. Koene, 20051021

/* Documentation:
   This implements base and derived classes for models of axon direction during
   network development that are based on curvature statistics (rather than
   causal hypotheses about environmental interactions).
   The base model provides a reference to the history of an axon's growth.
 */

#include <math.h>
#include "axon_direction_model.hh"
#include "neuron.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"
#include "environment_physics.hh"
#include "Spatial_Presentation.hh"

/* [***INCOMPLETE] Here I can include functions that parse the command line for
   general parameters to be used in specific direction models if no neuron or
   terminal segment specific parameters are given.
 */

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
direction_model default_direction_model = radial_dm; // identifier, changed when a different universal model is set
String universal_direction_model_root; // Used only to report if chaining at the universal set level.

general_direction_model_parameters_interface general_direction_model_parameters;

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_direction_model_base_ptr * direction_model_region_subset_schemas = NULL;


void general_direction_model_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("dirhistory_selection"))>=0) {
    if (downcase(clp.ParValue(n))==String("random_truncation")) direction_history_selection_model = random_truncation_dhs;
    else direction_history_selection_model = no_dhs;
  }
  if ((n=clp.Specifies_Parameter("general_rtdhs_probability"))>=0) random_truncation_dh_probability = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("general_rtdhs_minfraction"))>=0) {
    random_truncation_dh_minfraction = atof(clp.ParValue(n));
    if (random_truncation_dh_minfraction<0.0) random_truncation_dh_minfraction = 0.0;
    if (random_truncation_dh_minfraction>1.0) random_truncation_dh_minfraction = 1.0;
  }
  if ((n=clp.Specifies_Parameter("general_rtdhs_maxfraction"))>=0) {
    random_truncation_dh_maxfraction = atof(clp.ParValue(n));
    if (random_truncation_dh_maxfraction<random_truncation_dh_minfraction) random_truncation_dh_maxfraction = random_truncation_dh_minfraction;
    if (random_truncation_dh_maxfraction>1.0) random_truncation_dh_maxfraction = 1.0;
  }
}

String general_direction_model_parameters_interface::report_parameters() {
  String res("Fiber growth direction hypotheses:\n");
  switch (direction_history_selection_model) {
  case random_truncation_dhs:
    res += "  random fiber segment history truncation.\n    probability = ";
    res += String(random_truncation_dh_probability,"%.3f, truncated to a fraction between ");
    res += String(random_truncation_dh_minfraction,"%.3f and ");
    res += String(random_truncation_dh_maxfraction,"%.3f.\n");
    break;;
  default: res += "  fiber segment history is not modified.\n";
  }
  return res;
}

axon_direction_history::axon_direction_history(): dhs(NULL), effective_length(0.0) {
  // [***INCOMPLETE] This currently uses only the general default, even though
  // axon_direction_history objects are allocated one per terminal segment, so
  // that even terminal segment specific choices are possible.
  switch (general_direction_model_parameters.direction_history_selection_model) {
  case random_truncation_dhs: dhs = new random_truncation_dh_selection(); break;;
  default: break;; // to avoid a compilation warning message about no_dhs
  }
}

axon_direction_element::axon_direction_element(fibre_segment * _fs) {
  u = _fs->P1;
  u -= _fs->P0;
  l = _fs->Length();
  if (l!=0.0) u /= l;
  //r = l/2.0; // [***INCOMPLETE] May depend on integration method and r dependence function.
}

void pll_axon_direction_history::clear() {
  axon_direction_history::clear();
  adr.clear();
  if (dhs) dhs->history_reset();
}

void pll_axon_direction_history::update_effective_history() {
  if (dhs) dhs->effective_history(this);
}

fibre_segment * parse_axon_direction_history::seek(int idx) {
  fibre_segment * fs = firstafterbranch;
  for (; idx>0; idx--) {
    if (fs->Branch1()) fs = fs->Branch1();
    else if (fs->Branch2()) fs = fs->Branch2();
    else return NULL;
  }
  return fs;
}

void parse_axon_direction_history::clear(fibre_segment & fab) {
  axon_direction_history::clear();
  firstafterbranch = &fab;
  skipafterbranch = 0;
  segmentsafterbranch = 0;
  if (dhs) dhs->history_reset();
}

void parse_axon_direction_history::next_further() {
  if (segidx>=0) remaining -= segidxfs->Length();
  if ((remaining>0.0) && (segidx>skipafterbranch)) {
    segidx--;
    segidxfs = seek(segidx);
  } else segidx = -1;
}

double parse_axon_direction_history::segment_length() {
  if (segidx<0) return -1.0;
  double l = segidxfs->Length();
  if (l<remaining) return l;
  return remaining;
}

spatial & parse_axon_direction_history::segment_unit_direction() {
  if (segidx<0) u.set_all(0.0,0.0);
  else {
    u = segidxfs->P1;
    u -= segidxfs->P0;
    u /= segidxfs->Length();
  }
  return u;
}

void parse_axon_direction_history::update_effective_history() {
  if (dhs) dhs->effective_history(this);
}

void random_truncation_dh_selection::effective_history(axon_direction_history * adh) {
  // [***INCOMPLETE] This uses the general parameter values, not any specific
  // parameter values set for this history selection object.
  if (X_turn.get_rand_real1()<general_direction_model_parameters.random_truncation_dh_probability) {
    adh->effective_length *= X_turn.get_rand_range_real1(general_direction_model_parameters.random_truncation_dh_minfraction,general_direction_model_parameters.random_truncation_dh_maxfraction);
  }
}

direction_model_base::direction_model_base(direction_model_base * dmcontrib, double & dmweight, direction_model_base & schema):
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  admpdf(schema.admpdf),
#endif
  contributing(NULL), contributingweight(1.0), base_parameters(schema.base_parameters) {
  // [***INCOMPLETE] More specific choices for parameters, e.g. the veer
  // limits of the probability distribution, can be detected here. Those may
  // be specified for the network, the type of neuron, the neuron, axons or
  // dendrites, even for a specific arbor or growth cone (terminal segment).
  // The default values are set to match global parameters defined for this
  // particular or a related model (e.g. veeranglemin/max=turnanglemin/max).
  // Can I suggest inlining this constructor somehow (without moving it and
  // the include files that make turnanglemin/max available to the .hh file?
  if (dmcontrib) { // Clone chained models
    contributingweight = dmweight;
    contributing = dmcontrib->clone();
  }
}

direction_model_base::direction_model_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  admpdf(X_turn),
  contributing(NULL), contributingweight(1.0), base_parameters(turnanglemin,turnanglemax) {
  // [***INCOMPLETE] More specific choices for parameters, e.g. the veer
  // limits of the probability distribution, can be detected here. Those may
  // be specified for the network, the type of neuron, the neuron, axons or
  // dendrites, even for a specific arbor or growth cone (terminal segment).
  // The default values are set to match global parameters defined for this
  // particular or a related model (e.g. veeranglemin/max=turnanglemin/max).
  // Can I suggest inlining this constructor somehow (without moving it and
  // the include files that make turnanglemin/max available to the .hh file?
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  int n;
  double vmin = base_parameters.veeranglemin;
  double vmax = base_parameters.veeranglemax;
  if ((n=clp.Specifies_Parameter(thislabel+"veeranglemin"))>=0) {
    vmin = atof(clp.ParValue(n));
    if (vmin<0.0) {
      warning("Warning: "+thislabel+"veeranglemin="+clp.ParValue(n)+" not permitted, value must be >=0.0, not modified\n");
      vmin = base_parameters.veeranglemin;
    }
  }
  if ((n=clp.Specifies_Parameter(thislabel+"veeranglemax"))>=0) {
    vmax = atof(clp.ParValue(n));
    if (vmax<=0.0) {
      warning("Warning: "+thislabel+"veeranglemax="+clp.ParValue(n)+" not permitted, value must be >0.0, not modified\n");
      vmax = base_parameters.veeranglemax;
    }
  }
  if (vmax<=vmin) warning("Warning: "+thislabel+"veeranglemax must be > "+thislabel+"veeranglemin, not modified\n");
  else {
    base_parameters.veeranglemin = vmin;
    base_parameters.veeranglemax = vmax;
  }
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  admpdf.setup(base_parameters.veeranglemax-base_parameters.veeranglemin);
#endif
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"dm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!direction_model_selection(label,clp,contributing,contributing)) error("Error in direction_model_base(): No contributing direction model found for "+label+".\n");
}

String direction_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting tocontributing models.
  String res(report_parameters_specific());
  res += String(base_parameters.veeranglemin,"\n  [veeranglemin=%.2f] ");
  res += String(base_parameters.veeranglemax,"\n  [veeranglemax=%.2f] ");
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

void direction_model_base::turn(terminal_segment * ts, spatial & predicted, double initlen) {
  // Using the resulting preferred direction found in direction(), apply
  // random deviation.
  // ts points to the terminal segment before the turn
  // [***INCOMPLETE] Create a class as in DIL#20051129041848.1 and use it here.
#ifndef ADM_PRIORITIZE_SPEED_OVER_SPACE
  linear_pdf admpdf(base_parameters.veeranglemax - base_parameters.veeranglemin); // [***INCOMPLETE] The pdf is not yet an option.
#endif
#ifdef VECTOR3D
  // See TL#200511261312.2.
  //  double y = max_y*random_double();
  // double veerangle = (negadm_b + sqrt(SQadm_b + (twoa*y)))/a;
  double veerangle = admpdf.random_positive() + base_parameters.veeranglemin;
  spatial S1xyz;
  double rollangle = 2.0*M_PI*X_turn.get_rand_real1();
  spatial S1(1.0,rollangle,veerangle);
  prepare_parent_coordinate_system(predicted);
  daughter_to_network_coordinate_system(S1,S1xyz);
#endif
#ifdef VECTOR2D
  // See TL#200511261312.2.
  //double r = twomax_y*random_double();
  //double y = r;
  //if (r>max_y) y -= max_y;
  //double veerangle = (negadm_b + sqrt(SQadm_b + (twoa*y)))/a;
  double veerangle = admpdf.random_selection();
  if (veerangle>=0.0) veerangle += base_parameters.veeranglemin;
  else veerangle -= base_parameters.veeranglemin;
  spatial S1xyz(0.0,predicted.Y()+veerangle);
#endif
  // The following are the morphological updates applied:
  // 1. specify a turn that also results in a new terminal segment object
  S1xyz.set_X(initlen);
#ifdef VECTOR3D
  extra_fibre_data efd(rollangle);
  ts->turn_without_update_of_segment_vector(S1xyz,efd);
#endif
#ifdef VECTOR2D
  ts->turn_without_update_of_segment_vector(S1xyz);
#endif
  // 2. update turn (curvature) statistics
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
#ifdef INCLUDE_TURN_ANGLE_SAMPLING
  const stat_labels angledatalabel[2] = {axons_turn_angles,dendrites_turn_angles};
  if (sampled_output->netstats[angledatalabel[growing_dendrites]].Collect()) {
    STANDARDIZE_ANGLE_DATA(randomangle);
    // [***INCOMPLETE] The veerangle here is not the actual turn angle!
    sampled_output->netstats[angledatalabel[growing_dendrites]].Collect_Data(veerangle);
  }
#endif
#endif
}

tension_direction_model::tension_direction_model(axon_direction_history & _adh, direction_model_base * dmbcontrib, double & dmbweight,  tension_direction_model & schema): direction_model_base(dmbcontrib,dmbweight,schema), adh(&_adh), tension_parameters(schema.tension_parameters) {
}

tension_direction_model::tension_direction_model(axon_direction_history & _adh, String & thislabel, String & label, Command_Line_Parameters & clp): direction_model_base(thislabel,label,clp), adh(&_adh), tension_parameters(2.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"history_power"))>=0) tension_parameters.history_power =  atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
}

String tension_direction_model::report_parameters_specific() {
  String res(" segment history tension (history_power = ");
  res += String(tension_parameters.history_power,"%.3f)");
  return res;
}

direction_model_base * tension_direction_model::clone() {
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  return new tension_direction_model(*(new pll_axon_direction_history()),contributing,contributingweight,*this);
#else
  // [***INCOMPLETE] Clearly, this does not work as implemented here,
  // because ts2 is not available. Consider if I should keep the space
  // conserving version at all.
  return new tension_direction_model(*(new parse_axon_direction_history(ts2.TerminalSegment())),contributing,contributingweight,*this)
#endif
}

void tension_direction_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  // After branching the terminal segments on each branch need a
  // history. The first branch receives this history object after
  // clearing. The second branch receives a history object of
  // the same type.
  // [***INCOMPLETE] Using a #define to select the history class
  // defeats the purpose of compiling the code for multiple.
  ts2.set_direction_model(clone());
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  adh->clear();
#else
  adh->clear(ts1.TerminalSegment());
#endif
  ts1.set_direction_model(this);
}

void tension_direction_model::handle_turning(terminal_segment & ts) {
  // ts is the new terminal_segment.
  // Here, any computations that must be made on data, e.g. for the maintenance
  // of valid values in the history can be done.
  adh->update_effective_history();
  ts.set_direction_model(this);
}

void tension_direction_model::predict_direction(spatial & predicted) {
  // Calculate historical segment effects, multiply with unit direction
  // vectors and sum.
  // The historical segment effect is modeled in a manner analogous to
  // gravitation in this class.
  // It is expected that predicted is initialized to zero, unless the
  // following calculations should add to it.
  adh->restart_at_nearest();
  double distancebyfibre = 0.0, L; 
  for (; (L=adh->segment_length())>=0.0; adh->next_further()) if (L>0.0) {
    double mass = L * M_PI * adh->segment_radius_squared();
    double distance = distancebyfibre + (L/2.0);
    double relativegravitation = mass / pow(distance,tension_parameters.history_power);
    predicted += relativegravitation*adh->segment_unit_direction();
    distancebyfibre += L;
  }
}

spatial tension_direction_model::predict(double weight, terminal_segment * ts, neuron * e) {
  // Called in a chain of direction models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  spatial predicted;
  predict_direction(predicted);
  predicted.normalize();
  predicted *= weight;
  if (contributing) predicted += contributing->predict(contributingweight,ts,e);
  return predicted;
}

void tension_direction_model::direction(terminal_segment * ts, neuron * e, double initlen) {
  // ts points to the terminal segment before the turn
  // e points to the neuron to which the arbor belongs
  // The following are the preparations of this axon direction model:
  // 1. update the spatial parameters of the terminal fibre segment
  ts->update_segment_vector(); // hence, use ts->turn_free()
  // 2. include the terminal fibre segment in the history
  adh->add(ts->TerminalSegment(),ts->Length());
  // The following are the functions of this axon direction model:
  // 1. determine "pull" strengths of the preceding segments
  // NOTE: This version uses a gravitation analogy for the calculation of the
  // amplitude of the effect that each segment has on the predicted direction.
  // See TL#200511180719.1 and TL#200511181155.1.
  // 2. using pull strength as length sum direction vectors
  spatial predicted;
  predict_direction(predicted);
  if (contributing) {
    predicted.normalize(); // only necessary here
    predicted += contributing->predict(contributingweight,ts,e);
  }
  predicted.convert_to_spherical();
  if (predicted.X()==0.0) { // [***INCOMPLETE] See TL.
#ifdef DEBUGGING_DIRECTION
    debugging.tension_direction_model_no_direction++;
#endif
    predicted = ts->AngularCoords();
  }
  // [*** NOTE:] The probability distribution of the veer angle is currently
  // not dependent on the history size, which currently has an effect only
  // on the predicted direction (i.e. mean).
  // 3. make the turn with probabilistic veering around the predicted direction
  turn(ts,predicted,initlen);
}

radial_direction_model::radial_direction_model(direction_model_base * dmbcontrib, double & dmbweight,  radial_direction_model & schema): direction_model_base(dmbcontrib,dmbweight,schema), radial_parameters(schema.radial_parameters) {
}

radial_direction_model::radial_direction_model(String & thislabel, String & label, Command_Line_Parameters & clp): direction_model_base(thislabel,label,clp), radial_parameters(false) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"Samsonovich_hypothesis"))>=0) radial_parameters.Samsonovich_hypothesis = (downcase(clp.ParValue(n))==String("true"));
  // Chaining of models is taken care of in the base constructor.
}

String radial_direction_model::report_parameters_specific() {
  String res(" radial");
  if (radial_parameters.Samsonovich_hypothesis) res += " (with Samsonovich hypothesis)";
  else res += " (no Samsonovich hypothesis)";
  return res;
}

direction_model_base * radial_direction_model::clone() {
  return new radial_direction_model(contributing,contributingweight,*this);
}

void radial_direction_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_direction_model(clone());
  ts1.set_direction_model(this);
}

void radial_direction_model::handle_turning(terminal_segment & ts) {
  // ts is the new terminal_segment.
  ts.set_direction_model(this);
}

void radial_direction_model::predict_direction(spatial & predicted, neuron * n, terminal_segment * ts) {
  // This model function predicts a direction that is the same as the parent
  // direction, or a radial direction if the Samsonovich hypothesis is applied.
  // (This model is the replacement for the legacy direction code used in older
  // versions of NETMORPH, see TL#200808041538.1.)
  if (radial_parameters.Samsonovich_hypothesis) {
    fibre_segment * fs = ts->TerminalSegment();
    predicted = fs->P1;
    predicted -= fs->N()->Pos();
  } else {
    predicted = ts->AngularCoords();
    predicted.convert_from_spherical(); // to make summation of normalized contributions and then conversion back to spherical possible in direction()
  }
}

spatial radial_direction_model::predict(double weight, terminal_segment * ts, neuron * e) {
  // Called in a chain of direction models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  spatial predicted;
  predict_direction(predicted,e,ts);
  predicted.normalize();
  predicted *= weight;
  if (contributing) predicted += contributing->predict(contributingweight,ts,e);
  return predicted;
}

void radial_direction_model::direction(terminal_segment * ts, neuron * e, double initlen) {
  // ts points to the terminal segment before the turn
  // e points to the neuron to which the arbor belongs
  ts->update_segment_vector(); // hence, use ts->turn_free()
  spatial predicted;
  predict_direction(predicted,e,ts);
  if (contributing) {
    predicted.normalize(); // only necessary here
    predicted += contributing->predict(contributingweight,ts,e);
  }
  predicted.convert_to_spherical();
  if (predicted.X()==0.0) { // [***INCOMPLETE] See TL.
    predicted = ts->AngularCoords();
  }
  // make the turn with probabilistic veering around the predicted direction
  turn(ts,predicted,initlen);
}

vector_direction_model::vector_direction_model(direction_model_base * dmbcontrib, double & dmbweight,  vector_direction_model & schema): direction_model_base(dmbcontrib,dmbweight,schema), vector_parameters(schema.vector_parameters) {
}

vector_direction_model::vector_direction_model(String & thislabel, String & label, Command_Line_Parameters & clp): direction_model_base(thislabel,label,clp), vector_parameters() {
#ifdef VECTOR3D
  get_spatial(clp,thislabel+"direction",vector_parameters.direction);
#endif
#ifdef VECTOR2D
  warning("Warning: In 2D, the expected direction for the vector direction model is always (0,0).\n");
#endif
  // Chaining of models is taken care of in the base constructor.
}

String vector_direction_model::report_parameters_specific() {
  String res(" vector");
  double x,y,z;
  vector_parameters.direction.get_all(x,y,z);
  res += " (direction = [" + String(x,"%.3f,")+String(y,"%.3f,")+String(z,"%.3f])");
  return res;
}

direction_model_base * vector_direction_model::clone() {
  return new vector_direction_model(contributing,contributingweight,*this);
}

void vector_direction_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_direction_model(clone());
  ts1.set_direction_model(this);
}

void vector_direction_model::handle_turning(terminal_segment & ts) {
  // ts is the new terminal_segment.
  ts.set_direction_model(this);
}

void vector_direction_model::predict_direction(spatial & predicted, neuron * n, terminal_segment * ts) {
  // This model function predicts a direction along a specified vector.
  predicted = vector_parameters.direction;
}

spatial vector_direction_model::predict(double weight, terminal_segment * ts, neuron * e) {
  // Called in a chain of direction models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  spatial predicted;
  predict_direction(predicted,e,ts);
  predicted.normalize();
  predicted *= weight;
  if (contributing) predicted += contributing->predict(contributingweight,ts,e);
  return predicted;
}

void vector_direction_model::direction(terminal_segment * ts, neuron * e, double initlen) {
  // ts points to the terminal segment before the turn
  // e points to the neuron to which the arbor belongs
  ts->update_segment_vector(); // hence, use ts->turn_free()
  spatial predicted;
  predict_direction(predicted,e,ts);
  if (contributing) {
    predicted.normalize(); // only necessary here
    predicted += contributing->predict(contributingweight,ts,e);
  }
  predicted.convert_to_spherical();
  if (predicted.X()==0.0) { // [***INCOMPLETE] See TL.
    predicted = ts->AngularCoords();
  }
  // make the turn with probabilistic veering around the predicted direction
  turn(ts,predicted,initlen);
}

cell_attraction_direction_model::cell_attraction_direction_model(direction_model_base * dmbcontrib, double & dmbweight, cell_attraction_direction_model & schema): direction_model_base(dmbcontrib,dmbweight,schema) {
}

String cell_attraction_direction_model::report_parameters_specific() {
  String res(" cell attraction");
  return res;
}

direction_model_base * cell_attraction_direction_model::clone() {
  return new cell_attraction_direction_model(contributing,contributingweight,*this);
}

void cell_attraction_direction_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_direction_model(clone());
  ts1.set_direction_model(this);
}

void cell_attraction_direction_model::handle_turning(terminal_segment & ts) {
  ts.set_direction_model(this);
}

void cell_attraction_direction_model::predict_direction(spatial & predicted, neuron * n, spatial & growthcone) {
#ifdef TESTING_SIMPLE_ATTRACTION
  PLL_LOOP_FORWARD(neuron,n->Root()->head(),1) if (e!=n) if (n->attractedto==e->attracts) {
    spatial d(e->Pos());
    d -= growthcone;
    double SQdistance = d.len2();
    d /= SQdistance; // square distance gravitational analogy
    d /= sqrt(SQdistance); // working with weighted unit vectors
    predicted += d;
  }
#endif
}

spatial cell_attraction_direction_model::predict(double weight, terminal_segment * ts, neuron * e) {
  // Called in a chain of direction models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  spatial predicted;
  predict_direction(predicted,e,ts->TerminalSegment()->P1);
  predicted.normalize();
  predicted *= weight;
  if (contributing) predicted += contributing->predict(contributingweight,ts,e);
  return predicted;
}

void cell_attraction_direction_model::direction(terminal_segment * ts, neuron * e, double initlen) {
  // ts points to the terminal segment before the turn
  // e points to the neuron to which the arbor belongs
  // The following are the preparations of this direction model:
  // 1. update the spatial parameters of the terminal fibre segment
  ts->update_segment_vector(); // hence, use ts->turn_free()
  // The following are the functions of this direction model:
  // 1. predict the direction based on the sum of (chemical) attractions
  //    in a gravitational analogy
  spatial predicted;
  predict_direction(predicted,e,ts->TerminalSegment()->P1);
  if (contributing) {
    predicted.normalize(); // only necessary here
    predicted += contributing->predict(contributingweight,ts,e);
  }
  predicted.convert_to_spherical();
  if (predicted.X()==0.0) {
    warning("Warning: NO DIRECTION PREDICTED in cell_attraction_direction_model::direction()\n"); // [***INCOMPLETE] See TL.
    predicted = ts->AngularCoords();
  }
  // 3. make the turn with probabilistic veering around the predicted direction
  turn(ts,predicted,initlen);
}

String direction_model_idstr[NUM_dm] = {
  "radial",
  "segment_history_tension",
  "cell_attraction",
  "vector"
};

direction_model_base * direction_model_selection(String & label, Command_Line_Parameters & clp, direction_model_base * superior_set, direction_model_base_ptr & dmbptr) {
  // Note: Providing the result pointer in dmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a direction model is not defined explicitly inherit the direction model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  direction_model dm = NUM_dm;
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"direction_model"))>=0) {
    String dmstr(downcase(clp.ParValue(n)));
    for (direction_model i = direction_model(0); i < NUM_dm; i = direction_model(i+1)) if (dmstr==direction_model_idstr[i]) { dm = i; break; }
    if (dm==NUM_dm) {
      if (dmstr==String("legacy")) dm = radial_dm;
      else return NULL; // catch branch_angle_model syntax errors
    }
    if (label.empty()) default_direction_model = dm; // modify default if setting universal
    if (dmstr==String("segment_history_tension")) dm = segment_history_tension_dm;
  } else if (label.empty()) dm = default_direction_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"direction_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_direction_model_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"dm_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_direction_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (dmbptr) dmbptr->delete_shared();
  // Set and return specified model schema or universal model schema if at universal set level
  axon_direction_history * adh = NULL;
  switch (dm) {
  case segment_history_tension_dm: 
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
    adh = new pll_axon_direction_history();
#else
    adh = new parse_axon_direction_history(*this);
#endif
    dmbptr = new tension_direction_model(*adh,label,nextlabel,clp);
    break;;
  case cell_attraction_dm:
    dmbptr = new cell_attraction_direction_model(label,nextlabel,clp);
    break;;
  case radial_dm:
    dmbptr = new radial_direction_model(label,nextlabel,clp);
    break;;
  case vector_dm:
    dmbptr = new vector_direction_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    dmbptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through dmbptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (direction_model_region_subset_schemas[0][universal_sps]!=dmbptr)) {
    if (direction_model_region_subset_schemas[0][universal_sps]) direction_model_region_subset_schemas[0][universal_sps]->delete_shared();
    direction_model_region_subset_schemas[0][universal_sps] = dmbptr->clone();
  }
  return dmbptr;
}
