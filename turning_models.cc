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
// turning_models.cc
// Randal A. Koene, 20070509

/* Documentation:
   This implements base and derived classes for models of turning during
   network development that are based on bifurcation statistics (rather than
   causal hypotheses about environmental interactions).
 */

#include <math.h>
#include "turning_models.hh"
#include "neuron.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"

/* [***INCOMPLETE] Here I can include functions that parse the command line for
   general parameters to be used in specific turning models if no neuron or
   terminal segment specific parameters are given.
 */

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
turning_model default_turning_model = independent_rate_tm; // identifier, changed when a different universal model is set
String universal_turning_model_root; // Used only to report if chaining at the universal set level.

TSTM default_TSTM = tm_linear_rate; // identifier, changed when a different universal model is set
String universal_TSTM_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_TSTM_base_ptr * TSTM_region_subset_schemas = NULL;
 
turning_model_base::turning_model_base(turning_model_base * tmcontrib, double & tmweight, turning_model_base & schema):
  contributing(NULL), contributingweight(1.0) { //, base_parameters(schema.base_parameters) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  if (tmcontrib) { // Clone chained models
    contributingweight = tmweight;
    contributing = tmcontrib->clone();
  }
}

turning_model_base::turning_model_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  contributing(NULL), contributingweight(1.0) { //, base_parameters(turnanglemin,turnanglemax) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  // [***BEWARRE] OUTDATED - this is done in the TSTM models! The arbor level models should do something different!
  int n;
  /* double vmin = base_parameters.veeranglemin;
  double vmax = base_parameters.veeranglemax;
  if ((n=clp.Specifies_Parameter(thislabel+"veeranglemin"))>=0) {
    vmin = atof(clp.ParValue(n));
    if (vmin<0.0) {
      cout << "Warning: " << thislabel << "veeranglemin=" << clp.ParValue(n) << " not permitted, value must be >=0.0, not modified\n";
      vmin = base_parameters.veeranglemin;
    }
  }
  if ((n=clp.Specifies_Parameter(thislabel+"veeranglemax"))>=0) {
    vmax = atof(clp.ParValue(n));
    if (vmax<=0.0) {
      cout << "Warning: " << thislabel << "veeranglemax=" << clp.ParValue(n) << " not permitted, value must be >0.0, not modified\n";
      vmax = base_parameters.veeranglemax;
    }
  }
  if (vmax<=vmin) cout << "Warning: " << thislabel << "veeranglemax must be > " << thislabel << "veeranglemin, not modified\n";
  else {
    base_parameters.veeranglemin = vmin;
    base_parameters.veeranglemax = vmax;
    } */
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"tm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!turning_model_selection(label,clp,contributing)) error("Error in turning_model_base(): No contributing turning model found for "+label+".\n");
}

String turning_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting to contributing models.
  String res(report_parameters_specific());
  /*res += String(base_parameters.veeranglemin,"\n  [veeranglemin=%.2f] ");
    res += String(base_parameters.veeranglemax,"\n  [veeranglemax=%.2f] "); */
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

TSTM_base::TSTM_base(TSTM_base * tstmcontrib, double & tstmweight, TSTM_base & schema): cache_ts(NULL), cache_p_i(0.0),
  contributing(NULL), contributingweight(1.0) { //, base_parameters(schema.base_parameters)
  if (tstmcontrib) { // Clone chained models
    contributingweight = tstmweight;
    contributing = tstmcontrib->clone();
  }
}

TSTM_base::TSTM_base(String & thislabel, String & label, Command_Line_Parameters & clp): cache_ts(NULL), cache_p_i(0.0),
  contributing(NULL), contributingweight(1.0) { //, base_parameters(0.0)
  int n;
  // local parameters here
  // ...
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"tstm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!TSTM_selection(label,clp,contributing,contributing)) error("Error in TSTM_base(): No contributing TSTM found for "+label+".\n");
}

void TSTM_base::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_TSTM_model(clone());
  ts1.set_TSTM_model(this);
}

void TSTM_base::handle_turning(terminal_segment & ts) {
  ts.set_TSTM_model(this);
}

double TSTM_base::predict(double weight, double & probability) {
  // Called in a chain of TSTMs to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  // The actual weighting of contributions is done in the function predict_branching().
  // The sum of all weights is returned, so that division by that sum can produce
  // a properly scaled combined branching probability.
  probability += weight*predict_turn(); //weight);
  if (contributing) weight += contributing->predict(contributingweight,probability);
  return weight;
}

bool TSTM_base::turn_event(double r_selected, terminal_segment & ts, double p_i) {
  // PROTOCOL: This is the function that is called during iterative growth.
  cache_ts = &ts;
  cache_p_i = p_i;
  // 1. Combined probability
  double probability = predict_turn(); //1.0);
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,probability);
    probability /= summedcontributingweights;
  }
  // 2. Multiply with the arbor probability
  // [***INCOMPLETE] Arbor turning models are not yet in use.
  //probability *= cache_ts->Arbor()->BranchingModel()->branching_probability();
  return (r_selected<=probability);
}

String TSTM_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting to contributing models.
  String res(report_parameters_specific());
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

each_dt_TSTM::each_dt_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, each_dt_TSTM & schema): TSTM_base(tstmbcontrib,tstmbweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

each_dt_TSTM::each_dt_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp): TSTM_base(thislabel,label,clp), parameters(eq->T()) {
  // Chaining of models is taken care of in the base constructor.
  // Note: By initializing t_prev to eq->T() there is not immediately a turn at the moment
  // when a new TSTM object is constructed.
  //  int n;
  //  if ((n=clp.Specifies_Parameter(thislabel+"S"))>=0) parameters.S = atof(clp.ParValue(n));
}

TSTM_base * each_dt_TSTM::clone() { //terminal_segment * _ts) {
  return new each_dt_TSTM(contributing,contributingweight,*this);
}

double each_dt_TSTM::predict_turn() { //double weight) {
  if (eq->T()<=parameters.t_prev) return 0.0;
  parameters.t_prev = eq->T();
  return 1.0;
}

String each_dt_TSTM::report_parameters_specific() {
  String res(" each_dt");
  return res;
}

none_TSTM::none_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, none_TSTM & schema): TSTM_base(tstmbcontrib,tstmbweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

none_TSTM::none_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp): TSTM_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
  // Note: By initializing t_prev to eq->T() there is not immediately a turn at the moment
  // when a new TSTM object is constructed.
  //  int n;
  //  if ((n=clp.Specifies_Parameter(thislabel+"S"))>=0) parameters.S = atof(clp.ParValue(n));
}

TSTM_base * none_TSTM::clone() { //terminal_segment * _ts) {
  return new none_TSTM(contributing,contributingweight,*this);
}

double none_TSTM::predict_turn() { //double weight) {
  return -1.0; // any random value r between 0 and 1 can never be smaller than this, so that this predict that no turn events should ever occur.
}

String none_TSTM::report_parameters_specific() {
  String res(" none");
  return res;
}

linear_rate_TSTM::linear_rate_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, linear_rate_TSTM & schema): TSTM_base(tstmbcontrib,tstmbweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

linear_rate_TSTM::linear_rate_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp): TSTM_base(thislabel,label,clp), parameters(0.1) {
  // Chaining of models is taken care of in the base constructor.
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"turn_rate"))>=0) parameters.turn_rate = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"turn_separation"))>=0) parameters.turn_rate = 1.0/atof(clp.ParValue(n)); // see explanation in predict_turn()
}

TSTM_base * linear_rate_TSTM::clone() { //terminal_segment * _ts) {
  return new linear_rate_TSTM(contributing,contributingweight,*this);
}

double linear_rate_TSTM::predict_turn() { //double weight) {
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
  double P_turn = parameters.turn_rate*cache_ts->elongated_length();
#else
  double P_turn = parameters.turn_rate*_dt;
#endif
  if (P_turn>=1.0) warning("Warning: One or more turning points are expected (dendrite P_turn="+String(P_turn,"%.2f")+"), perhaps development time steps should be smaller\n");
return P_turn;
}

String linear_rate_TSTM::report_parameters_specific() {
  String res(" linear rate (turn_rate="+String(parameters.turn_rate,"%.3f"));
#ifdef ENABLE_FIXED_STEP_SIMULATION
  res += "turns/micrometer)";
#else
  res += "turns/second)";
#endif
  return res;
}

branch_coupled_TSTM::branch_coupled_TSTM(TSTM_base * tstmbcontrib, double & tstmbweight, branch_coupled_TSTM & schema): TSTM_base(tstmbcontrib,tstmbweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

branch_coupled_TSTM::branch_coupled_TSTM(String & thislabel, String & label, Command_Line_Parameters & clp): TSTM_base(thislabel,label,clp), parameters(9.25) {
  // Chaining of models is taken care of in the base constructor.
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tf"))>=0) parameters.turnfactor = atof(clp.ParValue(n));
}

TSTM_base * branch_coupled_TSTM::clone() { //terminal_segment * _ts) {
  return new branch_coupled_TSTM(contributing,contributingweight,*this);
}

double branch_coupled_TSTM::predict_turn() { //double weight) {
  return (1.0+parameters.turnfactor)*cache_p_i;
}

String branch_coupled_TSTM::report_parameters_specific() {
  String res(" branch_coupled (tf="+String(parameters.turnfactor,"%.3f)"));
  return res;
}

turning_model_base * turning_model_selection(String & label, Command_Line_Parameters & clp, turning_model_base_ptr & tmbptr) {
  // Note: Providing the result pointer in tmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // [***BEWARRE] OUTDATED - this is done in the TSTM models! The arbor level models should do something different!
  turning_model tm = default_turning_model;
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"turning_model"))>=0) {
    String tmstr(downcase(clp.ParValue(n)));
    if (tmstr==String("Independent_Rate")) tm = independent_rate_tm;
    else if (tmstr==String("Elongation_Locked")) tm = elongation_locked_tm;
    else if (tmstr==String("Bifurcation_Coupled")) tm = bifurcation_coupled_tm;
    else return NULL; // catch turning_model syntax errors
    if (label.empty()) default_turning_model = tm; // modify default if setting universal
  } else {
    if (!label.empty()) return NULL; // non-universal requires explicit turning model
  }
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"turning_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_turning_model_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"tm_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_turning_model_root = nextlabel; // to report at universal level
  }
  // Set and return specified model schema or universal model schema if at universal set level
  if (tmbptr) tmbptr->delete_shared(); // Clean up any previously defined schema
  switch (tm) {
  case independent_rate_tm: 
    //[***INCOMPLETE] tmbptr = new independent_rate_turning_model(label,nextlabel,clp);
    break;;
  case elongation_locked_tm:
    //[***INCOMPLETE] tmbptr = new elongation_locked_turning_model(label,nextlabel,clp);
    break;;
  case bifurcation_coupled_tm:
    //[***INCOMPLETE] tmbptr = new bifurcation_coupled_turning_model(label,nextlabel,clp);
    break;;
  default: break;; // to avoid a compilation warning message about the default
  }
  /*  if ((label.empty()) && (turning_model_subset_schemas[universal_sps]!=tmbptr)) {
    // This is in case a secondary schema reference is provided for the
    // universal set level for some reason, i.e. create a schema object for
    // turning_model_subset_schemas[universal_sps], as well as the recipient tmbptr.
    if (turning_model_subset_schemas[universal_sps]) turning_model_subset_schemas[universal_sps]->delete_shared();
    turning_model_subset_schemas[universal_sps] = tmbptr->clone();
    }*/
  return tmbptr;
}

String TSTM_idstr[tm_NUM] = {
  "none",
  "each_dt",
  "linear_rate",
  "branch_coupled"
};

TSTM_base * TSTM_selection(String & label, Command_Line_Parameters & clp, TSTM_base * superior_set, TSTM_base_ptr & tstmbptr) {
  // Note: Providing the result pointer in tstmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a turning model is not defined explicitly inherit the turning model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  TSTM tstm = tm_NUM; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"TSTM"))>=0) {
    String tstmstr(downcase(clp.ParValue(n)));
    for (TSTM i = TSTM(0); i < tm_NUM; i = TSTM(i+1)) if (tstmstr==TSTM_idstr[i]) { tstm = i; break; }
    if (tstm==tm_NUM) return NULL; // catch TSTM syntax errors
    if (label.empty()) default_TSTM = tstm; // modify default if setting universal
  } else if (label.empty()) tstm = default_TSTM; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"TSTM_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_TSTM_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"tstm_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_TSTM_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (tstmbptr) tstmbptr->delete_shared();
  // If specified, set and return specified model schema or universal model schema if at universal set level
  switch (tstm) {
  case tm_none:
    tstmbptr = new none_TSTM(label,nextlabel,clp);
    break;;
  case tm_each_dt:
    tstmbptr = new each_dt_TSTM(label,nextlabel,clp);
    break;;
  case tm_linear_rate:
    tstmbptr = new linear_rate_TSTM(label,nextlabel,clp);
    break;;
  case tm_branch_coupled:
    tstmbptr = new branch_coupled_TSTM(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    tstmbptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through tstmbptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (TSTM_region_subset_schemas[0][universal_sps]!=tstmbptr)) {
    if (TSTM_region_subset_schemas[0][universal_sps]) TSTM_region_subset_schemas[0][universal_sps]->delete_shared();
    TSTM_region_subset_schemas[0][universal_sps] = tstmbptr->clone();
  }
  return tstmbptr;
}

