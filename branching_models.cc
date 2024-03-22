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
// branching_models.cc
// Randal A. Koene, 20070509

/* Documentation:
   This implements base and derived classes for models of branching during
   network development that are based on bifurcation statistics (rather than
   causal hypotheses about environmental interactions).
 */

#include <math.h>
#include "branching_models.hh"
#include "neuron.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"

/* [***INCOMPLETE] Here I can include functions that parse the command line for
   general parameters to be used in specific branching models if no neuron or
   terminal segment specific parameters are given.
 */

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
branching_model default_branching_model = van_Pelt_bm; // identifier, changed when a different universal model is set
String universal_branching_model_root; // Used only to report if chaining at the universal set level.

TSBM default_TSBM = van_Pelt_tsbm; // identifier, changed when a different universal model is set
String universal_TSBM_root; // Used only to report if chaining at the universal set level.

branch_angle_model default_branch_angle_model = balanced_forces_bam; // identifier, changed when a different universal model is set
String universal_branch_angle_model_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_branching_model_base_ptr * branching_model_region_subset_schemas = NULL;

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_TSBM_base_ptr * TSBM_region_subset_schemas = NULL;
 
// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_branch_angle_model_base_ptr * branch_angle_model_region_subset_schemas = NULL;

branching_model_base::branching_model_base(branching_model_base * bmcontrib, double & bmweight, branching_model_base & schema):
  contributing(NULL), contributingweight(1.0), base_parameters(schema.base_parameters), cache1(0.0) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  if (bmcontrib) { // Clone chained models
    contributingweight = bmweight;
    contributing = bmcontrib->clone();
  }
}

branching_model_base::branching_model_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  contributing(NULL), contributingweight(1.0), base_parameters(NULL,0.0), cache1(0.0) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"min_node_interval"))>=0) {
    base_parameters.min_node_interval = atof(clp.ParValue(n));
    if (base_parameters.min_node_interval<0.0) {
      warning("Warning: "+thislabel+"min_node_interval="+clp.ParValue(n)+" interpreted as 0.0\n");
      base_parameters.min_node_interval = 0.0;
    }
  }
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"bm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!branching_model_selection(label,clp,contributing,contributing)) error("Error in branching_model_base(): No contributing branching model found for "+label+".\n");
}

double branching_model_base::predict(double weight, double & probability) {
  // Called in a chain of branching models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  // The actual weighting of contributions is done in the function predict_branching().
  // The sum of all weights is returned, so that division by that sum can produce
  // a properly scaled combined branching probability.
  probability += predict_branching(weight);
  if (contributing) weight += contributing->predict(contributingweight,probability);
  return weight;
}

double branching_model_base::branching_probability() {
  // PROTOCOL: This is the function that is called during iterative growth.
  //double predicted1=0.0, predicted2=0.0;
  double probability = predict_branching(1.0);
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,probability);
    probability /= summedcontributingweights;
  }
  return probability;
}

String branching_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting to contributing models.
  String res(" min_node_interval="+String(base_parameters.min_node_interval,"%.3f")+report_parameters_specific());
  /*res += String(base_parameters.veeranglemin,"\n  [veeranglemin=%.2f] ");
    res += String(base_parameters.veeranglemax,"\n  [veeranglemax=%.2f] "); */
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

van_Pelt_branching_model::van_Pelt_branching_model(branching_model_base * bmbcontrib, double & bmbweight, van_Pelt_branching_model & schema): branching_model_base(bmbcontrib,bmbweight,schema), probability(0.0), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
  creation_t = eq->T();
  last_t = creation_t;
}

van_Pelt_branching_model::van_Pelt_branching_model(String & thislabel, String & label, Command_Line_Parameters & clp): branching_model_base(thislabel,label,clp), probability(0.0), parameters(4.75,319680,0.5,gcc_whole_neuron) {
  // Chaining of models is taken care of in the base constructor.
  creation_t = eq->T();
  last_t = creation_t;
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"B_inf"))>=0) {
    double b = atof(clp.ParValue(n));
    if (b<0.0) {
      warning("Warning: Negative "+thislabel+"B_inf="+clp.ParValue(n)+" ignored\n");
    } else parameters.B_inf = b;
  }
  if ((n=clp.Specifies_Parameter(thislabel+"tau"))>=0) {
    double _tau = atof(clp.ParValue(n));
    if (_tau==0.0) {
      warning("Warning: Zero "+thislabel+"tau="+clp.ParValue(n)+" ignored\n");
    } else parameters.tau = _tau;
  }
  if ((n=clp.Specifies_Parameter(thislabel+"E"))>=0) parameters.E = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"E_competes_with"))>=0) {
    if (downcase(clp.ParValue(n))==String("whole_neuron")) parameters.competition_with_all_trees = gcc_whole_neuron;
    else if (downcase(clp.ParValue(n))==String("all_dendrites")) parameters.competition_with_all_trees = gcc_within_all_AorD;
    else if (downcase(clp.ParValue(n))==String("all_axons")) parameters.competition_with_all_trees = gcc_within_all_AorD;
    else if (downcase(clp.ParValue(n))==String("same_arbor")) parameters.competition_with_all_trees = gcc_within_tree;
    else warning("Unknown E competition value:"+clp.ParValue(n)+'\n');
  }
}

branching_model_base * van_Pelt_branching_model::clone(fibre_structure * _arbor) {
  branching_model_base * bmb = new van_Pelt_branching_model(contributing,contributingweight,*this);
  if (_arbor!=NULL) bmb->set_arbor(_arbor);
  return bmb;
}

//void van_Pelt_branching_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
// [***NOTE] Could set a new temp.variable in here that depends on N, so that its calculation is only done when needed.
//}

//void van_Pelt_branching_model::handle_turning(terminal_segment & ts) {
//}

double Btestarray[30] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
// 0 = t
// 1 = probability (Bt)
// 2 = TSBM probability
// 3 = TSBM probability * Bt
// 6 = number of probability calculations
// 11 = weight
// 12 = B_inf
// 13 = -t
// 14 = tau
// 15 = dt
// 16 = E
// 21 = sum for Cn
// 22 = divisor for Cn
// 
// 24 = binSgamma during init
// 25 = binSgamma before prob. comp.

double van_Pelt_branching_model::predict_branching(double weight) {
  // Computes the arbor-specific part of the van Pelt model branching probability,
  // That is the part of Eq.6 in [KOE:NETMORPH] before the parts involving S and C.
  // When S is not zero, those latter parts can be included by using the appropriate
  // TSBM.
  // Note that this model only calculates a value once per time step and uses a cache
  // of the last time for which a calculation was done to detect this. Also, since
  // only one calculation is needed for the entire arbor, we use the difference between
  // the last and current time for dt, instead of caching a fixed dt time step size.
  // That makes this model more versatile, and usable with variable time steps.
  // Parameters are:
  //   B_inf, tau, E
  if (last_t >= eq->T()) return probability;
  double dt = eq->T() - last_t;
  last_t = eq->T();
  double t = last_t - creation_t;
#ifdef USE_BTESTARRAY
  Btestarray[0] = t;
  Btestarray[11] = weight;
  Btestarray[12] = parameters.B_inf;
  Btestarray[13] = -t;
  Btestarray[14] = parameters.tau;
  Btestarray[15] = dt;
  Btestarray[16] = parameters.E;
#endif
  probability = weight * parameters.B_inf*exp(-t/parameters.tau)*(exp(dt/parameters.tau)-1.0);
#ifdef USE_BTESTARRAY
  Btestarray[10] = probability;
#endif
  switch (parameters.competition_with_all_trees) {
  case gcc_within_all_AorD:
    probability *= exp(-parameters.E*log((double) base_parameters.arbor->count_terminal_segments_in_PLL()));;
#ifdef USE_BTESTARRAY
    Btestarray[1] = probability;
#endif
    return probability;
    break;;
  case gcc_whole_neuron:
    probability *= exp(-parameters.E*log((double) (base_parameters.arbor->N()->total_input_terminal_segments()+base_parameters.arbor->N()->total_output_terminal_segments())));
#ifdef USE_BTESTARRAY
    Btestarray[1] = probability;
#endif
    return probability;
    break;;
  default: // gcc_within_tree
    probability *= exp(-parameters.E*log((double) base_parameters.arbor->count_terminal_segments()));
  }
#ifdef USE_BTESTARRAY
  Btestarray[6] += 1.0;
  Btestarray[1] = probability;
#endif
  return probability;
}

String van_Pelt_branching_model::report_parameters_specific() {
  String res(" van_Pelt (B_inf="+String(parameters.B_inf,"%.2f, tau=")+String(parameters.tau,"%.2f, E=")+String(parameters.E,"%.2f)"));
  return res;
}

TSBM_base::TSBM_base(TSBM_base * tsbmcontrib, double & tsbmweight, TSBM_base & schema):
  cache_ts(NULL), contributing(NULL), contributingweight(1.0) { //, base_parameters(schema.base_parameters) {
  if (tsbmcontrib) { // Clone chained models
    contributingweight = tsbmweight;
    contributing = tsbmcontrib->clone();
  }
}

TSBM_base::TSBM_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  cache_ts(NULL), contributing(NULL), contributingweight(1.0) { //, base_parameters(0.0) {
  int n;
  // local parameters here
  // ...
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"tsbm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!TSBM_selection(label,clp,contributing,contributing)) error("Error in TSBM_base(): No contributing TSBM found for "+label+".\n");
}

void TSBM_base::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_TSBM_model(clone());
  ts1.set_TSBM_model(this);
}

void TSBM_base::handle_turning(terminal_segment & ts) {
  ts.set_TSBM_model(this);
}

double TSBM_base::predict(double weight, double & probability) {
  // Called in a chain of TSBMs to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  // The actual weighting of contributions is done in the function predict_branching().
  // The sum of all weights is returned, so that division by that sum can produce
  // a properly scaled combined branching probability.
  probability += weight*predict_branching(); //weight);
  if (contributing) weight += contributing->predict(contributingweight,probability);
  return weight;
}

double TSBM_base::branching_probability(terminal_segment & ts) {
  // PROTOCOL: This is the function that is called during iterative growth.
  cache_ts = &ts;
  // 1. Local multiplicative adjustment (e.g. for S and C)
  double probability = predict_branching();
#ifdef USE_BTESTARRAY
  Btestarray[23] = probability;
#endif
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,probability);
    probability /= summedcontributingweights;
  }
  // 2. Multiply with the arbor probability
#ifdef USE_BTESTARRAY
  Btestarray[2] = probability;
#endif
  probability *= cache_ts->Arbor()->BranchingModel()->branching_probability();
#ifdef USE_BTESTARRAY
  Btestarray[3] = probability;
#endif
  return probability;
}

String TSBM_base::report_parameters() {
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

van_Pelt_TSBM::van_Pelt_TSBM(TSBM_base * tsbmbcontrib, double & tsbmbweight, van_Pelt_TSBM & schema): TSBM_base(tsbmbcontrib,tsbmbweight,schema), binSgamma(1.0), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

van_Pelt_TSBM::van_Pelt_TSBM(String & thislabel, String & label, Command_Line_Parameters & clp): TSBM_base(thislabel,label,clp), binSgamma(1.0), parameters(0.0) {
  // Chaining of models is taken care of in the base constructor.
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"S"))>=0) parameters.S = atof(clp.ParValue(n));
}

//***TEMP:
#define BINSGAMMA(_ts) exp(-(parameters.S*_ts->CentrifugalOrder())*LOG2)

#define LOG2 0.693147
void van_Pelt_TSBM::init_cache(terminal_segment * _ts) {
  binSgamma = exp(-(parameters.S*_ts->CentrifugalOrder())*LOG2);
#ifdef USE_BTESTARRAY
  Btestarray[24] = binSgamma;
#endif
  _ts->Arbor()->BranchingModel()->cache1 += binSgamma;
}

TSBM_base * van_Pelt_TSBM::clone(terminal_segment * _ts) {
  van_Pelt_TSBM * tsbmb = new van_Pelt_TSBM(contributing,contributingweight,*this);
  if (_ts!=NULL) tsbmb->init_cache(_ts);
  return tsbmb;
}

void van_Pelt_TSBM::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_TSBM_model(clone());
  ts1.set_TSBM_model(this);
  ts1.Arbor()->BranchingModel()->cache1 -= binSgamma; // update the sum for Cn
  init_cache(&ts1);
}

double van_Pelt_TSBM::predict_branching() { //double weight) {
  // Computes the S and C parts of Eq.6 in [KOE:NETMORPH] for the van Pelt model.
  // Note that binSgamma is computed and cached in the init_cache() function above.
#ifdef USE_BTESTARRAY
  Btestarray[21] = (double) cache_ts->Arbor()->count_terminal_segments();
  Btestarray[22] = cache_ts->Arbor()->BranchingModel()->cache1;
#endif
  //***double Cn = ((double) cache_ts->Arbor()->count_terminal_segments()) / cache_ts->Arbor()->BranchingModel()->cache1;
  double Cnsum = 0.0;
  PLL_LOOP_FORWARD(terminal_segment,cache_ts->Arbor()->TerminalSegments()->head(),1) Cnsum += BINSGAMMA(e);
  double Cn = ((double) cache_ts->Arbor()->count_terminal_segments()) / Cnsum;
  // [***NOTE] There may be cause to consider Cn over the whole neuron or over all axons or all dendrites (with or without an apical dendrite),
  // in which case I can add a switch here, as in the arbor model above. For now, we will assume that the order dependence plays a role only
  // within the same arbor.
#ifdef USE_BTESTARRAY
  Btestarray[25] = binSgamma;
#endif
  return BINSGAMMA(cache_ts)*Cn;
  //***return binSgamma*Cn;
}

String van_Pelt_TSBM::report_parameters_specific() {
  String res(" van_Pelt (S="+String(parameters.S,"%.2f)"));
  return res;
}

van_Pelt_specBM_TSBM::van_Pelt_specBM_TSBM(TSBM_base * tsbmbcontrib, double & tsbmbweight, van_Pelt_specBM_TSBM & schema): van_Pelt_TSBM(tsbmbcontrib,tsbmbweight,schema), specparameters(schema.specparameters) {
  // Chaining of models is taken care of in the base constructor.
}

van_Pelt_specBM_TSBM::van_Pelt_specBM_TSBM(String & thislabel, String & label, Command_Line_Parameters & clp, terminal_segment * ts): van_Pelt_TSBM(thislabel,label,clp), specparameters(NULL) {
  // Chaining of models is taken care of in the base constructor.
  // A specific branching model is identified here, which will then be shared with daughter branches.
  // [***NOTE] Should we expand thislabel with a tag such as "TSBM." to insure differentiation from other branching_model specifications?
  branching_model_selection(thislabel,clp,NULL,specparameters.bm);
  if (specparameters.bm) if (ts) specparameters.bm->set_arbor(ts->Arbor()); // The arbor has to be set either here or during cloning!
  //int n;
  //if ((n=clp.Specifies_Parameter(thislabel+"S"))>=0) parameters.S = atof(clp.ParValue(n));
}

TSBM_base * van_Pelt_specBM_TSBM::clone(terminal_segment * _ts) {
  //cout << "_CLONING_TSBM_"; cout.flush();
  van_Pelt_specBM_TSBM * tsbmb = new van_Pelt_specBM_TSBM(contributing,contributingweight,*this);
  if (_ts!=NULL) {
    tsbmb->init_cache(_ts);
    if (tsbmb->specparameters.bm) if (!(tsbmb->specparameters.bm->Arbor())) tsbmb->specparameters.bm->set_arbor(_ts->Arbor());
  }
  return tsbmb;
}

/* double van_Pelt_specBM_TSBM::predict_branching(double weight) {
  // Computes the S and C parts of Eq.6 in [KOE:NETMORPH] for the van Pelt model.
  // Uses a specific branching model, see TL#200807210612.5.
  // Note that binSgamma is computed and cached in the van_Pelt_TSBM::init_cache() function above.
  // In this function, the Arbor() is referenced to compute S and C.
  //double Cn = ((double) cache_ts->Arbor()->count_terminal_segments()) / cache_ts->Arbor()->BranchingModel()->cache1;
  double Cnsum = 0.0;
  PLL_LOOP_FORWARD(terminal_segment,cache_ts->Arbor()->TerminalSegments()->head(),1) Cnsum += BINSGAMMA(e);
  double Cn = ((double) cache_ts->Arbor()->count_terminal_segments()) / Cnsum;
  // [***NOTE] There may be cause to consider Cn over the whole neuron or over all axons or all dendrites (with or without an apical dendrite),
  // in which case I can add a switch here, as in the arbor model above. For now, we will assume that the order dependence plays a role only
  // within the same arbor.
  return BINSGAMMA(cache_ts)*Cn;
  //return binSgamma*Cn;
  } */

double van_Pelt_specBM_TSBM::branching_probability(terminal_segment & ts) {
  // PROTOCOL: This is the function that is called during iterative growth.
  // Uses a specific branching model, see TL#200807210612.5.
  cache_ts = &ts;
  // 1. Local multiplicative adjustment (e.g. for S and C)
  double probability = predict_branching();
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,probability);
    probability /= summedcontributingweights;
  }
  // 2. Multiply with the arbor probability
  if (!specparameters.bm) probability *= cache_ts->Arbor()->BranchingModel()->branching_probability(); // by default exactly as in van_Pelt_TSBM
  else probability *= specparameters.bm->branching_probability();
  return probability;
}

String van_Pelt_specBM_TSBM::report_parameters_specific() {
  String res(" specific BM "+van_Pelt_TSBM::report_parameters_specific());
  return res;
}

/* Use this function information to create polynomial O1 and O3 branching TSBMs:

Polynomial_O1_dendritic_growth_model::Polynomial_O1_dendritic_growth_model(double E, double c, double tf): van_Pelt_dendritic_growth_plus_turns_model(E,1.0,1.0,E,1.0,1.0,E,1.0,1.0,tm_each_dt,tf,0.05) {
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

Polynomial_O3_dendritic_growth_model::Polynomial_O3_dendritic_growth_model(double E, double c, double c1, double c2, double tf): Polynomial_O1_dendritic_growth_model(E,c,tf), _c1(c1), _c2(c2) {
  _L0_min = 0.0; // set by L(0) cubic approximation of dendrite length for mirrored data
  _L0_max = 0.0;
}

Polynomial_O3_dendritic_growth_model::Polynomial_O3_dendritic_growth_model(Command_Line_Parameters & clp): Polynomial_O1_dendritic_growth_model(clp), _c1(0.0/DAYSECONDSSQR), _c2(0.000039/DAYSECONDSCUB) {
  _c = 0.145511/DAYSECONDS;
  local_parse_CLP(clp);
}

double Polynomial_O1_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  //double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
  return _c*(t2-t1);//vnt_E; // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c =  0.145156 / (24.0*60.0*60.0))
}

double Polynomial_O3_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1, double t2) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  double t1sq = t1*t1, t2sq = t2*t2;
  return (_c*(t2-t1)) + (_c1*(t2sq-t1sq)) + (_c2*((t2sq*t2)-(t1sq*t1))); // set the c parameters by dividing daily fits to Ger Ramaker's data by (24.0*60.0*60.0)
}

void Polynomial_O1_dendritic_growth_model::set_fixed_time_step_size(double dt) {
  _dt=dt;
  _c_dt = _c*dt;
}

double Polynomial_O1_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1) {
  // [*** NOTE] The O3 polynomial fit to bifurcation in Ger Ramakers' data
  // does not use the t1 parameter.
  //double nt_E = exp(-_E*log((double) n->total_input_terminal_segments()));
  return _c_dt;// nt_E; // An O1 polynomial fit to bifurcation in Ger Ramakers' data (_c =  0.145156 / (24.0*60.0*60.0))
}

double Polynomial_O3_dendritic_growth_model::expected_number_of_branches(neuron * n, double t1) {
  // The result will be multiplied by nt_E to obtain a probability per
  // terminal segment.
  // This function does not use the parameter n.
  return expected_number_of_branches(n,t1,t1+_dt);
}

*/

branch_angle_model_base::branch_angle_model_base(branch_angle_model_base * bmcontrib, double & bmweight, branch_angle_model_base & schema):
  contributing(NULL), contributingweight(1.0) { //, base_parameters(schema.base_parameters) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  if (bmcontrib) { // Clone chained models
    contributingweight = bmweight;
    contributing = bmcontrib->clone();
  }
  pdf = schema.pdf->clone();
}

branch_angle_model_base::branch_angle_model_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  contributing(NULL), contributingweight(1.0), pdf(NULL) { //, base_parameters(turnanglemin,turnanglemax) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  int n;
  pdf = pdfselection(thislabel+"bam.",clp,X_branch);
  if (!pdf) pdf = new normal_pdf(X_branch,0.0,0.3,1.0); // a normal PDF is the default for BAMs
  else if (pdf->amplitude()>1.0) warning("Warning: Branch Angle Model has perturbation with amplitude greater than 1.0, unpredictable angles!\n");
  /* double vmin = base_parameters.veeranglemin;
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
      warning("Warning: "+thislabel+"veeranglemax="+clp.ParValue(n)+" not prmitted, value must be >0.0, not modified\n");
      vmax = base_parameters.veeranglemax;
    }
  }
  if (vmax<=vmin) warning("Warning: "+thislabel+"veeranglemax must be > "+thislabel+"veeranglemin, not modified\n");
  else {
    base_parameters.veeranglemin = vmin;
    base_parameters.veeranglemax = vmax;
    } */
  if (label.empty()) return;
  if ((n=clp.Specifies_Parameter(label+"bm_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  if (!branch_angle_model_selection(label,clp,contributing,contributing)) error("Error in branch_angle_model_base(): No contributing branch angle model found for "+label+".\n");
}

String branch_angle_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting to contributing models.
  String res(report_parameters_specific());
  res += " p_BAM=" + pdf->report_parameters();
  /*res += String(base_parameters.veeranglemin,"\n  [veeranglemin=%.2f] ");
    res += String(base_parameters.veeranglemax,"\n  [veeranglemax=%.2f] "); */
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

void branch_angle_model_base::bifurcate(terminal_segment * ts1, double branchangle1, terminal_segment * ts2, double branchangle2) {
  // This function receives expected branch angles, computed with branch_angle() functions.
  // Those branch angles may be modified to proportionally reduce or increase the two branch angles.
  // Such modification can be done in the form of a probabilistic perturbation or by scaling to the
  // sum of initial lengths (or elongation rates). Presently, perturbation takes the form of
  // normal distributed perturbation, as described in the gfneuroinformatics paper.
  // This function (and the Branch Angle Models as a whole) are only concerned with the PITCH
  // angle relative to the parent fiber direction. As such, they do not consider the ROLL in
  // 3D, or the transformation to absolute angular vectors for ts1 and ts2.
  // For more about the ensuing computations, please see the terminal_segment::branch() function.
  branchangle1 *= (1.0+pdf->random_selection());
  branchangle2 *= (1.0+pdf->random_selection());
  ts1->AngularCoords().set_Y(branchangle1);
  ts2->AngularCoords().set_Y(branchangle2);
}

balanced_forces_branch_angle_model::balanced_forces_branch_angle_model(branch_angle_model_base * bambcontrib, double & bambweight, balanced_forces_branch_angle_model & schema): branch_angle_model_base(bambcontrib,bambweight,schema) {
  parallelogramanglepdf = schema.parallelogramanglepdf->clone();
}

balanced_forces_branch_angle_model::balanced_forces_branch_angle_model(String & thislabel, String & label, Command_Line_Parameters & clp): branch_angle_model_base(thislabel,label,clp), parallelogramanglepdf(NULL) {
  parallelogramanglepdf = pdfselection(thislabel+"bam.bfbam.",clp,X_branch);
  if (!parallelogramanglepdf) parallelogramanglepdf = new normal_pdf(X_branch,0.5*M_PI,0.5,M_PI-0.1);
  else if (parallelogramanglepdf->amplitude()>=M_PI) warning("Warning: The parallelogram angle PDF of the Balanced Forces Branch Angle Model has an amplitude of 180 degrees or more!\n");
}

String balanced_forces_branch_angle_model::report_parameters_specific() {
  String res(" balanced forces (parallelogram angle PDF: ");
  res += parallelogramanglepdf->report_parameters() + ')';
  return res;
}

branch_angle_model_base * balanced_forces_branch_angle_model::clone() {
  return new balanced_forces_branch_angle_model(contributing,contributingweight,*this);
}

void balanced_forces_branch_angle_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_branch_angle_model(clone());
  ts1.set_branch_angle_model(this);
}

void balanced_forces_branch_angle_model::handle_turning(terminal_segment & ts) {
  ts.set_branch_angle_model(this);
}

void balanced_forces_branch_angle_model::predict_branch_angles(double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2, double weight) {
  // The results are added to the values in predicted1 and predicted2 weighted by weight.
  // Here forces are derived from the initial elongation rates of the new terminal segments.
  // Those forces create a parallelogram with a combined force vector along the parent
  // direction. The resulting angles of the force vectors with the projected parent direction
  // are the predicted branch angles.
  // In addition to using the ratio of the forces to determine relative angles, the angles
  // can be narrowed or widened in proportion, by lengthening or shortening the vector along
  // the parent direction, with which the parallelogram is made. An equivalent alternative
  // is to use a right angled parallelogram (which is easy) to find the relative angles
  // first, and to then proportionally reduce or enlarge both angles.
  // The reduction or enlargement could be (a) probabilistic or (b) determined by the
  // absolute elongation rate of the sum of the two elongation rates, since we can imagine
  // that more rapid forward progression of a growth cone could lead to greater continuing
  // forward motion along both branches.
  // The probabilistic option fits in well with the "veering" approach usually done in the
  // base function, such as the bifurcation() function.
  // Either (a) or (b), or both (a) and (b) could be applied.
  // Note that using the distributed elongations or the random values here does not matter,
  // since the total distributed cancels out.
  double parallelogramangle = parallelogramanglepdf->random_positive();
  double force1 = ts1->Length(), force2 = ts2->Length();
  if ((force1==0.0) && (force2==0.0)) { // The distributed elongation was zero, 
    force1 = X_branch.get_rand_real1();
    force2 = X_branch.get_rand_real1();
  }
  if (force1==0.0) {
    /*if (force2==0.0) {
      warning("Warning: Both forces in balanced_forces_branch_angle_model::predict_branch_angles() are zero, branches do not diverge!\n");
      return;
      }*/
    // angle1 is parallelogramangle, angle2 is zero
    predicted1 += (weight*parallelogramangle);
    return;
  } else if (force2==0.0) {
    // angle1 is zero, angle2 is parallelogramangle
    predicted2 += (weight*parallelogramangle);
    return;
  }
  double angle1 = atan(force2/force1);
  double angle2 = parallelogramangle - angle1;
  predicted1 += weight*angle1; predicted2 += weight*angle2;
}

double balanced_forces_branch_angle_model::predict(double weight, double & predicted1, double & predicted2, terminal_segment * ts1, terminal_segment * ts2) {
  // Called in a chain of branch angle models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  // The actual weighting of contributions is done in the function predict_branch_angles().
  // The sum of all weights is returned, so that division by that sum can produce
  // properly scaled combined branch angles.
  predict_branch_angles(predicted1,predicted2,ts1,ts2,weight);
  if (contributing) weight += contributing->predict(contributingweight,predicted1,predicted2,ts1,ts2);
  return weight;
}

void balanced_forces_branch_angle_model::branch_angle(terminal_segment * ts1, terminal_segment * ts2) { // fibre_segment * parent
  // PROTOCOL: This is the function that is called during iterative growth.
  // The following are the preparations of this branch angle model:
  // [...]
  // The following are the functions of this branch angle model:
  // 1. predict the branch angles based on a balance of forces analogy
  double predicted1=0.0, predicted2=0.0;
  predict_branch_angles(predicted1,predicted2,ts1,ts2);
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,predicted1,predicted2,ts1,ts2);
    predicted1 /= summedcontributingweights; predicted2 /= summedcontributingweights;
  }
  // [***BEWARE] Is it possible for both angles to be zero? 
  // 2. bifurcate with a roll around the parent direction, with possibly perturbed predicted branch angles, and spherical coordinates for the new terminal segments
  bifurcate(ts1,predicted1,ts2,predicted2);
}

String branching_model_idstr[NUM_bm] = {
  "van_pelt",
  "polynomial_o1",
  "polynomial_o3"
};

branching_model_base * branching_model_selection(String & label, Command_Line_Parameters & clp, branching_model_base * superior_set, branching_model_base_ptr & bmbptr) {
  // Note: Providing the result pointer in bmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a branching model is not defined explicitly inherit the branching model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  branching_model bm = NUM_bm; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"branching_model"))>=0) {
    String bmstr(downcase(clp.ParValue(n)));
    for (branching_model i = branching_model(0); i < NUM_bm; i = branching_model(i+1)) if (bmstr==branching_model_idstr[i]) { bm = i; break; }
    if (bm==NUM_bm) return NULL; // catch branching_model syntax errors
    if (label.empty()) default_branching_model = bm; // modify default if setting universal
  } else if (label.empty()) bm = default_branching_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"branching_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_branching_model_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"bm_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_branching_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (bmbptr) bmbptr->delete_shared();
  // If specified, set and return specified model schema or universal model schema if at universal set level
  switch (bm) {
  case van_Pelt_bm:
    bmbptr = new van_Pelt_branching_model(label,nextlabel,clp);
    break;;
  case polynomial_O1_bm:
    //[***INCOMPLETE]    bmbptr = new polynomial_O1_branching_model(label,nextlabel,clp);
    break;;
  case polynomial_O3_bm:
    //[***INCOMPLETE]    bmbptr = new polynomial_O3_branching_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    //cout << "Branching model... defaulting\n";
    if (!superior_set) return NULL;
    bmbptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through bmbptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (branching_model_region_subset_schemas[0][universal_sps]!=bmbptr)) {
    if (branching_model_region_subset_schemas[0][universal_sps]) branching_model_region_subset_schemas[0][universal_sps]->delete_shared();
    branching_model_region_subset_schemas[0][universal_sps] = bmbptr->clone();
  }
  return bmbptr;
}

String TSBM_idstr[NUM_tsbm] = {
  "van_pelt",
  "van_pelt_specbm"
};

TSBM_base * TSBM_selection(String & label, Command_Line_Parameters & clp, TSBM_base * superior_set, TSBM_base_ptr & tsbmbptr, terminal_segment * ts) {
  // Note: Providing the result pointer in tsbmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a branching model is not defined explicitly inherit the branching model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  TSBM tsbm = NUM_tsbm; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"TSBM"))>=0) {
    String tsbmstr(downcase(clp.ParValue(n)));
    for (TSBM i = TSBM(0); i < NUM_tsbm; i = TSBM(i+1)) if (tsbmstr==TSBM_idstr[i]) { tsbm = i; break; }
    if (tsbm==NUM_tsbm) return NULL; // catch TSBM syntax errors
    if (label.empty()) default_TSBM = tsbm; // modify default if setting universal
  } else if (label.empty()) tsbm = default_TSBM; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"TSBM_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_TSBM_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"tsbm_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_TSBM_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (tsbmbptr) tsbmbptr->delete_shared();
  // If specified, set and return specified model schema or universal model schema if at universal set level
  switch (tsbm) {
  case van_Pelt_tsbm:
    tsbmbptr = new van_Pelt_TSBM(label,nextlabel,clp);
    break;;
  case van_Pelt_specBM_tsbm:
    tsbmbptr = new van_Pelt_specBM_TSBM(label,nextlabel,clp,ts);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    tsbmbptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through tsbmbptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (TSBM_region_subset_schemas[0][universal_sps]!=tsbmbptr)) {
    if (TSBM_region_subset_schemas[0][universal_sps]) TSBM_region_subset_schemas[0][universal_sps]->delete_shared();
    TSBM_region_subset_schemas[0][universal_sps] = tsbmbptr->clone();
  }
  return tsbmbptr;
}

String branch_angle_model_idstr[NUM_bam] = {
  "balanced_forces",
  "fork_main_daughter" // downcased!
};

branch_angle_model_base * branch_angle_model_selection(String & label, Command_Line_Parameters & clp, branch_angle_model_base * superior_set, branch_angle_model_base_ptr & bambptr) {
  // Note: Providing the result pointer in bambptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a branch angle model is not defined explicitly inherit the branch angle model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  branch_angle_model bam = NUM_bam; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"branch_angle_model"))>=0) {
    String bamstr(downcase(clp.ParValue(n)));
    for (branch_angle_model i = branch_angle_model(0); i < NUM_bam; i = branch_angle_model(i+1)) if (bamstr==branch_angle_model_idstr[i]) { bam = i; break; }
    if (bam==NUM_bam) return NULL; // catch branch_angle_model syntax errors
    if (label.empty()) default_branch_angle_model = bam; // modify default if setting universal
  } else if (label.empty()) bam = default_branch_angle_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"branch_angle_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_branch_angle_model_root = nextlabel; // to report at universal level
  } else if ((n=clp.Specifies_Parameter(label+"bam_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_branch_angle_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (bambptr) bambptr->delete_shared();
  // If specified, set and return specified model schema or universal model schema if at universal set level
  switch (bam) {
  case balanced_forces_bam: 
    bambptr = new balanced_forces_branch_angle_model(label,nextlabel,clp);
    break;;
  case fork_main_daughter_bam:
    //[***INCOMPLETE]    bambptr = new fork_main_daughter_branch_angle_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    bambptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through bambptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (branch_angle_model_region_subset_schemas[0][universal_sps]!=bambptr)) {
    if (branch_angle_model_region_subset_schemas[0][universal_sps]) branch_angle_model_region_subset_schemas[0][universal_sps]->delete_shared();
    branch_angle_model_region_subset_schemas[0][universal_sps] = bambptr->clone();
  }
  return bambptr;
}
