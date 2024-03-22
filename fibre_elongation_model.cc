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
// fibre_elongation_model.cc
// Randal A. Koene, 20051221

#include "fibre_elongation_model.hh"
#include "fibre_structure.hh"
#include "dendritic_growth_model.hh"
#include "axon_direction_model.hh"
#include "neuron.hh"
#include "Txt_Object.hh"
#include "diagnostic.hh"

#define DAYSECONDS 86400.0
#define DAYSECONDSSQR (86400.0*86400.0)
#define DAYSECONDSCUB (86400.0*86400.0*86400.0)

#ifdef TESTING_ELONGATION_TOTAL
double usedarborelongationfraction = 0.0;
#endif
// Elongation model selection indicators at the general level
arbor_elongation_model general_arbor_elongation_model = van_Pelt_aem;
String general_arbor_elongation_model_root;

terminal_segment_elongation_model general_terminal_segment_elongation_model = simple_tsem;
String general_terminal_segment_elongation_model_root;

//arbor_elongation_model_base * general_arbor_elongation_model_schema = NULL;
//terminal_segment_elongation_model_base * general_terminal_segment_elongation_model_schema = NULL;

// Schema object pointer arrays
// [***NOTE] See TL#200603120618.2.

// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
elongation_rate_initialization_model default_elongation_rate_initialization_model = length_distribution_eri; // identifier, changed when a different universal model is set
String universal_elongation_rate_initialization_model_root; // Used only to report if chaining at the universal set level.

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_arbor_elongation_model_base_ptr * arbor_elongation_model_region_subset_schemas = NULL;

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_terminal_segment_elongation_model_base_ptr * terminal_segment_elongation_model_region_subset_schemas = NULL;

// Region and natural subset specific schemas (initialized in neuron.cc:general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net))
region_elongation_rate_initialization_model_base_ptr * elongation_rate_initialization_model_region_subset_schemas = NULL;

arbor_elongation_model_base::arbor_elongation_model_base(arbor_elongation_model_base * aemcontrib, double & aemweight, arbor_elongation_model_base & schema): contributing(NULL), contributingweight(1.0), prev_t(eq->T()), base_parameters(schema.base_parameters) {
  if (aemcontrib) { // Clone chained models
    contributingweight = aemweight;
    contributing = aemcontrib->clone();
  }
  pdf = schema.pdf->clone();
}

arbor_elongation_model_base::arbor_elongation_model_base(String & thislabel, String & label, Command_Line_Parameters & clp): contributing(NULL), contributingweight(1.0), prev_t(eq->T()), pdf(NULL), base_parameters(0.0) {
  // This constructor is intended for use by schemas, which can set up model
  // chains, as well as look for general (thislabel=="") or subset parameter
  // settings that differ from the compile-time defaults.
  // (See TL#200603120618.5.)
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  pdf = pdfselection(thislabel+"aem.",clp,X_elongate);
  if (!pdf) pdf = new normal_pdf(X_elongate,0.0,0.3,1.0); // a normal PDF is the default for AEMs
  else if (pdf->amplitude()>1.0) warning("Warning: Arbor Elongation Model has perturbation with amplitude greater than 1.0, retraction may occur!\n");
  if (label.empty()) return;
  int n;
  if ((n=clp.Specifies_Parameter(label+"aem_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  arbor_elongation_model_selection(label,clp,contributing,contributing);
}

void arbor_elongation_model_base::reset(double t) {
  prev_t = t;
  base_parameters.dL = 0.0;
  if (contributing) contributing->reset(t);
}

int arbor_elongation_model_base::total_terminal_segments(fibre_structure * arbor, growth_cone_competition_type scope) {
  // This function counts all the terminal segments in all the arbors that are linked in the list that contains "arbor".
  int res = 0;
  switch (scope) {
  case gcc_within_tree:
    res = arbor->count_terminal_segments();
    break;;
  case gcc_within_all_AorD:
    PLL_LOOP_FORWARD(fibre_structure,arbor->head(),1) res += e->count_terminal_segments();
    break;;
  default: // gcc_whole_neuron
    res = arbor->N()->total_input_terminal_segments() + arbor->N()->total_output_terminal_segments();
  };
  return res;
}

double arbor_elongation_model_base::predict(double weight, fibre_structure * arbor) {
  // This chain is called by elongation(), and therefore only once for
  // each t>prev_t.
  // [***INCOMPLETE] The manner of combination of chained arbor elongation
  // models is not yet based on a specific rationale (see TL#200603131455.1).
  double dL = weight*predict_elongation(arbor); // resources for arbor elongation
  if (contributing) dL += contributing->predict(contributingweight,arbor);
  return dL;
}

double arbor_elongation_model_base::elongation(fibre_structure * arbor) {
  // (TL#200603021224.1 documents the start of considerations about this
  // implementation.)
  if (eq->T()<=prev_t) return base_parameters.dL;
  // 1. predict elongation with this arbor elongation model
  base_parameters.dL = predict_elongation(arbor);
  // 2. combine with predictions by chained models
  if (contributing) {
    base_parameters.dL += contributing->predict(contributingweight,arbor);
  }
  // 3. probabilistic deviation
  // [***NOTE] Below is the case when NOT LEGACY_ELONGATION, perturbation
  // is normally added here!
#ifndef LEGACY_ELONGATION
#ifdef TESTING_AEM_RANDOM_DISTRIBUTION
  double r = pdf->random_selection();
  aemrandomdist[((unsigned int) (r*10.0))+10]++;
  base_parameters.dL *= (1.0+r);
#else
  base_parameters.dL *= (1.0+pdf->random_selection());
#endif
#endif
  // return the new cached value
  return base_parameters.dL; // resources for arbor elongation
}

String arbor_elongation_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting tocontributing models.
  String res(report_parameters_specific());
  res += " p_AEM=" + pdf->report_parameters();
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

//van_Pelt_arbor_elongation_model::van_Pelt_arbor_elongation_model(): _nu0(15.5/DAYSECONDS), _F(0.74) {
//}

van_Pelt_arbor_elongation_model::van_Pelt_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, van_Pelt_arbor_elongation_model & schema): arbor_elongation_model_base(aemcontrib,aemweight,schema), parameters(schema.parameters) {
  // This constructor is specifically meant to be used by clone() calls.
  // See TL#200603120707.2.
  //schema.copy_parameters(parameters);
  //_nu0 = schema.get_nu0();
  //_F = schema.get_F();
  // Chaining of models is taken care of in the base constructor.
}

van_Pelt_arbor_elongation_model::van_Pelt_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): arbor_elongation_model_base(thislabel,label,clp), parameters(15.5/DAYSECONDS,0.74,gcc_whole_neuron) {
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  // Detect schema specific parameter values
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"growth_nu0"))>=0) parameters._nu0 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"growth_F"))>=0) parameters._F = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"F_competes_with"))>=0) {
    if (downcase(clp.ParValue(n))==String("whole_neuron")) parameters._competition_with_all_trees = gcc_whole_neuron;
    else if (downcase(clp.ParValue(n))==String("all_dendrites")) parameters._competition_with_all_trees = gcc_within_all_AorD;
    else if (downcase(clp.ParValue(n))==String("all_axons")) parameters._competition_with_all_trees = gcc_within_all_AorD;
    else if (downcase(clp.ParValue(n))==String("same_arbor")) parameters._competition_with_all_trees = gcc_within_tree;
    else warning("Unknown F competition value:"+clp.ParValue(n)+'\n');
  }
  // Chaining of models is taken care of in the base constructor.
}

arbor_elongation_model_base * van_Pelt_arbor_elongation_model::clone() {
  // See tension_direction_model::clone() for a more complex example of a
  // cloning function.
  return new van_Pelt_arbor_elongation_model(contributing,contributingweight,*this);
}

double van_Pelt_arbor_elongation_model::predict_elongation(fibre_structure * arbor) {
  // prev_t should be set in here, since prev_t is cached separately in
  // each chained arbor elongation model, and this function is called both
  // by elongation() and by predict().
  // This model function is called directly or via predict() only once for
  // each update time. (See update time test above.)
  //#define TEST_VAN_PELT_ARBOR_ELONGATION_MODEL
#ifdef TEST_VAN_PELT_ARBOR_ELONGATION_MODEL
  int nofF = total_terminal_segments(arbor,parameters._competition_with_all_trees);
  int nofarbor = arbor->count_terminal_segments();
  double nt_F = exp(-parameters._F*log((double) nofF));
  nt_F *= (eq->T()-prev_t)*parameters._nu0; // van Pelt segment elongation scaled by N^(-F)
  // [***NOTE] fixed step alternative: nt_F *= _dt_nu0;
  prev_t = eq->T();
  double res = nt_F*((double) nofarbor);
  //if ((eq->T()<200.0) || (eq->T()>1.8143e+06)) cout << "At t=" << eq->T() << ": n(t)*dt*nu0*N^(-F) = (" << nofarbor << ")*(" << (eq->T()-prev_t) << ")*(" << parameters._nu0 << ")*(" << nofF << ")^(-" << parameters._F << ") = " << res << '\n';
  cout << " N=" << nofF << " n=" << nofarbor << " E[dL]=" << res;
  return res;
#else
  double nt_F = exp(-parameters._F*log((double) total_terminal_segments(arbor,parameters._competition_with_all_trees)));
  nt_F *= (eq->T()-prev_t)*parameters._nu0; // van Pelt segment elongation scaled by N^(-F)
  // [***NOTE] fixed step alternative: nt_F *= _dt_nu0;
  prev_t = eq->T();
  return nt_F*((double) arbor->count_terminal_segments());
#endif
}

String van_Pelt_arbor_elongation_model::report_parameters_specific() {
  String res(" van Pelt, nu0=");
  res += String(parameters._nu0*DAYSECONDS,"%.2f/day F=");
  res += String(parameters._F,"%.2f");
  res += " F competes ";
  switch (parameters._competition_with_all_trees) {
  case gcc_whole_neuron: res += "with whole neuron"; break;;
  case gcc_within_all_AorD: res += "with all axons or with all dendrites"; break;;
  default: res += "within the axon/dendrite arbor"; // gcc_within_tree
  };
  return res;
}

Polynomial_O1_arbor_elongation_model::Polynomial_O1_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, Polynomial_O1_arbor_elongation_model & schema): arbor_elongation_model_base(aemcontrib,aemweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

Polynomial_O1_arbor_elongation_model::Polynomial_O1_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): arbor_elongation_model_base(thislabel,label,clp), parameters(5.631771/DAYSECONDS,1.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"growth_nu0"))>=0) parameters._nu0 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"growth_F"))>=0) parameters._F = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
}

arbor_elongation_model_base * Polynomial_O1_arbor_elongation_model::clone() {
  return new Polynomial_O1_arbor_elongation_model(contributing,contributingweight,*this);
}

double Polynomial_O1_arbor_elongation_model::predict_elongation(fibre_structure * arbor) {
  // This model function is called directly or via predict() only once for
  // each update time. (See update time test above.)
  double dL = (eq->T()-prev_t)*parameters._nu0*exp((1.0 - parameters._F)*log((double) arbor->count_terminal_segments())); // resources for arbor elongation
  // [***NOTE] fixed step alternative: _dt_nu0*exp((1.0 - _F)*log((double) arbor->count_terminal_segments()));
  prev_t = eq->T();
  return dL;
}

String Polynomial_O1_arbor_elongation_model::report_parameters_specific() {
  String res(" Polynomial_01, nu0=");
  res += String(parameters._nu0*DAYSECONDS,"%.2f/day F=");
  res += String(parameters._F,"%.2f");
  return res;
}

Polynomial_O3_arbor_elongation_model::Polynomial_O3_arbor_elongation_model(arbor_elongation_model_base * aemcontrib, double & aemweight, Polynomial_O3_arbor_elongation_model & schema): arbor_elongation_model_base(aemcontrib,aemweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

Polynomial_O3_arbor_elongation_model::Polynomial_O3_arbor_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): arbor_elongation_model_base(thislabel,label,clp), parameters(0.0/DAYSECONDSSQR,0.006564/DAYSECONDSCUB,3.839628/DAYSECONDS,1.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"growth_nu0"))>=0) parameters._nu0 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"growth_nu1"))>=0) parameters._nu1 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"growth_nu2"))>=0) parameters._nu2 = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"growth_F"))>=0) parameters._F = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
}

arbor_elongation_model_base * Polynomial_O3_arbor_elongation_model::clone() {
  return new Polynomial_O3_arbor_elongation_model(contributing,contributingweight,*this);
}

double Polynomial_O3_arbor_elongation_model::predict_elongation(fibre_structure * arbor) {
  // L(t2) - L(t1), where each is approximated by a cubic polynomial function
  // This model function is called directly or via predict() only once for
  // each update time. (See update time test above.)
  double t2 = eq->T();
  double t1sq = prev_t*prev_t, t2sq = t2*t2;
  double dL = (parameters._nu0*(t2-prev_t)) + (parameters._nu1*(t2sq-t1sq)) + (parameters._nu2*((t2sq*t2)-(t1sq*prev_t)));
  prev_t = t2;
  return dL*exp((1.0 - parameters._F)*log((double) arbor->count_terminal_segments()));
}

String Polynomial_O3_arbor_elongation_model::report_parameters_specific() {
  String res(" Polynomial_03, nu0=");
  res += String(parameters._nu0*DAYSECONDSCUB,"%.2f/day^3 nu1=");
  res += String(parameters._nu1*DAYSECONDSCUB,"%.2f/day^3 nu2=");
  res += String(parameters._nu2*DAYSECONDSCUB,"%.2f/day^3 F=");
  res += String(parameters._F,"%.2f");
  return res;
}

terminal_segment_elongation_model_base::terminal_segment_elongation_model_base(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, terminal_segment_elongation_model_base & schema): contributing(NULL), contributingweight(1.0), prev_t(eq->T()), l_i_cache_t(eq->T()), base_parameters(schema.base_parameters) {
  // [***NOTE:] In the current implementation, any TSEM that requires ts during initialization (e.g. nonnorm_BESTL) must be at the head of a chain.
  if (tsemcontrib) { // Clone chained models
    contributingweight = tsemweight;
    contributing = tsemcontrib->clone(NULL);
  }
  pdf = schema.pdf->clone();
  // [***INCOMPLETE] The REDUCED_MEMORY_PTSEM option is in need of a place for initialization, where (a) the arbor can be reached in order to set tsem_pdf, where (b) the applicable prefix is known (e.g. "all_axons."), and (c) needs to take care of deleting the allocated PDF object as well. Also, I need to add more compiler directives around references to pdf when compiling with REDUCED_MEMORY_PTSEM.
}

terminal_segment_elongation_model_base::terminal_segment_elongation_model_base(String & thislabel, String & label, Command_Line_Parameters & clp): contributing(NULL), contributingweight(1.0), prev_t(eq->T()), l_i_cache_t(eq->T()), base_parameters(0.0) {
  // This constructor is intended for use by schemas, which can set up model
  // chains, as well as look for general (thislabel=="") or subset parameter
  // settings that differ from the compile-time defaults.
  // (See TL#200603120618.5.)
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  pdf = pdfselection(thislabel+"tsem.",clp,X_elongate);
  if (!pdf) pdf = new normal_pdf(X_elongate,0.0,0.3,1.0); // a normal PDF is the default for TSEMs
  else {
    if (pdf->amplitude()>1.0) warning("Warning: Terminal Segment Elongation Model has perturbation with amplitude greater than 1.0, retraction may occur!\n");
    else if ((pdf->amplitude()==0.0) && (pdf->mean_X()>0.0)) warning("Warning: TSEM perturbation is added to 1.0 before multiplication with the computed elongation speed. When using a delta PDF this means that values greater than zero can cause ramping up to 'infinite' elongation speed and resulting invalid computations.\n");
  }
  if (label.empty()) return;
  int n;
  if ((n=clp.Specifies_Parameter(label+"tsem_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  terminal_segment_elongation_model_selection(label,clp,contributing,contributing);
}

void terminal_segment_elongation_model_base::reset(double t) {
  // [***NOTE] This only appears to be used in a "glimpse ahead" in network.cc.
  prev_t = t;
  l_i_cache_t = t;
  base_parameters.l_i_cache = 0.0;
  if (contributing) contributing->reset(t);
}

void terminal_segment_elongation_model_base::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_elongation_model(clone(&ts2));
  ts1.set_elongation_model(this);
}

void terminal_segment_elongation_model_base::handle_turning(terminal_segment & ts) {
  ts.set_elongation_model(this);
}

double terminal_segment_elongation_model_base::predict(double weight) {
  // This chain is called by perturbed_expected_elongation(), and therefore
  // only once for each t>l_i_cache_t.
  // This function is called by:
  //   perturbed_expected_elongation()
  // This function calls:
  //   predict_elongate()
  //   predict()
  // [***INCOMPLETE] The manner of combination of chained arbor elongation
  // models is not yet based on a specific rationale (see TL#200603140534.1).
  double l_i = weight*predict_elongate();
  if (contributing) l_i += contributing->predict(contributingweight);
  prev_t=eq->T(); // store for used by calculations in chained model
  return l_i;
}

#ifdef REDUCED_MEMORY_PTSEM
probability_distribution_function * tsem_pdf_cache = NULL;
#endif

double terminal_segment_elongation_model_base::perturbed_expected_elongation() {
  // This function is called by:
  //   elongate()
  // This function calls:
  //   predict_elongate()
  //   predict()
  // See the comment about the base_parameters.l_i_cache in the header file.
  // [***INCOMPLETE] The manner of combination of chained arbor elongation
  // models is not yet based on a specific rationale (see TL#200603140534.1).
  if (eq->T()<=l_i_cache_t) return base_parameters.l_i_cache;
  // 1. predict expected elongation
  // Using different temporary variable in case contributing models need
  // the original cached value.
  double l_i = predict_elongate();
  // 2. combine with predictions by chained models
  if (contributing) {
    l_i += contributing->predict(contributingweight);
  }
  // 3. probabilistic deviation (i.e. perturbed expected elongation)
#ifdef LEGACY_ELONGATION
  base_parameters.l_i_cache = l_i;
#else
#ifdef REDUCED_MEMORY_PTSEM
  base_parameters.l_i_cache = l_i*(1.0+tsem_pdf_cache->random_selection());
#else
  base_parameters.l_i_cache = l_i*(1.0+pdf->random_selection());
#endif
#endif
  // [***INCOMPLETE] remember to provide a uniform random pdf option for
  // legacy testing purposes (e.g. 1.0*random_double(0.2,1.8))
  l_i_cache_t = eq->T();
  return base_parameters.l_i_cache;
}

void terminal_segment_elongation_model_base::elongate(terminal_segment * ts) {
  // PROTOCOL: This is called during growth.
  // This function calls:
  //   perturbed_expected_elongation()
  //   sum_of_perturbed_expected_elongations()
  //   ElongationModel()->elongation()
  //   ts->add_length()
  //   ts->update_segment_vector()
  // (TL#200603021224.1 documents the start of considerations about this
  // implementation.)
  if (eq->T()<=prev_t) return;
  // 1., 2. & 3 compute perturbed expected elongation
  // The three regular steps of the model call have been turned into
  // a separate function here, so that the calculations can be invoked
  // for all terminal segments as needed, while the sum of those results
  // is needed for the elongation by proportional allocation of resources.
  // [***NOTE] If the elongation calculations for specific terminal segments
  // need to access other information about the terminal segment, or about
  // other objects in the ts->Arbor(), then we can pass the ts parameter to
  // other member functions, such as perturbed_expected_elongation().
#ifdef REDUCED_MEMORY_PTSEM
  tsem_pdf_cache = ts->Arbor()->TSEM_Pdf();
#endif
  double l_i = perturbed_expected_elongation();
  // 4. normalize with regard to the sum of all perturbed expected elongations
  //    of all terminal segments on the same arbor
  register fibre_structure * arbor = ts->Arbor();
  double sum_l_i = arbor->sum_of_perturbed_expected_elongations();
  if (sum_l_i!=0.0) l_i /= sum_l_i; // otherwise all are "delayed"
#ifdef LEGACY_ELONGATION
  l_i *= random_double(0.2,1.8);
#endif
  // 5. elongate by proportion of available resources
  double L_i = l_i * arbor->ElongationModel()->elongation(arbor);
#ifdef TESTING_ELONGATION_TOTAL
  usedarborelongationfraction += l_i;
#endif
#ifdef DEBUGGING_ELONGATION
  if (L_i<0.0) debugging.tsemb_elongate_decreasing++;
  else if (L_i==0.0) debugging.tsemb_elongate_stopped++;
#endif
  ts->add_length(L_i);
  if (figattr_update_terminal_segments_visibly) ts->update_segment_vector();
  prev_t = eq->T();
  return;
}

String terminal_segment_elongation_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting tocontributing models.
  String res(report_parameters_specific());
  res += " p_TSEM=" + pdf->report_parameters();
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

// As described in TL#200603071348.1 and in the methods and future
// directions in the framework-generation paper, the elongation function
// (1) checks how much time has passed since the last elongation update
// for this terminal segment, (2) requests resources available to the
// arbor, and (3) determines this terminal segment's updated proportion
// of resource consumption. Then the terminal segment length can be updated.

simple_terminal_segment_elongation_model::simple_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, simple_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

simple_terminal_segment_elongation_model::simple_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp) {
  // This simple model has no local parameters to specify.
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * simple_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new simple_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double simple_terminal_segment_elongation_model::predict_elongate() {
  // This model function is called directly by perturbed_expected_elongation()
  // or via predict() (also due to perturbed_expected_elongation()) only once
  // for each update time. (See cache time and update time tests above.)
  return 1.0; // expected_L_i
}

String simple_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" simple");
  return res;
}

// *** inertia may want to use a cached result of the previous NORMALIZED
//     perturbed expected elongation, and then store its calculated result
//     in the l_i_cache for access by other segment updates with comparable
//     scale. I have to read the paper very carefully, to see exactly what
//     values are to be used in which steps.

inertia_terminal_segment_elongation_model::inertia_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, inertia_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

inertia_terminal_segment_elongation_model::inertia_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * inertia_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new inertia_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double inertia_terminal_segment_elongation_model::predict_elongate() {
  // This model function is called directly by perturbed_expected_elongation()
  // or via predict() (also due to perturbed_expected_elongation()) only once
  // for each update time. (See cache time and update time tests above.)
  // [***NOTE] The paper describes this as being normalized, but since the
  // normalization is applied equally to all terminal segments, it is not
  // necessary, unless we wish to express the predicted elongation of a
  // new terminal segment in relation to a proportions between 0.0 and 1.0.
  // Note: At the head of a chain, the combined result provides inertia,
  // elsewhere in a chain, a local inertial component is contributed.
  if (base_parameters.l_i_cache==0.0) base_parameters.l_i_cache = 1.0+pdf->random_selection(); // updated, see TL#200709150740.1
  return base_parameters.l_i_cache;
}

String inertia_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" inertia");
  return res;
}

second_order_terminal_segment_elongation_model::second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, second_order_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
}

second_order_terminal_segment_elongation_model::second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp), parameters(0.5,86400.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tsem_inertia"))>=0) parameters.inertia = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * second_order_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new second_order_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double second_order_terminal_segment_elongation_model::predict_elongate() {
  // This model function is called directly by perturbed_expected_elongation()
  // or via predict() (also due to perturbed_expected_elongation()) only once
  // for each update time. (See cache time and update time tests above.)
  // [***NOTE] The paper describes this as being normalized, but since the
  // normalization is applied equally to all terminal segments, it is not
  // necessary, unless we wish to express the predicted elongation of a
  // new terminal segment in relation to a proportions between 0.0 and 1.0.
  // The acceleration in the first terminal segment of an arbor is initialized
  // to 0.5, so that all terminal segments initially tend to increase their
  // relative speed of elongation. The acceleration is constrained between
  // hard-coded limits of -1.0 and 1.0, between which perturbations are
  // scaled by inertia. The resulting elongation is always >= 0.0.
  // Note: At the head of the chain, the second order model is applied to the
  // perturbed result of chained models. Elsewhere in the chain, a local
  // variable is maintained in base_parameters.l_i_cache, which is updated
  // without further perturbation and contributes to the chain response.
  double dt = eq->T()-prev_t;
  parameters.l_i_acceleration += (dt/parameters.inertia)*pdf->random_selection(); // Note: "l_i_acceleration" is a bit of a misnomer, it is called the velocity in the paper.
  if (parameters.l_i_acceleration>1.0) parameters.l_i_acceleration=1.0;
  else if (parameters.l_i_acceleration<-1.0) parameters.l_i_acceleration=-1.0;
  // modeled analogous to v = u + at
#define expected_L_i base_parameters.l_i_cache
  expected_L_i = base_parameters.l_i_cache + parameters.l_i_acceleration*dt; // updated, see TL#200709150740.1, was separate "double expected_L_i" local variable
  if (expected_L_i<0.0) return 0.0;
  return expected_L_i;
#undef expected_L_i
}

String second_order_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" second order, inertia=");
  res += String(parameters.inertia,"%.2f");
  return res;
}

constrained_second_order_terminal_segment_elongation_model::constrained_second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, constrained_second_order_terminal_segment_elongation_model & schema): second_order_terminal_segment_elongation_model(tsemcontrib,tsemweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

constrained_second_order_terminal_segment_elongation_model::constrained_second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): second_order_terminal_segment_elongation_model(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * constrained_second_order_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new constrained_second_order_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double constrained_second_order_terminal_segment_elongation_model::predict_elongate() {
  double expected_L_i = second_order_terminal_segment_elongation_model::predict_elongate();
  // With an initial acceleration of 0.5 and assuming mean perturbation around
  // zero, some terminal segments may reach a relative speed of DAYSECONDS in
  // seven days.
#define TSEM_SPEED_CONSTRAINT 3.5*DAYSECONDS
  if (expected_L_i>TSEM_SPEED_CONSTRAINT) return TSEM_SPEED_CONSTRAINT;
  return expected_L_i;
}

initialized_CSO_terminal_segment_elongation_model::initialized_CSO_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, initialized_CSO_terminal_segment_elongation_model & schema): constrained_second_order_terminal_segment_elongation_model(tsemcontrib,tsemweight,schema), icsoparameters(schema.icsoparameters) {
  // Beware: That which I am calling "l_i_acceleration" in the source code
  // is $\dot{\TSEMquote}(t)$ in ~/doc/tex/generation-framework/generation-framework.tex.
  // This stems from the former emphasis of the update of fiber length over
  // time as an elongation rate, and there for the rate of change in that
  // was an acceleration.
  // Chaining of models is taken care of in the base constructor.
  if (icsoparameters.initial_l_i_acceleration>=0.0) parameters.l_i_acceleration = icsoparameters.initial_l_i_acceleration;
  // Otherwise the cloned value from parent or schema parameters is used.
}

#define ICSO_CLONE_PARENT_OR_SCHEMA -1.0
initialized_CSO_terminal_segment_elongation_model::initialized_CSO_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): constrained_second_order_terminal_segment_elongation_model(thislabel,label,clp), icsoparameters(ICSO_CLONE_PARENT_OR_SCHEMA) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tsem_initial_acceleration"))>=0) icsoparameters.initial_l_i_acceleration = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
  if (icsoparameters.initial_l_i_acceleration>=0.0) parameters.l_i_acceleration = icsoparameters.initial_l_i_acceleration;
  // Through this constructor, no parent or schema information is received
  // (this object likely IS a schema), so that the first initial value of the
  // acceleration is set by the base class (probable default = 0.5).
}

terminal_segment_elongation_model_base * initialized_CSO_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new initialized_CSO_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

String initialized_CSO_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" initialized "+constrained_second_order_terminal_segment_elongation_model::report_parameters());
  res += String(icsoparameters.initial_l_i_acceleration," tsem_initial_acceleration=%.2f (0.0 simulates delayed branching)");
  return res;
}

decaying_second_order_terminal_segment_elongation_model::decaying_second_order_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, decaying_second_order_terminal_segment_elongation_model & schema): second_order_terminal_segment_elongation_model(tsemcontrib,tsemweight,schema), dsoparameters(schema.dsoparameters) {
  // Chaining of models is taken care of in the base constructor.
}

decaying_second_order_terminal_segment_elongation_model::decaying_second_order_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): second_order_terminal_segment_elongation_model(thislabel,label,clp), dsoparameters(DAYSECONDS/4.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tsem_decay"))>=0) dsoparameters.t_decay = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * decaying_second_order_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new decaying_second_order_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double decaying_second_order_terminal_segment_elongation_model::predict_elongate() {
  // As in the second order model, except that acceleration decays toward zero
  // unless perturbed.
  // Note: At the head of the chain, the decaying second order model is applied to
  // the perturbed result of chained models. Elsewhere in the chain, a local
  // variable is maintained in base_parameters.l_i_cache, which is updated
  // without further perturbation and contributes to the chain response.
  double dt = eq->T()-prev_t;
  parameters.l_i_acceleration *= exp(-dt/dsoparameters.t_decay);
  parameters.l_i_acceleration += (dt/parameters.inertia)*pdf->random_selection();
  if (parameters.l_i_acceleration>1.0) parameters.l_i_acceleration=1.0;
  else if (parameters.l_i_acceleration<-1.0) parameters.l_i_acceleration=-1.0;
  // modeled analogous to v = u + at
#define expected_L_i base_parameters.l_i_cache
  expected_L_i = base_parameters.l_i_cache + parameters.l_i_acceleration*dt; // updated, see TL#200709150740.1, was separate "double expected_L_i" local variable
  if (expected_L_i<0.0) return 0.0;
  return expected_L_i;
#undef expected_L_i
}

String decaying_second_order_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" decaying "+second_order_terminal_segment_elongation_model::report_parameters());
  res += String(dsoparameters.t_decay," t_decay=%.2f");
  return res;
}

delayed_branch_terminal_segment_elongation_model::delayed_branch_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, delayed_branch_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema), parameters(schema.parameters) {
  // Chaining of models is taken care of in the base constructor.
  // [***INCOMPLETE] This requires a random_selection function with a
  // decaying probability, perhaps just a random_positive of a normal
  // function, using the delay time constant.
  elongation_commencement = eq->T() + 0.0;
}

delayed_branch_terminal_segment_elongation_model::delayed_branch_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp), parameters(3600.0) {
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tsem_tdelay"))>=0) parameters.delay_time_constant = atof(clp.ParValue(n));
  // Chaining of models is taken care of in the base constructor.
  // [***INCOMPLETE] This requires a random_selection function with a
  // decaying probability, perhaps just a random_positive of a normal
  // function, using the delay time constant.
  // [***INCOMPLETE] We might want to avoid a delay for the first segment,
  // but perhaps not, since that could allow different fiber structures to
  // begin growing at slightly different times.
  // The paper specifies that this should be a model option. Implement this
  // DBSTEM_at_first_segment parameter.
  elongation_commencement = eq->T() + 0.0;
}

terminal_segment_elongation_model_base * delayed_branch_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new delayed_branch_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double delayed_branch_terminal_segment_elongation_model::predict_elongate() {
  // An exponential function with delay_time_constant is used to determine
  // the likelihood that a delayed branch begins to elongate. When that
  // happens, the terminal segment continues to elongate according to a
  // different elongation function. Until then, the proportion of elongation
  // resources requested is 0.0.
  if (eq->T()<elongation_commencement) return 0.0;
  // *** old method: double min_dt = parameters.creation_time - eq->T();
  // *** old method: double threshold = exp(min_dt/parameters.delay_time_constant);
  // *** old method: if (random_double()<=threshold) return 0.0; // no elongation yet
  // [***INCOMPLETE] The case where elongation begins.
  // This specific model should always return 0.0 (see paper), but it should
  // activating the chain of models at this point if it is not already
  // active.
  // There are two ways to control the application of chained models and
  // consequently two ways to remember the chain during the delay period:
  // 1) Detach the chain and attach it to a pointer cache variable.
  // 2) Define a replacement for the perturbed_expected_elongation()
  //    function that can take the delay into account by testing
  //    if ((eq->T()>=elongation_commencement) && contributing) rather than
  //    just if (contributing).
  return 0.0;
}

String delayed_branch_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" delayed, tdelay=");
  res += String(parameters.delay_time_constant,"%.2f");
  return res;
}

BESTL_terminal_segment_elongation_model::BESTL_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, BESTL_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema), branchpdf(NULL) {
  // Chaining of models is taken care of in the base constructor.
  branchpdf = schema.branchpdf->clone();
}

BESTL_terminal_segment_elongation_model::BESTL_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
  base_parameters.l_i_cache = 1.0; // this is the assumed mean of the BESTL quota distribution
  // Note: The mean can shift if other models contribute to the elongation quota determination.
  branchpdf = pdfselection(thislabel+"tsem.branch.",clp,X_elongate);
  if (!branchpdf) branchpdf = new normal_pdf(X_elongate,2.0,1.0); // a normal PDF is the default for TSEM branch length initialization
}

terminal_segment_elongation_model_base * BESTL_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  return new BESTL_terminal_segment_elongation_model(contributing,contributingweight,*this);
}

double BESTL_terminal_segment_elongation_model::predict_elongate() {
  // This model expects no change of elongation rate until the next bifurcation.
  // The base_parameters.l_i_cache value is initialized in the constructor and
  // by an Elongation_Rate_Initialization model.
  // Note: If the BESTL_TSEM is at the head of a chain then the use of the
  // base_parameters.l_i_cache means that it attempts to maintain the rate of
  // elongation computed by the chain (unless further modified by chained
  // contributions). Elsewhere in a chain, the BESTL_TSEM contributes a stable
  // expected elongation quota.
  // For an unperturbed BESTL_TSEM, use the perturbation PDF delta_PDF (with
  // the default mean value 0.0).
  return base_parameters.l_i_cache; // see TL#200709150740.1
}

void BESTL_terminal_segment_elongation_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_elongation_model(clone(&ts2));
  ts1.set_elongation_model(this);
  if (branchpdf) {
    ts1.AngularCoords().add_X(branchpdf->random_positive());
    ts2.AngularCoords().add_X(branchpdf->random_positive());
  }
}

String BESTL_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" BESTL: branchPDF: ");
  if (branchpdf) res += branchpdf->report_parameters();
  return res;
}

BESTL_nonnormalizing_terminal_segment_elongation_model::BESTL_nonnormalizing_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, BESTL_nonnormalizing_terminal_segment_elongation_model & schema): BESTL_terminal_segment_elongation_model(tsemcontrib,tsemweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

BESTL_nonnormalizing_terminal_segment_elongation_model::BESTL_nonnormalizing_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): BESTL_terminal_segment_elongation_model(thislabel,label,clp) {
  // [*** Note:] See TL#200807020156.5 - WITHDRAWN! See TL#200807020156.10.
  //probability_distribution_function * rootratepdf = pdfselection(thislabel+"tsem.rootrate.",clp,X_elongate);
  //if (rootratepdf) base_parameters.l_i_cache = rootratepdf->random_positive();
  //else base_parameters.l_i_cache = 0.00013889;
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * BESTL_nonnormalizing_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  terminal_segment_elongation_model_base * tsemb = new BESTL_nonnormalizing_terminal_segment_elongation_model(contributing,contributingweight,*this);
  if (ts) {
    if ((!ts->Prev()) && (!ts->Next())) { // We still use the regular ERI->initialize() call from the terminal_segment::branch() function when branching
      if (!(ts->ElongationRateInitializationModel())) error("ElongationRateInitializationModel() is NULL in BESTL_nonnormalizing_terminal_segment_elongation_model::clone().\n");
      else ts->ElongationRateInitializationModel()->root_initialize(tsemb->base_parameters.l_i_cache); // Insures that initial segments also draw from the ERI distribution (see TL#200807020156.10)
    }
  }
  return tsemb;
}

void BESTL_nonnormalizing_terminal_segment_elongation_model::elongate(terminal_segment * ts) {
  // PROTOCOL: This is called during growth.
  // This function calls:
  //   perturbed_expected_elongation()
  //   ElongationModel()->elongation()
  //   ts->add_length()
  //   ts->update_segment_vector()
  if (eq->T()<=prev_t) return;
  // 1., 2. & 3 compute perturbed expected elongation
  // The three regular steps of the model call have been turned into
  // a separate function here, so that the calculations can be invoked
  // for all terminal segments as needed, while the sum of those results
  // is needed for the elongation by proportional allocation of resources.
  // [***NOTE] If the elongation calculations for specific terminal segments
  // need to access other information about the terminal segment, or about
  // other objects in the ts->Arbor(), then we can pass the ts parameter to
  // other member functions, such as perturbed_expected_elongation().
#ifdef REDUCED_MEMORY_PTSEM
  tsem_pdf_cache = ts->Arbor()->TSEM_Pdf();
#endif
  double l_i = perturbed_expected_elongation();
  // Note: There is no normalization to other growth cones here!
#ifdef LEGACY_ELONGATION
  l_i *= random_double(0.2,1.8);
#endif
  // 5. elongate by proportion of available resources
  // double L_i = l_i * ts->Arbor()->ElongationModel()->elongation(ts->Arbor()); // [***NOTE] See TL#200807020156.5.
  double L_i = l_i * (eq->T()-prev_t); // All is done locally for these absolute elongation rates! (See TL#200807020156.5.)
#ifdef TESTING_ELONGATION_TOTAL
  usedarborelongationfraction += l_i;
#endif
#ifdef DEBUGGING_ELONGATION
  if (L_i<0.0) debugging.tsemb_elongate_decreasing++;
  else if (L_i==0.0) debugging.tsemb_elongate_stopped++;
#endif
  ts->add_length(L_i);
  // Test actual elongation pieces created per dt:
  //cout << String(L_i,"_%.8f") << '&' << String(l_i,"%.8f");
  if (figattr_update_terminal_segments_visibly) ts->update_segment_vector();
  prev_t = eq->T();
  return;
}

String BESTL_nonnormalizing_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" non-normalizing ");
  res += BESTL_terminal_segment_elongation_model::report_parameters_specific();
  return res;
}

BESTLNN_pyramidal_AD_terminal_segment_elongation_model::BESTLNN_pyramidal_AD_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, BESTLNN_pyramidal_AD_terminal_segment_elongation_model & schema): BESTL_nonnormalizing_terminal_segment_elongation_model(tsemcontrib,tsemweight,schema), parameters(schema.parameters), trunklengthpdf(NULL), obliquespdf(NULL), obliqueanglepdf(NULL) {
  // Chaining of models is taken care of in the base constructor.
  trunklengthpdf = schema.trunklengthpdf->clone();
  // Note that we do not generate a new parameters.trunk_length_SQ value here, so that the expected trunk length is not changed when passing branch nodes that connect to obliques.
  obliquespdf = schema.obliquespdf->clone();
  // Note that we do not generate a new number of obliques here, so that the expected number of oblique branches does not change as we pass branch nodes that connect to them.
  if (schema.obliqueanglepdf) obliqueanglepdf = schema.obliqueanglepdf->clone();
}

void BESTLNN_pyramidal_AD_terminal_segment_elongation_model::initialize_parameters() {
  // This is used to initialize trunk length and number of obliques parameters when new
  // PDFs may be selected by a constructor or when the TSEM of a root fiber terminal
  // segment is allocated (see TL#200808050218.2).
  parameters.trunk_length_SQ = trunklengthpdf->random_positive();
  parameters.trunk_length_SQ *= parameters.trunk_length_SQ;
  parameters.num_obliques = (int) obliquespdf->random_positive();
  parameters.dist2nextoblique = 0;
  set_next_oblique_distance();
}

BESTLNN_pyramidal_AD_terminal_segment_elongation_model::BESTLNN_pyramidal_AD_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): BESTL_nonnormalizing_terminal_segment_elongation_model(thislabel,label,clp), parameters(700.0*700.0,0,DBL_MAX,"pyrAP"), trunklengthpdf(NULL), obliquespdf(NULL), obliqueanglepdf(NULL) {
  // [*** Note:] See TL#200807020156.5 - WITHDRAWN! See TL#200807020156.10.
  trunklengthpdf = pdfselection(thislabel+"tsem.trunklength.",clp,X_elongate);
  if (!trunklengthpdf) trunklengthpdf = new normal_pdf(X_elongate,700.0,100.0);
  obliquespdf = pdfselection(thislabel+"tsem.obliques.",clp,X_elongate);
  if (!obliquespdf) obliquespdf = new normal_pdf(X_elongate,7.0,3.0);
  obliqueanglepdf = pdfselection(thislabel+"tsem.obliqueangle.",clp,X_elongate); // if NULL then there is no perturbation
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"tsem.prefix"))>=0) parameters.prefix = clp.ParValue(n);
  initialize_parameters();// Here we generate new values, since the PDFs may be new (see TL#200808050218.2).
  // Chaining of models is taken care of in the base constructor.
}

terminal_segment_elongation_model_base * BESTLNN_pyramidal_AD_terminal_segment_elongation_model::clone(terminal_segment * ts) {
  BESTLNN_pyramidal_AD_terminal_segment_elongation_model * tsemb = new BESTLNN_pyramidal_AD_terminal_segment_elongation_model(contributing,contributingweight,*this);
  if (ts) {
    if ((!ts->Prev()) && (!ts->Next())) { // We still use the regular ERI->initialize() call from the terminal_segment::branch() function when branching
      tsemb->initialize_parameters(); // Note: Fixed according to TL#200808050218.2.
      if (!(ts->ElongationRateInitializationModel())) error("ElongationRateInitializationModel() is NULL in BESTL_nonnormalizing_terminal_segment_elongation_model::clone().\n");
      else ts->ElongationRateInitializationModel()->root_initialize(tsemb->base_parameters.l_i_cache); // Insures that initial segments also draw from the ERI distribution (see TL#200807020156.10)
    }
  }
  return tsemb;
}

void BESTLNN_pyramidal_AD_terminal_segment_elongation_model::set_next_oblique_distance() {
  if (parameters.num_obliques>0) {
    parameters.dist2nextoblique = sqrt(parameters.dist2nextoblique) + (sqrt(parameters.trunk_length_SQ)/((double) parameters.num_obliques));
    parameters.dist2nextoblique *= parameters.dist2nextoblique;
    //cout << "NUM OBLIQUES: " << parameters.num_obliques << " NEXT OBLIQUE: " << parameters.dist2nextoblique << '\n';
  } else parameters.dist2nextoblique = DBL_MAX;
}

void BESTLNN_pyramidal_AD_terminal_segment_elongation_model::set_tuft_models(terminal_segment * ts) {
  // BEWARE: This function causes a self-destruct of the present terminal segment. Do not try to
  // modify any variables of the terminal segment object or its dependents after calling this
  // function!
  //cout << "TUFTING!"; cout.flush();
  if (Txt_tuftrootbranchnodelist) (*Txt_tuftrootbranchnodelist) += String((long) ts->TerminalSegment()) + '\n';
  ts->branch(NULL);
  String prefixstr(parameters.prefix);
  if (!prefixstr.empty()) prefixstr += ".tuft.";
  else prefixstr = "tuft.";
  // TSEMs
  terminal_segment_elongation_model_base * tsembptr1 = NULL;
  terminal_segment_elongation_model_selection(prefixstr,*main_clp,terminal_segment_elongation_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],tsembptr1);
  most_recent_branch1_ts->set_elongation_model(tsembptr1);
  terminal_segment_elongation_model_base * tsembptr2 = tsembptr1->clone(most_recent_branch2_ts);
  most_recent_branch2_ts->set_elongation_model(tsembptr2);
  // ERI models
  elongation_rate_initialization_model_base * eribptr = NULL;
  elongation_rate_initialization_model_selection(prefixstr,*main_clp,elongation_rate_initialization_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],eribptr);
  most_recent_branch1_ts->set_ERI_model(eribptr);
  eribptr->root_initialize(tsembptr1->base_parameters.l_i_cache); // reinitialize according to the new ERI model
  eribptr = eribptr->clone();
  most_recent_branch2_ts->set_ERI_model(eribptr);
  eribptr->root_initialize(tsembptr2->base_parameters.l_i_cache); // reinitialize according to the new ERI model
  // TSB models
  TSBM_base * tsbmbptr = NULL;
  TSBM_selection(prefixstr,*main_clp,TSBM_region_subset_schemas[0][all_pyramidal_dendrites_sps],tsbmbptr,most_recent_branch1_ts);
  most_recent_branch1_ts->set_TSBM_model(tsbmbptr);
  most_recent_branch2_ts->set_TSBM_model(tsbmbptr->clone(most_recent_branch2_ts));
  // Branch angle models
  branch_angle_model_base * bambptr = NULL;
  branch_angle_model_selection(prefixstr,*main_clp,branch_angle_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],bambptr);
  most_recent_branch1_ts->set_branch_angle_model(bambptr);
  most_recent_branch2_ts->set_branch_angle_model(bambptr->clone());
  // [***NOTE] Possibly recompute branch angles here in a manner similar to that used in terminal_segment::branch().
  // Direction models
  direction_model_base * dmbptr = NULL;
  direction_model_selection(prefixstr,*main_clp,direction_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],dmbptr);
  most_recent_branch1_ts->set_direction_model(dmbptr);
  most_recent_branch2_ts->set_direction_model(dmbptr->clone());
}

void BESTLNN_pyramidal_AD_terminal_segment_elongation_model::set_oblique_models(terminal_segment * ts) {
  // BEWARE: This function causes a self-destruct of the present terminal segment. Do not try to
  // modify any variables of the terminal segment object or its dependents after calling this
  // function!
  //cout << "OBLIQUE!"; cout.flush();
  if (Txt_obliquerootbranchnodelist) (*Txt_obliquerootbranchnodelist) += String((long) ts->TerminalSegment()) + '\n';
  // 1. For obliques, preserve the growth angle of the main branch and prepare an oblique angle
  spatial mainbranchacoords(ts->AngularCoords());
  mainbranchacoords.set_X(0.0);
  double veerangle = M_PI/2.0;
  if (obliqueanglepdf) veerangle += obliqueanglepdf->random_selection();
#ifdef VECTOR3D
  spatial obliqueacoords;
  veer(ts,veerangle,obliqueacoords,-1.0); // for the oblique branch
#endif
#ifdef VECTOR2D
  spatial obliqueacoords(ts->AngularCoords());
  if (X_turn.get_rand_real1()<0.5) obliqueacoords.add_Y(-veerangle);
  else obliqueacoords.add_Y(veerangle);
#endif
  obliqueacoords.set_X(0.0);
  ts->branch(NULL);
  // 2. Modify the growth angles of oblique branches to preserve the angle of the main branch and give obliques orthogonal angles
  most_recent_branch1_ts->set_angularcoords(mainbranchacoords); // Here, environmental pressures are also applied.
  most_recent_branch1_ts->update_segment_vector();
  most_recent_branch2_ts->set_angularcoords(obliqueacoords); // Here, environmental pressures are also applied.
  most_recent_branch2_ts->update_segment_vector();
  // 3. Set models for the oblique branch
  String prefixstr(parameters.prefix);
  if (!prefixstr.empty()) prefixstr += ".oblique.";
  else prefixstr = "oblique.";
  // TSEMs
  terminal_segment_elongation_model_base * tsembptr = NULL;
  terminal_segment_elongation_model_selection(prefixstr,*main_clp,terminal_segment_elongation_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],tsembptr);
  // most_recent_branch1_ts retains the current apical trunk settings
  most_recent_branch2_ts->set_elongation_model(tsembptr);
  // ERI models
  elongation_rate_initialization_model_base * eribptr = NULL;
  elongation_rate_initialization_model_selection(prefixstr,*main_clp,elongation_rate_initialization_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],eribptr);
  // most_recent_branch1_ts retains the current apical trunk settings
  most_recent_branch2_ts->set_ERI_model(eribptr);
  eribptr->root_initialize(tsembptr->base_parameters.l_i_cache); // reinitialize according to the new ERI model
  // TSB models
  TSBM_base * tsbmbptr = NULL;
  TSBM_selection(prefixstr,*main_clp,TSBM_region_subset_schemas[0][all_pyramidal_dendrites_sps],tsbmbptr,most_recent_branch2_ts);
  most_recent_branch2_ts->set_TSBM_model(tsbmbptr);
  // Branch angle models
  branch_angle_model_base * bambptr = NULL;
  branch_angle_model_selection(prefixstr,*main_clp,branch_angle_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],bambptr);
  most_recent_branch2_ts->set_branch_angle_model(bambptr);
  // [***NOTE] Possibly recompute branch angles here in a manner similar to that used in terminal_segment::branch().
  // Direction models
  direction_model_base * dmbptr = NULL;
  direction_model_selection(prefixstr,*main_clp,direction_model_region_subset_schemas[0][all_pyramidal_dendrites_sps],dmbptr);
  most_recent_branch2_ts->set_direction_model(dmbptr);
}

void BESTLNN_pyramidal_AD_terminal_segment_elongation_model::elongate(terminal_segment * ts) {
  // PROTOCOL: This is called during growth.
  // This function calls:
  //   perturbed_expected_elongation()
  //   ElongationModel()->elongation()
  //   ts->add_length()
  //   ts->update_segment_vector()
  if (eq->T()<=prev_t) return;
  // 1., 2. & 3 compute perturbed expected elongation
  // The three regular steps of the model call have been turned into
  // a separate function here, so that the calculations can be invoked
  // for all terminal segments as needed, while the sum of those results
  // is needed for the elongation by proportional allocation of resources.
  // [***NOTE] If the elongation calculations for specific terminal segments
  // need to access other information about the terminal segment, or about
  // other objects in the ts->Arbor(), then we can pass the ts parameter to
  // other member functions, such as perturbed_expected_elongation().
#ifdef REDUCED_MEMORY_PTSEM
  tsem_pdf_cache = ts->Arbor()->TSEM_Pdf();
#endif
  double l_i = perturbed_expected_elongation();
  // Note: There is no normalization to other growth cones here!
#ifdef LEGACY_ELONGATION
  l_i *= random_double(0.2,1.8);
#endif
  // 5. elongate by proportion of available resources
  // double L_i = l_i * ts->Arbor()->ElongationModel()->elongation(ts->Arbor()); // [***NOTE] See TL#200807020156.5.
  double L_i = l_i * (eq->T()-prev_t); // All is done locally for these absolute elongation rates! (See TL#200807020156.5.)
#ifdef TESTING_ELONGATION_TOTAL
  usedarborelongationfraction += l_i;
#endif
#ifdef DEBUGGING_ELONGATION
  if (L_i<0.0) debugging.tsemb_elongate_decreasing++;
  else if (L_i==0.0) debugging.tsemb_elongate_stopped++;
#endif
  ts->add_length(L_i);
  // Test actual elongation pieces created per dt:
  //cout << String(L_i,"_%.8f");
  // Is it time to tuft?
  ts->update_segment_vector();
  double dist2tosoma = ts->TerminalSegment()->cartesian_soma_center_distance2();
  prev_t = eq->T(); // do this first to avoid self-destruct problems (see set_tuft_models() and set_oblique_models())
  if (dist2tosoma>=parameters.trunk_length_SQ) set_tuft_models(ts);
  else { // Is it time for an oblique branch?
    if (dist2tosoma>=parameters.dist2nextoblique) {
      set_next_oblique_distance(); // do this first to avoid self-destruct problems (see set_oblique_models())
      set_oblique_models(ts);
    }
  }
  return;
}

String BESTLNN_pyramidal_AD_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" pyramidal apical dendrite prototyping ");
  res += BESTL_nonnormalizing_terminal_segment_elongation_model::report_parameters_specific();
  res += " trunk length PDF: " + trunklengthpdf->report_parameters();
  res += " obliques PDF: " + obliquespdf->report_parameters();
  res += " prefix: " + parameters.prefix;
  return res;
}

elongation_rate_initialization_model_base::elongation_rate_initialization_model_base(elongation_rate_initialization_model_base * ericontrib, double & eriweight, elongation_rate_initialization_model_base & schema): contributing(NULL), contributingweight(1.0) { //, base_parameters(schema.base_parameters) {
  if (ericontrib) { // Clone chained models
    contributingweight = eriweight;
    contributing = ericontrib->clone();
  }
  //pdf = schema.pdf->clone();
}

elongation_rate_initialization_model_base::elongation_rate_initialization_model_base(String & thislabel, String & label, Command_Line_Parameters & clp): contributing(NULL), contributingweight(1.0) { //, base_parameters(0.0) {
  // This constructor is intended for use by schemas, which can set up model
  // chains, as well as look for general (thislabel=="") or subset parameter
  // settings that differ from the compile-time defaults.
  // (See TL#200603120618.5.)
  // [***NOTE] This is not a CLP_Modifiable, because I do not wish to have to
  // store the label in order to be able to parse_CLP at any time.
  //pdf = pdfselection(thislabel+"tsem.",clp,X_elongate);
  //if (!pdf) pdf = new normal_pdf(X_elongate,0.0,0.3,1.0); // a normal PDF is the default for TSEMs
  //else if (pdf->amplitude()>1.0) warning("Warning: Elongation Rate Initialization Model has perturbation with amplitude greater than 1.0!\n");
  if (label.empty()) return;
  int n;
  if ((n=clp.Specifies_Parameter(label+"eri_weight"))>=0) contributingweight = atof(clp.ParValue(n));
  elongation_rate_initialization_model_selection(label,clp,contributing,contributing);
}

void elongation_rate_initialization_model_base::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
  ts2.set_ERI_model(clone());
  ts1.set_ERI_model(this);
}

void elongation_rate_initialization_model_base::handle_turning(terminal_segment & ts) {
  ts.set_ERI_model(this);
}

#define SWAPDOUBLES(a,b) { \
double swaptmp = a; \
a = b; \
b = swaptmp; \
}

double elongation_rate_initialization_model_base::Get_Mean_and_Max_Quotas(terminal_segment * t, double & maxquota) {
  if (!t) {
    maxquota = 1.0;
    return 1.0;
  }
  double meanquota = 0.0;
  int quotacnt = 0;
  PLL_LOOP_FORWARD(terminal_segment,t,1) {
    double quota = e->ElongationModel()->base_parameters.l_i_cache;
    meanquota += quota;
    if (quota>maxquota) maxquota = quota;
    quotacnt++;
  }
  if (quotacnt<=1) return meanquota;
  return meanquota / ((double) quotacnt);
}

double elongation_rate_initialization_model_base::predict(double weight, double & predictedquota1, double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2) {
  // Called in a chain of ERI models to apply, where the first model
  // in the chain is automatically assigned the weight 1.0 and the weighting
  // of other direction models is proportional to that (i.e. it can be > 1.0).
  // The actual weighting of contributions is done in the function predict_initial_quota().
  // The sum of all weights is returned, so that division by that sum can produce
  // properly scaled combined branch angles.
  predict_initial_quota(predictedquota1,predictedquota2,ts1,ts2,weight);
  if (contributing) weight += contributing->predict(contributingweight,predictedquota1,predictedquota2,ts1,ts2);
  return weight;
}

void elongation_rate_initialization_model_base::initialize(terminal_segment * ts1, terminal_segment * ts2) {
  // PROTOCOL: This is called at bifurcation during growth.
  // This function calls:
  //   predict_initial_quota()
  //   predict()
  double predictedquota1 = 0.0, predictedquota2 = 0.0;
  predict_initial_quota(predictedquota1,predictedquota2,ts1,ts2);
  if (contributing) {
    double summedcontributingweights = 1.0 + contributing->predict(contributingweight,predictedquota1,predictedquota2,ts1,ts2);
    predictedquota1 /= summedcontributingweights; predictedquota2 /= summedcontributingweights;
  }
  // [***NOTE] Do we need to do anything about the initial elongation quotas stored in
  //           local l_i_cache variables of contributing models, or should those remain
  //           explicitly unaffected, in case they convey a local normalized contribution?
  ts1->ElongationModel()->base_parameters.l_i_cache = predictedquota1;
  ts2->ElongationModel()->base_parameters.l_i_cache = predictedquota2;
}

String elongation_rate_initialization_model_base::report_parameters() {
  // This calls report_parameters_specific() in any derived class and
  // propagates parameter reporting tocontributing models.
  String res(report_parameters_specific());
  //res += " p_TSEM=" + pdf->report_parameters();
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

length_distribution_eri_model::length_distribution_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, length_distribution_eri_model & schema): elongation_rate_initialization_model_base(ericontrib,eriweight,schema) {
  // Chaining of models is taken care of in the base constructor.
  pdf = schema.pdf->clone();
}

length_distribution_eri_model::length_distribution_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): elongation_rate_initialization_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
  //cout << "SETTING UP LENGTH DISTR. model\n";
  pdf = pdfselection(thislabel+"eri.",clp,X_elongate);
  //if (!pdf) cout << "NO SPECIFIC PDF FOUND\n";
  if (!pdf) pdf = new normal_pdf(X_elongate,0.0,1.0,3.0); // a normal PDF is the default for TSEMs
  //  else if (pdf->amplitude()>1.0) warning("Warning: Elongation Rate Initialization Model has perturbation with amplitude greater than 1.0!\n");
}

elongation_rate_initialization_model_base * length_distribution_eri_model::clone() {
  return new length_distribution_eri_model(contributing,contributingweight,*this);
}

void length_distribution_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  // Note that this model relates to considerations described in
  // TL#200709152200.1 and linked task log entries.
  // (See source code dated prior to 20080425 for the older implementation
  // that was changed in response to problems dealt with in TL#200804182330.1.)
  // We compute a (normal) random number. The standard mean value is the mean
  // of existing quotas. Negative values are compressed through the sigmoid
  // function 2*(1/1+e(-t)), so that a value of -inf equates to a zero quota
  // and a value of 0 equates to 1.0*meanquota of previously existing growth
  // cones. Positive values are converted as (0.5*t)+1, so that 0 is again
  // equal to 1.0*meanquota and 2 becomes 2*meanquota. Clearly now, it is
  // possible change both the possible relative slowdown or speedup by changing
  // the standard deviation of the normal function. The std 1.0 will be set as
  // the compiled default. Changing the mean value away from 0.0 can do more
  // drastic things by biasing new initial elongation rates to be slower or
  // faster than the elongation rates of existing growth cones.
  //cout << "ENTERED LENGTH-DISTR. quota prediction" << pdf->report_parameters() << '\n';
  double len1 = ts1->Length();
  double len2 = ts2->Length();
  double Xa = pdf->random_selection();
  double Xb = pdf->random_selection();
  if (Xa<0.0) Xa = 2.0/(1.0+exp(-Xa));
  else Xa = (0.5*Xa)+1.0;
  if (Xb<0.0) Xb = 2.0/(1.0+exp(-Xb));
  else Xb = (0.5*Xb)+1.0;
  // Initial quotas take into account the ratio len1:len2, or at least their order
  if (len1>len2) {
    if (Xa<Xb) SWAPDOUBLES(Xa,Xb);
  } else {
    if (Xa>Xb) SWAPDOUBLES(Xa,Xb);
  }
  // If the terminal segment elongation models of ts1 and ts2 do not have any
  // contributing models (contributing==NULL), and they are of a BESTL type,
  // then there is no need to compute a mean to compete with, and all will
  // be initialized around the same mean 1.0. Otherwise, a mean is computed
  // by inspecting all terminal segments. For simplicity, we always obtain
  // the current mean and maximum here.
  // [***NOTE] This depends on ts1->Next() being the old first terminal segment
  // in the list.
  double maxquota = 1.0;
  double meanquota = Get_Mean_and_Max_Quotas(ts1->Next(),maxquota);
  // What to do about unbounded probability distributions (amplitude == DBL_MAX - mean_X)?
  // Watch out for amplitude == 0
  // Suggest the use of a truncation value
  // Scaling
  //double maxamp = 2*(maxquota-meanquota);
  Xa *= meanquota;
  Xb *= meanquota;
  predictedquota1 += (weight*Xa);
  predictedquota2 += (weight*Xb);
#ifdef TESTQUOTAS
  quotas += String(eq->T(),"t=%.2f, ") + String(meanquota,"MQ=%.5f, ") + String(predictedquota1,"b1Q=%.5f, ") + String(predictedquota2,"b2Q=%.2f\n");
#endif
}

String length_distribution_eri_model::report_parameters_specific() {
  String res(" length distribution");
  res += pdf->report_parameters();
  return res;
}

nonnorm_BESTL_length_distribution_eri_model::nonnorm_BESTL_length_distribution_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, nonnorm_BESTL_length_distribution_eri_model & schema): length_distribution_eri_model(ericontrib,eriweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

nonnorm_BESTL_length_distribution_eri_model::nonnorm_BESTL_length_distribution_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): length_distribution_eri_model(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
}

elongation_rate_initialization_model_base * nonnorm_BESTL_length_distribution_eri_model::clone() {
  return new nonnorm_BESTL_length_distribution_eri_model(contributing,contributingweight,*this);
}

void nonnorm_BESTL_length_distribution_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  // This model applies the same correlation between initial branch lengths and
  // initial branch elongation rates as the length_distribution_eri_model, but
  // considers elongation rates absolute in accordance with notes in TL#200807020156.5.
  double len1 = ts1->Length();
  double len2 = ts2->Length();
  double Xa = pdf->random_positive();
  double Xb = pdf->random_positive();
  // Initial quotas take into account the ratio len1:len2, or at least their order
  if (len1>len2) {
    if (Xa<Xb) SWAPDOUBLES(Xa,Xb);
  } else {
    if (Xa>Xb) SWAPDOUBLES(Xa,Xb);
  }
  predictedquota1 += (weight*Xa);
  predictedquota2 += (weight*Xb);
}

void nonnorm_BESTL_length_distribution_eri_model::root_initialize(double & l_i_cache) {
  // [***INCOMPLETE] We may add "contributing" parts here, just as in the initialize() function, resulting in a weighted sum.
  l_i_cache = pdf->random_positive();
}

String nonnorm_BESTL_length_distribution_eri_model::report_parameters_specific() {
  String res(" non-normalizing BESTL length distribution");
  res += pdf->report_parameters();
  return res;
}

pure_stochastic_eri_model::pure_stochastic_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, pure_stochastic_eri_model & schema): elongation_rate_initialization_model_base(ericontrib,eriweight,schema) {
  // Chaining of models is taken care of in the base constructor.
  pdf = schema.pdf->clone();
}

pure_stochastic_eri_model::pure_stochastic_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): elongation_rate_initialization_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
  pdf = pdfselection(thislabel+"eri.",clp,X_elongate);
  if (!pdf) pdf = new normal_pdf(X_elongate,1.0,0.3,2.0); // a normal PDF is the default for TSEMs
  else if (pdf->amplitude()>1.0) warning("Warning: Elongation Rate Initialization Model has perturbation with amplitude greater than 1.0!\n");
}

elongation_rate_initialization_model_base * pure_stochastic_eri_model::clone() {
  return new pure_stochastic_eri_model(contributing,contributingweight,*this);
}

void pure_stochastic_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  // Note that this model relates to considerations described in
  // TL#200709152200.1 and linked task log entries.
  // (See source code dated prior to 20080425 for the older implementation
  // that was changed in response to problems dealt with in TL#200804182330.1.)
  // We compute a (normal) random number. The standard mean value is the mean
  // of existing quotas. Negative values are compressed through the sigmoid
  // function 2*(1/1+e(-t)), so that a value of -inf equates to a zero quota
  // and a value of 0 equates to 1.0*meanquota of previously existing growth
  // cones. Positive values are converted as (0.5*t)+1, so that 0 is again
  // equal to 1.0*meanquota and 2 becomes 2*meanquota. Clearly now, it is
  // possible change both the possible relative slowdown or speedup by changing
  // the standard deviation of the normal function. The std 1.0 will be set as
  // the compiled default. Changing the mean value away from 0.0 can do more
  // drastic things by biasing new initial elongation rates to be slower or
  // faster than the elongation rates of existing growth cones.
  double Xa = pdf->random_selection();
  double Xb = pdf->random_selection();
  if (Xa<0.0) Xa = 2.0/(1.0+exp(-Xa));
  else Xa = (0.5*Xa)+1.0;
  if (Xb<0.0) Xb = 2.0/(1.0+exp(-Xb));
  else Xb = (0.5*Xb)+1.0;
  // If the terminal segment elongation models of ts1 and ts2 do not have any
  // contributing models (contributing==NULL), and they are of a BESTL type,
  // then there is no need to compute a mean to compete with, and all will
  // be initialized around the same mean 1.0. Otherwise, a mean is computed
  // by inspecting all terminal segments. For simplicity, we always obtain
  // the current mean and maximum here.
  // [***NOTE] This depends on ts1->Next() being the old first terminal segment
  // in the list.
  double maxquota = 1.0;
  double meanquota = Get_Mean_and_Max_Quotas(ts1->Next(),maxquota);
  // What to do about unbounded probability distributions (amplitude == DBL_MAX - mean_X)?
  // Watch out for amplitude == 0
  // Suggest the use of a truncation value
  // Scaling
  //double maxamp = 2*(maxquota-meanquota);
  Xa *= meanquota;
  Xb *= meanquota;
  predictedquota1 += (weight*Xa);
  predictedquota2 += (weight*Xb);
}

String pure_stochastic_eri_model::report_parameters_specific() {
  String res(" pure stochastic");
  return res;
}

zero_eri_model::zero_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, zero_eri_model & schema): elongation_rate_initialization_model_base(ericontrib,eriweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

zero_eri_model::zero_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): elongation_rate_initialization_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
}

elongation_rate_initialization_model_base * zero_eri_model::clone() {
  return new zero_eri_model(contributing,contributingweight,*this);
}

void zero_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
//predictedquota1 += 0.0; This code is omitted, as it does not actually do anything.
//predictedquota2 += 0.0;
}

String zero_eri_model::report_parameters_specific() {
  String res(" zero");
  return res;
}

unitary_eri_model::unitary_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, unitary_eri_model & schema): elongation_rate_initialization_model_base(ericontrib,eriweight,schema), unitary_parameters(schema.unitary_parameters) {
  // Chaining of models is taken care of in the base constructor.
}

unitary_eri_model::unitary_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): elongation_rate_initialization_model_base(thislabel,label,clp), unitary_parameters(1.0) {
  // Chaining of models is taken care of in the base constructor
  int n;
  if ((n=clp.Specifies_Parameter(label+"eri_initialquota"))>=0) unitary_parameters.initialquota = atof(clp.ParValue(n));
}

elongation_rate_initialization_model_base * unitary_eri_model::clone() {
  return new unitary_eri_model(contributing,contributingweight,*this);
}

void unitary_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  predictedquota1 += unitary_parameters.initialquota;
  predictedquota2 += unitary_parameters.initialquota;
}

String unitary_eri_model::report_parameters_specific() {
  String res(" unitary");
  return res;
}

continue_defaults_eri_model::continue_defaults_eri_model(elongation_rate_initialization_model_base * ericontrib, double & eriweight, continue_defaults_eri_model & schema): elongation_rate_initialization_model_base(ericontrib,eriweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

continue_defaults_eri_model::continue_defaults_eri_model(String & thislabel, String & label, Command_Line_Parameters & clp): elongation_rate_initialization_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor
}

elongation_rate_initialization_model_base * continue_defaults_eri_model::clone() {
  return new continue_defaults_eri_model(contributing,contributingweight,*this);
}

void continue_defaults_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  predictedquota1 += ts1->ElongationModel()->base_parameters.l_i_cache;
  predictedquota2 += ts2->ElongationModel()->base_parameters.l_i_cache;
}

String continue_defaults_eri_model::report_parameters_specific() {
  String res(" continue defaults");
  return res;
}

void terminal_segment_elongation_event::event() {
  // [***NOTE] I can log this event (see IFspike): devres->elongation(ts,t); 
  ts->ElongationModel()->elongate(ts);
}

String arbor_elongation_model_idstr[NUM_aem] = {
  "van_pelt", // downcased!
  "polynomial_o1",
  "polynomial_o3"
};

arbor_elongation_model_base * arbor_elongation_model_selection(String & label, Command_Line_Parameters & clp, arbor_elongation_model_base * superior_set, arbor_elongation_model_base_ptr & aembptr) {
  // Note: Providing the result pointer in aembptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a arbor elongation model is not defined explicitly inherit the arbor elongation model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  arbor_elongation_model aem = NUM_aem; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"arbor_elongation_model"))>=0) {
    String aemstr(downcase(clp.ParValue(n)));
    for (arbor_elongation_model i = arbor_elongation_model(0); i < NUM_aem; i = arbor_elongation_model(i+1)) if (aemstr==arbor_elongation_model_idstr[i]) { aem = i; break; }
    if (aem==NUM_aem) return NULL; // catch arbor_elongation_model syntax errors
    if (label.empty()) general_arbor_elongation_model = aem; // modify default if setting universal
  } else if (label.empty()) aem = general_arbor_elongation_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"arbor_elongation_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) general_arbor_elongation_model_root = nextlabel; // if at most general level
  } else if ((n=clp.Specifies_Parameter(label+"aem_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) general_arbor_elongation_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (aembptr) aembptr->delete_shared();
  // Set and return specified model schema or general model schema if at most general level
  switch (aem) {
  case van_Pelt_aem:
    aembptr = new van_Pelt_arbor_elongation_model(label,nextlabel,clp);
    break;;
  case polynomial_O1_aem:
    aembptr = new Polynomial_O1_arbor_elongation_model(label,nextlabel,clp);
    break;;
  case polynomial_O3_aem:
    aembptr = new Polynomial_O3_arbor_elongation_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    aembptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through aembptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (arbor_elongation_model_region_subset_schemas[0][universal_sps]!=aembptr)) {
    if (arbor_elongation_model_region_subset_schemas[0][universal_sps]) arbor_elongation_model_region_subset_schemas[0][universal_sps]->delete_shared();
    arbor_elongation_model_region_subset_schemas[0][universal_sps] = aembptr->clone();
  }
  return aembptr;
}

String terminal_segment_elongation_model_idstr[NUM_tsem] = {
  "simple",
  "inertia",
  "second_order",
  "constrained_second_order",
  "initialized_cso",
  "decaying_second_order",
  "bestl",
  "nonnorm_bestl", // downcased!
  "pyrad_bestlnn" // downcased!
};

terminal_segment_elongation_model_base * terminal_segment_elongation_model_selection(String & label, Command_Line_Parameters & clp, terminal_segment_elongation_model_base * superior_set, terminal_segment_elongation_model_base_ptr & tsembptr) {
  // Note: Providing the result pointer in tsembptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which a terminal segment elongation model is not defined explicitly inherit the terminal segment elongation model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  terminal_segment_elongation_model tsem = NUM_tsem; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"terminal_segment_elongation_model"))>=0) {
    String tsemstr(downcase(clp.ParValue(n)));
    for (terminal_segment_elongation_model i = terminal_segment_elongation_model(0); i < NUM_tsem; i = terminal_segment_elongation_model(i+1)) if (tsemstr==terminal_segment_elongation_model_idstr[i]) { tsem = i; break; }
    if (tsem==NUM_tsem) return NULL; // catch terminal_segment_elongation_model syntax errors
    if (label.empty()) general_terminal_segment_elongation_model = tsem; // modify default if setting universal
  } else if (label.empty()) tsem = general_terminal_segment_elongation_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"terminal_segment_elongation_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) general_terminal_segment_elongation_model_root = nextlabel; // if at most general level
  } else if ((n=clp.Specifies_Parameter(label+"tsem_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) general_terminal_segment_elongation_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (tsembptr) tsembptr->delete_shared();
  // Set and return specified model schema or general model schema if at most general level
  switch (tsem) {
  case simple_tsem:
    tsembptr = new simple_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case inertia_tsem:
    tsembptr = new inertia_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case second_order_tsem:
    tsembptr = new second_order_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case constrained_second_order_tsem:
    tsembptr = new constrained_second_order_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case initialized_CSO_tsem:
    tsembptr = new initialized_CSO_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case decaying_second_order_tsem:
    tsembptr = new decaying_second_order_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case BESTL_tsem:
    tsembptr = new BESTL_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case nonnorm_BESTL_tsem:
    tsembptr = new BESTL_nonnormalizing_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  case pyrAD_BESTLNN_tsem:
    tsembptr = new BESTLNN_pyramidal_AD_terminal_segment_elongation_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    tsembptr = superior_set->clone(NULL);
    break;;
  }
  // In case a secondary schema reference is provided through aembptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (terminal_segment_elongation_model_region_subset_schemas[0][universal_sps]!=tsembptr)) {
    if (terminal_segment_elongation_model_region_subset_schemas[0][universal_sps]) terminal_segment_elongation_model_region_subset_schemas[0][universal_sps]->delete_shared();
    terminal_segment_elongation_model_region_subset_schemas[0][universal_sps] = tsembptr->clone(NULL);
  }
  return tsembptr;
}

String elongation_rate_initialization_model_idstr[NUM_eri] = {
  "length_distribution",
  "nonnorm_bestl_length_distribution",
  "pure_stochastic",
  "zero",
  "unitary",
  "continue_defaults"
};

elongation_rate_initialization_model_base * elongation_rate_initialization_model_selection(String & label, Command_Line_Parameters & clp, elongation_rate_initialization_model_base * superior_set, elongation_rate_initialization_model_base_ptr & eribptr) {
  // Note: Providing the result pointer in eribptr enables this function to
  // delete a previously defined schema before linking to a new one.
  // PROTOCOL: [See 200709130802.] Sets for which an elongation rate initialization model is not defined explicitly inherit the elongation rate initialization model of
  //           their specific direct superior by cloning it, thereby also inheriting parameters specified therein.
  elongation_rate_initialization_model eri = NUM_eri; // implies undefined
  String nextlabel;
  int n;
  // Determine if a specific model is selected for the set with the given label
  if ((n=clp.Specifies_Parameter(label+"elongation_rate_initialization_model"))>=0) {
    String eristr(downcase(clp.ParValue(n)));
    for (elongation_rate_initialization_model i = elongation_rate_initialization_model(0); i < NUM_eri; i = elongation_rate_initialization_model(i+1)) if (eristr==elongation_rate_initialization_model_idstr[i]) { eri = i; break; }
    if (eri==NUM_eri) return NULL; // catch elongation_rate_initialization_model syntax errors
    if (label.empty()) default_elongation_rate_initialization_model = eri; // modify default if setting universal
  } else if (label.empty()) eri = default_elongation_rate_initialization_model; // the universal set must have a model, since it cannot inherit
  // Determine if a following model selection should include chaining detection
  if ((n=clp.Specifies_Parameter(label+"elongation_rate_initialization_model_label"))>=0) {
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_elongation_rate_initialization_model_root = nextlabel; // if at most general level
  } else if ((n=clp.Specifies_Parameter(label+"eri_label"))>=0) { // alternate syntax
    nextlabel = label+clp.URI_unescape_ParValue(n)+'.';
    if (label.empty()) universal_elongation_rate_initialization_model_root = nextlabel; // to report at universal level
  }
  // Clean up any previously defined schema
  if (eribptr) eribptr->delete_shared();
  // Set and return specified model schema or universal model schema if at universal level
  switch (eri) {
  case length_distribution_eri:
    eribptr = new length_distribution_eri_model(label,nextlabel,clp);
    break;;
  case nonnorm_BESTL_length_distribution_eri:
    eribptr = new nonnorm_BESTL_length_distribution_eri_model(label,nextlabel,clp);
    break;;
  case pure_stochastic_eri:
    eribptr = new pure_stochastic_eri_model(label,nextlabel,clp);
    break;;
  case zero_eri:
    eribptr = new zero_eri_model(label,nextlabel,clp);
    break;;
  case unitary_eri:
    eribptr = new unitary_eri_model(label,nextlabel,clp);
    break;;
  case continue_defaults_eri:
    eribptr = new continue_defaults_eri_model(label,nextlabel,clp);
    break;;
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    eribptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through aembptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (elongation_rate_initialization_model_region_subset_schemas[0][universal_sps]!=eribptr)) {
    if (elongation_rate_initialization_model_region_subset_schemas[0][universal_sps]) elongation_rate_initialization_model_region_subset_schemas[0][universal_sps]->delete_shared();
    elongation_rate_initialization_model_region_subset_schemas[0][universal_sps] = eribptr->clone();
  }
  return eribptr;
}
