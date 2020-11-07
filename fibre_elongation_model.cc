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
#include "neuron.hh"
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
arbor_elongation_model_base_ptr general_arbor_elongation_model_base_subset_schemas[NUM_NATURAL_SPS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
terminal_segment_elongation_model_base_ptr general_terminal_segment_elongation_model_base_subset_schemas[NUM_NATURAL_SPS] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
elongation_rate_initialization_model default_elongation_rate_initialization_model = length_distribution_eri; // identifier, changed when a different universal model is set
String universal_elongation_rate_initialization_model_root; // Used only to report if chaining at the universal set level.
elongation_rate_initialization_model_base_ptr elongation_rate_initialization_model_subset_schemas[NUM_NATURAL_SPS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

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
  if (tsemcontrib) { // Clone chained models
    contributingweight = tsemweight;
    contributing = tsemcontrib->clone();
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
  else if (pdf->amplitude()>1.0) warning("Warning: Terminal Segment Elongation Model has perturbation with amplitude greater than 1.0, retraction may occur!\n");
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
  ts2.set_elongation_model(clone());
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

terminal_segment_elongation_model_base * simple_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * inertia_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * second_order_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * constrained_second_order_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * initialized_CSO_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * decaying_second_order_terminal_segment_elongation_model::clone() {
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

terminal_segment_elongation_model_base * delayed_branch_terminal_segment_elongation_model::clone() {
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

BESTL_terminal_segment_elongation_model::BESTL_terminal_segment_elongation_model(terminal_segment_elongation_model_base * tsemcontrib, double & tsemweight, BESTL_terminal_segment_elongation_model & schema): terminal_segment_elongation_model_base(tsemcontrib,tsemweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

BESTL_terminal_segment_elongation_model::BESTL_terminal_segment_elongation_model(String & thislabel, String & label, Command_Line_Parameters & clp): terminal_segment_elongation_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
  base_parameters.l_i_cache = 1.0; // this is the assumed mean of the BESTL quota distribution
  // Note: The mean can shift if other models contribute to the elongation quota determination.
}

terminal_segment_elongation_model_base * BESTL_terminal_segment_elongation_model::clone() {
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

String BESTL_terminal_segment_elongation_model::report_parameters_specific() {
  String res(" BESTL");
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
  pdf = pdfselection(thislabel+"eri.",clp,X_elongate);
  if (!pdf) pdf = new normal_pdf(X_elongate,1.0,0.3,2.0); // a normal PDF is the default for TSEMs
  else if (pdf->amplitude()>1.0) warning("Warning: Elongation Rate Initialization Model has perturbation with amplitude greater than 1.0!\n");
}

elongation_rate_initialization_model_base * length_distribution_eri_model::clone() {
  return new length_distribution_eri_model(contributing,contributingweight,*this);
}

void length_distribution_eri_model::predict_initial_quota(double & predictedquota1,double & predictedquota2, terminal_segment * ts1, terminal_segment * ts2,double weight) {
  double len1 = ts1->Length();
  double len2 = ts2->Length();
  // A distribution about 1.0 is used by default.
  double Xa = pdf->random_positive();
  double Xb = pdf->random_positive();
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
  // If necessary, initial quotas are scaled, so that 1.0 maps to the mean and
  // 2.0 is twice the maximum. Notice that specifying different mean, std and
  // truncation parameters for the PDF can shift the way in which new branch
  // elongation quotas compare to the quotas of the existing tree.
  // The difference between 2.0 and 1.0 is 1.0, which must be mapped to max-mean
  double maxamp = 2*(maxquota-meanquota);
  Xa = ((Xa - 1.0) * maxamp) + meanquota;
  Xb = ((Xb - 1.0) * maxamp) + meanquota;
  predictedquota1 += (weight*Xa);
  predictedquota2 += (weight*Xb);
}

String length_distribution_eri_model::report_parameters_specific() {
  String res(" length distribution");
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
  // A distribution about 1.0 is used by default.
  double Xa = pdf->random_positive();
  double Xb = pdf->random_positive();
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
  // If necessary, initial quotas are scaled, so that 1.0 maps to the mean and
  // 2.0 is twice the maximum. Notice that specifying different mean, std and
  // truncation parameters for the PDF can shift the way in which new branch
  // elongation quotas compare to the quotas of the existing tree.
  // The difference between 2.0 and 1.0 is 1.0, which must be mapped to max-mean
  double maxamp = 2*(maxquota-meanquota);
  Xa = ((Xa - 1.0) * maxamp) + meanquota;
  Xb = ((Xb - 1.0) * maxamp) + meanquota;
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
  if ((label.empty()) && (general_arbor_elongation_model_base_subset_schemas[universal_sps]!=aembptr)) {
    if (general_arbor_elongation_model_base_subset_schemas[universal_sps]) general_arbor_elongation_model_base_subset_schemas[universal_sps]->delete_shared();
    general_arbor_elongation_model_base_subset_schemas[universal_sps] = aembptr->clone();
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
  "bestl" // downcased!
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
  default: // If not specified, inherit from immediate superior set
    if (!superior_set) return NULL;
    tsembptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through aembptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (general_terminal_segment_elongation_model_base_subset_schemas[universal_sps]!=tsembptr)) {
    if (general_terminal_segment_elongation_model_base_subset_schemas[universal_sps]) general_terminal_segment_elongation_model_base_subset_schemas[universal_sps]->delete_shared();
    general_terminal_segment_elongation_model_base_subset_schemas[universal_sps] = tsembptr->clone();
  }
  return tsembptr;
}

String elongation_rate_initialization_model_idstr[NUM_eri] = {
  "length_distribution",
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
  if ((label.empty()) && (elongation_rate_initialization_model_subset_schemas[universal_sps]!=eribptr)) {
    if (elongation_rate_initialization_model_subset_schemas[universal_sps]) elongation_rate_initialization_model_subset_schemas[universal_sps]->delete_shared();
    elongation_rate_initialization_model_subset_schemas[universal_sps] = eribptr->clone();
  }
  return eribptr;
}
