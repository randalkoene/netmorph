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
turning_model_base_ptr turning_model_subset_schemas[NUM_NATURAL_SPS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

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

turning_model_base * turning_model_selection(String & label, Command_Line_Parameters & clp, turning_model_base_ptr & tmbptr) {
  // Note: Providing the result pointer in tmbptr enables this function to
  // delete a previously defined schema before linking to a new one.
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
  if ((label.empty()) && (turning_model_subset_schemas[universal_sps]!=tmbptr)) {
    // This is in case a secondary schema reference is provided for the
    // universal set level for some reason, i.e. create a schema object for
    // turning_model_subset_schemas[universal_sps], as well as the recipient tmbptr.
    if (turning_model_subset_schemas[universal_sps]) turning_model_subset_schemas[universal_sps]->delete_shared();
    turning_model_subset_schemas[universal_sps] = tmbptr->clone();
  }
  return tmbptr;
}
