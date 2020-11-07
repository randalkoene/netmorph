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
branching_model_base_ptr branching_model_subset_schemas[NUM_NATURAL_SPS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

branch_angle_model default_branch_angle_model = balanced_forces_bam; // identifier, changed when a different universal model is set
String universal_branch_angle_model_root; // Used only to report if chaining at the universal set level.
branch_angle_model_base_ptr branch_angle_model_subset_schemas[NUM_NATURAL_SPS] = {
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL,
  NULL
};

branching_model_base::branching_model_base(branching_model_base * bmcontrib, double & bmweight, branching_model_base & schema):
  contributing(NULL), contributingweight(1.0), base_parameters(schema.base_parameters) {
  // [***INCOMPLETE] Can I suggest inlining this constructor somehow?
  if (bmcontrib) { // Clone chained models
    contributingweight = bmweight;
    contributing = bmcontrib->clone();
  }
}

branching_model_base::branching_model_base(String & thislabel, String & label, Command_Line_Parameters & clp):
  contributing(NULL), contributingweight(1.0), base_parameters(0.0) {
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
  String res(report_parameters_specific());
  /*res += String(base_parameters.veeranglemin,"\n  [veeranglemin=%.2f] ");
    res += String(base_parameters.veeranglemax,"\n  [veeranglemax=%.2f] "); */
  if (contributing) {
    res += String(contributingweight,"\n  [contrib.w=%.2f] ");
    res += contributing->report_parameters();
  }
  return res;
}

van_Pelt_branching_model::van_Pelt_branching_model(branching_model_base * bmbcontrib, double & bmbweight, van_Pelt_branching_model & schema): branching_model_base(bmbcontrib,bmbweight,schema) {
  // Chaining of models is taken care of in the base constructor.
}

van_Pelt_branching_model::van_Pelt_branching_model(String & thislabel, String & label, Command_Line_Parameters & clp): branching_model_base(thislabel,label,clp) {
  // Chaining of models is taken care of in the base constructor.
}

branching_model_base * van_Pelt_branching_model::clone() {
  return new van_Pelt_branching_model(contributing,contributingweight,*this);
}

//void van_Pelt_branching_model::handle_branching(terminal_segment & ts1, terminal_segment & ts2) {
// [***NOTE] Could set a new temp.variable in here that depends on N, so that its calculation is only done when needed.
//}

//void van_Pelt_branching_model::handle_turning(terminal_segment & ts) {
//}

double van_Pelt_branching_model::predict_branching(double weight) {
  double probability = 0.0;
  // [***INCOMPLETE] This is a stub. Move the actual branching model code here from dendritic_growth_model.cc.
  probability *= weight;
  return probability;
}

String van_Pelt_branching_model::report_parameters_specific() {
  String res(" van_Pelt");
  return res;
}

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
  if (!pdf) pdf = new normal_pdf(X_elongate,0.0,0.3,1.0); // a normal PDF is the default for BAMs
  else if (pdf->amplitude()>1.0) warning("Warning: Branch Angle Model has perturbation with amplitude greater than 1.0, unpredictable angles!\n");
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
}

String balanced_forces_branch_angle_model::report_parameters_specific() {
  String res(" balanced forces");
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
  const double M_PI_HALF = M_PI/2.0;
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
    // angle1 is M_PI_HALF, angle2 is zero
    predicted1 += (weight*M_PI_HALF);
    return;
  } else if (force2==0.0) {
    // angle1 is zero, angle2 is M_PI_HALF
    predicted2 += (weight*M_PI_HALF);
    return;
  }
  double angle1 = atan(force2/force1);
  double angle2 = M_PI_HALF - angle1;
#ifdef TEST_ANGLES
  if (angle1<0.0) cout << "Angle1 is smaller than zero!\n";
  if (angle2>M_PI_HALF) cout << "Angle2 is greater than M_PI_HALF!\n";
#endif
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
    if (!superior_set) return NULL;
    bmbptr = superior_set->clone();
    break;;
  }
  // In case a secondary schema reference is provided through bmbptr for the universal set level, create a clone for the universal set
  if ((label.empty()) && (branching_model_subset_schemas[universal_sps]!=bmbptr)) {
    if (branching_model_subset_schemas[universal_sps]) branching_model_subset_schemas[universal_sps]->delete_shared();
    branching_model_subset_schemas[universal_sps] = bmbptr->clone();
  }
  return bmbptr;
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
  if ((label.empty()) && (branch_angle_model_subset_schemas[universal_sps]!=bambptr)) {
    if (branch_angle_model_subset_schemas[universal_sps]) branch_angle_model_subset_schemas[universal_sps]->delete_shared();
    branch_angle_model_subset_schemas[universal_sps] = bambptr->clone();
  }
  return bambptr;
}
