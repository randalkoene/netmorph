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
// neurite_diameter_model.cc
// Randal A. Koene, 20070709

#include <math.h>
#include "fibre_structure.hh"
#include "network.hh"
#include "neurite_diameter_model.hh"

// constants

const char neurite_diameter_model_label[NUM_ndm][12] = {
  "none",
  "asymptotic",
  "rall"
};

// variables and objects

normal_pdf defaultrallepowerpdf(X_branch,1.47,0.3);
normal_pdf defaultralldtermpdf(X_branch,0.7,0.3);

// classes

void asymptotic_soma_proportion_ndm::diameters() {
}

void asymptotic_soma_proportion_ndm::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("ndm.asymptotic_proportion"))>=0) proportion = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("ndm.asymptotic_lambda"))>=0) lambda = atof(clp.ParValue(n));
}

String asymptotic_soma_proportion_ndm::report_parameters() {
  String res("Asymptotic Soma Proportion neurite diameter model: proportion=");
  res += String(proportion,"%.3f");
  res += " lambda=";
  res += String(lambda,"%.3f");
  return res;
}

double rall_power_law_ndm::fibre_diameter(fibre_segment * fs) {
  // This function recursively reaches down to terminal segments, then retracts,
  // while computing new diameters at bifurcations and setting the diameters of
  // all fiber segments.
  if (fs->branch1) {
    if (fs->branch2) { // this is the parent of a bifurcation
      double d_1 = fibre_diameter(fs->branch1);
      double d_2 = fibre_diameter(fs->branch2);
      double e_power = PDF_e_power->random_positive();
      double d_epower = exp(e_power*log(d_1)) + exp(e_power*log(d_2));
      double d = exp(log(d_epower)/e_power);
      fs->set_diameter(d);
      return d;
    } else { // this is between bifurcations
      double d = fibre_diameter(fs->branch1);
      fs->set_diameter(d);
      return d;
    }
  } else if (fs->branch2) { // this is between bifurcations
      double d = fibre_diameter(fs->branch2);
      fs->set_diameter(d);
      return d;
  } else { // this is a terminal segment
    double d = PDF_d_term->random_positive();
    fs->set_diameter(d);
    return d;
  }
}

void rall_power_law_ndm::diameters() {
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) { // for each neuron in the network
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->InputStructure()->head(),1,fs) { // for dendrites
      fibre_diameter(fs);
    }
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,fs) { // for axons
      fibre_diameter(fs);
    }
  }
}

void rall_power_law_ndm::parse_CLP(Command_Line_Parameters & clp) {
  //  int n;
  PDF_e_power = pdfselection("ndm.e_power",clp,X_branch);
  if (!PDF_e_power) PDF_e_power = defaultrallepowerpdf.clone();
  PDF_d_term = pdfselection("ndm.d_term",clp,X_branch);
  if (!PDF_d_term) PDF_d_term = defaultralldtermpdf.clone();
}

String rall_power_law_ndm::report_parameters() {
  String res("Rall Power Law neurite diameter model: e_power ");
  if (PDF_e_power) res += PDF_e_power->report_parameters();
  res += " d_term ";
  if (PDF_d_term) res += PDF_d_term->report_parameters();
  res += "\n===============\n|| DEVELOPER NOTE: Neurite diameter models are not yet set-schema module compliant. Selection affects the entire network!\n===============\n";
  return res;
}

// functions

Neurite_Diameter_Model * select_neurite_diameter_model(network & net, Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("neurite_diameter_model"))>=0) {
    String downcaselabel(downcase(clp.ParValue(n)));
    for (neurite_diameter_model_id ndmi = no_ndm; ndmi < NUM_ndm; ndmi = (neurite_diameter_model_id) (ndmi + 1)) {
      if (downcaselabel.matches(neurite_diameter_model_label[ndmi])) {
	switch (ndmi) {
	case no_ndm: return NULL; break;;
	case asymptotic_ndm: return new asymptotic_soma_proportion_ndm(net,clp); break;;
	case rall_ndm: return new rall_power_law_ndm(net,clp); break;;
	case NUM_ndm: break;;
	  //default: ;
	}
      }
    }
    warning("Warning: Unknown neurite diameter model with label '"+downcaselabel+"', no model applied\n");
  }
  return NULL;
}
