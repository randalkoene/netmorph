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
// neurite_diameter_model.hh
// Randal A. Koene, 20070709

/* Documentation:
   This implements base and derived classes for models of neurite diameter
   determination.

   At this time, these classes are applied equally to all neurons in a
   network and therefore do not participate in the NETMORPH interface
   protocol for model and parameter selection according to the smallest
   matching set.

   See task log entries around 20070627.
 */

#ifndef __NEURITE_DIAMETER_MODEL_HH
#define __NEURITE_DIAMETER_MODEL_HH

#include "templates.hh"
//#include "fibre_structure.hh"
//#include "spatial.hh"
#include "Command_Line_Parameters.hh"
//#include "global.hh"

#ifndef __NETWORK_HH
class network;
#endif
#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
#endif

enum neurite_diameter_model_id { no_ndm, asymptotic_ndm, rall_ndm, NUM_ndm };

extern const char neurite_diameter_model_label[NUM_ndm][12];

class Neurite_Diameter_Model {
protected:
  network * net;
public:
  Neurite_Diameter_Model(network & _net): net(&_net) {}
  virtual ~Neurite_Diameter_Model() {}
  virtual void diameters() = 0;
  virtual void parse_CLP(Command_Line_Parameters & clp) = 0;
  virtual String report_parameters() = 0;
};

class asymptotic_soma_proportion_ndm: public Neurite_Diameter_Model {
  // This model is based on the simple premise that the maximum neurite diameter
  // is some proportion of the soma diameter. According to a lambda parameter,
  // the length of the arbor determines how closely that asymptotic proportion
  // is approached exponentially.
protected:
  double proportion;
  double lambda;
public:
  asymptotic_soma_proportion_ndm(network & _net): Neurite_Diameter_Model(_net), proportion(0.3), lambda(50.0) {}
  asymptotic_soma_proportion_ndm(network & _net, Command_Line_Parameters & clp): Neurite_Diameter_Model(_net), proportion(0.3), lambda(50.0) { parse_CLP(clp); }
  //virtual ~asymptotic_soma_proportion_ndm() {}
  virtual void diameters();
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class rall_power_law_ndm: public Neurite_Diameter_Model {
  // This model implements the Rall model of dendritic segment diameters, based
  // on the equation d_p^e = d_1^e + d_2^e, where d_1 and d_2 are the daughter
  // branch diameters of a parent with diameter d_p, and where e is 3/2.
  // Based on observations by Larkman (1991) of terminal segments with a narrow
  // range of diameters, we assume that all terminal segments have an equal
  // diameter d_term, which is set by default to 1.25 micrometers.
  // Also see [PELT:VARIABILITY].
protected:
  probability_distribution_function * PDF_e_power;
  probability_distribution_function * PDF_d_term;
  double fibre_diameter(fibre_segment * fs);
public:
  rall_power_law_ndm(network & _net): Neurite_Diameter_Model(_net), PDF_e_power(NULL), PDF_d_term(NULL) {
    // cloning should take place here from schema
    error("The default constructor of rall_power_law_ndm is not yet a usable set-modular component!\n");
  }
  rall_power_law_ndm(network & _net, Command_Line_Parameters & clp): Neurite_Diameter_Model(_net), PDF_e_power(NULL), PDF_d_term(NULL) { parse_CLP(clp); }
  virtual ~rall_power_law_ndm() { if (PDF_e_power) delete PDF_e_power; if (PDF_d_term) delete PDF_d_term; }
  virtual void diameters();
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

Neurite_Diameter_Model * select_neurite_diameter_model(network & net, Command_Line_Parameters & clp);

#endif
