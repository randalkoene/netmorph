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
// synapse_formation_model.hh
// Randal A. Koene, 20050107
//
// Classes that contain synapse formation models and functions that apply
// those models.

#ifndef __SYNAPSE_FORMATION_MODEL_HH
#define __SYNAPSE_FORMATION_MODEL_HH

#include "fibre_structure.hh"
#include "Connection_Statistics.hh"
#include "neuron.hh"
#include "spatial.hh"
#include "Command_Line_Parameters.hh"

#ifndef __NETWORK_HH
class network;
#endif

/* [***INCOMPLETE] This set of classes needs to be cast into the standardized NETMOPRH protocol:
   (1) Separate the models into models that determine candidate synaptic sites and models
       that place actual synapses at or near a subset of the candidate synaptic sites.
   (2) Provide model and parameter selection functions, enable chaining.
 */

class general_synapse_formation_model_parameters_interface: public CLP_Modifiable {
  // This class serves to create the general_synapse_formation_model_parameters
  // object for parameters that should apply to all synapses in generated
  // networks. A global reference object is initialized in nibr.cc.
public:
  general_synapse_formation_model_parameters_interface() {}
  virtual ~general_synapse_formation_model_parameters_interface() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

extern general_synapse_formation_model_parameters_interface general_synapse_formation_model_parameters;

// The base model of synapse formation establishes candidate synapses during
// the growth of fibre structure and then establishes actual synapses, a
// subset of the candidates that takes the likelihood of synapse formation
// into account and seeks to achieve the connection statistics specified.
//
// Two functions of classes derived from synapse_formation_model are
// essential for the generation of synaptic connections:
// 1. evaluate_possible_connection(): This uses any appropriate measure to
//    determine if a candidate synapse is in order (e.g. distance).
// 2. establish_connections(): This is a response to any relevant conditions
//    (e.g. electrical activity, distance, pre- and postsynaptic neuron types,
//    age of neurons or of their connection, etc.) that determines how many
//    receptors of which type exist at a connection. 
//
// An alternative model of synapse formation that establishes actual synapses
// sooner and thereby enables early termination of network development is
// described as a consideration for future work (see
// <A HREF="../../../doc/html/minduploading.org/rak/nibr/nibr-future-work.html#consideration-connection-evaluation">Considerations of the connection evaluation</A>).

class synapse_formation_model {
protected:
  network * net;
  Connection_Statistics_Root * cstats;
  double proximitythreshold;
  // [***INCOMPLETE] Currently in the general parameters below: probability_distribution_function * pdf; // PDF used for selection of actual synapses at candidate sites
public:
  synapse_formation_model(network * n, Connection_Statistics_Root * cr): net(n), cstats(cr), proximitythreshold(0.0) {}
  virtual ~synapse_formation_model() {}
  virtual double evaluate_possible_connection(fibre_segment * axonsegment, fibre_segment * dendritesegment);
  virtual void establish_connections();
  double Proximity_Threshold() { return proximitythreshold; }
  void set_proximity_threshold(double pthres) { proximitythreshold = pthres; }
};

extern synapse_formation_model * sfm; // *** this could be passed as a parameter argument instead of used as a global parameter

extern double possible_syntypes[UNTYPED_NEURON+1][UNTYPED_NEURON+1][syntype_candidate];

extern double maxfibredistance[UNTYPED_NEURON+1][UNTYPED_NEURON+1];

extern probability_distribution_function * sfmpdf;

extern bool no_autapses;

/*
Selecting actual synapses (and their receptors) for each connection:
(1) Initially, connections contains one candidate synapse each.
(2) A synapse formation model, such as the one accessed through sfm knows
    the network to which it applies. Therefore, establish_connections() can
    handle the connections of each neuron in turn.
(3) The function possible_syntypes can be used to obtain an initial probable
    distribution of specific synaptic receptor types for a given pair of
    pre- and postsynaptic neurons.
(4) A synapse is created for each type of receptor at a connection.
(5) The distribution ratios are used as initial parameters for the individual
    synapses and affect the values computed by their transmission functions.
(6) If a connection ends up having no synaptic receptors then no synapses are
    created. This can happen for instance if connections between a specific
    pair of pre- and postsynaptic neurons are not possible.
(7) The candidate synapse is removed.
*/

#endif
