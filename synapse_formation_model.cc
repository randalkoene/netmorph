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
// synapse_formation_model.cc
// Randal A. Koene, 20050107

#include <math.h>
#include "diagnostic.hh"
#include "synapse_formation_model.hh"
#include "synapse_structure.hh"
#include "synapse.hh"
#include "network.hh"

//#define DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS

synapse_formation_model defaultsfm(NULL,NULL);
synapse_formation_model * sfm = &defaultsfm;

general_synapse_formation_model_parameters_interface general_synapse_formation_model_parameters;

// An initial probability distribution of synapse types:
// Each value indicates the likelihood that a synaptic connection between
// the pre- and postsynaptic neurons includes a specific receptor type.
// A value less than 0 indicates that no such receptor can exist at that
// synaptic connection.
double possible_syntypes[UNTYPED_NEURON+1][UNTYPED_NEURON+1][syntype_candidate] = {
  { // pre==PRINCIPAL_NEURON
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.4, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==INTERNEURON
    { // post==PRINCIPAL_NEURON
      -9.9, -9.9, 1.0
    },
    { // post==INTERNEURON
      -9.9, -9.9, 1.0
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      -9.9, -9.9, 1.0
    },
    { // post==PYRAMIDAL
      -9.9, -9.9, 1.0
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==MULTIPOLAR_NONPYRAMIDAL
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.4, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==PYRAMIDAL
    { // post==PRINCIPAL_NEURON
      1.0, -9.9, -9.9
    },
    { // post==INTERNEURON
      1.0, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      1.0, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      1.0, 0.9, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  },
  { // pre==UNTYPED_NEURON
    { // post==PRINCIPAL_NEURON
      -9.9, -9.9, -9.9
    },
    { // post==INTERNEURON
      -9.9, -9.9, -9.9
    },
    { // post==MULTIPOLAR_NONPYRAMIDAL
      -9.9, -9.9, -9.9
    },
    { // post==PYRAMIDAL
      -9.9, -9.9, -9.9
    },
    { // post==UNTYPED_NEURON
      -9.9, -9.9, -9.9
    }
  }
};

// The following array provides pre- and postsynaptic neuron type specific
// information about the maximum distance between fibres at which a synapse
// can be established. Typically, a combination that commonly exhibits
// dendritic spikes can bridge a distance of up to 1.0 microns, while those
// without spines can bridge a distance of up to 0.1 microns.
double maxfibredistance[UNTYPED_NEURON+1][UNTYPED_NEURON+1] = {
  { // pre==PRINCIPAL_NEURON
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==INTERNEURON
      0.1, 0.1, 0.1, 0.1, 0.0
  },
  { // pre==MULTIPOLAR_NONPYRAMIDAL
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==PYRAMIDAL
      1.0, 0.1, 1.0, 1.0, 0.0
  },
  { // pre==UNTYPED_NEURON
      0.0, 0.0, 0.0, 0.0, 0.0
  }
};

uniform_pdf defaultsfmpdf(X_synapses);
probability_distribution_function * sfmpdf = &defaultsfmpdf;

void general_synapse_formation_model_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  // [***NOTE] This does not yet include a facility for setting these
  // parameters differently for specific synapse_formation_model objects.
  int n;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1))
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1))
      for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1))
	if ((n=clp.Specifies_Parameter("P_receptor."+String(neuron_short_name[n_pre])+'.'+String(neuron_short_name[n_post])+'.'+String(synapse_type_name[receptortype])))>=0) possible_syntypes[n_pre][n_post][receptortype] = atof(clp.ParValue(n));
  sfmpdf = pdfselection("synapse_formation",clp,X_synapses);
  if (!sfmpdf) sfmpdf = &defaultsfmpdf;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1))
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1))
      if ((n=clp.Specifies_Parameter("D_synmax."+String(neuron_short_name[n_pre])+'.'+String(neuron_short_name[n_post])))>=0) maxfibredistance[n_pre][n_post] = atof(clp.ParValue(n));
}

String general_synapse_formation_model_parameters_interface::report_parameters() {
  String res("Synapse formation model hypotheses: ");
  res += "\n  Receptor type probabilities (P_receptor):\n                  POST\n";
  String setofreceptorsstr;
  for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1)) {
    String spacerstr("       ");
    for (int i = 0; ((i<6) && (synapse_type_name[receptortype][i]!='\0')); i++) spacerstr[i] = synapse_type_name[receptortype][i];
    setofreceptorsstr += spacerstr;
  }
  String setofneuronsstr;
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) {
#ifdef BIGSTRING_REPLICATE_PROBLEM
    int spacebufsize = (syntype_candidate*7)+2;
    char * spacebuf = new char[spacebufsize+1];
    for (int i = 0; i<spacebufsize; i++) spacebuf[i] = ' ';
    spacebuf[spacebufsize] = '\0';
    String spacerstr(spacebuf);
    delete[] spacebuf;
#else
    String spacerstr(replicate(' ',(syntype_candidate*7)+2));
#endif
    for (unsigned int i = 0; ((i<spacerstr.length()) && (neuron_short_name[n_pre][i]!='\0')); i++) spacerstr[i] = neuron_short_name[n_pre][i];
    setofneuronsstr += spacerstr;
  }
  res += "                  "+setofneuronsstr+'\n';
  res += "    PRE         ";
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) res += "  "+setofreceptorsstr;
  res += '\n';
  for (neuron_type n_pre = (neuron_type) 0; n_pre <= UNTYPED_NEURON; n_pre = (neuron_type) (n_pre + 1)) {
    String spacerstr("                 ");
    for (int i = 0; ((i<12) && (neuron_short_name[n_pre][i]!='\0')); i++) spacerstr[i+4] = neuron_short_name[n_pre][i];
    res += spacerstr;
    for (neuron_type n_post = (neuron_type) 0; n_post <= UNTYPED_NEURON; n_post = (neuron_type) (n_post + 1)) {
      for (synapse_type receptortype = (synapse_type) 0; receptortype < syntype_candidate; receptortype = (synapse_type) (receptortype + 1)) res += String(possible_syntypes[n_pre][n_post][receptortype],"%6.2f ");
      if (n_post<UNTYPED_NEURON) res += "  ";
      else res += '\n';
    }
  }
  res += "  " + sfmpdf->report_parameters() + '\n';
  return res;
}

double synapse_formation_model::evaluate_possible_connection(fibre_segment * axonsegment, fibre_segment * dendritesegment) {
  // A model of synapse formation is applied to determine the likelihood that
  // a synapse forms between the axon segment and the dedrite segment.
  // In the base model, a non-zero likelihood leads to the creation of a
  // candidate synapse, i.e. a synapse with a development state set to the
  // virtual/quantum theoretical value of "candidate".
  // Returns the likelihood, which is also stored as a parameter of the
  // candidate synapse if one is created.
  // The calling function can decide to test connprocid in the other fibre
  // segment if it wishes to avoid evaluating the same pair of segments
  // more than once. Such multiple processing may otherwise occur in
  // specific cases, such as if segments that are allocated to multiple
  // spatial partitions appear in the same or neighboring paritions more
  // than once. In the current implementation of the growth model, this
  // does not otherwise occur, since each axon or dendrite segment is
  // evaluated only once for its possible connections with fibre segments
  // that developed earlier.
  // *** If multiple synapses between the same pair of segments are a
  // possibility then those can be created here or during a call to
  // establish_connections().
  DIAGNOSTIC_EVALUATE_POSSIBLE_CONNECTION_PROFILING;
  axonsegment->connprocid = dendritesegment;
  dendritesegment->connprocid = axonsegment;
  // Determine the likelihood
  if (axonsegment->N()==dendritesegment->N()) return 0.0;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_synapse_structure);
  synapse_structure * s = new synapse_structure(*axonsegment,*dendritesegment); // ***WHY DO I CREATE THIS EVEN IF IT IS SUBSEQUENTLY NOT USED?
  DIAGNOSTIC_AFTER_ALLOCATION(new_synapse_structure,1,sizeof(*s));
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
  double likelihood = 1.0 - (dist3D_Segment_to_Segment(*axonsegment,*dendritesegment,s)/proximitythreshold);
#else
  double likelihood = 1.0 - (dist3D_Segment_to_Segment(*axonsegment,*dendritesegment,s)/maxfibredistance[axonsegment->N()->TypeID()][dendritesegment->N()->TypeID()]);
#endif
  ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) {
  ///  cout << "\nSYNEVAL: T=" << eq->T() << " axonsegment=" << (unsigned long) axonsegment << " dendritesegment=" << (unsigned long) dendritesegment;
  ///}
  if (likelihood>0.0) {
    ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) cout << " CANDIDATE\n"; cout.flush();
    // Candidate synapses are created in the fibre structure, but not listed
    // in connections between neurons until the actual synapses are
    // established.
    connectionptr c = axonsegment->N()->connect_to(dendritesegment->N());
    new candidate_synapse(*c,likelihood,s);
    // ***synapse_inventory[syntype_iGluR][1]++;
    return likelihood;
  }
  ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) cout << '\n'; cout.flush();
  // plug a potential memory leak!
  DIAGNOSTIC_BEFORE_ALLOCATION(new_synapse_structure);
  delete s;
  DIAGNOSTIC_AFTER_ALLOCATION(new_synapse_structure,-1,sizeof(synapse_structure));
  return likelihood;
}

void synapse_formation_model::establish_connections() {
  // Selects which candidate synapses become actual synapses.
  // The type of actual neuron created depends on the location of the
  // synapse and on the types of neurons involved.
  if (!net) {
    cout << "Warning: net==NULL in synapse_formation_model::establish_connections()\n";
    return;
  }
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    // ***synapse_inventory[syntype_GluR][1]++;
    PLL_LOOP_FORWARD_NESTED(connection,e->OutputConnections()->head(),1,c) {
      // ***synapse_inventory[syntype_iGluR][1]++;
      synapse * s_next;
      for (synapse * s = c->Synapses()->head(); (s); s = s_next) {
	s_next = s->Next(); // cached, since s may be removed
	// ***synapse_inventory[syntype_mGluR][1]++;
	if (s->type_ID()==syntype_candidate) {
	  // ***synapse_inventory[syntype_mGluR][0]++;
	  candidate_synapse * cs = static_cast<candidate_synapse *>(s);
	  // For this synapse formation model, the rules are outlined in
	  // Task Log entry TL#200508060747.1.
	  double likelihood = cs->Likelihood(); // == ratio of proximitythreshold microns
	  neuron_type pretype = cs->Presynaptic_Neuron()->TypeID();
	  neuron_type posttype = cs->Postsynaptic_Neuron()->TypeID();
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
	  double distance = likelihood * proximitythreshold;
	  likelihood = 1.0 - (distance/maxfibredistance[pretype][posttype]); // can be less than 0
	  if (likelihood>0.0) {
#endif
	    // ***synapse_inventory[syntype_GluR][0]++;
	    // Do the following for each possible synapse type (fully defined
	    // derived classes). This creates synapses (receptors) with initvalue
	    // ratio of initial amplitude parameter.
	    // [*** NOTE] The following method is based on an unsupported
	    // hypothesis with many assumptions that ignore the effects of
	    // early network activity. The method is described in
	    // TL#200508080920.5.
	    if (sfmpdf->random_positive()<=likelihood) { // a synapse is formed
	      candidates_conversion_attempted++;
	      double * initvalueptr = possible_syntypes[pretype][posttype];
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_AMPAR]) new AMPA_receptors(*c,likelihood*INITIAL_STRENGTH_AMPA,cs->Structure());
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_NMDAR]) new NMDA_receptors(*c,likelihood*INITIAL_STRENGTH_NMDA,cs->Structure());
	      if (X_synapses.get_rand_real1()<=initvalueptr[syntype_GABAR]) new GABA_receptors(*c,likelihood*INITIAL_STRENGTH_GABA,cs->Structure());
	    }
#ifdef DELAYED_SPECIFIC_PROXIMITY_THRESHOLDS
	  }
#endif
	  // ***synapse_inventory[syntype_mGluR][0]++;
	  cs->remove(); // clean up candidate synapses
	}
      }
    }
  }
}
