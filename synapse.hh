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
// synapse.hh
// Randal A. Koene, 20041118
//
// Classes the define synapse types.

#ifndef __SYNAPSE_HH
#define __SYNAPSE_HH

// Required class headers

#include "templates.hh"
#include "Fig_Object.hh"
#include "state_storable.hh"
#include "Txt_Object.hh"
#include "event.hh"

// Class pointer forward declarations (to avoid recursive inclusion)
#ifndef __SPATIAL_HH
class spatial;
#endif
#ifndef __CONNECTION_HH
class connection;
#endif
#ifndef __SYNAPSE_STRUCTURE_HH
class synapse_structure;
#endif
#ifndef __NEURON_HH
class neuron;
#endif

enum synapse_type { syntype_AMPAR, syntype_NMDAR, syntype_GABAR, syntype_candidate, syntype_GluR, syntype_iGluR, syntype_mGluR, syntype_IDs }; // To avoid base classes that are not used directly, syntype_candidate acts as the maximum type index, otherwise syntype_IDs indicates the number of valid IDs

extern const char synapse_type_name[syntype_IDs][10];

extern bool figattr_receptor[syntype_candidate];

#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
#define SYNGENESIS 0
#define SYNLOSS 1
extern long unsigned int synapse_inventory[syntype_IDs][2];
extern long unsigned int candidates_conversion_attempted;
#endif

// [***NOTE] The values defined here were picked without a physiological basis.
#define INITIAL_STRENGTH_AMPA 10.0
#define INITIAL_STRENGTH_NMDA  2.0
#define INITIAL_STRENGTH_GABA 10.0

class synapse: public PLLHandle<synapse>, public state_storable {
protected:
  synapse_type type_id;
  connection * c;
  synapse_structure * s; // defines spatial information about the synapse
  // (stages of the details may also consist of nodes at branch points)
  // presynaptic detail
  // postsynaptic detail
  synapse(connection & conn, synapse_type id, synapse_structure * synstruc = NULL); // Only derived classes may be instantiated.
public:
  ~synapse();
  synapse_type type_ID() { return type_id; }
  connection * Connection() { return c; }
  synapse_structure * Structure() { return s; }
  neuron * Presynaptic_Neuron();
  neuron * Postsynaptic_Neuron();
  //virtual double Strength() { return 0.0; }
  virtual double Peak_Response() { return 0.0; }
  virtual double Reversal_Potential() { return 0.0; }
  void move_add(spatial & addvec);
  void fanin_rot();
  Fig_Object * net_Fig();
  Txt_Object * net_Txt();
};

typedef synapse * synapseptr;

class synaptogenesis_data {
  // This data structure keeps track of synaptogenesis information, such as the time at which a specific synapse was created.
protected:
  double * data;
  synapseptr * synlist;
  unsigned int datalen;
  unsigned int dataidx;
public:
  synaptogenesis_data(): datalen(10240), dataidx(0) { data = new double[10240]; synlist = new synapseptr[10240]; }
  ~synaptogenesis_data() { delete[] synlist; delete[] data; }
  void add(synapse * s, double t);
  synapseptr get_synapse(unsigned int i) { if (i<dataidx) return synlist[i]; else return NULL; }
  double get_t_genesis(unsigned int i) { if (i<dataidx) return data[i]; else return -1.0; }
  double find_t_genesis(synapse * s);
};

extern synaptogenesis_data * SynaptoGenesis_Data;

class candidate_synapse: public synapse {
protected:
  double likelihood;
public:
  candidate_synapse(connection & conn, double l, synapse_structure * synstruc = NULL): synapse(conn,syntype_candidate,synstruc), likelihood(l) {}
  double Likelihood() { return likelihood; }
};

class glutamate_receptors: public synapse {
protected:
  double strength; // (probably) measured here in nano Siemens
public:
  glutamate_receptors(connection & conn, double initvalue = 0.0, synapse_structure * synstruc = NULL, synapse_type id = syntype_GluR): synapse(conn,id,synstruc), strength(initvalue) {
    if (SynaptoGenesis_Data) SynaptoGenesis_Data->add(this,eq->T());
  }
  virtual double Strength() { return strength; }
  virtual double Peak_Response() { return strength; } // See eq.3 in the STM-FIFO paper.
  //virtual double Reversal_Potential() { return 0.0; }
};

class ionotropic_glutamate_receptors: public glutamate_receptors {
public:
  ionotropic_glutamate_receptors(connection & conn, double initvalue = 0.0, synapse_structure * synstruc = NULL, synapse_type id = syntype_iGluR): glutamate_receptors(conn,initvalue,synstruc,id) {}
  //virtual double Reversal_Potential() { return 0.0; }
};

class metabotropic_glutamate_receptors: public glutamate_receptors {
public:
  metabotropic_glutamate_receptors(connection & conn, double initvalue = 0.0, synapse_structure * synstruc = NULL, synapse_type id = syntype_mGluR): glutamate_receptors(conn,initvalue,synstruc,id) {}
};

class AMPA_receptors: public ionotropic_glutamate_receptors {
public:
  AMPA_receptors(connection & conn, double initvalue = INITIAL_STRENGTH_AMPA, synapse_structure * synstruc = NULL, synapse_type id = syntype_AMPAR): ionotropic_glutamate_receptors(conn,initvalue,synstruc,id) {}
};

class NMDA_receptors: public ionotropic_glutamate_receptors {
public:
  NMDA_receptors(connection & conn, double initvalue = INITIAL_STRENGTH_NMDA, synapse_structure * synstruc = NULL, synapse_type id = syntype_NMDAR): ionotropic_glutamate_receptors(conn,initvalue,synstruc,id) {}
};

class GABA_receptors: public metabotropic_glutamate_receptors {
public:
  GABA_receptors(connection & conn, double initvalue = INITIAL_STRENGTH_GABA, synapse_structure * synstruc = NULL, synapse_type id = syntype_GABAR): metabotropic_glutamate_receptors(conn,initvalue,synstruc,id) {}
  virtual double Reversal_Potential() { return -65.0; }
};

#endif
