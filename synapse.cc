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
// synapse.cc
// Randal A. Koene, 20041118

#include "synapse.hh"
#include "connection.hh"
#include "synapse_structure.hh"
#include "Color_Table.hh"
#include "global.hh"

synaptogenesis_data * SynaptoGenesis_Data = NULL;

bool figattr_receptor[syntype_candidate] = {
  true,
  true,
  true,
};

const char synapse_type_name[syntype_IDs][10] = {
  "AMPAR",
  "NMDAR",
  "GABAR",
  "candidate",
  "GluR",
  "iGluR",
  "mGluR",
};

#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
long unsigned int synapse_inventory[syntype_IDs][2] = {
  { 0, 0},
  { 0, 0},
  { 0, 0},
  { 0, 0},
  { 0, 0},
  { 0, 0},
  { 0, 0}
};
long unsigned int candidates_conversion_attempted = 0;
#endif

synapse::synapse(connection & conn, synapse_type id, synapse_structure * synstruc): type_id(id), c(&conn), s(synstruc) {
  // Here the synapse is added to a connection and depending on its derived
  // type synapse building functions may create details such as neuroanatomy
  // and neurophysiology, including morphology and arborization for a
  // axon-dendritic or axon-somatic synapse.
  c->Synapses()->link_before(this);
  if (s) s->require();
#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
  synapse_inventory[type_id][SYNGENESIS]++;
#endif
}

synapse::~synapse() {
  if (s) s->release();
#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
  synapse_inventory[type_id][SYNLOSS]++;
#endif
  //cout << "DS" << String((long) this) << '\n'; cout.flush();
}

neuron * synapse::Presynaptic_Neuron() {
  return c->PreSynaptic();
}

neuron * synapse::Postsynaptic_Neuron() {
  return c->PostSynaptic();
}

Fig_Object * synapse::net_Fig() {
  synapse_type stype = this->type_ID();
  if (stype<syntype_candidate) {
    if (!figattr_receptor[stype]) return 0;
    _figvar_connection_color=colortable->colnum(static_cast<color_table_name>(CT_AMPAR+stype));
  }
  Fig_Circle * fc = 0;
  if (s) return s->net_Fig();
  return fc;
  // [***INCOMPLETE] I may wish to add the following:
  // a) radio button selection of _show_synapse_structure and _show_synapse_dot
  // b) derived type net_Fig() that call this net_Fig() first, thereby setting
  //    up the color, which can then draw a dot of a size that relates to some
  //    parameter value (e.g. strength)
}

Txt_Object * synapse::net_Txt() {
  (*Txt_synapselist) += String(Txt_synapseindex);
  (*Txt_synapselist) += ',';
  (*Txt_synapselist) += synapse_type_name[type_id];
  ///(*Txt_synapselist) += ','+String((long) this);
  ///(*Txt_synapselist) += ','+String((long) s);
  double x=0.0, y=0.0, z=0.0;
  ///(s->P0).get_all(x,y,z);
  ///(*Txt_synapselist) += String(x,",%f");
  ///(*Txt_synapselist) += String(y,",%f");
  ///(*Txt_synapselist) += String(z,",%f");
  (s->P1).get_all(x,y,z);
  (*Txt_synapselist) += String(x,",%f");
  (*Txt_synapselist) += String(y,",%f");
  (*Txt_synapselist) += String(z,",%f,");
  ///(*Txt_synapselist) += String((long) (s->AxonSegment())) + ',';
  ///(*Txt_synapselist) += String((long) (s->DendriteSegment())) + ',';
  (*Txt_synapselist) += String((long) Presynaptic_Neuron());
  (*Txt_synapselist) += ',';
  (*Txt_synapselist) += String((long) Postsynaptic_Neuron());
  if (SynaptoGenesis_Data) (*Txt_synapselist) += String(SynaptoGenesis_Data->find_t_genesis(this),",%f\n");
  Txt_synapseindex++;
  return NULL;
}

void synaptogenesis_data::add(synapse * s, double t) {
  if (dataidx>=datalen) {
    unsigned int newdatalen = datalen + 10240;
    synapseptr * sptr = new synapseptr[newdatalen];
    double * d = new double[newdatalen];
    for (unsigned int i = 0; i<dataidx; i++) {
      sptr[i] = synlist[i];
      d[i] = data[i];
    }
    datalen = newdatalen;
    delete[] synlist;
    delete[] data;
    synlist = sptr;
    data = d;
  }
  synlist[dataidx] = s;
  data[dataidx] = t;
  dataidx++;
}

double synaptogenesis_data::find_t_genesis(synapse * s) {
  for (unsigned int i = 0; i<dataidx; i++) if (synlist[i] == s) return data[i];
  return -1.0;
}
