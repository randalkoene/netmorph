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
// connection.cc
// Randal A. Koene, 20041118

#include "connection.hh"
#include "neuron.hh"
#include "Sampled_Output.hh"

// variables

sum_collapse_synapses defaultAbstoConn;
abstract_to_connection * abstoconn = &defaultAbstoConn;

double max_abstract_strength = 0.0;
double * figattr_connections_thresholdlist = NULL;
int figattr_connections_thresholdlistlength = 0;
double figattr_connections_threshold = 0.0;

// classes

void sum_collapse_synapses::abstract() {
  strength = 0.0;
  double Vm = conn->PostSynaptic()->Membrane_Potential();
  PLL_LOOP_FORWARD(synapse,conn->Synapses()->head(),1) strength += (e->Reversal_Potential()-Vm)*e->Peak_Response();
  if (strength>max_abstract_strength) max_abstract_strength=strength;
  // [***INCOMPLETE] There is currently no accounting for propagation delays.
}

posttopre_connection::~posttopre_connection() {
  //  cout << "DPC" << (int) this << "(c" << (int) c << ")\n"; cout.flush();
  if (c->signal()!=PLL_MULTI_LIST_DELEGATOR_DESTRUCTING) {
    c->ptpc = NULL; // this dimension is already destructing
    c->remove();
  }
}

connection::~connection() {
  //  cout << "DC" << (int) this << '\n'; cout.flush();
  PLLHandle_Signal pllsignalbuf = pllsignal;
  pllsignal = PLL_MULTI_LIST_DELEGATOR_DESTRUCTING;
  if (ptpc) ptpc->remove(); // remove this dimension from its list
  pllsignal = pllsignalbuf;
}

connection::connection(neuron * presyn, neuron * postsyn): presynaptic(presyn), postsynaptic(postsyn) {
  ptpc = new posttopre_connection(*this);
  presyn->outputconnections.link_before(this);
  postsyn->inputconnections.link_before(ptpc);
  atc = abstoconn->create(*this);
}

connection * connection::create(neuron * presyn, neuron * postsyn) {
  connection * c = connection_exists(presyn,postsyn);
  if (c) return c;
  return new connection(presyn,postsyn);
}

Fig_Object * connection::net_Fig() {
  // *** I could move the figattr_show_abstract_connections test here
  if (!fig_in_zoom(postsynaptic->Pos())) return NULL;
  if (Abstract_Strength()<figattr_connections_threshold) return NULL;
  double postx, posty;
  postsynaptic->Pos().plane_mapped(postx,posty);
  long X,Y;
  X=FIGSCALED(postx); Y=FIGSCALED(posty);
  return new Fig_Line(0,1,_figvar_connection_color,7,_figvar_connection_depth,-1,0.0,0,0,_figvar_connection_axon_x,_figvar_connection_axon_y,X,Y);
}

// functions

connection * connection_exists(neuron * presyn, neuron * postsyn) {
  if ((!presyn) || (!postsyn)) return NULL;
  PLL_LOOP_FORWARD(connection,presyn->OutputConnections()->head(),1) if (e->PostSynaptic()==postsyn) return e;
  return NULL;
}

