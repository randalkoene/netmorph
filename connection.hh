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
// connection.hh
// Randal A. Koene, 20041118
//
// Classes the define connections between neuron objects.

#ifndef __CONNECTION_HH
#define __CONNECTION_HH

// forward declarations for pointers here and in included headers
class connection;
class posttopre_connection;

class neuron; // forward declaration, resolved in neuron.hh

//#include "neuron.hh" // include in .cc instead to avoid problems
#include <stddef.h>
#include "synapse.hh"
#include "templates.hh"
#include "Fig_Object.hh"
#include "state_storable.hh"

// variables

extern double max_abstract_strength;
extern double * figattr_connections_thresholdlist;
extern int figattr_connections_thresholdlistlength;
extern double figattr_connections_threshold;

// classes

class abstract_to_connection {
  // These are classes that convert the effects of neurophysiological
  // detail on the postsynaptic response to input activity into an
  // abstract form that is more easily computed during simulation.
protected:
  abstract_to_connection() {}
public:
  virtual ~abstract_to_connection() {}
  virtual abstract_to_connection * create(connection & c) = 0;
  virtual void abstract() = 0;
  virtual double Strength() = 0; // Some abstract measure useful for comparisons, even if not for activity calculations (depending on the abstract representation).
  virtual double PropDelay() = 0;
};

extern abstract_to_connection * abstoconn;

class sum_collapse_synapses: public abstract_to_connection {
  // This method is the most simplistic: It sums the peak
  // response amplitudes of synapses to obtain a combined
  // strength of input through a connection.
  // In doing so, it assumes a specific postsynaptic membrane
  // potential Vm and simultaneous arrival of input through
  // all synapses of this connection. The differences between
  // synaptic reversal potentials Erev and Vm are used together
  // with the peak responses of the synaptic conductances to
  // obtain approximated peak currents. Those are summed and
  // become the abstract strength value.
  // If abstract() is called only once, instead of during
  // each iteration of activity processing, then Vm is Vrest
  // so that constant strength contributions are made and
  // so that summed currents through multiple connections
  // take a constant amount of time to lead to a spike (if
  // threshold can be reached).
  // A more detailed spiking neuron network representation
  // is achieved by calling abstract() during each iteration
  // so that synaptic contributions elicited at different
  // times or changes of the synapses can be accounted for.
  // Note: A different approach is to produce analytical
  // or sampled numerical combined responses. Those can then
  // be used exactly as in a simplified spiking neuronal
  // network simulation.
  // Note: Another approach is to convert time dependent
  // processes to an event representation with summable
  // spline pieces, the event-predictive emulation method.
protected:
  connection * conn;
  double strength;
  double delay;
public:
  sum_collapse_synapses(): conn(NULL), strength(0.0), delay(1.0) {}
  sum_collapse_synapses(connection & c): conn(&c), strength(0.0), delay(1.0) {}
  virtual abstract_to_connection * create(connection & c) {
    return new sum_collapse_synapses(c);
  }
  virtual void abstract();
  virtual double Strength() { return strength; }
  virtual double PropDelay() { return delay; }
};

class posttopre_connection: public PLLHandle<posttopre_connection> {
  // this is a generic way to connect an object in many
  // different linked lists (as copied from splineevent/continuous_node.hh)
private:
  connection * c;
public:
  posttopre_connection(connection & conn) { c = &conn; }
  ~posttopre_connection(); // (see nibr.cc)
};

class connection: public PLLHandle<connection>, public state_storable {
  // in list of outputconnections of presynaptic neurons
  friend class posttopre_connection;
  //friend synapse;
private:
  posttopre_connection * ptpc; // in list of inputconnections of postsynaptic neuron
  connection(): presynaptic(NULL), postsynaptic(NULL) {
    ptpc = new posttopre_connection(*this);
    atc = abstoconn->create(*this);
  }
  connection(neuron * presyn, neuron * postsyn); // (see nibr.cc)
protected:
  neuron * presynaptic;
  neuron * postsynaptic;
  PLLRoot<synapse> synapses; // this connection represents a set of realistic synapses, the neuroanatomy and neurophysiology of their morphology and of the arborization that leads to them (axon-dendritic or axon-somatic)
  abstract_to_connection * atc;
public:
  static connection * create() { return new connection(); }
  static connection * create(neuron * presyn, neuron * postsyn); // (see nibr.cc)
  ~connection();
  neuron * PreSynaptic() { return presynaptic; }
  neuron * PostSynaptic() { return postsynaptic; }
  PLLRoot<synapse> * Synapses() { return &synapses; }
  double Abstract_Strength() { return atc->Strength(); }
  double Abstract_PropDelay() { return atc->PropDelay(); }
  void abstract() { atc->abstract(); }
  virtual Fig_Object * net_Fig(); // (see nibr.cc)
};

typedef connection * connectionptr;

// functions

connection * connection_exists(neuron * presyn, neuron * postsyn);

#endif
