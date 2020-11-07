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
// IFactivity.cc
// Randal A. Koene, 20051202

#include "IFactivity.hh"
#include "neuron.hh"
#include "connection.hh"
#include "Sampled_Output.hh"

// [***NOTE] For more versatility in terms of the chosen type of Activity and
// Event handlers, each addressed component, e.g. connection objects and
// presynaptic_structure objects, can establish their own Activity objects.
// Instead of creating new objects in a loop within IFspike::event(), I can
// then just call e->queue() in that loop and let the Activity handler create
// a specific type of Event.

/* Specifically (20070509):
   - let connection or neuron objects maintain pointers to the handlers to use for
     such things as spikes, psps, etc.
     Event * psphandler, PSPHandler() { return psphandler; }
   - call something such as connection->PSPHandler()->queue(t+e->Abstract_PropDelay(),e), which creates the event
   - for example a fixed_time_step_PSP could be used to either keep queueing itself at successive time steps,
     or to queue a fixed_time_step_Membrane_Potential event that keeps queueing itself at successive time steps.
     (At least for the fixed step approach, I can use the integrate and fire neuron functions that are
     given in my KOE:STMFIFO paper. I can use a dynamic version of those, or functions given in chapter 4 of
     Gerstner and Kistler's book.)
 */

void IFspike::event() {
  actres->spike(n,t); // data output spike onset time
  if (n->has_abstracted_connections()) {
    PLL_LOOP_FORWARD(connection,n->OutputConnections()->head(),1) eq->add(new IFabstractedpsp(t+e->Abstract_PropDelay(),e));
  } else {
    // [***NOTE] Here I assume that a single spike is generated in the soma and
    // that it spreads into each available axon.
    PLL_LOOP_FORWARD(fibre_structure,n->OutputStructure()->head(),1) eq->add(new IFinteractiveaxonspike(t,e));
  }
  // [***INCOMPLETE] I must implement propagation for
  // non-abstracted connections.
  // [***INCOMPLETE] Other stuff, such as resetting and predicting after
  // reset, since new conditions and intrinsic dynamics apply.
}

void IFabstractedpsp::event() {
  actres->psp(c->PreSynaptic(),c->PostSynaptic(),t); // data output PSP onset time
  // [***INCOMPLETE] Other stuff, such as predicting for the target, given
  // the new input conditions.
}

void IFinteractiveaxonspike::event() {
  // Any data for actres?
  // [***NOTE] Is there anything that can interfere with propagation, i.e. is
  // there any reason to make new predictions?
  // [***NOTE] I try to keep activity related calculations separate from
  // structural definitions, so I obtain the radius from the fiber, but
  // calculate the propagation delay here.
  double spikespeed = (f->Average_Radius()/DEFAULTFIBRERADIUS)*DEFAULTUNMYELINATEDSPIKESPEED;
  // [***INCOMPLETE] I can compute the speed once and put it in a cache, e.g.
  // the f->cache.d variable.
  double propdelay = f->Length()/spikespeed;
  if (f->Branch1()) eq->add(new IFinteractiveaxonspike(t+propdelay,f->Branch1()));
  if (f->Branch2()) eq->add(new IFinteractiveaxonspike(t+propdelay,f->Branch2()));
  // [***INCOMPLETE] I must encounter synapses!
  /*
    For this non-electrical model, do the following:
    - When actual synapses are formed from candidate synapses (synapse
      formation model), set "leads to synapse" in all its parent fibers.
      That may require parent pointers in fiber, or it may require doing the
      setting as a recursive call returns.
      In synapse_formation_model::establish_connections(), candidate synapses
      are dealt with as listed in connections, and their structures are passed
      to the actual synapses that are created. Therefore, we cannot simply set
      a flag recursively. We must actively seek out parent fiber.
    - As fiber update() functions are called, cache parameters for activity.
      Give fibers an Activity link, which can do the spike scheduling and
      hold parameters.
   */
}

/* [***INCOMPLETE] Neuron specific activity events, e.g. for dentate gyrus neurons.
   See DIL#20010201093512.1 and DIL#20050308134422.1.

   IPSCs of inhibitory GABAergic synapses between DG basket cells (interneurons) and
   DG granule cells have:
   - fast rise
   - biexponential decay with mean time constants of 2 ms and 9 ms
   - mean quantal conductance change was 1.7 nS
   - synaptic transmission involved 3 to 7 release sites
   - each release site released transmitter with a probability of about 0.5
     (The chance of zero release, i.e. transmission failure, is therefore between
     0.125 and 0.0039.)
   - paired pulse depression occurred (max. 37 percent in a 10 ms interval) with
     a recovery time constant of 2 seconds (presynaptic, but independent of
     previous release)
   - frequency dependent depression of transmission was half-max at 5 Hz after
     1000 presynaptic spikes

   More details are provided in the paper [KRAUS:DG].

   This is also a good example for the type of specific parameter data to look for
   about other neurons and the synapses between them.
 */
