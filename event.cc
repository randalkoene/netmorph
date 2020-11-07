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
// event.cc
// Randal A. Koene, 20050330

#include "event.hh"
#include "network.hh"
#include "global.hh"

Event_Queue * eq = NULL; // globally accessible
Null_Activity NullActivity;

double max_random_spike_interval = 100.0;

Event::Event(): t(eq->T()) {}

void Event_Queue::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("events_seconds"))>=0) T_end = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("events_append_seconds"))>=0) T_end = max_growth_time + atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("maxrandomspikeinterval"))>=0) max_random_spike_interval = atof(clp.ParValue(n));
}

String Event_Queue::report_parameters() {
  String res("General events parameters:\n");
  res += "  Duration of events simulation = "+String(T_end,"%.3f seconds.\n");
  return res;
}

// [***NOTE] If the bool result is never used then I could make this function
// type void for a small potential speed-up.
bool Event_Queue::add(Event * e) {
  if (!e) return false;
  // [***INCOMPLETE] This is a terribly slow way to sort!
  if (!head()) return link_before(e); // is only
  PLL_LOOP_FORWARD_NESTED(Event,head(),1,q) if (q->T()>=e->T()) return insert_before(q,e);
  return link_before(e); // at end
}

void Event_Queue::variablesteps() {
  Event * e;
  while ((e=head())!=NULL) {
    if ((t=e->T())>=T_end) break; // set global t to event t
    e->event();
    e->remove();
  }
}

void Event_Queue::fixedsteps(double dt) {
  
}

void Random_Spikes::event() {
  // The event handlers within objects that do scheduling, predicting
  // and the allocation of specific event elements can deal with any
  // number of different types of events. For that reason, the ones
  // that deal with activity in the network are addressed through
  // activity() functions.
  // [***NOTE] I may implement chains of event handlers in the same
  // fashion as I did for direction models.
  network * n = eq->Net();
  int N = X_misc.get_rand_range(0,n->PLLRoot<neuron>::length()-1); // [***NOTE] There should probably be a X_event.
  n->PLLRoot<neuron>::el(N)->activity()->queue(); // specific handler module schedules specific type of spike
  eq->add(new Random_Spikes(t+(maxinterval*X_misc.get_rand_real1())));
}

Null_Activity::~Null_Activity() {
  warning("Warning: Attempting to destruct Null_Activity object.\n");
}

void Null_Activity::queue(int eventtype) {
  warning("Null_Activity\n");
}

void Null_Activity::queue(double t, int eventtype) {
  warning("Null_Activity\n");
}
