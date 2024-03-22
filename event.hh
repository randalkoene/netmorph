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
// event.hh
// Randal A. Koene, 20050330
//
// Classes and definitions for basic events.
// (See <A HREF="../../../doc/html/lists/task-log.20050325.html#200503301409">TL#200503301409</A>, and of course all event-predictive emulation resources.)
// Note: I am making this 2D/3D dependent so that I don't need to make
// IFactivity independent of spatial considerations.

#ifndef __EVENT_HH
#define __EVENT_HH

#include "templates.hh"
#include "state_storable.hh"
#include "Command_Line_Parameters.hh"

class Event_Queue;
class Event;
#ifndef __NETWORK_HH
class network;
#endif

extern Event_Queue * eq; // globally accessible
extern double max_random_spike_interval;

class Event: public PLLHandle<Event> {
protected:
  double t;
public:
  Event();
  // [***NOTE] Should I make scheduling with t<Root()->T() impossible?
  Event(double _t): t(_t) {}
  double T() { return t; }
  virtual void event() = 0;
};

class Event_Queue: public PLLRoot<Event>, public CLP_Modifiable, public state_storable {
protected:
  network * net;
  double t; // most recent update/modification time
  double T_end; // this is set through parse_CLP()
public:
  // [*** NOTE] See ~/src/nnmodels/splineevent/ for function examples.
  // The default event end time is set to the default max_growth_time plus
  // ten seconds.
  //#define DEFAULT_EVENT_END_TIME (20.0*24.0*60.0*60.0)+10.0
  // [***UPDATE 20080728] We set the default event end time to zero,
  // so that no activity events are processed after simulated development.
#define DEFAULT_EVENT_END_TIME 0.0
  Event_Queue(network * n): net(n), t(0.0), T_end(DEFAULT_EVENT_END_TIME) {}
  network * Net() { return net; }
  double T() { return t; }
  double T_End() { return T_end; }
  bool add(Event * e);
  void variablesteps();
  void fixedsteps(double dt = 0.003);
  virtual void parse_CLP(Command_Line_Parameters & clp) ;
  virtual String report_parameters();
  // [***INCOMPLETE] Must provide state_storable functions.
};

class Random_Spikes: public Event {
protected:
  double maxinterval;
public:
  Random_Spikes(): maxinterval(max_random_spike_interval) {}
  Random_Spikes(double _t): Event(_t), maxinterval(max_random_spike_interval) {}
  virtual void event();
};

class Activity {
public:
  Activity() {}
  virtual ~Activity() {}
  virtual void queue(int eventtype = 0) = 0;
  virtual void queue(double t, int eventtype = 0) = 0;
};

class Null_Activity: public Activity {
public:
  Null_Activity() {}
  ~Null_Activity();
  virtual void queue(int eventtype = 0);
  virtual void queue(double t, int eventtype = 0);
};

extern Null_Activity NullActivity;

#endif
