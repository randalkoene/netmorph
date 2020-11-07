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
// IFactivity.hh
// Randal A. Koene, 20051202
//
// Classes and definitions for basic IF neural network activity.

// Note: I am declaring this file dependent on the spatial compilation
// for 2D or 3D conditions, i.e. producing separate object files. Although
// this may not seem to be an obvious requirement for IFactivity, it may
// be useful if I decide to add some distance dependent activity effect.
// (And it saves me the trouble of having to define a neuron_base class
// that contains non-spatial parts of the neuron class.)

#ifndef __IFACTIVITY_HH
#define __IFACTIVITY_HH

#include "event.hh"

#ifndef __NEURON_HH
class neuron;
#endif
#ifndef __CONNECTION_HH
class connection;
#endif
#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
#endif

// The following default value is in micrometers per second
#define DEFAULTUNMYELINATEDSPIKESPEED 1000000.0

class IFspike: public Event {
protected:
  neuron * n;
public:
  IFspike(neuron * _n): n(_n) {}
  IFspike(double _t, neuron * _n): Event(_t), n(_n) {}
  virtual void event();
};

class IFactivity: public Activity {
  // This is an Activity handler for neurons. If this handler is chosen and
  // linked to a neuron through neuron::set_activity(), then activity events
  // of that neuron are treated as activity in an integrate-and-fire (IF)
  // model of a neuron. When Activity is then queued for that neuron, the
  // Event that is added to the Event_Queue is an IFspike, and consequent
  // events that result from the IFspike::event() also involve model choices
  // that suit an integrate-and-fire neuron model.
protected:
  neuron * n;
public:
  IFactivity(neuron * _n): n(_n) {}
  virtual void queue(int eventtype = 0) { eq->add(new IFspike(n)); }
  virtual void queue(double t, int eventtype = 0) { eq->add(new IFspike(t,n)); }
};

class IFabstractedpsp: public Event {
protected:
  connection * c;
public:
  IFabstractedpsp(connection * _c): c(_c) {}
  IFabstractedpsp(double _t, connection * _c): Event(_t), c(_c) {}
  virtual void event();
};

class IFinteractiveaxonspike: public Event {
protected:
  fibre_segment * f;
public:
  IFinteractiveaxonspike(fibre_segment * _f): f(_f) {}
  IFinteractiveaxonspike(double _t, fibre_segment * _f): Event(_t), f(_f) {}
  virtual void event();
};

#endif
