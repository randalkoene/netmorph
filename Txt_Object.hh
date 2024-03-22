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
// Txt_Object.hh
// Randal A. Koene, 20070201
//
// Classes that define Txt_Object and objects that inherit it.

/* This format was developed in response to output requested by Jaap van Pelt on 20070126 (see jaap folder on askja).

"Wellicht is een lijst met synapse gegevens het snelst te produceren, daar de structuur ervan simpel is.

synapsnr , x , y , z , presyncellnr , postsyncellnr , time of origin

Tevens zal Peter een neuron lijst nodig hebben die voor ieder neuron de coordinaten van het celcentrum aangeven dus:

cellnr , xcenter , ycenter , zcenter"
*/

#ifndef __TXT_OBJECT_HH
#define __TXT_OBJECT_HH

#include <limits.h>
#include "BigString.hh"
#include "templates.hh"

// definitions

#define COLUMN_LABELS_NEURONS "#  neuron index number, neuron label, neuron type, region ID, soma center x coordinate, soma center y coordinate, soma center z coordinate\n"
#define COLUMN_LABELS_SYNAPSES "#  synapse index number, postsynaptic receptor type, presynaptic x coordinate, presynaptic y coordinate, presynaptic z coordinate, postsynaptic x coordinate, postsynaptic y coordinate, postsynaptic z coordinate, presynaptic axon fiber piece ID, postsynaptic dendrite fiber piece ID, presynaptic neuron label, postsynaptic neuron label, simulated time of synaptogenesis\n"
#define COLUMN_LABELS_ROOT_NODES "#  node index number, fiber piece label, fiber structure type, x coordinate, y coordinate, z coordinate, soma neuron label, simulated time of node genesis\n"
#define COLUMN_LABELS_FIBER_NODES "#  node index number, fiber piece label, fiber structure type, x coordinate, y coordinate, z coordinate, soma neuron label, parent label, simulated time of node genesis\n"
#define COLUMN_LABELS_TUFT_NODES "#  fiber piece label\n"
#define COLUMN_LABELS_OBLIQUE_NODES "#  fiber piece label\n"

// functions

// classes

class Txt_Object: public PLLHandle<Txt_Object> {
public:
  Txt_Object() {}
  virtual String str() { return String(""); }
};

class Txt_Header: public Txt_Object {
protected:
  String text;
public:
  Txt_Header();
  virtual String str() { return text; }
};

class Txt_Group: public Txt_Object {
protected:
  String comment;
  PLLRoot<Txt_Object> txtobjects;
public:
  Txt_Group(String c): comment(c) {}
  Txt_Group() {}
  virtual String str();
  void Add_Txt_Object(Txt_Object * txtobject);
};

extern String * Txt_neuronlist;
extern String * Txt_synapselist;
extern String * Txt_fiberrootlist;
extern String * Txt_continuationnodelist;
extern String * Txt_bifurcationnodelist;
extern String * Txt_terminalgrowthconelist;
extern String * Txt_tuftrootbranchnodelist;
extern String * Txt_obliquerootbranchnodelist;
extern long Txt_neuronindex;
extern long Txt_synapseindex;
extern long Txt_nodeindex;

#endif
