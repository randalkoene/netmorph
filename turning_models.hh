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
// turning_models.hh
// Randal A. Koene, 20070509

/* Documentation:
   This implements base and derived classes for models of bifurcation selection
   at the arbor level and for branch angle determination at the level of
   growth cones or terminal segments.

   These classes stive to adhere to the NETMORPH interface protocol for
   model and parameter selection according to the smallest matching set.
 */

#ifndef __TURNING_MODELS_HH
#define __TURNING_MODELS_HH

#include "templates.hh"
//#include "fibre_structure.hh"
#include "spatial.hh"
#include "Command_Line_Parameters.hh"
#include "global.hh"

#ifndef __FIBRE_STRUCTURE_HH
class fibre_segment;
class terminal_segment;
#endif
#ifndef __NEURON_HH
class neuron;
#endif

class turning_model_base {
  // Each axon/dendrite arbor is designated a turning_model object, which
  // is responsible for determining the momentary turning rate, and
  // therefore the probability that a change of direction occurs at a growth cone.
  // In an event-driven implementation of the simulation, the turning_model
  // objects can also select growth cones at which turning occurs and
  // send those events to the queue.
  // Selection of a turning_model is done by finding a smallest matching
  // set, where the universal set encompasses all.
  // [***NOTE] In a future implementation, this base model can be a template
  // to insure that all turning_model objects define the right set of
  // constructors, and that procedures are carried out correctly in all
  // turning_model objects.
protected:
  turning_model_base * contributing;
  double contributingweight;
  /*  struct turning_model_base_parameters {
    sometype something;
    turning_model_base_parameters(sometime _something): something(_something) {}
    } base_parameters; */
  virtual ~turning_model_base() { if (contributing) contributing->delete_shared(); }
public:
  turning_model_base(turning_model_base * tmcontrib, double & tmweight, turning_model_base & schema);
  turning_model_base(String & thislabel, String & label, Command_Line_Parameters & clp);
  virtual turning_model_base * clone() = 0;
  virtual void delete_shared() { delete this; }
  virtual void handle_turning(terminal_segment & ts1, terminal_segment & ts2) = 0;
  virtual void handle_turning(terminal_segment & ts) = 0;
  //  virtual spatial predict(double weight, terminal_segment * ts, neuron * e) = 0;
  //  virtual void turning(terminal_segment * ts, neuron * e) = 0;
  virtual String report_parameters_specific() = 0;
  String report_parameters();
};

typedef turning_model_base * turning_model_base_ptr;

// The turn angle is determined by a Direction Models in the axon_direction_model.hh/cc files.

// global variables
// See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.
extern turning_model default_turning_model; // identifier, changed when a different universal model is set
extern String universal_turning_model_root; // Used only to report if chaining at the universal set level.
extern turning_model_base_ptr turning_model_subset_schemas[];

// functions
turning_model_base * turning_model_selection(String & label, Command_Line_Parameters & clp, turning_model_base_ptr & tmbptr);

#endif
// __TURNING_MODELS_HH

