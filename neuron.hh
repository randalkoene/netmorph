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
// neuron.hh
// Randal A. Koene, 20041118
//
// Classes the define neuron types.

#ifndef __NEURON_HH
#define __NEURON_HH

// forward declarations for pointers here and in included headers
class neuron;
class multipolar_nonpyramidal;
class pyramidal;
class bipolar;
class interneuron;
class electrode;
class neuronsetel;
class neuronset;

#ifndef __NETWORK_HH
class network;
#endif

#ifndef __NETWORK_GENERATED_STATISTICS_HH
class network_statistics_base;
#endif

#include <set>

#include <include/templates.hh>
#include "global.hh"
#include "Color_Table.hh"
#include "prepost_structure.hh"
#include "connection.hh"
#include "slice.hh"
#include "spatial.hh"
#include "event.hh"
#include "Network_Generated_Statistics.hh"
#include "Command_Line_Parameters.hh"
#include "state_storable.hh"
#include "Txt_Object.hh"
#include "VRML_Object.hh"
#include "Catacomb_Object.hh"
#include "Include/Embeddable.h"

enum neuron_type { PRINCIPAL_NEURON, INTERNEURON, MULTIPOLAR_NONPYRAMIDAL, BIPOLAR, PYRAMIDAL, UNTYPED_NEURON };
// Note: Untyped is the last one in the enumeration list so that I can do
// such things as int a[UNTYPED_NEURON+1] easily.

enum basal_force_model { unrestricted_bfm, surface_division_bfm, forced_drift_bfm, NUM_bfm };

#define VmREST -60.0

extern const char neuron_type_name[][25];
extern const char neuron_short_name[][12];

extern const char natural_subset_idstr[][40];

extern int min_basal[];
extern int max_basal[];
extern double min_angle[];
extern double max_angle[];
extern int max_axons[];

extern bool pia_attraction_repulsion_hypothesis;

class general_neuron_parameters_interface: public CLP_Modifiable {
  // This class serves to create the general_neuron_parameters object
  // for parameters that should apply to all neurons in generated networks.
  // [*** NOTE] This may be incorporated into network or neuron classes
  // if preferable.
protected:
  network * cached_net_ptr;
public:
  general_neuron_parameters_interface(): cached_net_ptr(NULL) {}
  virtual ~general_neuron_parameters_interface() {}
  network * Cached_Net() { return cached_net_ptr; }
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual void parse_CLP(Command_Line_Parameters & clp, network & net);
  virtual String report_parameters();
};

extern general_neuron_parameters_interface general_neuron_parameters;

extern long num_neurons_created;

/**
 * Base class for operators that should be applied to the whole set of
 * neurons.
 * See for example how this is used by BuildFromNetmorphNetwork when
 * Netmorph is embedded in NES.
 * (RK 20240713)
 */
class neuron_list_op {
public:
  // This is run whenever entering a new neuron.
  virtual void op(neuron* n) = 0;
};

/**
 * This data is used when cell attraction is identical for all
 * growthcones of a specific neuron. This is used by cell_attraction
 * in axon_direction_model.
 */
struct neuron_chemical_factors {
  std::set<int> attractedto; // Indices of chemical factors attracted to.
  std::set<int> repelledby; // Indices of chemical factors repelled by.

  std::set<int> attractor; // Indices of chemical factors this attracts with.
  std::set<int> repellor; // Indices of chemical factors this repels with.
};

// <A NAME="one-clear-on-multiple-linked-lists">Only one clear() must be called on a multiply linked list</A>
// See <A HREF="../../../doc/html/lists/task-log.20040915.html#200410060901">TL#200410060901</A>.
/**
 * Beware: We've had some issues using the this pointer in functions
 * of the virtual base class. It might be necessary to use the
 * this_neuron() function instead.
 */
class neuron: public PLLHandle<neuron>, public state_storable {
  friend class connection;

protected:
  long numberID; // numerical ID that is unique for one neuron that is labeled by its address, and which is assigned automatically

  spatial P; // spatial position of center of soma
  double radius; // half diameter in micrometers

  // For the following, also see TL#200603101357.2.
  PLLRoot<fibre_structure> outputstructure; // axonal morphology (allocate presynaptic_structure)
  PLLRoot<fibre_structure> inputstructure; // dendritic morphology (allocate postsynaptic_structure)
  PLLRoot<posttopre_connection> inputconnections; // abstract connection reference post from pre
  PLLRoot<connection> outputconnections; // abstract connection reference pre to post
  bool abstracted_connections;

  double Vm;
  Activity * a;

  int chemdata_index = 0;
  int RotateToNextAttractor(int reuse_attractor);

public:
  neuron(): radius(5.0), Vm(VmREST), abstracted_connections(false), a(&NullActivity), figneuroncolor(colortable->colnum(CT_neuron_untyped)), figconnectionscolor(colortable->colnum(CT_connection_excitatory)), figsynapsescolor(colortable->colnum(CT_synapses)), figdendritescolor(colortable->colnum(CT_dendrites)), figaxonscolor(colortable->colnum(CT_axon_excitatory)), figattr(0) { numberID = num_neurons_created; num_neurons_created++; }
  neuron(spatial & npos): P(npos), radius(5.0), Vm(VmREST), abstracted_connections(false), a(&NullActivity), figneuroncolor(colortable->colnum(CT_neuron_untyped)), figconnectionscolor(colortable->colnum(CT_connection_excitatory)), figdendritescolor(colortable->colnum(CT_dendrites)), figaxonscolor(colortable->colnum(CT_axon_excitatory)), figattr(0) { numberID = num_neurons_created; num_neurons_created++; }
  ~neuron() { if (a!=&NullActivity) delete a; }

  // self references and identifiers
  virtual void* this_neuron() = 0; // See Beware note above.
  String label();
  long numerical_ID() { return numberID; }
  String numerical_IDStr() { return String(numberID); }

  virtual neuron_type TypeID() { return UNTYPED_NEURON; }
  const char * TypeChars() { return neuron_short_name[TypeID()]; }
  String TypeStr() { return String(TypeChars()); }

  // configuration  
  virtual void parse_CLP(Command_Line_Parameters & clp) = 0;
  virtual void neuron_specific_configurator(String config_clp, String valuestr);
  String neuron_specific_reports();

  // architecture, morphology
  void set_position(spatial & npos) { P = npos; }
  void set_position_in_Z_plane(double xpos, double ypos) { P.set_all(xpos,ypos); }
  spatial & Pos() { return P; }
  void set_radius(double r) { radius=r; }
  double Radius() { return radius; }

  int total_input_segments();
  int total_output_segments();
  int total_input_terminal_segments();
  int total_output_terminal_segments();

  void move(spatial & pos); // Re-center a neuron on pos
  void fanin_rot(); // in spherical coordinates, rotate all points to theta=0

  // chemical attraction
  neuron_chemical_factors chemdata;
  bool is_attractor() const { return !chemdata.attractor.empty(); }
  bool is_repellor() const { return !chemdata.repellor.empty(); }
  int NextAttractor(int reuse_attractor, bool possibly_rotate_attractors);
#ifdef TESTING_SIMPLE_ATTRACTION
  int attracts = 0;
  int attractedto = 0;
#endif

  // dynamics
  void set_activity(Activity * _a) { a = _a; }
  Activity * activity() { return a; }
  virtual double Membrane_Potential() { return Vm; }

  // connections
  virtual void initialize_output_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void initialize_input_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  PLLRoot<connection> * OutputConnections() { return &outputconnections; }
  PLLRoot<fibre_structure> * InputStructure() { return &inputstructure; }
  PLLRoot<fibre_structure> * OutputStructure() { return &outputstructure; }
  void clear_output_connections() { outputconnections.clear(); }
  void clear_input_connections() { inputconnections.clear(); }
  connection * connect_to(neuron * postsyn); // (see nibr.cc)
  connection * connect_from(neuron * presyn); // (see nibr.cc)
  bool is_connected_to(neuron* postsyn) const;
  bool is_connected_from(neuron* presyn) const;
  void abstract_connections();
  bool has_abstracted_connections() { return abstracted_connections; }
  int num_input_connections() const;
  int num_output_connections() const;
  int number_of_synapses();

  // information output
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(network_statistics_base & nsb);
#endif

  // development output
  virtual Fig_Object * net_Fig();
  virtual Txt_Object * net_Txt();
  virtual bool net_NES(NetData& netdata);
  virtual VRML_Object * net_VRML();
  //virtual Catacomb_Object * net_Catacomb();
#ifdef VECTOR3D
  virtual void net_Slice(Slice * slice);
#endif

  long figneuroncolor;
  long figconnectionscolor;
  long figsynapsescolor;
  long figdendritescolor;
  long figaxonscolor;
  int figattr; // attributes used during .fig file output
  void set_fig_colors(long n, long c, long d, long a) { figneuroncolor=n; figconnectionscolor=c; figdendritescolor=d; figaxonscolor=a; }
  inline bool fig_presyn_is_visible() { return (!(figattr & FIGMASK_SHOWPRESYN)); }
  inline bool fig_postsyn_is_visible() { return (!(figattr & FIGMASK_SHOWPOSTSYN)); }
  inline bool fig_cell_is_visible() { return (!(figattr & FIGMASK_SHOWCELL)); }
  inline bool fig_report_depth() { return (figattr & FIGMASK_REPORTDEPTH); }

  // traversing operations
  void tree_op(fibre_tree_op& op); // Fiber tree traversal operation, applied over all fiber structures.
  void neuron_op(neuron_list_op& op);
  void synapse_op(synapse_tree_op& op); // Synapse tree traversal operation, applied over all connections.

  // cache variables
  common_cache cache;
};

typedef neuron * neuronptr;

class neuronptrlist: public PLLHandle<neuronptrlist> {
protected:
  neuronptr nptr;
public:
  neuronptrlist(neuronptr _nptr): nptr(_nptr) {}
  operator const neuronptr() const { return nptr; }
  neuronptr N() { return nptr; }
};

class principal: public neuron {
public:
  principal() { set_fig_colors(colortable->colnum(CT_neuron_principal),colortable->colnum(CT_connection_excitatory),colortable->colnum(CT_dendrites),colortable->colnum(CT_axon_excitatory)); } //parse_CLP(*main_clp); }
  principal(spatial npos): neuron(npos) { figneuroncolor=colortable->colnum(CT_neuron_principal); figconnectionscolor=colortable->colnum(CT_connection_excitatory); } //parse_CLP(*main_clp); }
  virtual void* this_neuron() { return (void*) this; };
  virtual neuron_type TypeID() { return PRINCIPAL_NEURON; }
  virtual void parse_CLP(Command_Line_Parameters & clp);
};

typedef principal * principalptr;

class multipolar_nonpyramidal: public principal {
protected:
  void initialize_multiple_poles_common(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev, bool isinputstructure, int _typeid = -1); // was isinputstructure = false
  void initialize_multipolar_output_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle);
  void initialize_multipolar_input_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, int _typeid = -1);
  void initialize_multipolar_output_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev);
  void initialize_multipolar_input_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev, int _typeid = -1);
  void surface_division_basal_force_model(int numpoles, double mintotlength, double maxtotlength, double minangle, double maxangle, const spatial & rc, bool isinputstructure, int _typeid);
  void forced_drift_basal_force_model();
public:
  multipolar_nonpyramidal() {}
  multipolar_nonpyramidal(spatial & npos): principal(npos) {}
  virtual void* this_neuron() { return (void*) this; };
  virtual neuron_type TypeID() { return MULTIPOLAR_NONPYRAMIDAL; }
  virtual void initialize_output_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void initialize_input_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void parse_CLP(Command_Line_Parameters & clp);
};

class bipolar: public multipolar_nonpyramidal {
public:
  bipolar() {}
  bipolar(spatial & npos): multipolar_nonpyramidal(npos) {}
  virtual void* this_neuron() { return (void*) this; };
  virtual neuron_type TypeID() { return BIPOLAR; }
  virtual void initialize_output_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void initialize_input_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
};

class pyramidal: public multipolar_nonpyramidal {
public:
  pyramidal() {}
  pyramidal(spatial & npos): multipolar_nonpyramidal(npos) {}
  virtual void* this_neuron() { return (void*) this; };
  virtual neuron_type TypeID() { return PYRAMIDAL; }
  virtual void initialize_output_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void initialize_input_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
};

class interneuron: public multipolar_nonpyramidal {
public:
  interneuron() { set_fig_colors(colortable->colnum(CT_neuron_interneuron),colortable->colnum(CT_connection_inhibitory),colortable->colnum(CT_dendrites),colortable->colnum(CT_axon_inhibitory)); }
  interneuron(spatial & npos): multipolar_nonpyramidal(npos) { figneuroncolor=colortable->colnum(CT_neuron_interneuron); figconnectionscolor=colortable->colnum(CT_connection_inhibitory); }
  virtual void* this_neuron() { return (void*) this; };
  virtual neuron_type TypeID() { return INTERNEURON; }
  virtual void initialize_output_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
  virtual void initialize_input_structure(double mintotlength, double maxtotlength); // (see nibr.cc)
};

class electrode: public neuron {
public:
  electrode() { radius = 12.0/2.0; } // diameter from [PELT:LONGTERMDYN]
  electrode(double r) { radius = r; }
  virtual void* this_neuron() { return (void*) this; };
  virtual void parse_CLP(Command_Line_Parameters & clp) {}
  virtual Fig_Object * net_Fig();
};

class neuronsetel: public PLLHandle<neuronsetel> {
protected:
  neuron * n;
public:
  neuronsetel(neuron & nref): n(&nref) {}
  neuron * Neuron() { return n; }
};

class neuronset: public PLLRoot<neuronsetel> {
public:
  neuronset() {}
};

//extern neuron * processing_neuron; // *** this cannot be global when parallel processing

#endif
