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
// network.hh
// Randal A. Koene, 20041118
//
// Classes that define neuronal networks and related functions.

#ifndef __NETWORK_HH
#include "templates.hh"
#include "Network_Statistics.hh"
#include "Spatial_Segment_Subset.hh"
#include "Command_Line_Parameters.hh"
#include "Results.hh"
#include "Connection_Statistics.hh"
#include "dendritic_growth_model.hh"
#include "slice.hh"
#include "spatial.hh"
#include "event.hh"
#include "state_storable.hh"
#include "Txt_Object.hh"
#include "VRML_Object.hh"
#include "Catacomb_Object.hh"
#define __NETWORK_HH

// variables

extern bool sampling; // sampling flag to distinguish callers (should be set and reset by caller)
extern bool numneurons_is_default; // flag to indicate that numneurons was set by a fallback to default 

// classes

class network_statistics_data;

#ifndef __NETWORK_GENERATED_STATISTICS_HH
class network_statistics_base;
#endif
#ifndef __NEURITE_DIAMETER_MODEL_HH
class Neurite_Diameter_Model;
#endif

class region: public PLLHandle<region> {
  // A region object is an arbitrary collection of neurons. The network can
  // assign any number of regions, in accordance with groupings established,
  // for example during neuron placement.
protected:
  String name;
  PLLRoot<neuronptrlist> nlist;
public:
  region(String _name): name(_name) {}
  region(String _name, neuron & n): name(_name) { nlist.link_before(new neuronptrlist(&n)); }
  String Name() { return name; }
  PLLRoot<neuronptrlist> & Nlist() { return nlist; }
  region & append(neuron & n) { nlist.link_before(new neuronptrlist(&n)); return *this; }
  region & operator+=(neuron & n) { return append(n); }
  region & append(String _name, neuron & n);
  region & conditional_create(String _name, neuron & n);
  int find(neuron & n);
  spatial center();
  Catacomb_Group * net_Catacomb();
};

class regionslist: public PLLRoot<region> {
public:
  regionslist() {}
  regionslist(String _name, neuron & n) {
    link_before(new region(_name,n));
  }
  region & append(String _name, neuron & n) {
    if (head()) return head()->conditional_create(_name,n);
    link_before(new region(_name,n));
    return *head();
  }
  int find(neuron & n);
};

class network: public PLLRoot<neuron>, public Event_Queue {
protected:
  String netinfo;
  bool edges; // if false, functional edge effects are avoided by looping connections around, e.g. turning a rectangle into a torus
  bool candidate_synapses;
  bool synapses_during_development;
  spatial center;
  neuron * paramspreadmin, * paramspreadmax; // samples used to determine the spread
  Network_Statistics_Root * nsr;
  Spatial_Segment_Subset * sss; // used within develop_connection_structure()
  Neurite_Diameter_Model * ndm;
  regionslist regions;
  void add_typed_neuron(neuron_type nt, neuron * psmin, neuron * psmax);
  void remove_abstract_connections_without_synapses();
public:
  network(): Event_Queue(this), edges(true), candidate_synapses(true), synapses_during_development(true), sss(NULL), ndm(NULL) {}
  network(int numneurons, bool e = true); // (see nibr.cc)
  network(int numneurons, double principalprobability, bool e = true); // (see nibr.cc)
  network(int numneurons, neuron * psmin, neuron * psmax, double principalprobability, bool e = true); // (see nibr.cc)
  network(int numneurons, neuron * psmin, neuron * psmax, Network_Statistics_Root & netstats, bool _edges = true); // (see nibr.cc)
  virtual ~network() { delete sss; }
  regionslist & Regions() { return regions; }
  spatial & Center() { return center; }
  double Time() { return t; }
  bool Seek_Candidate_Synapses() { return candidate_synapses; }
  Shape_Hexagon_Result shape_hexagon(Command_Line_Parameters & clp);
  Shape_Rectangle_Result shape_rectangle(Command_Line_Parameters & clp);
  Shape_Circle_Result shape_circle(Command_Line_Parameters & clp);
#ifdef VECTOR3D
  Shape_Box_Result shape_box(Command_Line_Parameters & clp);
  Shape_Regions_Result shape_regions(Command_Line_Parameters & clp, neuron * psmin = NULL, neuron * psmax = NULL);
#endif
  Shape_Result shape_network(Command_Line_Parameters & clp, neuron * psmin = NULL, neuron * psmax = NULL);
  neuron * nearest(spatial & p);
  neuronset * inrange(neuron * n, double radius, neuron_type ntype = UNTYPED_NEURON); // (see nibr.cc)
  void uniform_random_connectivity(double range, int minconn, int maxconn, neuron_type sourcetype = UNTYPED_NEURON, neuron_type targettype = UNTYPED_NEURON, bool clearexisting = true); // (see nibr.cc)
  void develop_connection_structure(Connection_Statistics_Root & cstats, dendritic_growth_model & dgm, dendritic_growth_model & agm); // (see nibr.cc)
  Spatial_Segment_Subset * spatial_segment_subset() { return sss; }
  void spatial_extents(spatial & minpoint, spatial & maxpoint);
  void initialize_spatial_segment_subset(Connection_Statistics_Root & cstats, dendritic_growth_model & dgm, dendritic_growth_model & agm);
  int number_of_connections(); // (see nibr.cc)
  int number_of_synapses(); // (see nibr.cc)
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(network_statistics_base & nsb);
#endif
  double mean_number_of_output_terminal_segments(network_statistics_data * nsd = NULL);
  double mean_number_of_input_terminal_segments(network_statistics_data * nsd = NULL);
  double mean_presynaptic_structure_length(network_statistics_data & termnsd, network_statistics_data & nsd, network_statistics_data * internsd = NULL);
  double mean_postsynaptic_structure_length(network_statistics_data & termnsd, network_statistics_data & nsd, network_statistics_data * internsd = NULL);
  void move(spatial & pos); // move all neurons in the network to center their soma on pos
  void fanin_rot(); // in spherical coordinates, rotate all points to theta=0
  void abstract_connections();
  Fig_Group * abstract_connections_Fig();
  Fig_Group * Fig_Abstract_Connections(String figname, double width = 0.0, bool overwrite = false);
  Fig_Group * net_Fig(); // (see nibr.cc)
  Fig_Text * time_Fig(Fig_Object & rfo);
  Fig_Group * progress_bar_Fig(); // (see nibr.cc)
  static Fig_Group * scale_bar_Fig(Fig_Object & rfo);
  static Fig_Group * scale_bar_Fig(long xoffset, long yoffset);
  static Fig_Group * axis_arrows_Fig(Fig_Object & rfo);
  Fig_Group * Fig_Output(String figname, double width = 0.0, bool overwrite = false); // (see nibr.cc)
  Fig_Group * Fig_Neurons(String figname, double width = 0.0);
  Txt_Group * net_Txt();
  bool Txt_Output(String txtname);
  bool Data_Output_Synapse_Distance(String dataname);
  bool Data_Output_Connection_Distance(String dataname);
  VRML_Group * net_VRML();
  bool VRML_Output(String ac3dname);
  Catacomb_Group * net_Catacomb();
  bool Catacomb_Output(String ccmname);
#ifdef VECTOR3D
  void net_Slice(Slice & slice);
  bool Slice_Output(String slicename, Slice & slice);
  Fig_Group * Slice_Outlines_Output(String figname, Slice & slice, double width = 0.0, bool overwrite = false);
#endif
  void set_figattr(int fattr); // (see nibr.cc)
  void visible_pre_and_target_post(neuron & n); // (see nibr.cc)
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class electrodearray: public network {
public:
  electrodearray() {}
  electrodearray(int numelectrodes) {
    netinfo = "initialized with "+String((long) numelectrodes)+" electrodes\n";
    for (int i=0; i<numelectrodes; i++) PLLRoot<neuron>::link_before(new electrode());
  }
  void set_center(const spatial & c) { center = c; }
};

class connectivity_graph {
protected:
  network * net;
public:
  connectivity_graph(network & _n);
  ~connectivity_graph() { delete[] X; delete[] Y; delete[] R; delete[] n; }
  int N;
  long *X, *Y, *R;
  long graph_R;
  neuronptr * n; // [***NOTE] Or I could use the neuronset class.
  double maxstrength;
  int find(neuron * m);
};

#ifdef VECTOR3D

class shaped_volume: public CLP_Modifiable {
protected:
  String * label;
public:
  shaped_volume(String * l = NULL): label(l) {}
  virtual ~shaped_volume() {}
  virtual String str() = 0;
  virtual double displaywidth() = 0;
  virtual spatial random_location() = 0;
};

class disc_shaped_volume: public shaped_volume {
protected:
  spatial center;
  double radius;
  double thickness;
public:
  disc_shaped_volume(String * l = NULL): shaped_volume(l), radius(1.0), thickness(1.0) {}
  virtual String str() { return String("disc"); };
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual double displaywidth() { return 2.0*radius; }
  virtual spatial random_location();
};

class box_shaped_volume: public shaped_volume {
protected:
  spatial center;
  double height;
  double width;
  double depth;
public:
  box_shaped_volume(String * l = NULL): shaped_volume(l), height(1.0), width(1.0), depth(1.0) {}
  virtual String str() { return String("box"); };
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual double displaywidth() {
    if (height>width) {
      if (height>depth) return height;
      return depth;
    } else {
      if (width>depth) return width;
      return depth;
    }
  }
  virtual spatial random_location();
};

class sphere_shaped_volume: public shaped_volume {
protected:
  spatial center;
  double radius;
public:
  sphere_shaped_volume(String * l = NULL): shaped_volume(l), radius(1.0) {}
  virtual String str() { return String("sphere"); };
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual double displaywidth() { return 2.0*radius; }
  virtual spatial random_location();
};

class region_parameters: public PLLHandle<region_parameters> {
protected:
  String label;
  shaped_volume * sv;
  neuronptr * memberneurons;
  int nummemberneurons; // only those from the general pool
  int generalandspecificneurons; // including those specifically defined for the region
  double minneuronseparation;
  int specificneurons[UNTYPED_NEURON];
public:
  region_parameters(const char * l, Command_Line_Parameters & clp);
  ~region_parameters() { delete sv; delete[] memberneurons; }
  String Label() { return label; }
  shaped_volume & shape() { return *sv; }
  virtual String report_parameters();
  neuron * add_neurons(PLLRoot<neuron> * all, neuron * n, neuron * psmin = NULL, neuron * psmax = NULL);
  neuronptr * neurons() { return memberneurons; }
  int N() { return generalandspecificneurons; }
  int general_pool_N() { return nummemberneurons; }
  int specific_N() { return generalandspecificneurons-nummemberneurons; }
  void set_general_pool(unsigned int num) {
    generalandspecificneurons -= nummemberneurons;
    nummemberneurons = num;
    generalandspecificneurons += nummemberneurons;
  }
  double average_distance_to_nearest_neighbor();
};

class region_parameters_root: public PLLRoot<region_parameters>, public CLP_Modifiable {
protected:
  int numregions;
public:
  region_parameters_root(): numregions(0) {}
  virtual ~region_parameters_root() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  int total_N() {
    int num = 0;
    PLL_LOOP_FORWARD(region_parameters,head(),1) num += e->N();
    return num;
  }
  int total_general_pool_N() {
    int num = 0;
    PLL_LOOP_FORWARD(region_parameters,head(),1) num += e->general_pool_N();
    return num;
  }
  int total_specific_N() {
    int num = 0;
    PLL_LOOP_FORWARD(region_parameters,head(),1) num += e->specific_N();
    return num;
  }
};

#endif

#endif
