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
// fibre_structure.hh
// Randal A. Koene, 20041119, 20060311
//
// Classes that define fibre structure and functions that apply both to
// presynaptic (axonal) and postsynaptic (dendritic) fibre structure..

// [Update AC 201110322:]  Added APICAL flag to fibre_segment

#ifndef __FIBRE_STRUCTURE_HH
#include <unistd.h>
#include <math.h>
// Required class headers
#include "global.hh"
#include "templates.hh"
#include "spatial.hh"
#include "state_storable.hh"
#include "Command_Line_Parameters.hh"
#include "fibre_elongation_model.hh"
//#include "dendritic_growth_model.hh"
#include "branching_models.hh"
#include "turning_models.hh"
#include "environment_physics.hh"
#include "event.hh"
#include "Color_Table.hh"
#include "diagnostic.hh"
#include "Txt_Object.hh"
#include "VRML_Object.hh"
#define __FIBRE_STRUCTURE_HH

// Class pointer forward declarations (to avoid recursive inclusion)
class fibre_segment;
class fibre_structure;

#ifndef __DENDRITIC_GROWTH_MODEL_HH
class Delayed_Branching_Model;
#endif
#ifndef __NEURON_HH
class neuron;
#endif
#ifndef __FIG_OBJECT_HH
class Fig_Object;
#endif
#ifndef __NETWORK_GENERATED_STATISTICS_HH
class network_statistics_base;
class network_statistics_data;
#endif
#ifndef __AXON_DIRECTION_MODEL_HH
class direction_model_base;
#endif
#ifndef __BRANCHING_MODELS_HH
class branching_model_base;
class branch_angle_model_base;
#endif

typedef fibre_segment * fibre_segment_ptr;

// The following is in micrometers
#define DEFAULTFIBRERADIUS 1.0

#ifdef VECTOR3D
// Currently, the only extra fibre data that is needed for some work on
// arbor fibre segments is a double, namely the theta angle relative to
// a coordinate system based on the parent segment.
#define extra_fibre_data double

extern extra_fibre_data root_efd;
#endif

union common_cache {
  char   ch;
  int    i;
  long   l;
  float  f;
  double d;
};

extern fibre_structure_class_id parsing_fs_type;

class general_fibre_structure_parameters_interface: public CLP_Modifiable {
  // This is used to set and provide general (default) parameter settings
  // used in functions of fibre_segment and related classes.
public:
  general_fibre_structure_parameters_interface(): fibrediameter(false) {}
  virtual ~general_fibre_structure_parameters_interface() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  bool fibrediameter;
};

extern general_fibre_structure_parameters_interface general_fibre_structure_parameters;

struct fibre_collected_data {
  // This structure simplifies the passing of data during Collect_Data() calls.
public:
  network_statistics_base * nsb; // The object that collects statistics.
  double arborlength; // Total sum of fiber segment lengths in arbor.
  fibre_segment * prevbifurcation; // Points to fiber before last bifurcation.
  double lensincesoma; // Fiber length since the soma (along this path).
  double lensincebifurcation; // Fiber length since the last bifurcation.
  int turnssincebifurcation; // Turns since the last bifurcation.
  bool ispresynaptic; // Collecting data in pre- or postsynaptic structure.
  fibre_collected_data(network_statistics_base & _nsb, bool _ispresynaptic): nsb(&_nsb), arborlength(0.0), prevbifurcation(NULL), lensincesoma(0.0), lensincebifurcation(0.0), turnssincebifurcation(0), ispresynaptic(_ispresynaptic) {}
  inline void branch_at_end(fibre_segment & fs) {
    prevbifurcation = &fs;
    lensincebifurcation = 0.0;
    turnssincebifurcation = 0;
  }
};

class nodegenesis_data {
  // This data structure keeps track of nodegenesis information, such as the time at which a specific node was created.
protected:
  double * data;
  fibre_segment_ptr * nodelist;
  unsigned int datalen;
  unsigned int dataidx;
public:
  nodegenesis_data(): datalen(10240), dataidx(0) { data = new double[10240]; nodelist = new fibre_segment_ptr[10240]; }
  ~nodegenesis_data() { delete[] nodelist; delete[] data; }
  void add(fibre_segment * s, double t);
  fibre_segment_ptr get_node(unsigned int i) { if (i<dataidx) return nodelist[i]; else return NULL; }
  double get_t_genesis(unsigned int i) { if (i<dataidx) return data[i]; else return -1.0; }
  double find_t_genesis(fibre_segment * s);
};

extern nodegenesis_data * NodeGenesis_Data;

class fibre_segment: public Segment, public state_storable {
protected:
  neuron * n;
#ifdef VECTOR3D
  extra_fibre_data d;
#endif
  fibre_segment * parent, * branch1, * branch2;
  // *** I can make this able to have more than one branch, or I can
  // create different classes for that and change fibre_structure from
  // inheriting fibre_segment to linking to a type of fibre_segment.
  int APICAL;
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  double diameter; // [***NOTE] I can (a) replace this with a pointer to a bunch of specialized data, or (b) used derived classes such as the "solid" classes defined below.
#endif
protected:
  void fibre_segment_initialization() {
    if (NodeGenesis_Data) NodeGenesis_Data->add(this,eq->T());
  }
public:
#ifdef VECTOR3D
  fibre_segment(neuron & _n, fibre_segment * p, extra_fibre_data & efd): n(&_n), d(efd), parent(p), branch1(NULL), branch2(NULL), APICAL(parent->APICAL)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  { fibre_segment_initialization(); }
  fibre_segment(neuron & _n, fibre_segment * p, extra_fibre_data & efd, Segment & s, int AP): Segment(s), n(&_n), d(efd), parent(p), branch1(NULL), branch2(NULL), APICAL(AP)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  { fibre_segment_initialization(); }
  fibre_segment(neuron & _n, fibre_segment * p, extra_fibre_data & efd, spatial & segstart, spatial & segend): Segment(segstart,segend), n(&_n), d(efd), parent(p), branch1(NULL), branch2(NULL), APICAL(parent->APICAL)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  { fibre_segment_initialization(); }
  //  extra_fibre_data & Extra_Fibre_Data() { return d; }
  extra_fibre_data Extra_Fibre_Data() { return d; }
  double parent_theta() { return d; }
  bool branch(extra_fibre_data & efd1, extra_fibre_data & efd2); // create a bifurcation point
  bool turn(extra_fibre_data & efd); // create a continuation point
  fibre_segment * continuation_node_to_branch(extra_fibre_data & efd);
#endif
#ifdef VECTOR2D
  fibre_segment(neuron & _n, fibre_segment * p): n(&_n), parent(p), branch1(NULL), branch2(NULL), APICAL(parent->APICAL)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  {}
  fibre_segment(neuron & _n, fibre_segment * p, Segment & s, int AP): Segment(s), n(&_n), parent(p), branch1(NULL), branch2(NULL), APICAL(AP)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  {}
  fibre_segment(neuron & _n, fibre_segment * p, spatial & segstart, spatial & segend): Segment(segstart,segend), n(&_n), parent(p), branch1(NULL), branch2(NULL), APICAL(parent->APICAL)
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
						      , diameter(0.0)
#endif
  {}
  bool branch(); // create a bifurcation point
  bool turn(); // create a continuation point
  fibre_segment * continuation_node_to_branch();
#endif
  //***fibre_segment(double a, double l): branch1(NULL), branch2(NULL), angle(a), seglength(l) {}
  virtual ~fibre_segment() { delete branch1; delete branch2; }
  neuron * N() { return n; }
  fibre_segment * Branch1() { return branch1; }
  fibre_segment * Branch2() { return branch2; }
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  double Diameter() { return diameter; }
  void set_diameter(double d) { diameter = d; }
#endif
  bool has_daughters() { return ((branch1) || (branch2)); }
  void branches_at_continuation_nodes(Delayed_Branching_Model * dbm);
  virtual double Average_Radius() { return DEFAULTFIBRERADIUS; }
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(fibre_collected_data & fcd);
#endif
  virtual int count_segments();
  virtual int count_terminal_segments();
  virtual void count_segment_types(unsigned long & continuation, unsigned long & bifurcation, unsigned long terminal);
  virtual double sum_of_intermediate_segment_lengths(network_statistics_data * internsd = NULL);
  double cartesian_soma_center_distance2();
  void fiber_location(double dist, spatial & loc);
  double cartesian_to_fiber_ratio(fibre_collected_data & fcd);
  double soma_to_term_cartesian_to_fiber_ratio(double len);
  void move_add(spatial & addvec);
  void fanin_rot();
  virtual Fig_Object * net_Fig();
  virtual Txt_Object * bifurcationnode_net_Txt(long parentnodeindex);
  virtual Txt_Object * continuationnode_net_Txt(long parentnodeindex);
  virtual Txt_Object * net_Txt(long parentnodeindex);
  virtual VRML_Object * net_VRML();
  fibre_segment * connprocid; // flag used during connection evaluation
  common_cache cache;
};

#ifdef INCLUDE_SOLID_FIBRE_SEGMENT
class solid_fibre_segment: public fibre_segment {
  // Includes a specific fibre radius (thickness) at each end.
protected:
  double R1, R2; // in micrometers
public:
#ifdef VECTOR3D
  solid_fibre_segment(neuron & _n, extra_fibre_data & efd): fibre_segment(_n,efd), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
  solid_fibre_segment(neuron & _n, extra_fibre_data & efd, Segment & s): fibre_segment(_n,efd,s), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
  solid_fibre_segment(neuron & _n, extra_fibre_data & efd, spatial & segstart, spatial & segend): fibre_segment(_n,efd,segstart,segend), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
#endif
#ifdef VECTOR2D
  solid_fibre_segment(neuron & _n): fibre_segment(_n), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
  solid_fibre_segment(neuron & _n, Segment & s): fibre_segment(_n,s), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
  solid_fibre_segment(neuron & _n, spatial & segstart, spatial & segend): fibre_segment(_n,segstart,segend), R1(DEFAULTFIBRERADIUS), R2(DEFAULTFIBRERADIUS) {}
#endif
  virtual double Average_Radius() { return 0.5*(R1+R2); }
};
#endif

class terminal_segment: public PLLHandle<terminal_segment>, public state_storable {
  // this class is simply a more rapid way to locate terminal segments
  // Environmental Pressures are applied to newly created terminal segments and
  // when terminal segment length or direction are changed, as noted below.
protected:
  fibre_structure * arbor;
  fibre_segment * terminalsegment;
  int centrifugalorder; // Gamma parameter in Eq.6 of [KOE:NETMORPH]
  spatial angularcoords; // first element is length, others are angles
  branch_angle_model_base * bamodel;
  direction_model_base * dirmodel;
  terminal_segment_elongation_model_base * elmodel;
  elongation_rate_initialization_model_base * erimodel;
  TSBM_base * tsbmodel;
  TSTM_base * tstmodel;
#ifdef ENABLE_FIXED_STEP_SIMULATION
  double fixedstepelongation; // cached value of last elongation in a fixed step simulation
#endif
public:
  terminal_segment(fibre_structure & a, fibre_segment & ts, int corder = 0): arbor(&a), terminalsegment(&ts), centrifugalorder(corder), dirmodel(NULL), elmodel(NULL), erimodel(NULL), tsbmodel(NULL)
#ifdef ENABLE_FIXED_STEP_SIMULATION
					, fixedstepelongation(0.0)
#endif
  {}
  terminal_segment(fibre_structure & a, fibre_segment & ts, spatial & acoords, int corder = 0): arbor(&a), terminalsegment(&ts), centrifugalorder(corder), angularcoords(acoords), dirmodel(NULL), elmodel(NULL), erimodel(NULL), tsbmodel(NULL)
#ifdef ENABLE_FIXED_STEP_SIMULATION
					, fixedstepelongation(0.0)
#endif
  {
    pb.direction_boundary_effect(terminalsegment->P0,angularcoords); // [***NOTE] Environmental Pressures are applied here.
    update_segment_vector();
  }
  // [***INCOMPLETE] There should be a destructor that deletes unshared model objects using their delete_shared() function!
  fibre_structure * Arbor() { return arbor; }
  fibre_segment * TerminalSegment() { return terminalsegment; }
  int CentrifugalOrder() { return centrifugalorder; }
  branch_angle_model_base * BranchAngleModel() { return bamodel; }
  direction_model_base * DirectionModel() { return dirmodel; }
  terminal_segment_elongation_model_base * ElongationModel() { return elmodel; }
  elongation_rate_initialization_model_base * ElongationRateInitializationModel() { return erimodel; }
  TSBM_base * TSBModel() { return tsbmodel; }
  TSTM_base * TSTModel() { return tstmodel; }
  void set_branch_angle_model(branch_angle_model_base * bam) { bamodel = bam; }
  void set_direction_model(direction_model_base * dm) { dirmodel = dm; }
  void set_elongation_model(terminal_segment_elongation_model_base * em) { elmodel = em; }
  void set_ERI_model(elongation_rate_initialization_model_base * eri) { erimodel = eri; }
  void set_TSBM_model(TSBM_base * tsbm) { tsbmodel = tsbm; }
  void set_TSTM_model(TSTM_base * tstm) { tstmodel = tstm; }
  virtual void set_angularcoords(spatial & acoords) {
    angularcoords = acoords;
    pb.direction_boundary_effect(terminalsegment->P0,angularcoords); // [***NOTE] Environmental Pressures are applied here.
  }
  virtual void set_length(double l) {
    angularcoords.set_X(l);
    pb.direction_boundary_effect(terminalsegment->P0,angularcoords); // [***NOTE] Environmental Pressures are applied here.
  }
  virtual void add_length(double l) {
    angularcoords.add_X(l);
#ifdef ENABLE_FIXED_STEP_SIMULATION
    fixedstepelongation = l;
#endif
    pb.direction_boundary_effect(terminalsegment->P0,angularcoords);  // [***NOTE] Environmental Pressures are applied here.
  }
#ifdef ENABLE_FIXED_STEP_SIMULATION
  double elongated_length() { return fixedstepelongation; }
  double randomize_last_segment_elongation(double mininterval);
#endif
  spatial & AngularCoords() { return angularcoords; }
  double Length() { return angularcoords.X(); }
  double Theta() { return angularcoords.Y(); }
#ifdef VECTOR3D
  double Phi() { return angularcoords.Z(); }
#endif
  void update_segment_vector() { terminalsegment->set_by_angular_coordinates(angularcoords); } // updates the fibre_segment end point
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  void Collect_Data(fibre_collected_data & fcd);
#endif
#ifdef VECTOR3D
  bool branch(PLLRoot<terminal_segment> * prevts);
  bool turn(neuron * n);
  bool turn_without_update_of_segment_vector(spatial & acoords, extra_fibre_data & efd);
#endif
#ifdef VECTOR2D
  bool branch(PLLRoot<terminal_segment> * prevts);
  bool turn(neuron * n);
  bool turn_without_update_of_segment_vector(spatial & acoords);
#endif
};

typedef terminal_segment * terminal_segment_ptr;

class fibre_structure: public PLLHandle<fibre_structure>, public fibre_segment {
  // This object is the root of an arbor of fibre_segment branches.
protected:
  PLLRoot<terminal_segment> terminalsegments;
  terminal_segment_ptr * terminalsegmentsarray;
  int terminalarraylength;
  arbor_elongation_model_base * elmodel;
  branching_model_base * bmmodel;
#ifdef REDUCED_MEMORY_PTSEM
  probability_distribution_function * tsem_pdf; // When not REDUCED_MEMORY_PTSEM this is in terminal_segment_elongation_model
#endif
  double sum_l_i_cache;
  double sum_l_i_cache_t;
  long colnum; // numerical color reference ID in color table
public:
  // Note: Creating fibre_structure as presynaptic_structure or as
  // postsynaptic_structure assures that elmodel is valid.
#ifdef VECTOR3D
  fibre_structure(neuron & _n): fibre_segment(_n,NULL,root_efd), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, Segment & s, int AP): fibre_segment(_n,NULL,root_efd,s,AP), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, spatial & segstart, spatial & segend): fibre_segment(_n,NULL,root_efd,segstart,segend), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, Segment & s, spatial & acoords, int AP): fibre_segment(_n,NULL,root_efd,s,AP), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this,acoords)); }
#endif
#ifdef VECTOR2D
  fibre_structure(neuron & _n): fibre_segment(_n,NULL), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, Segment & s, int AP): fibre_segment(_n,NULL,s,AP), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, spatial & segstart, spatial & segend): fibre_segment(_n,NULL,segstart,segend), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }

  fibre_structure(neuron & _n, Segment & s, spatial & acoords, int AP): fibre_segment(_n,NULL,s,AP), terminalsegmentsarray(NULL), terminalarraylength(0), elmodel(NULL), bmmodel(NULL), sum_l_i_cache(0.0), sum_l_i_cache_t(eq->T()), colnum(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this,acoords)); }
#endif
  //fibre_structure(double a, double l): fibre_segment(a,l), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,*this)); }
  virtual ~fibre_structure() { delete[] terminalsegmentsarray; if (elmodel) delete elmodel; if (bmmodel) delete bmmodel; }
  virtual int count_terminal_segments() { return terminalsegments.length(); }
  virtual int count_terminal_segments_in_PLL(); // E.g. all terminal segments of all axons, or all terminal segments of all dendrites
  PLLRoot<terminal_segment> * TerminalSegments() { return &terminalsegments; }
  virtual terminal_segment_ptr * array_of_terminal_segments();
  int TerminalArrayLength() { return terminalarraylength; }
  arbor_elongation_model_base * ElongationModel() { return elmodel; }
  branching_model_base * BranchingModel() { return bmmodel; }
  virtual double total_arbor_length(double * totterminalseglength = NULL, network_statistics_data * internsd = NULL);
  void center_of_gravity(spatial & cg);
#ifndef ELONGATION_EVENT_PREDICTION
  virtual void batch_elongate();
#endif
  double sum_of_perturbed_expected_elongations();
  void move_add(spatial & addvec);
  void fanin_rot();
  virtual Fig_Object * net_Fig();
  virtual Txt_Object * net_Txt();
  double cache; // cache used during statistics calculation
#ifdef REDUCED_MEMORY_PTSEM
  probability_distribution_function * TSEM_Pdf() { return tsem_pdf; }
#endif
};

extern terminal_segment * most_recent_branch1_ts;
extern terminal_segment * most_recent_branch2_ts;

struct structure_initialization {
  double L0_min;
  double L0_max;
  long colnum;
  structure_initialization(): L0_min(-1.0), L0_max(-1.0), colnum(0) {} // these values must identify "unspecified" (used in structure_initialization_selection())
  void init(natural_schema_parent_set n, double _L0_min = 9.0, double _L0_max = 11.0);
  String report_parameters();
  void initialize(fibre_structure * fs_root);
};

typedef structure_initialization * region_structure_initialization_ptr;
extern region_structure_initialization_ptr * structure_initialization_region_subset_schemas;

structure_initialization * structure_initialization_selection(String & label, Command_Line_Parameters & clp, structure_initialization & superior_set, structure_initialization & strucinit);

/* HOW DO YOU IMPLEMENT THIS?
   I need to be able to attach to a common pointer on pre/post synaptic
   structure of a neuron, yet I also need to be able to operate as
   either fibre_segment or solid_fibre_segment.

class solid_fibre_structure: public solid_fibre_segment {
  // This object is the root of an arbor of solid_fibre_segment branches.
protected:
  PLLRoot<terminal_segment> terminalsegments;
  terminal_segment_ptr * terminalsegmentsarray;
  int terminalarraylength;
public:
#ifdef VECTOR3D
  fibre_structure(neuron & _n): fibre_segment(_n,root_efd), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, Segment & s): fibre_segment(_n,root_efd,s), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, spatial & segstart, spatial & segend): fibre_segment(_n,root_efd,segstart,segend), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, Segment & s, spatial & acoords): fibre_segment(_n,root_efd,s), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,acoords)); }
#endif
#ifdef VECTOR2D
  fibre_structure(neuron & _n): fibre_segment(_n), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, Segment & s): fibre_segment(_n,s), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, spatial & segstart, spatial & segend): fibre_segment(_n,segstart,segend), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this)); }
  fibre_structure(neuron & _n, Segment & s, spatial & acoords): fibre_segment(_n,s), terminalsegmentsarray(NULL), terminalarraylength(0), cache(0.0) { terminalsegments.link_before(new terminal_segment(*this,acoords)); }
#endif
};*/

#endif
