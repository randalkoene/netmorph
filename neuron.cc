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
// neuron.cc
// Randal A. Koene, 20041118

#include <math.h>
#include "neuron.hh"
#include "fibre_elongation_model.hh"
#include "axon_direction_model.hh"
#include "Sampled_Output.hh"
#include "BigString.hh"

//neuron * processing_neuron = NULL; // *** this cannot be global when parallel processing

const char neuron_type_name[][25] =
  {"Principal Neuron",
   "Interneuron",
   "Multipolar NonPyramidal",
   "Pyramidal",
   "Untyped Neuron"};

const char neuron_short_name[][12] =
  {"principal",
   "interneuron",
   "multipolar",
   "pyramidal",
   "untyped"};

#ifdef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
const char natural_subset_idstr[NUM_NATURAL_SPS][40] = {
  "",
  "all_axons.",
  "all_dendrites.",
  "all_pyramidal_axons.",
  "all_pyramidal_dendrites.",
  "all_interneuron_axons.",
  "all_interneuron_dendrites.",
  "all_apical_pyramidal_dendrites."};
#else
const char natural_subset_idstr[NUM_NATURAL_SUBSETS_AEM][40] = {
  "all_axons.",
  "all_dendrites.",
  "all_pyramidal_axons.",
  "all_pyramidal_dendrites.",
  "all_interneuron_axons.",
  "all_interneuron_dendrites.",
  "all_apical_pyramidal_dendrites."};
#endif

// Whether to use the hypothesis that apical attraction to chemical at pia
// and axonal repulsion from that same chemical not only affects growth
// direction but also initial segment positioning, i.e. place axon and
// apical dendrites opposite each other rather than relating them only to
// the center of gravity of basal dendritic trees.
bool pia_attraction_repulsion_hypothesis = true;

/* Parameter specification following model selection:
   -------------------------------------------------
   The general_neuron_parameters_interface demonstrates a clear standardized
   method for the specification of parameter values for successively more
   constrained subsets of components, after a model has been selected.
   This process involves the allocation of "schema" objects with parameter
   values, as specified at a given resolution (e.g. the default settings at
   the global or general level), of which the schema object subsequently
   modifies some values that are specific to the subset that the schema
   represents. The schema objects are "cloned" in the components that
   belong to the corresponding subset.

   See http://rak.minduploading.org:8080/caspan/Members/randalk/model-specification-implementation/.

   For the most advanced example that follows guidelines (as of 20060313)
   see fibre_elongation_model.hh and fibre_elongation_model.cc, and the
   comments and TL references there.

   For comments about parameter copying and for the implementation rule
   regarding t_0 of a model instance (the prev_t parameter), see
   TL#200603140355.2.
*/

// [***INCOMPLETE] Make this settable through the CLP interface.
// [***INCOMPLETE] The model selection and parameter setting interface should
// allow selection and setting at various resolutions, and a propagation
// method (e.g. cloning) should facilitate it. See both DIL#20051105072005.1
// and TL#200603111137.1 and the thread that precedes it.
#ifndef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  direction_model general_axon_direction_model = segment_history_tension_dm;
  direction_model general_dendrite_direction_model = legacy_dm;
  String general_axon_direction_model_root;
  String general_dendrite_direction_model_root;
#endif

general_neuron_parameters_interface general_neuron_parameters;

// Local temporary variables used during initialization:
#define PYRAMIDAL_MIN_BASAL 2
#define PYRAMIDAL_MAX_BASAL 5
int pyramidal_min_basal = PYRAMIDAL_MIN_BASAL;
int pyramidal_max_basal = PYRAMIDAL_MAX_BASAL;
#define PYRAMIDAL_BASAL_MINANGLE (0.1*M_PI)
#define PYRAMIDAL_BASAL_MAXANGLE (0.5*M_PI)
double pyramidal_basal_minangle = PYRAMIDAL_BASAL_MINANGLE;
double pyramidal_basal_maxangle = PYRAMIDAL_BASAL_MAXANGLE;

basal_force_model bfm = unrestricted_bfm;
String bfmstr[NUM_bfm] = {
  "unrestricted",
  "surface_division",
  "forced_drift",
};

long num_neurons_created = 0;

void natural_set_model_error(natural_schema_parent_set natural_set_id, String modelstr) {
  // PROTOCOL: Do not allow sets with undefined model schemas.
  String setstr(natural_subset_idstr[natural_set_id]);
  if (natural_set_id==universal_sps) setstr = "universal";
  error("Error: Natural set '"+setstr+"' has no associated "+modelstr+" model schema.\n");
}

void general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  // PROTOCOL: Defined sets, including natural sets, receive schema objects, either due
  // to explicit specification, or by inheriting the model and parameter choices from
  // an immediate superior set. [See task log entries of DIL#20070612153815.1.]
  int n;
  if ((n=clp.Specifies_Parameter("pia_attraction_repulsion_hypothesis"))>=0) pia_attraction_repulsion_hypothesis = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("use_specified_basal_direction"))>=0) {
    use_specified_basal_direction = (downcase(clp.ParValue(n))==String("true"));
    if (use_specified_basal_direction) {
      if ((n=clp.Specifies_Parameter("specified_basal_direction_theta"))>=0) specified_basal_direction_theta = atof(clp.ParValue(n));
#ifdef VECTOR3D
      if ((n=clp.Specifies_Parameter("specified_basal_direction_phi"))>=0) specified_basal_direction_phi = atof(clp.ParValue(n));
#endif
      bool compute_direction_from_euler = false;
      double x = 0.0, y = 0.0, z = 0.0;
      if ((n=clp.Specifies_Parameter("specified_basal_direction_relx"))>=0) { x = atof(clp.ParValue(n)); compute_direction_from_euler = true; }
      if ((n=clp.Specifies_Parameter("specified_basal_direction_rely"))>=0) { y = atof(clp.ParValue(n)); compute_direction_from_euler = true; }
#ifdef VECTOR3D
      if ((n=clp.Specifies_Parameter("specified_basal_direction_relz"))>=0) { z = atof(clp.ParValue(n)); compute_direction_from_euler = true; }
#endif
      if (compute_direction_from_euler) {
	spatial apdir(x,y,z);
	apdir.convert_to_spherical();
	specified_basal_direction_theta = apdir.Y();
#ifdef VECTOR3D
	specified_basal_direction_phi = apdir.Z();
#endif
      }
    }
  }

  // General model selection (see TL#200603120422.1), and general parameters
  String genmodel; // genmodel="" is the base model identifier
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    genmodel = natural_subset_idstr[i];
#ifdef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
    if (!direction_model_selection(genmodel,clp,direction_model_subset_schemas[superior_natural_set[i]],direction_model_subset_schemas[i])) natural_set_model_error(i,"direction");
#endif
    if (!branching_model_selection(genmodel,clp,branching_model_subset_schemas[superior_natural_set[i]],branching_model_subset_schemas[i])) natural_set_model_error(i,"branching");
    if (!branch_angle_model_selection(genmodel,clp,branch_angle_model_subset_schemas[superior_natural_set[i]],branch_angle_model_subset_schemas[i])) natural_set_model_error(i,"branch angle");
    if (!arbor_elongation_model_selection(genmodel,clp,general_arbor_elongation_model_base_subset_schemas[superior_natural_set[i]],general_arbor_elongation_model_base_subset_schemas[i])) natural_set_model_error(i,"arbor elongation");
    if (!terminal_segment_elongation_model_selection(genmodel,clp,general_terminal_segment_elongation_model_base_subset_schemas[superior_natural_set[i]],general_terminal_segment_elongation_model_base_subset_schemas[i])) natural_set_model_error(i,"terminal segment elongation");
    if (!elongation_rate_initialization_model_selection(genmodel,clp,elongation_rate_initialization_model_subset_schemas[superior_natural_set[i]],elongation_rate_initialization_model_subset_schemas[i])) natural_set_model_error(i,"elongation rate initialization");
  }
#ifndef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  if ((n=clp.Specifies_Parameter("axon_direction_model"))>=0) {
    if (downcase(clp.ParValue(n))==String("segment_history_tension")) general_axon_direction_model = segment_history_tension_dm;
    else if (downcase(clp.ParValue(n))==String("cell_attraction")) general_axon_direction_model = cell_attraction_dm;
    else general_axon_direction_model = legacy_dm;
  }
  if ((n=clp.Specifies_Parameter("axon_direction_model_label"))>=0) general_axon_direction_model_root = clp.URI_unescape_ParValue(n);
  if ((n=clp.Specifies_Parameter("dendrite_direction_model"))>=0) {
    if (downcase(clp.ParValue(n))==String("segment_history_tension")) general_dendrite_direction_model = segment_history_tension_dm;
    else if (downcase(clp.ParValue(n))==String("cell_attraction")) general_dendrite_direction_model = cell_attraction_dm;
    else general_dendrite_direction_model = legacy_dm;
  }
  if ((n=clp.Specifies_Parameter("dendrite_direction_model_label"))>=0) general_dendrite_direction_model_root = clp.URI_unescape_ParValue(n);
#endif
  if ((n=clp.Specifies_Parameter("figattr_connections_threshold"))>=0) {
    if (figattr_connections_thresholdlist) delete[] figattr_connections_thresholdlist;
    figattr_connections_thresholdlist = str2vec(clp.ParValue(n),&figattr_connections_thresholdlistlength);
    figattr_connections_threshold = figattr_connections_thresholdlist[0];
    for (int i = 1; i<figattr_connections_thresholdlistlength; i++) if (figattr_connections_thresholdlist[i]<figattr_connections_threshold) figattr_connections_threshold = figattr_connections_thresholdlist[i];
  }
}

String general_neuron_parameters_interface::report_parameters() {
  String res("Pyramidal neuron morphology hypothesis:\n");
  if (pia_attraction_repulsion_hypothesis) res += "  Opposing chemical attraction and repulsion place axon and apical dendritic arbor at opposite ends.\n";
  else res += "  Axon location and location of apical dendritic arbor are determined by the center of gravity of basal dendritic arbors.\n";
  if (use_specified_basal_direction) res += "  Use the specified basal direction with theta=" + String(specified_basal_direction_theta,"%.3f and phi=") + String(specified_basal_direction_phi,"%.3f\n");
  else res += "  Compute basal directions\n";

  // axon direction model

#ifndef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  res += "General axon direction model: ";
  switch (general_axon_direction_model) {
  case segment_history_tension_dm: res += "segment history tension."; break;;
  case cell_attraction_dm: res += "cell attraction."; break;;
  default: res += "non-specific legacy functions.";
  }
  if (!general_axon_direction_model_root.empty()) res += " (contributing models with root "+general_axon_direction_model_root+".)";
  // dendrite direction model
  res += "\nGeneral dendrite direction model: ";
  switch (general_dendrite_direction_model) {
  case segment_history_tension_dm: res += "segment history tension."; break;;
  case cell_attraction_dm: res += "cell attraction."; break;;
  default: res += "non-specific legacy functions.";
  }
  if (!general_dendrite_direction_model_root.empty()) res += " (contributing models with root "+general_dendrite_direction_model_root+".)";
#endif

  // general arbor elongation model
  res += "\nGeneral arbor elongation model: ";
  switch (general_arbor_elongation_model) {
  case polynomial_O1_aem: res += "first order polynomial approximation."; break;;
  case polynomial_O3_aem: res += "third order polynomial approximation."; break;;
  default: res += "van Pelt and Uylings.";
  }
  if (!general_arbor_elongation_model_root.empty()) res += " (contributing models with root "+general_arbor_elongation_model_root+".)";
  //res += general_arbor_elongation_model_schema->report_parameters();
  // general terminal segment elongation model
  res += "\nGeneral terminal segment elongation model: ";
  switch (general_terminal_segment_elongation_model) {
  case inertia_tsem: res += "inertia."; break;;
  default: res += "simple (or not).";
  }
  if (!general_terminal_segment_elongation_model_root.empty()) res += " (contributing models with root "+general_terminal_segment_elongation_model_root+".)";
  //res += general_terminal_segment_elongation_model_schema->report_parameters();
  // natural subsets
  // [***INCOMPLETE] Does this report the model class, or any information
  // about contributing models or their parameters?
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (general_arbor_elongation_model_base_subset_schemas[i]) {
      res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' arbor elongation model: ";
      res += general_arbor_elongation_model_base_subset_schemas[i]->report_parameters();
    }
  }
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (general_terminal_segment_elongation_model_base_subset_schemas[i]) {
      res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' terminal segment elongation model: ";
      res += general_terminal_segment_elongation_model_base_subset_schemas[i]->report_parameters();
    }
  }
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (elongation_rate_initialization_model_subset_schemas[i]) {
      if (i==universal_sps) res += "\nUniversal elongation rate initialization model: ";
      else { res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' elongation rate initialization model: "; }
      res += elongation_rate_initialization_model_subset_schemas[i]->report_parameters();
      if (i==universal_sps) if (!universal_elongation_rate_initialization_model_root.empty()) res += " (contributing models with root "+universal_elongation_rate_initialization_model_root+".)";
    }
  }
  if (!elongation_rate_initialization_model_subset_schemas[universal_sps]) res += "\nUniversal elongation_rate_initialization model: length_distribution";
#ifdef INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (direction_model_subset_schemas[i]) {
      if (i==universal_sps) res += "\nUniversal direction model: ";
      else { res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' direction model: "; }
      res += direction_model_subset_schemas[i]->report_parameters();
      if (i==universal_sps) if (!universal_direction_model_root.empty()) res += " (contributing models with root "+universal_direction_model_root+".)";
    }
  }
  if (!direction_model_subset_schemas[universal_sps]) res += "\nUniversal direction model: legacy";
#endif
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (branching_model_subset_schemas[i]) {
      if (i==universal_sps) res += "\nUniversal branching model: ";
      else { res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' branching model: "; }
      res += branching_model_subset_schemas[i]->report_parameters();
      if (i==universal_sps) if (!universal_branching_model_root.empty()) res += " (contributing models with root "+universal_branching_model_root+".)";
    }
  }
  if (!branching_model_subset_schemas[universal_sps]) res += "\nUniversal branching model: van_Pelt";
  for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
    if (branch_angle_model_subset_schemas[i]) {
      if (i==universal_sps) res += "\nUniversal branch angle model: ";
      else { res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' branch angle model: "; }
      res += branch_angle_model_subset_schemas[i]->report_parameters();
      if (i==universal_sps) if (!universal_branch_angle_model_root.empty()) res += " (contributing models with root "+universal_branch_angle_model_root+".)";
    }
  }
  if (!branch_angle_model_subset_schemas[universal_sps]) res += "\nUniversal branch angle model: balanced_forces";
  res += '\n';
  res += "Lowest threshold strength at which abstract connections are shown (when drawn): " + String(figattr_connections_threshold,"%.3f\n");
  if (figattr_connections_thresholdlist) {
    res += "Using list of thresholds in abstract connection figures: ";
    for (int i = 0; i<figattr_connections_thresholdlistlength; i++) res += String(figattr_connections_thresholdlist[i],"%.3f ");
    res += '\n';
  }
  return res;
}

connection * neuron::connect_to(neuron * postsyn) {
  // make a conection from this presynaptic neuron to a postsynaptic neuron
  return connection::create(this,postsyn);
}

connection * neuron::connect_from(neuron * presyn) {
  // make a conection to this postsynaptic neuron from a presynaptic neuron
  return connection::create(presyn,this);
}

void neuron::initialize_output_structure(double mintotlength, double maxtotlength) {
  // For untyped neurons, presynaptic structure is assumed to consist of a
  // single axonal tree.
  // initial segment orientation
  spatial acoords;
  acoords.RandomAngular(mintotlength,maxtotlength,0.0,2.0*M_PI);
  // initial segment P0
  Segment initseg(P,P);
  double len = acoords.X();
  acoords.set_X(radius);
  initseg.set_by_angular_coordinates(acoords);
  initseg.P0 = initseg.P1;
  acoords.set_X(len);
  outputstructure.link_before(new presynaptic_structure(*this,initseg,acoords));
}

void neuron::initialize_input_structure(double mintotlength, double maxtotlength) {
  // For untyped neurons, postsynaptic structure is assumed to consist of a
  // single dendritic tree.
  // initial segment orientation
  spatial acoords;
  acoords.RandomAngular(mintotlength,maxtotlength,0.0,2.0*M_PI);
  // initial segment P0
  Segment initseg(P,P);
  double len = acoords.X();
  acoords.set_X(radius);
  initseg.set_by_angular_coordinates(acoords);
  initseg.P0 = initseg.P1;
  acoords.set_X(len);
  // this call sets up terminal segment data and updates the fibre segment
  inputstructure.link_before(new postsynaptic_structure(*this,initseg,acoords));
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void neuron::Collect_Data(network_statistics_base & nsb) {
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) static_cast<presynaptic_structure *>(e)->Collect_Data(nsb);
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) static_cast<postsynaptic_structure *>(e)->Collect_Data(nsb);
}
#endif

int neuron::number_of_synapses() {
  int sum = 0;
  PLL_LOOP_FORWARD(connection,outputconnections.head(),1) sum += e->Synapses()->length();
  return sum;
}

int neuron::total_input_segments() {
  int res = 0;
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) res += e->count_segments();
  return res;
}

int neuron::total_output_segments() {
  int res = 0;
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) res += e->count_segments();
  return res;
}

int neuron::total_input_terminal_segments() {
  int res = 0;
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) res += e->count_terminal_segments();
  return res;
}

int neuron::total_output_terminal_segments() {
  int res = 0;
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) res += e->count_terminal_segments();
  return res;
}

void neuron::abstract_connections() {
  PLL_LOOP_FORWARD(connection,outputconnections.head(),1) e->abstract();
  // [***INCOMPLETE] The following should probably be un-set if arbor
  // development or some such process again affects this neuron.
  abstracted_connections = true;
}

Fig_Object * neuron::net_Fig() {
  /* Documentation:
     The parameter "figattr_box_fibre_independently" allows fibre that is
     within a zoom box to be drawn, even if the neuron that the fibre belongs
     to is not within the box. For sanity reasons, setting this parameter
     automatically also sets figattr_fibres_nobox=false.
   */
  double x, y;
  static const int somafillstyle[2] = { -1, 20 };
  P.plane_mapped(x,y);
  //  cout << "Nx=" << x << ", Ny=" << y << '\n';
  //cout << '<'; cout.flush();
  Fig_Circle * neuron_fig = NULL;
  if (fig_in_zoom(P)) {
    if ((figattr_show_neurons) && (fig_cell_is_visible())) {
      long X,Y,R;
      long color=0; if (figattr_use_color) color=figneuroncolor;
      X=(long) (x*figscale); Y=(long) (y*figscale);
      R=(long) (figscale*radius); if (R==0) R = 1;
      neuron_fig = new Fig_Circle(0,1,color,color,998,somafillstyle[figattr_fill_somas],0.0,0.0,X,Y,R,X,Y,X+150,Y+75); // ,9,20 [pre20070709, fillcolor used to be 7]
    }
  } else if (!figattr_box_fibre_independently) return NULL;
  Fig_Group * connections_fig = NULL;
  if ((figattr_show_abstract_connections) && (outputconnections.head())) {
    connections_fig = new Fig_Group(); // *** can add info in comment here
    long x1,y1;
    if (figattr_use_color) _figvar_connection_color=figconnectionscolor;
    double axonoffset=0.70711*radius;
    x1=FIGSCALED(x+axonoffset); y1=FIGSCALED(y+axonoffset); 
    _figvar_connection_axon_x=FIGSCALED(x+axonoffset+axonoffset);
    _figvar_connection_axon_y=FIGSCALED(y+axonoffset+axonoffset); 
    // A small axon segment is drawn at 45 degrees anti-clockwise from the
    // downward direction. From there, lines are drawn to connected neuron
    // centers.
    if (fig_cell_is_visible()) connections_fig->Add_Fig_Object(new Fig_Line(0,1,_figvar_connection_color,7,_figvar_connection_depth,-1,0.0,0,0, x1,y1,_figvar_connection_axon_x,_figvar_connection_axon_y));
    if (fig_presyn_is_visible()) {
      PLL_LOOP_FORWARD(connection,outputconnections.head(),1) connections_fig->Add_Fig_Object(e->net_Fig());
    } else {
      PLL_LOOP_FORWARD(connection,outputconnections.head(),1) if (e->PostSynaptic()->fig_postsyn_is_visible()) connections_fig->Add_Fig_Object(e->net_Fig());
    }
    if (fig_report_depth()) cout << "Reported Connections Depth = " << _figvar_connection_depth << '\n';
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
  //cout << '1'; cout.flush();
  Fig_Group * postsynaptic_structures_fig = NULL;
  if ((figattr_show_postsynaptic_structure) && (inputstructure.head())) {
    postsynaptic_structures_fig = new Fig_Group(); // *** can add info in comment here
    if (figattr_use_color) _figvar_connection_color=figdendritescolor;
    PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) postsynaptic_structures_fig->Add_Fig_Object(e->net_Fig());
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
  Fig_Group * presynaptic_structures_fig = NULL;
  if ((figattr_show_presynaptic_structure) && (inputstructure.head())) {
    presynaptic_structures_fig = new Fig_Group(); // *** can add info in comment here
    if (figattr_use_color) _figvar_connection_color=figaxonscolor;
    PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) presynaptic_structures_fig->Add_Fig_Object(e->net_Fig());
    /*    PLL_LOOP_FORWARD(presynaptic_structure,outputstructure.head(),1) {
      Fig_Object * fo = e->net_Fig();
presynaptic_structures_fig->Add_Fig_Object(fo);
 cout << "<PRE>\n";
 cout << fo->str();
 cout << "</PRE>\n";
 }*/
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
  Fig_Group * synapse_structures_fig = NULL;
  if ((figattr_show_synapse_structure) && (outputconnections.head())) {
    synapse_structures_fig = new Fig_Group(); // *** can add info in comment here
    if (figattr_use_color) _figvar_connection_color=figsynapsescolor; // [***NOTE] This is superceded by color allocation within synapse::net_Fig().
    // ***synapse_inventory[syntype_GluR][0]++;

    PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
      // ***synapse_inventory[syntype_iGluR][0]++;
      PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) {
	// ***synapse_inventory[syntype_mGluR][0]++;
	synapse_structures_fig->Add_Fig_Object(s->net_Fig());
      }
    }
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
  //cout << '>'; cout.flush();
  if ((!connections_fig) && (!postsynaptic_structures_fig) && (!presynaptic_structures_fig) && (!synapse_structures_fig)) return neuron_fig;
  Fig_Group * neuron_structure_fig = new Fig_Group();
  neuron_structure_fig->Add_Fig_Object(neuron_fig);
  neuron_structure_fig->Add_Fig_Object(connections_fig);
  neuron_structure_fig->Add_Fig_Object(postsynaptic_structures_fig);
  neuron_structure_fig->Add_Fig_Object(presynaptic_structures_fig);
  neuron_structure_fig->Add_Fig_Object(synapse_structures_fig);
  return neuron_structure_fig;
}

Txt_Object * neuron::net_Txt() {
  (*Txt_neuronlist) += String(Txt_neuronindex);
  (*Txt_neuronlist) += ',';
  (*Txt_neuronlist) += String((long) this);
  double x=0.0, y=0.0, z=0.0;
  P.get_all(x,y,z);
  (*Txt_neuronlist) += String(x,",%f");
  (*Txt_neuronlist) += String(y,",%f");
  (*Txt_neuronlist) += String(z,",%f\n");
  if (outputconnections.head()) {
    PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
      PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) s->net_Txt();
    }
  }
  if (outattr_track_nodegenesis) {
    parsing_fs_type = dendrite_fs;
    PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) e->net_Txt();
    parsing_fs_type = axon_fs;
    PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) e->net_Txt();
  }
  Txt_neuronindex++;
  return NULL;
}

void neuron::net_Slice(Slice * slice) { // [***INCOMPLETE] Change the return type.
  // Efficiency: Determine which spatial segments contain or intersect with a slice. Determine which fiber segments
  // belong to those spatial segments. For those fiber segments, call net_Slice().
  // [***INCOMPLETE] Within slices, draw neurons, call net_Slice() for fibers.
  // [***NOTE} I am implementing this in 2 steps: (1) I implement calling all fibers and having them test for intersection with slices,
  // (2) I make it more efficient using spatial segments.
  // [***NOTE] I will look for notes that I made about efficient implementations for the generation of voxels or slice images at a specific resolution.
  // [***NOTE] Keep in mind the imaging process used.
  // A. Collect neuron figure object.
#ifdef INCLUDE_INCOMPLETE_SLICE_CODE
  Fig_Soma * neuron_fig = NULL; // Fig_Soma draws a soma with cell walls between axon and dendrite segments in the visible surface.
  if (slice->contains(P,radius)) {
    long X,Y,R;
    long color=0; if (figattr_use_color) color=figneuroncolor;
    X=(long) (x*figscale); Y=(long) (y*figscale); // *** CORRECT THIS
    R=(long) (figscale*radius); if (R==0) R = 1;
    neuron_fig = new Fig_Soma(0,1,color,7,998,-1,0.0,0.0,X,Y,R,X,Y,X+150,Y+75); // ,9,20
  }
#endif
  // B. Collect presynaptic and postsynaptic fiber objects.
  // C. Collect synapse objects.
  // Combine all the objects and return them.
}

void multipolar_nonpyramidal::surface_division_basal_force_model(int numpoles, double mintotlength, double maxtotlength, double minangle, double maxangle, const spatial & rc, bool isinputstructure, int _typeid) {
  // 1. Compute the average surface area available to each.
    double dnumpoles = (double) numpoles;
    mintotlength /= dnumpoles; maxtotlength /= dnumpoles;
    // Approximating: Calculating angles on available circle rather than sphere, then comparing with angles between vectors.
    spatial * S1array = new spatial[numpoles];
  // 2. Place randomly, use surface distance to angle conversion to test distance to nearest neighbors for acceptance or renewed placement attempt. 
#ifdef VECTOR3D
    double minsepangle = (2.0*(maxangle-minangle)) / dnumpoles;
    double sinphi(sin(rc.Z())), cosphi(cos(rc.Z()));
    double sintheta(sin(rc.Y())), costheta(cos(rc.Y()));
    for (int i = 0; i<numpoles; i++) {
      bool acceptinit = false; int cnt = 0;
      do {
	cnt++;
	S1array[i].set_all(1.0,2.0*M_PI*X_misc.get_rand_real1(),X_misc.get_rand_range_real1(minangle,maxangle));
	S1array[i].convert_from_spherical();
	acceptinit = true;
	for (int j = 0; j<i; j++) {
	  double dotp = S1array[j]*S1array[i]; // and they all have length 1.0
	  if (acos(dotp)<minsepangle) { acceptinit = false; break; };
	}
	if ((!acceptinit) && (cnt>100)) {
	  cout << "Warning - unable to place arbors at desired minimum angular separation.\n";
	  break;
	}
      } while (!acceptinit);
    }
#endif
    for (int i = 0; i<numpoles; i++) {
#ifdef VECTOR3D
      // Double rotation transformation T S, for phi and theta of the parent
      spatial acoords;
      double cosphiSx(cosphi*S1array[i].X()), sinphiSz(sinphi*S1array[i].Z());
      //   Sx(in x,y,z) = cos(phi)cos(theta)Sx-sin(theta)Sy+sin(phi)cos(theta)Sz
      acoords.set_X((costheta*cosphiSx) - (sintheta*S1array[i].Y()) + (costheta*sinphiSz));
      //   Sy(in x,y,z) = cos(phi)sin(theta)Sx+cos(theta)Sy+sin(phi)sin(theta)Sz
      acoords.set_Y((sintheta*cosphiSx) + (costheta*S1array[i].Y()) + (sintheta*sinphiSz));
      //   Sz(in x,y,z) = -sin(phi)Sx+0Sy+cos(phi)Sz
      acoords.set_Z((cosphi*S1array[i].Z()) - (sinphi*S1array[i].X()));
      acoords.convert_to_spherical();
      acoords.set_X(radius);
#endif
#ifdef VECTOR2D
      cout << "Warning - surface division BFM is not yet implemented for 2D\n";
      double anglespread = maxangle - minangle;
      double randomangle = anglespread*((2.0*X_misc.get_rand_real1()) - 1.0); // +/- random angle
      if (randomangle<0.0) randomangle -= minangle; // negative random between minangle and maxangle
      else randomangle += minangle; // positive random between minangle and maxangle
      spatial acoords(radius,rc.Y()+randomangle); // valid relative location on membrane
#endif
      Segment initseg(P,P);
      initseg.set_by_angular_coordinates(acoords); // segment [P,P+Euler(acoords)]
      initseg.P0 = initseg.P1; // segment starts on membrane at P+Euler(acoords)
      acoords.set_X(X_misc.get_rand_range_real1(mintotlength,maxtotlength)); // initial segment length
      if (isinputstructure) inputstructure.link_before(new postsynaptic_structure(*this,initseg,acoords,_typeid));
      else outputstructure.link_before(new presynaptic_structure(*this,initseg,acoords,_typeid));
    }
    delete[] S1array;
}

void multipolar_nonpyramidal::forced_drift_basal_force_model() {
  // 1. Place all in one spot.
  // 2. Use direction of push to move apart, plus deviation (p + e).
  // 3. Stop when the total satisfaction variable does not improve for three cycles.
}

//#define DEBUG_INIT_SEGMENT_ANGLES
void multipolar_nonpyramidal::initialize_multiple_poles_common(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev, bool isinputstructure, int _typeid) {
  // mintotlength and maxtotlength constrain the length of the sum of the
  // initial segments.
  // numpoles is the number of initial segments (arbors)
  // minangle and maxangle determine the minimum and maximum angles from the
  // region center (regioncenter) between which poles are allocated
  // regioncenter is a RELATIVE vector from center of the neuron to the
  // approximate center of the distribution of initial segments
  // (e.g. regioncenter = cg - P) in SPHERICAL coordinates
  // maxdev is the maximum amplitude of the angular deviation of regioncenter
  // New initial segments are linked to the end of the list of initial
  // segments.
  spatial rc(regioncenter);
  if (numpoles<2) {
    rc.AddRandomAngularAllConstrained(mintotlength,maxtotlength,-maxdev,maxdev);
    double len = rc.X(); // length of possible initial segment at regioncenter
    rc.set_X(radius); // deviated relative vector for regioncenter to membrane
    Segment initseg(P,P);
    initseg.set_by_angular_coordinates(rc); // segment [P,P+Euler(rc)]
    initseg.P0 = initseg.P1; // segment starts on membrane at P+Euler(rc)
    rc.set_X(len); // initial segment length
    if (isinputstructure) inputstructure.link_before(new postsynaptic_structure(*this,initseg,rc,_typeid));

    else outputstructure.link_before(new presynaptic_structure(*this,initseg,rc,_typeid)); // initial segment fibre from P0 membrane location to length len in rc deviated regioncenter direction
  } else {
    if (bfm==surface_division_bfm) surface_division_basal_force_model(numpoles,mintotlength,maxtotlength,minangle,maxangle,rc,isinputstructure,_typeid);
    else {
      // [POSSIBLY INCOMPLETE***] Perhaps vary in 2D/3D around the cone circle
      // between minimum and maximum increases of angle:
      // anglespread = maxangle - minangle;
      // minangleadd = angularspread/dnumpoles;
      // maxangleadd = minangleadd;
      // minangleadd -= (multipolar_initangledevratio*minangleadd);
      // maxangleadd += (multipolar_initangledevratio*maxangleadd);
      double dnumpoles = (double) numpoles;
      mintotlength /= dnumpoles; maxtotlength /= dnumpoles;
#ifdef VECTOR3D
      double sinphi(sin(rc.Z())), cosphi(cos(rc.Z()));
      double sintheta(sin(rc.Y())), costheta(cos(rc.Y()));
#endif
      for (int i = 0; i<numpoles; i++) {
#ifdef VECTOR3D
	double Stheta = 2.0*M_PI*X_misc.get_rand_real1();
#ifdef DEBUG_INIT_SEGMENT_ANGLES
	cout << "initStheta=" << Stheta << '\n'; cout.flush();
#endif
	spatial S1(1.0,Stheta,X_misc.get_rand_range_real1(minangle,maxangle)), acoords;
	S1.convert_from_spherical();
	// Double rotation transformation T S, for phi and theta of the parent
	double cosphiSx(cosphi*S1.X()), sinphiSz(sinphi*S1.Z());
	//   Sx(in x,y,z) = cos(phi)cos(theta)Sx-sin(theta)Sy+sin(phi)cos(theta)Sz
	acoords.set_X((costheta*cosphiSx) - (sintheta*S1.Y()) + (costheta*sinphiSz));
	//   Sy(in x,y,z) = cos(phi)sin(theta)Sx+cos(theta)Sy+sin(phi)sin(theta)Sz
	acoords.set_Y((sintheta*cosphiSx) + (costheta*S1.Y()) + (sintheta*sinphiSz));
	//   Sz(in x,y,z) = -sin(phi)Sx+0Sy+cos(phi)Sz
	acoords.set_Z((cosphi*S1.Z()) - (sinphi*S1.X()));
	acoords.convert_to_spherical();
	acoords.set_X(radius);
#ifdef DEBUG_INIT_SEGMENT_ANGLES
	cout << "initTHETA=" << acoords.Y() << " initPHI=" << acoords.Z() << '\n'; cout.flush();
#endif
#endif
#ifdef VECTOR2D
	double anglespread = maxangle - minangle;
	double randomangle = anglespread*((2.0*X_misc.get_rand_real1()) - 1.0); // +/- random angle
	if (randomangle<0.0) randomangle -= minangle; // negative random between minangle and maxangle
	else randomangle += minangle; // positive random between minangle and maxangle
	spatial acoords(radius,rc.Y()+randomangle); // valid relative location on membrane
#ifdef DEBUG_INIT_SEGMENT_ANGLES
	cout << "initTHETA=" << acoords.Y() << '\n'; cout.flush();
#endif
#endif
	Segment initseg(P,P);
	initseg.set_by_angular_coordinates(acoords); // segment [P,P+Euler(acoords)]
	initseg.P0 = initseg.P1; // segment starts on membrane at P+Euler(acoords)
#ifdef DEBUG_INIT_SEGMENT_ANGLES
#ifdef VECTOR3D
	cout << "neuronCENTER=(" << P.X() << ',' << P.Y() << ',' << P.Z() << ')';
	cout << " segP0=(" << initseg.P0.X() << ',' << initseg.P0.Y() << ',' << initseg.P0.Z() << ")\n";
	cout.flush();
#endif
#endif
	acoords.set_X(X_misc.get_rand_range_real1(mintotlength,maxtotlength)); // initial segment length
	if (isinputstructure) inputstructure.link_before(new postsynaptic_structure(*this,initseg,acoords,_typeid));
	else outputstructure.link_before(new presynaptic_structure(*this,initseg,acoords,_typeid));
      }
    }
  }
}

void multipolar_nonpyramidal::initialize_multipolar_output_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle) {
  spatial rc(radius,0.0,0.0);
  if (random_orientation) rc.RandomAngular(radius,radius,0,2*M_PI);
  initialize_multiple_poles_common(mintotlength,maxtotlength,numpoles,minangle,maxangle,rc,0.0,false);
}

void multipolar_nonpyramidal::initialize_multipolar_input_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, int _typeid) {
  // _typeid may be used during postsynaptic_structure initialization.
  spatial rc(radius,0.0,0.0);
  if (!use_specified_basal_direction) {
    if (random_orientation) rc.RandomAngular(radius,radius,0,2*M_PI);
  } else {
    rc.set_X(1.0); // is changed in the called function anyway
    rc.set_Y(specified_basal_direction_theta);
#ifdef VECTOR3D
    rc.set_Z(specified_basal_direction_phi);
#endif
  }
  initialize_multiple_poles_common(mintotlength,maxtotlength,numpoles,minangle,maxangle,rc,0.0,true,_typeid);
}

void multipolar_nonpyramidal::initialize_multipolar_output_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev) {
  // regioncenter is NOT presumed to be relative or in spherical coordinates
  spatial rc(regioncenter);
  rc -= P; // relative vector from P to regioncenter
  rc.convert_to_spherical();
  initialize_multiple_poles_common(mintotlength,maxtotlength,numpoles,minangle,maxangle,rc,maxdev,false);
}

void multipolar_nonpyramidal::initialize_multipolar_input_structure(double mintotlength, double maxtotlength, int numpoles, double minangle, double maxangle, const spatial & regioncenter, double maxdev, int _typeid) {
  // _typeid may be used during postsynaptic_structure initialization.
  // regioncenter is NOT presumed to be relative or in spherical coordinates
  spatial rc(regioncenter);
  //if (!use_specified_basal_direction) {
  rc -= P; // relative vector from P to regioncenter
  rc.convert_to_spherical();
    /*  NOT HERE! A SPECIFIC REGIONCENTER IS GIVEN! } else {
    rc.set_X(1.0); // is changed in the called function anyway
    rc.set_Y(specified_basal_direction_theta);
    #ifdef VECTOR3D
    rc.set_Z(specified_basal_direction_phi);
    #endif
  }*/
  initialize_multiple_poles_common(mintotlength,maxtotlength,numpoles,minangle,maxangle,rc,maxdev,true,_typeid);
}

void multipolar_nonpyramidal::initialize_output_structure(double mintotlength, double maxtotlength) {
  int numpoles = X_misc.get_rand_range(multipolar_axon_mintrees,multipolar_axon_maxtrees);
  initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,2.0*M_PI);
}

void multipolar_nonpyramidal::initialize_input_structure(double mintotlength, double maxtotlength) {
  int numpoles = X_misc.get_rand_range(multipolar_mintrees,multipolar_maxtrees);
  initialize_multipolar_input_structure(mintotlength,maxtotlength,numpoles,0.0,2.0*M_PI);
}

void pyramidal::parse_CLP(Command_Line_Parameters & clp) {
  // Note: Facilities for setting the parameters of individual neurons, as well as of all pyramidal neurons can be implemented here.
  int n;
  // All pyramidal neurons:
  pyramidal_min_basal = PYRAMIDAL_MIN_BASAL;
  pyramidal_max_basal = PYRAMIDAL_MAX_BASAL;
  if ((n=clp.Specifies_Parameter("pyramidal.min_basal"))>=0) pyramidal_min_basal = atoi(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("pyramidal.max_basal"))>=0) pyramidal_max_basal = atoi(clp.ParValue(n));
  pyramidal_basal_minangle = PYRAMIDAL_BASAL_MINANGLE;
  pyramidal_basal_maxangle = PYRAMIDAL_BASAL_MAXANGLE;
  if ((n=clp.Specifies_Parameter("pyramidal.basal.minangle"))>=0) pyramidal_basal_minangle = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("pyramidal.basal.maxangle"))>=0) pyramidal_basal_maxangle = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("pyramidal.basal.force_model"))>=0) for (int i = 0; i < (int) NUM_bfm; i++) if (downcase(clp.ParValue(n))==bfmstr[i]) { bfm = (basal_force_model) i; break; }
  // Specific to this pyramidal neuron:
}

void pyramidal::initialize_output_structure(double mintotlength, double maxtotlength) {
  int numpoles = 1;
  double maxangle = M_PI/2.0;
  double maxdeviation = pyramidaltreesmaxdeviation;
  if (!random_orientation) {
    maxangle = random_orientation_amplitude;
    maxdeviation = random_orientation_amplitude;
  }
  if (inputstructure.tail()) {
    spatial cg;
    if (pia_attraction_repulsion_hypothesis) {
      cg = inputstructure.tail()->P0;
      cg -= P; // relative to neuron center
      cg.Negate();
      cg += P;
    } else {
      inputstructure.head()->center_of_gravity(cg);
    }
    initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,maxangle,cg,maxdeviation);
  } else initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,maxangle);
}

void pyramidal::initialize_input_structure(double mintotlength, double maxtotlength) {
  // The typical input structure should consist of a basal dendritic region
  // and an apical dendrite at the opposite side of the neuron. The axon
  // should project out of the center of the basal dendritic region.
  // Note that in network::develop_connection_structure this function is
  // called before the function for output_structure, since dgm is
  // initialized before agm.
  // The Apical Dendrite is initialized last. According to the comments and
  // code of initialize_multiple_poles_common(), this means that the Apical
  // Dendrite is the last dendritic arbor in the list of input structures
  // (this is accessed as the tail() of the list).
  int numbasal = X_misc.get_rand_range(pyramidal_min_basal,pyramidal_max_basal);
  double mintotapicallength = mintotlength / ((double) (numbasal+1));
  double mintotbasallength = mintotlength - mintotapicallength;
  double maxtotapicallength = maxtotlength / ((double) (numbasal+1));
  double maxtotbasallength = maxtotlength - maxtotapicallength;
  double maxangle = 0.1*M_PI;
  double maxdeviation = pyramidaltreesmaxdeviation;
  if (!random_orientation) {
    maxangle = random_orientation_amplitude;
    maxdeviation = random_orientation_amplitude;
  }
  if (outputstructure.tail()) { // axon was placed first
    spatial cg;
    outputstructure.head()->center_of_gravity(cg);
    // 1. initializing the basal dendritic region
    // note: mintotlength and maxtotlength are divided by the number of
    // basal dendrites to determine the initial length
    initialize_multipolar_input_structure(mintotbasallength,maxtotbasallength,numbasal,pyramidal_basal_minangle,pyramidal_basal_maxangle,cg,pyramidaltreesmaxdeviation);
    // 2. initializing the apical dendrite
    if (!pia_attraction_repulsion_hypothesis) inputstructure.head()->center_of_gravity(cg); // prioretize center of gravity of basal dendritic trees
    cg -= P;
    cg.Negate();
    cg += P;
    initialize_multipolar_input_structure(mintotapicallength,maxtotapicallength,1,0.0,maxangle,cg,maxdeviation,ALL_APICAL_PYRAMIDAL_DENDRITES_AEM);
  } else { // placing dendritic trees first
    // 1. initializing the basal dendritic region
    // note: mintotlength and maxtotlength are divided by the number of
    // basal dendrites to determine the initial length
    initialize_multipolar_input_structure(mintotbasallength,maxtotbasallength,numbasal,pyramidal_basal_minangle,pyramidal_basal_maxangle);
    // 2. initializing the apical dendrite
    spatial cg;
    inputstructure.head()->center_of_gravity(cg);
    cg -= P;
    cg.Negate();
    cg += P;
    initialize_multipolar_input_structure(mintotapicallength,maxtotapicallength,1,0.0,maxangle,cg,maxdeviation,ALL_APICAL_PYRAMIDAL_DENDRITES_AEM);
  }
}

void interneuron::initialize_output_structure(double mintotlength, double maxtotlength) {
  // 3. initializing the axon
  int numpoles = 1; // random_int(interneuron_axon_mintrees,interneuron_axon_maxtrees);
  if (inputstructure.tail()) {
    spatial cg;
    inputstructure.head()->center_of_gravity(cg);
    cg -= P;
    cg.Negate();
    cg += P;
    initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,0.1*M_PI,cg,pyramidaltreesmaxdeviation);
  } else initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,2.0*M_PI);
}

void interneuron::initialize_input_structure(double mintotlength, double maxtotlength) {
  int numpoles = X_misc.get_rand_range(2,4);
  initialize_multipolar_input_structure(mintotlength,maxtotlength,numpoles,0.0,2.0*M_PI);
}

Fig_Object * electrode::net_Fig() {
  if (!figattr_show_electrodes) return NULL;
  double x, y;
  P.plane_mapped(x,y);
  if (!fig_in_zoom(P)) return NULL;
  long X,Y,R;
  long fillstyle=-1; if (figattr_use_color) fillstyle=20;
  X=FIGSCALED(x); Y=FIGSCALED(y); R=FIGSCALED(radius);
  if (R==0) R = 1;
  return new Fig_Circle(0,1,0,3,999,fillstyle,0.0,0.0,X,Y,R,X,Y,X+150,Y+75);}
