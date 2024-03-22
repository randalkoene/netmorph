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
#include "network.hh"

//neuron * processing_neuron = NULL; // *** this cannot be global when parallel processing

const char neuron_type_name[][25] =
  {"Principal Neuron",
   "Interneuron",
   "Multipolar NonPyramidal",
   "Bipolar",
   "Pyramidal",
   "Untyped Neuron"};

const char neuron_short_name[][12] =
  {"principal",
   "interneuron",
   "multipolar",
   "bipolar",
   "pyramidal",
   "untyped"};

// [UPDATE 20080702:] The case in which we define INCLUDE_SCHEMA_PARENT_SET_PROTOCOL_DIRECTION_MODELS is now REQUIRED for compilation.

const char natural_subset_idstr[NUM_NATURAL_SPS][40] = {
  "",
  "all_axons.",
  "all_dendrites.",
  "all_multipolar_axons.",
  "all_multipolar_dendrites.",
  "all_bipolar_axons.",
  "all_bipolar_dendrites.",
  "all_pyramidal_axons.",
  "all_pyramidal_dendrites.",
  "all_interneuron_axons.",
  "all_interneuron_dendrites.",
  "all_apical_pyramidal_dendrites."};

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

general_neuron_parameters_interface general_neuron_parameters;

// Arrays containing compile-time default values
// principal, interneuron, multipolar, bipolar, pyramidal, untyped
int min_basal[UNTYPED_NEURON+1] = { 2, 2, 2, 1, 4, 2 };
int max_basal[UNTYPED_NEURON+1] = { 5, 4, 5, 1, 8, 5 };
double min_angle[UNTYPED_NEURON+1] = { 0.1*M_PI, 0.1*M_PI, 0.1*M_PI, 0.1*M_PI, 0.1*M_PI, 0.1*M_PI };
double max_angle[UNTYPED_NEURON+1] = { 0.5*M_PI, 2.0*M_PI, 2.0*M_PI, 0.25*M_PI, 0.5*M_PI, 0.5*M_PI };
int max_axons[UNTYPED_NEURON+1] = { 1, 1, 1, 1, 1, 1 };

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

void region_natural_set_model_error(region * r, natural_schema_parent_set natural_set_id, String modelstr) {
  // PROTOCOL: Do not allow sets with undefined model schemas.
  String setstr(natural_subset_idstr[natural_set_id]);
  if (natural_set_id==universal_sps) setstr = "universal";
  String regstr;
  if (!r) regstr = "universal network";
  else regstr = r->Name();
  error("Error: Region '"+regstr+"' natural set '"+setstr+"' has no associated "+modelstr+" model schema.\n");
}

int superior_schema_index(region * r, int schema_index) {
  if (!r) return superior_natural_set[schema_index];
  return schema_index; // same schema in universal network
}

void general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  error("Error: Please use the general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net) function in order to access network regions.\n");
}

void general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net) {
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
  if ((n=clp.Specifies_Parameter("apical_specified_overrides_centroid"))>=0) apical_specified_overrides_centroid = (downcase(clp.ParValue(n))==String("true"));

  // Population specific model selection (see TL#200807020156.12):
  cached_net_ptr = &net;
  String popid; // popid="" is the universal network model identifier
  // 1. We prepare schema arrays according to the number of regions defined
  int numschemasets = net.Regions().length()+1; // index 0 is for the universal network
  structure_initialization_region_subset_schemas = new region_structure_initialization_ptr[numschemasets];
  direction_model_region_subset_schemas = new region_direction_model_base_ptr[numschemasets];
  branching_model_region_subset_schemas = new region_branching_model_base_ptr[numschemasets];
  TSBM_region_subset_schemas = new region_TSBM_base_ptr[numschemasets];
  TSTM_region_subset_schemas = new region_TSTM_base_ptr[numschemasets];
  branch_angle_model_region_subset_schemas = new region_branch_angle_model_base_ptr[numschemasets];
  arbor_elongation_model_region_subset_schemas = new region_arbor_elongation_model_base_ptr[numschemasets];
  terminal_segment_elongation_model_region_subset_schemas = new region_terminal_segment_elongation_model_base_ptr[numschemasets];
  elongation_rate_initialization_model_region_subset_schemas = new region_elongation_rate_initialization_model_base_ptr[numschemasets];
  for (int i=0; i<numschemasets; i++) {
    structure_initialization_region_subset_schemas[i] = new structure_initialization[NUM_NATURAL_SPS]; // self-clearing constructor
    direction_model_region_subset_schemas[i] = new direction_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) direction_model_region_subset_schemas[i][j] = NULL;
    branching_model_region_subset_schemas[i] = new branching_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) branching_model_region_subset_schemas[i][j] = NULL;
    TSBM_region_subset_schemas[i] = new TSBM_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) TSBM_region_subset_schemas[i][j] = NULL;
    TSTM_region_subset_schemas[i] = new TSTM_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) TSTM_region_subset_schemas[i][j] = NULL;
    branch_angle_model_region_subset_schemas[i] = new branch_angle_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) branch_angle_model_region_subset_schemas[i][j] = NULL;
    arbor_elongation_model_region_subset_schemas[i] = new arbor_elongation_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) arbor_elongation_model_region_subset_schemas[i][j] = NULL;
    terminal_segment_elongation_model_region_subset_schemas[i] = new terminal_segment_elongation_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) terminal_segment_elongation_model_region_subset_schemas[i][j] = NULL;
    elongation_rate_initialization_model_region_subset_schemas[i] = new elongation_rate_initialization_model_base_ptr[NUM_NATURAL_SPS];
    for (int j=0; j<NUM_NATURAL_SPS; j++) elongation_rate_initialization_model_region_subset_schemas[i][j] = NULL;
  }
  //    Some initial values for structure initialization values in the universal set
  for (int j=0; j<NUM_NATURAL_SPS; j++) structure_initialization_region_subset_schemas[0][j].init(natural_schema_parent_set(j));
  // 2. We start with the universal network models, then work our way through existing populations
  region * r = NULL; int regid = 0;
  // [***INCOMPLETE] We may need to modify the way in which universal "contributing" models are specified and used in the selection functions.
  do {
    // General model selection (see TL#200603120422.1), and general parameters
    String genmodel; // genmodel="" is the base model identifier
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      genmodel = popid+natural_subset_idstr[i];
      if (!structure_initialization_selection(genmodel,clp,structure_initialization_region_subset_schemas[0][superior_schema_index(r,i)],structure_initialization_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"structure initialization");
      if (!direction_model_selection(genmodel,clp,direction_model_region_subset_schemas[0][superior_schema_index(r,i)],direction_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"direction");
      if (!branching_model_selection(genmodel,clp,branching_model_region_subset_schemas[0][superior_schema_index(r,i)],branching_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"branching");
      if (!TSBM_selection(genmodel,clp,TSBM_region_subset_schemas[0][superior_schema_index(r,i)],TSBM_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"TSBM");
      if (!TSTM_selection(genmodel,clp,TSTM_region_subset_schemas[0][superior_schema_index(r,i)],TSTM_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"TSTM");
      if (!branch_angle_model_selection(genmodel,clp,branch_angle_model_region_subset_schemas[0][superior_schema_index(r,i)],branch_angle_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"branch angle");
      if (!arbor_elongation_model_selection(genmodel,clp,arbor_elongation_model_region_subset_schemas[0][superior_schema_index(r,i)],arbor_elongation_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"arbor elongation");
      if (!terminal_segment_elongation_model_selection(genmodel,clp,terminal_segment_elongation_model_region_subset_schemas[0][superior_schema_index(r,i)],terminal_segment_elongation_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"terminal segment elongation");
      if (!elongation_rate_initialization_model_selection(genmodel,clp,elongation_rate_initialization_model_region_subset_schemas[0][superior_schema_index(r,i)],elongation_rate_initialization_model_region_subset_schemas[regid][i])) region_natural_set_model_error(r,i,"elongation rate initialization");
    }
    // Get the first region if we were working on the universal network, otherwise get the next region.
    if (popid.empty()) r = net.Regions().head();
    else r = r->Next();
    if (r) popid = r->Name()+'.';
    regid++;
  } while (r);

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

  if (!cached_net_ptr) error("Error: Please use the general_neuron_parameters_interface::parse_CLP(Command_Line_Parameters & clp, network & net) function in order to cache a pointer for network regions access.\n");
  int numschemasets = cached_net_ptr->Regions().length()+1;
  for (int regid = 0; regid<numschemasets; regid++) {
    res += "\nRegion '";
    if (regid<1) res += "universal network"; else res += cached_net_ptr->Regions().el(regid-1)->Name();
    res += "':";
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' structure initialization: ";
      res += structure_initialization_region_subset_schemas[regid][i].report_parameters();
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (arbor_elongation_model_region_subset_schemas[regid][i]) {
	res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' arbor elongation model: ";
	res += arbor_elongation_model_region_subset_schemas[regid][i]->report_parameters();
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (terminal_segment_elongation_model_region_subset_schemas[regid][i]) {
	res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' terminal segment elongation model: ";
	res += terminal_segment_elongation_model_region_subset_schemas[regid][i]->report_parameters();
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (elongation_rate_initialization_model_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal elongation rate initialization model: ";
	else { res += "\nSubset '"; res += natural_subset_idstr[i]; res += "' elongation rate initialization model: "; }
	res += elongation_rate_initialization_model_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_elongation_rate_initialization_model_root.empty()) res += " (contributing models with root "+universal_elongation_rate_initialization_model_root+".)";
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (direction_model_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal direction model: ";
	else { res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' direction model: "; }
	res += direction_model_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_direction_model_root.empty()) res += " (contributing models with root "+universal_direction_model_root+".)";
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (branching_model_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal branching model: ";
	else { res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' branching model: "; }
	res += branching_model_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_branching_model_root.empty()) res += " (contributing models with root "+universal_branching_model_root+".)";
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (TSBM_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal TSBM: ";
	else { res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' TSBM: "; }
	res += TSBM_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_TSBM_root.empty()) res += " (contributing models with root "+universal_TSBM_root+".)";
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (TSTM_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal TSTM: ";
	else { res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' TSTM: "; }
	res += TSTM_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_TSTM_root.empty()) res += " (contributing models with root "+universal_TSTM_root+".)";
      }
    }
    for (natural_schema_parent_set i = universal_sps; i<NUM_NATURAL_SPS; i = natural_schema_parent_set(i+1)) {
      if (branch_angle_model_region_subset_schemas[regid][i]) {
	if (i==universal_sps) res += "\nUniversal branch angle model: ";
	else { res += "\n  Subset '"; res += natural_subset_idstr[i]; res += "' branch angle model: "; }
	res += branch_angle_model_region_subset_schemas[regid][i]->report_parameters();
	if (i==universal_sps) if (!universal_branch_angle_model_root.empty()) res += " (contributing models with root "+universal_branch_angle_model_root+".)";
      }
    }
  }
  if (!elongation_rate_initialization_model_region_subset_schemas[0][universal_sps]) res += "\nUniversal elongation_rate_initialization model: length_distribution";
  if (!direction_model_region_subset_schemas[0][universal_sps]) res += "\nUniversal direction model: radial";
  if (!branching_model_region_subset_schemas[0][universal_sps]) res += "\nUniversal branching model: van_Pelt";
  if (!branch_angle_model_region_subset_schemas[0][universal_sps]) res += "\nUniversal branch angle model: balanced_forces";
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

void neuron::move(spatial & pos) {
  // Recenter a neuron on pos.
  spatial addvec(pos);
  addvec -= P;
  P = pos; // move the soma center
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) e->move_add(addvec);
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) e->move_add(addvec);
  PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
    PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) {
      s->move_add(addvec);
    }
  }
}

void neuron::fanin_rot() {
  net_fanin_histogram.collect=(statsattr_fan_in_analysis==dendrite_fs);
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) e->fanin_rot();
  net_fanin_histogram.collect=(statsattr_fan_in_analysis==axon_fs);
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) e->fanin_rot();
  PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
    PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) {
      s->fanin_rot();
    }
  }
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
    if (fig_report_depth()) progress("Reported Connections Depth = "+String((long) _figvar_connection_depth)+'\n');
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
  (*Txt_neuronlist) += ',';
  (*Txt_neuronlist) += neuron_short_name[TypeID()];
  (*Txt_neuronlist) += ',';
  int regnum = -1;
  if (eq) regnum = eq->Net()->Regions().find(*this);
  if (regnum>=0) (*Txt_neuronlist) += eq->Net()->Regions().el(regnum)->Name();
  else (*Txt_neuronlist) += "universal";
  double x=0.0, y=0.0, z=0.0;
  P.get_all(x,y,z);
  (*Txt_neuronlist) += String(x,",%f");
  (*Txt_neuronlist) += String(y,",%f");
  (*Txt_neuronlist) += String(z,",%f\n");
  progress(" synapses (");
  if (outputconnections.head()) {
    PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
      PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) s->net_Txt();
      progress("*");
    }
  }
  if (outattr_track_nodegenesis) {
    progress(") structure(");
    parsing_fs_type = dendrite_fs;
    PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) e->net_Txt();
    progress("dendrites,");
    parsing_fs_type = axon_fs;
    PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) e->net_Txt();
    progress("axons)\n");
  } else progress(")\n");
  Txt_neuronindex++;
  return NULL;
}

VRML_Object * neuron::net_VRML() {
  if (fig_in_zoom(P)) {
    (*VRML_neuronlist) += "<Transform translation='";
    double x=0.0, y=0.0, z=0.0;
    P.get_all(x,y,z);
    (*VRML_neuronlist) += String(x,"%.2f");
    (*VRML_neuronlist) += String(y," %.2f");
    (*VRML_neuronlist) += String(z," %.2f");
    (*VRML_neuronlist) += "'>\n<ProtoInstance name='Soma'>\n";
    (*VRML_neuronlist) += "</ProtoInstance>\n</Transform>\n";
  } else if (!figattr_box_fibre_independently) return NULL;
  // [***INCOMPLETE] Continue implementing below.
  //if (outputconnections.head()) {
  //  PLL_LOOP_FORWARD(connection,outputconnections.head(),1) {
  //    PLL_LOOP_FORWARD_NESTED(synapse,e->Synapses()->head(),1,s) s->net_VRML();
  //  }
  //}
  parsing_fs_type = dendrite_fs;
  PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) e->net_VRML();
  parsing_fs_type = axon_fs;
  PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) e->net_VRML();
  VRML_neuronindex++;
  return NULL;
}

#ifdef VECTOR3D
void neuron::net_Slice(Slice * slice) { // [***INCOMPLETE] Change the return type.
  // Efficiency: Determine which spatial segments are contained in or intersect with a slice. Determine which fiber segments
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
  // B. Collect presynaptic and postsynaptic fiber objects.
  Fig_Group * postsynaptic_structures_fig = NULL;
  if (inputstructure.head()) {
    postsynaptic_structures_fig = new Fig_Group();
    PLL_LOOP_FORWARD(fibre_structure,inputstructure.head(),1) postsynaptic_structures_fig->Add_Fig_Object(e->net_Slice());
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
  Fig_Group * presynaptic_structures_fig = NULL;
  if (inputstructure.head()) {
    presynaptic_structures_fig = new Fig_Group();
    PLL_LOOP_FORWARD(fibre_structure,outputstructure.head(),1) presynaptic_structures_fig->Add_Fig_Object(e->net_Slice());
    _figvar_connection_depth--; if (_figvar_connection_depth<10) _figvar_connection_depth = 997;
  }
#endif
  // C. Collect synapse objects.
  // Combine all the objects and return them.
}
#endif

void principal::parse_CLP(Command_Line_Parameters & clp) {
  // Note: Facilities for setting the parameters of individual neurons, as well as of all pyramidal neurons can be implemented here.
  // [***INCOMPLETE] Ideally, this should be settable per region label and per structure, i.e. basal dendrites, apical dendrites and axons.
  int n;
  neuron_type nt = this->TypeID();
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".min_basal"))>=0) min_basal[nt] = atoi(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".max_basal"))>=0) max_basal[nt] = atoi(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".basal.minangle"))>=0) min_angle[nt] = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".basal.maxangle"))>=0) max_angle[nt] = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".basal.force_model"))>=0) for (int i = 0; i < (int) NUM_bfm; i++) if (downcase(clp.ParValue(n))==bfmstr[i]) { bfm = (basal_force_model) i; break; }
  // Specific to this neuron:
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
	if ((!acceptinit) && (cnt>500)) {
	  warning("Warning: Unable to place arbors at desired minimum angular separation.\n  Affected neuron type: "+String(neuron_type_name[TypeID()])+", neuron ID:"+String((long) numerical_ID())+", angle constraints: ["+String(minangle,"%.3f,")+String(maxangle,"%.3f")+"], number of arbors: "+String((long) numpoles)+".\n");
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
      warning("Warning: Surface division BFM is not yet implemented for 2D\n");
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
#ifdef DIVIDE_INITLENGTHS_OVER_ARBORS
      double dnumpoles = (double) numpoles;
      mintotlength /= dnumpoles; maxtotlength /= dnumpoles;
#endif
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

/*! \brief Initialization of dendrite root nodes and root fiber segments.

    This function calls the general root node placement function initialize_multiple_poles_common(),
    which uses the minangle and maxangle parameters when allocating 2 or more fiber structure roots,
    but uses a single maxdev parameter when allocating one fiber structure root near a regioncenter
    vector rc. In this case, the maxdev parameter is 0.0 in the call, which means that any single
    dendrite (e.g. an apical dendrite) will be placed along a vector rc that is random or specified.
    In the case of an apical dendrite, this particular function is only called if the flags
    apical_specified_overrides_centroid and use_specified_basal_direction are both set. Then the
    vector rc is negated, so that the apical dendrite is placed exaclty along the opposite direction
    of the specified basal direction. This option removes all variability from the initial direction
    of apical root fiber segments.    
 */
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
    if (_typeid==ALL_APICAL_PYRAMIDAL_DENDRITES_AEM) { // you only end up here if apical_specified_overrides_centroid==true
      rc.convert_from_spherical();
      rc.Negate();
      rc.convert_to_spherical();
    }
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
  int numpoles = X_misc.get_rand_range(1,max_axons[TypeID()]);
  if (inputstructure.tail()) {
    spatial cg;
    inputstructure.head()->center_of_gravity(cg);
    cg -= P;
    cg.Negate();
    cg += P;
    initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,0.1*M_PI,cg,multipolartreesmaxdeviation);
  } else initialize_multipolar_output_structure(mintotlength,maxtotlength,numpoles,0.0,2.0*M_PI);
}

void multipolar_nonpyramidal::initialize_input_structure(double mintotlength, double maxtotlength) {
  int numpoles = X_misc.get_rand_range(min_basal[TypeID()],max_basal[TypeID()]);
  initialize_multipolar_input_structure(mintotlength,maxtotlength,numpoles,min_angle[TypeID()],max_angle[TypeID()]);
}

void multipolar_nonpyramidal::parse_CLP(Command_Line_Parameters & clp) {
  // Note: Facilities for setting the parameters of individual neurons, as well as of all pyramidal neurons can be implemented here.
  // [***INCOMPLETE] Ideally, this should be settable per region label and per structure, i.e. basal dendrites, apical dendrites and axons.
  principal::parse_CLP(clp);
  int n;
  neuron_type nt = TypeID();
  if ((n=clp.Specifies_Parameter(String(neuron_short_name[nt])+".max_axons"))>=0) max_axons[nt] = atoi(clp.ParValue(n));
  // Specific to this neuron:
}

void bipolar::initialize_output_structure(double mintotlength, double maxtotlength) {
  int numpoles = 1;
  double maxangle = M_PI/2.0;
  double maxdeviation = bipolartreesmaxdeviation;
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

void bipolar::initialize_input_structure(double mintotlength, double maxtotlength) {
  // A bipolar cell typically has one axon and one dendrite extending at opposite
  // ends.
  int numbasal = X_misc.get_rand_range(min_basal[TypeID()],max_basal[TypeID()]);
#ifdef DIVIDE_INITLENGTHS_OVER_ARBORS
  mintotlength /= (double) (numbasal+1);
  maxtotlength /= (double) (numbasal+1);
#endif
  double maxangle = 0.1*M_PI;
  double maxdeviation = bipolartreesmaxdeviation;
  if (!random_orientation) {
    maxangle = random_orientation_amplitude;
    maxdeviation = random_orientation_amplitude;
  }
  if (outputstructure.tail()) { // axon was placed first
    spatial cg;
    outputstructure.head()->center_of_gravity(cg);
    initialize_multipolar_input_structure(mintotlength,maxtotlength,numbasal,min_angle[TypeID()],max_angle[TypeID()],cg,bipolartreesmaxdeviation);
  } else { // placing dendritic trees first
    initialize_multipolar_input_structure(mintotlength,maxtotlength,numbasal,min_angle[TypeID()],max_angle[TypeID()]);
  }
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
  int numbasal = X_misc.get_rand_range(min_basal[TypeID()],max_basal[TypeID()]);
#ifdef DIVIDE_INITLENGTHS_OVER_ARBORS
  double mintotapicallength = mintotlength / ((double) (numbasal+1));
  double mintotbasallength = mintotlength - mintotapicallength;
  double maxtotapicallength = maxtotlength / ((double) (numbasal+1));
  double maxtotbasallength = maxtotlength - maxtotapicallength;
#else
  double mintotapicallength = mintotlength;
  double mintotbasallength = mintotlength;
  double maxtotapicallength = maxtotlength;
  double maxtotbasallength = maxtotlength;
#endif
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
    initialize_multipolar_input_structure(mintotbasallength,maxtotbasallength,numbasal,min_angle[TypeID()],max_angle[TypeID()],cg,pyramidaltreesmaxdeviation);
    // 2. initializing the apical dendrite
    if (!pia_attraction_repulsion_hypothesis) inputstructure.head()->center_of_gravity(cg); // prioretize center of gravity of basal dendritic trees
    cg -= P;
    cg.Negate();
    cg += P;
    initialize_multipolar_input_structure(mintotapicallength,maxtotapicallength,1,0.0,maxangle,cg,maxdeviation,ALL_APICAL_PYRAMIDAL_DENDRITES_AEM);
  } else { // placing dendritic trees first
    // 1. initializing the basal dendritic region
    initialize_multipolar_input_structure(mintotbasallength,maxtotbasallength,numbasal,min_angle[TypeID()],max_angle[TypeID()]);
    // 2. initializing the apical dendrite
    if ((use_specified_basal_direction) && (apical_specified_overrides_centroid)) initialize_multipolar_input_structure(mintotapicallength,maxtotapicallength,1,0.0,maxangle,ALL_APICAL_PYRAMIDAL_DENDRITES_AEM);
    else {
      spatial cg;
      inputstructure.head()->center_of_gravity(cg);
      cg -= P;
      cg.Negate();
      cg += P;
      initialize_multipolar_input_structure(mintotapicallength,maxtotapicallength,1,0.0,maxangle,cg,maxdeviation,ALL_APICAL_PYRAMIDAL_DENDRITES_AEM);
    }
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
  int numpoles = X_misc.get_rand_range(min_basal[TypeID()],max_basal[TypeID()]);
  initialize_multipolar_input_structure(mintotlength,maxtotlength,numpoles,min_angle[TypeID()],max_angle[TypeID()]);
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
