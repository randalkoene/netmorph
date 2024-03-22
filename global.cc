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
// global.cc
// Randal A. Koene, 20041118

#include <math.h>
#include <time.h>
#include "global.hh"
#include "file.hh"
#include "Fig_Object.hh"

#define STRVECBUFSIZE 256

#ifdef TESTING_AEM_RANDOM_DISTRIBUTION
unsigned long int aemrandomdist[21] = {0,0,0,0,0,0,0,0,0,0,0};
#endif

unsigned int random_seed = 0;
bool outattr_show_progress = false;
bool outattr_show_figure = true;
bool outattr_show_stats = true;
bool outattr_track_synaptogenesis = false;
bool outattr_track_nodegenesis = false;
bool outattr_make_full_Txt = false;
bool outattr_Txt_sequence = false;
bool outattr_Txt_separate_files = false;
bool outattr_make_full_X3D = false;
bool outattr_make_full_Catacomb = false;
fibre_structure_class_id statsattr_fan_in_analysis = NUM_fs; // default, implies "none"
bool outattr_synapse_distance_frequency = false;
bool outattr_connection_distance_frequency = false;
double outattr_distance_frequency_distbinsize = 50.0;

bool figattr_show_electrodes = false;
bool figattr_show_neurons = true;
bool figattr_show_abstract_connections = false;
bool figattr_show_presynaptic_structure = true;
bool figattr_show_postsynaptic_structure = true;
bool figattr_show_synapse_structure = true;
bool figattr_show_spatial_segment_subsets = false;
bool figattr_show_connection_evaluation = false;
bool figattr_use_color = false;
bool figattr_fill_somas = true;
bool figattr_update_terminal_segments_visibly = false;
bool figattr_zoom = false;
bool figattr_fibres_nobox = false;
bool figattr_box_fibre_independently = true;
bool figattr_make_full_Fig = false;
bool figattr_make_zoom_Fig = false;
bool figattr_make_connections_Fig = false;
bool figattr_make_abstract_Fig = false;
bool figattr_make_neurons_Figs = false;
int figattr_sample_show_progress = 0; // 1 = text, 2 = progress bar
bool figattr_show_scale = true;
bool figattr_show_axis_arrows = true;

double figattr_focus_box_x1 = 1.0;
double figattr_focus_box_y1 = 0.0;
double figattr_focus_box_x2 = -1.0;
double figattr_focus_box_y2 = 0.0;
double figattr_focus_box_z1 = 0.0;
double figattr_focus_box_z2 = 0.0;

int _figvar_connection_depth = 997;
long _figvar_connection_axon_x = 0;
long _figvar_connection_axon_y = 0;
long _figvar_connection_color = 0;

Fig_Group * _figgroup = NULL;

double max_growth_time = (20.0*24.0*60.0*60.0);
bool reduce_rand_calls = false;

bool random_orientation = true;
double random_orientation_amplitude = 0.2*M_PI;
bool apical_specified_overrides_centroid = false;
bool use_specified_basal_direction = false;
double specified_basal_direction_theta = 0.0;
double specified_basal_direction_phi = 0.0;
double multipolar_initangledevratio = 0.3;
double multipolar_axon_initangledevratio = 0.3;

double pyramidaltreesmaxdeviation = 0.1*M_PI; // = M_PI/2.0; // variation around designated positions for the axon and apical dendrite // *** this should be in Pyramidal_statistics or in the pspreadmin/pspreadmax initialization
double bipolartreesmaxdeviation = 0.1*M_PI;
double multipolartreesmaxdeviation = 0.1*M_PI;

double branchanglemin = M_PI/8.0; // *** this should be in Multipolar_Statistics
double branchanglemax = M_PI/4.0; // *** this should be in Multipolar_Statistics
double branchangleothermax = 3.0*(M_PI/16.0); // *** this should be in Multipolar_Statistics
double branchIIanglemin = M_PI/3.0; // *** this should be in Multipolar_Statistics
double branchIIanglemax = M_PI/2.0; // *** this should be in Multipolar_Statistics
double branchIIangleothermax = 0.2*M_PI; // *** this should be in Multipolar_Statistics
double branchIIprobability = 0.5; // *** this should be in Multipolar_Statistics
bool dendrite_two_branch_types = true; // *** this should be in Multipolar_Statistics
bool axon_two_branch_types = true; // *** this should be in Multipolar_Statistics
double turnanglemin = M_PI/16.0; // *** this should be in Multipolar_Statistics
double turnanglemax = M_PI/4.0; // *** this should be in Multipolar_Statistics

bool fibreswithturns = true;
bool initiallengthatactualfirstbranch = false;
bool branchsomewhereinsegment = true;

int spatialsegmentsubsetsizelimit = 65;
unsigned int maxnumberofspatialsegmentsubsets = INT_MAX;
double minimum_span_spatialsegmentsubset = 2.0; // (default 2 micron)
double specified_max_spatial_segment_coordinate = -1.0; 

bool evaluate_connectivity_immediately = true; // (see remark in nibr.hh)

unsigned long warning_fibre_segment_net_Fig_zero_length = 0;
unsigned long warning_Spatial_Segment_Subset_intersects_zero_length = 0;
int warnings_on = WARN_FILE;
int reports_on = WARN_FILE;
int progress_on = WARN_STDOUT;

#ifdef TESTQUOTAS
String quotas;
#endif

// A superior that is the same set means the set has no superior, which can only
// be true for the universal set.
natural_schema_parent_set superior_natural_set[NUM_NATURAL_SPS] = {
  universal_sps,
  universal_sps,
  universal_sps,
  all_axons_sps,
  all_dendrites_sps,
  all_axons_sps,
  all_dendrites_sps,
  all_axons_sps,
  all_dendrites_sps,
  all_axons_sps,
  all_dendrites_sps,
  all_pyramidal_dendrites_sps
};

// These are the random number generators.
// Their seed values depend on each other, so that simulation results
// are deterministically reproducible.
rng X_misc((unsigned int) time(NULL));
rng X_elongate(X_misc.get_rand());
rng X_branch(X_elongate.get_rand());
rng X_turn(X_branch.get_rand());
rng X_synapses(X_turn.get_rand());

String outputdirectory("");
String absoluteoutputURL("http://rak.minduploading.org/nibr/output/"); // default when called from a form

String warningfile("");
String reportfile("");
String progressfile("");

// functions

void error(String msg) {
  cerr << msg;
  exit(1);
}

void warning(String msg) {
  if (warnings_on<=WARN_OFF) return;
  if ((warnings_on<WARN_FILE) || (warningfile.empty())) { cerr << msg; cerr.flush(); }
  if ((warnings_on>WARN_STDOUT) && (!warningfile.empty())) append_file_from_String(warningfile,msg);
}

void report(String msg) {
  if (reports_on<=WARN_OFF) return;
  if ((reports_on<WARN_FILE) || (reportfile.empty())) { cout << msg; cout.flush(); }
  if ((reports_on>WARN_STDOUT) && (!reportfile.empty())) append_file_from_String(reportfile,msg);
}

void progress(String msg) {
  if (progress_on<=WARN_OFF) return;
  if ((progress_on<WARN_FILE) || (progressfile.empty())) { cout << msg; cout.flush(); }
  if ((progress_on>WARN_STDOUT) && (!progressfile.empty())) append_file_from_String(progressfile,msg);
}

void global_parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("random_orientation"))>=0) {
    random_orientation = (downcase(clp.ParValue(n))==String("true"));
    if (!random_orientation) random_orientation_amplitude = atof(clp.ParValue(n));
  }
  if ((n=clp.Specifies_Parameter("turnanglemin"))>=0) {
    turnanglemin = atof(clp.ParValue(n));
    if (turnanglemin<0.0) error("Error: (global_parse_CLP) turnanglemin must be >= 0.0\n");
  }
  if ((n=clp.Specifies_Parameter("turnanglemax"))>=0) {
    turnanglemax = atof(clp.ParValue(n));
    if (turnanglemax<0.0) error("Error: (global_parse_CLP) turnanglemax must be > 0.0\n");
  }
  if (turnanglemax<=turnanglemin) error("Error: (global_parse_CLP) turnanglemax must be > turnanglemin\n");
  if ((n=clp.Specifies_Parameter("dendrite_two_branch_types"))>=0) dendrite_two_branch_types = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("axon_two_branch_types"))>=0) axon_two_branch_types = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_Fig_rescale"))>=0) figattr_Fig_rescale = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_fibres_nobox"))>=0) figattr_fibres_nobox = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_box_fibre_independently"))>=0) figattr_box_fibre_independently = (downcase(clp.ParValue(n))==String("true"));
  if (figattr_box_fibre_independently) figattr_fibres_nobox=false; // for sanity
  if ((n=clp.Specifies_Parameter("figattr_fill_somas"))>=0) figattr_fill_somas = (downcase(clp.ParValue(n))==String("true"));  
  if ((n=clp.Specifies_Parameter("figattr_show_scale"))>=0) figattr_show_scale = (downcase(clp.ParValue(n))==String("true"));  
  if ((n=clp.Specifies_Parameter("figattr_show_axis_arrows"))>=0) figattr_show_axis_arrows = (downcase(clp.ParValue(n))==String("true"));  
  if ((n=clp.Specifies_Parameter("partition_to_synapses"))>=0) evaluate_connectivity_immediately = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("max_spatial_segment_subsets"))>=0) maxnumberofspatialsegmentsubsets = atoi(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("min_spatial_segment_span"))>=0) minimum_span_spatialsegmentsubset = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("max_spatial_segment_coordinate"))>=0) specified_max_spatial_segment_coordinate = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("outattr_synapse_distance_frequency"))>=0) outattr_synapse_distance_frequency = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_connection_distance_frequency"))>=0) outattr_connection_distance_frequency = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_distance_frequency_distbinsize"))>=0) outattr_distance_frequency_distbinsize = atof(clp.ParValue(n));
}

String global_report_parameters() {
  String res("Global parameters:\n  Initial somatic locations of dendrites and axons are ");
  if (random_orientation) res += "randomly oriented\n";
  else res += "not randomly oriented, but perturbed with an angular amplitude of "+String(random_orientation_amplitude,"%.3f\n");
  res += "  turnanglemin = " + String(turnanglemin,"%.3f turnanglemax = ") + String(turnanglemax,"%.3f\n") + "  dendrite branches are of";
  if (dendrite_two_branch_types) res += " two types\n  axon branches are of";
  else res += " one type\n  axon branches are of";
  if (axon_two_branch_types) res += " two types\n";
  else res += " one type\n  ";
  if (!figattr_Fig_rescale) res +="do note ";
  res += "rescale final Fig output coordinates to match specified figureTeXwidth\n";
  if (!figattr_box_fibre_independently) res += "do not ";
  res += "show fibres of neurons outside the focus box\n  ";
  if (!figattr_fibres_nobox) res += "do not ";
  res += "test if fibres are outside the focus box\n";
  if (!figattr_fill_somas) res += "do not ";
  res += "draw filled somata\n";
  if (!figattr_show_scale) res += "do not ";
  res += "draw a scale bar\n";
  if (!figattr_show_axis_arrows) res += "do not ";
  res += "draw x-, y-, and z- axis arrows\n";
  if (warnings_on) res += "Global NETMORPH parameters:\n  Messages classified as NETMORPH warnings are displayed.\n";
  else res += "Global NETMORPH parameters:\n  Messages classified as NETMORPH warnings are NOT displayed.\n";
  if (evaluate_connectivity_immediately) res += "Global synapse parameters:\n  Candidate synapses are sought as soon as fibers are allocated to spatial segments.\n";
  else res += "Global synapse parameters:\n  Candidate synapses are not sought as soon as fibers are allocated to spatial segments.\n";
  res += "Global Spatial Segment Subset parameters:\n  Maximum number of spatial segment subsets = " + String((long) maxnumberofspatialsegmentsubsets) + "\n";
  res += "  Minimum span of a Spatial Segment Subset (requested) = " + String(minimum_span_spatialsegmentsubset,"%.2f micrometers\n");
  if (specified_max_spatial_segment_coordinate>0.0) res += "  Spatial Segment Subset maximum coordinates explicitly set to " + String(specified_max_spatial_segment_coordinate,"%.2f micrometers\n");
  if (outattr_synapse_distance_frequency) res += "Produce frequency by distance data for synapses.\n";
  if (outattr_connection_distance_frequency) res += "Produce frequency by distance data for connections.\n";
  if (outattr_synapse_distance_frequency || outattr_connection_distance_frequency) res += "With distance bin size " + String(outattr_distance_frequency_distbinsize,"%.2f") + " micrometers.\n"; 
  return res;
}

void set_random_seed(unsigned int seed) {
  random_seed = seed;
  report("random seed set to "+String((long) random_seed)+'\n');
  X_misc.init(seed);
  X_elongate.init(X_misc.get_rand());
  X_branch.init(X_elongate.get_rand());
  X_turn.init(X_branch.get_rand());
  X_synapses.init(X_turn.get_rand());
  srand(seed); // [***INCOMPLETE] remove this
}

void set_time_random_seed() {
  set_random_seed((unsigned int) time(NULL));
}

/*double random_double(double rmin, double rmax) {
  return ((((double) rand())/((double) RAND_MAX))*(rmax-rmin)) + rmin;
}

double random_double() {
  // [***NOTE] This is always a positive number between 0.0 and 1.0.
  return ((double) rand())/((double) RAND_MAX);
}

int random_int(int rmin, int rmax) {
  return (rand() % ((rmax-rmin)+1)) + rmin;
}

int random_int(int rmax) {
  return rand() % (rmax+1);
  }*/

#ifdef INCLUDE_PDF_SAMPLING
#define PDF_SAMPLE_RESOLUTION 25
#define PDF_SAMPLE_ARRAYSIZE ((2*PDF_SAMPLE_RESOLUTION)+1)
unsigned long * probability_distribution_function::sample(unsigned long * samplearray, double s) {
  int i;
  if (!samplearray) {
    samplearray = new unsigned long[PDF_SAMPLE_ARRAYSIZE];
    for (i = 0; i<PDF_SAMPLE_ARRAYSIZE; i++) samplearray[i] = 0;
  }
  double granularity = samplemax/(double) PDF_SAMPLE_RESOLUTION;
  if (abs(s)>samplemax) return samplearray;
  if (s>=0.0) i = (int) ((s/granularity) + 0.5);
  else i = (int) ((s/granularity) - 0.5);
  samplearray[i+PDF_SAMPLE_RESOLUTION]++;
  return samplearray;
}
String probability_distribution_function::sample_histogram_data() {
  String res("x = [ ");
  double granularity = samplemax/(double) PDF_SAMPLE_RESOLUTION;
  double min_x = -(granularity * (double) PDF_SAMPLE_RESOLUTION);
  double x = min_x;
  for (int i = 0; i<PDF_SAMPLE_ARRAYSIZE; i++) { res += String(x,"%f "); x += granularity; }
  res += "];\n";
  if (!positivesamples) res += "PDF_pos = [];\n";
  else {
    res += "PDF_pos = [ ";
    for (int i = 0; i<PDF_SAMPLE_ARRAYSIZE; i++) res += String(positivesamples[i],"%f ");
    res += "];\n";
  }
  if (!selectionsamples) res += "PDF = [];\n";
  else {
    res += "PDF = [ ";
    for (int i = 0; i<PDF_SAMPLE_ARRAYSIZE; i++) res += String(selectionsamples[i],"%f ");
    res += "];\n";
  }
  return res;
}
#endif

delta_pdf::delta_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni), v(0.0) {
  int n;
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  if ((n=clp.Specifies_Parameter(thislabel+"value"))>=0) v = atof(clp.ParValue(n));
}

linear_pdf::linear_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni) {
  double max_x = 1.0, b = 1.0;
  int n;
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  if ((n=clp.Specifies_Parameter(thislabel+"max_x"))>=0) max_x = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"height_b"))>=0) b = atof(clp.ParValue(n));
  setup(max_x,b);
}

String linear_pdf::report_local_parameters() {
  String res(2.0*max_y/-negb,"max_x=%.3f");
  return res + String(-negb," b=%.3f");
}

String linear_pdf::report_parameters() {
  return String("linear PDF: "+report_local_parameters());
}

void spline_normal_pdf::setup(double max_x, double sigvalue, double sigproportion) {
  _max_x = max_x;
  _sigvalue = sigvalue;
  _sigproportion = sigproportion;
  double l2 = max_x-sigvalue, b1, b2;
  if (l2<=0.0) b2 = 0.0;
  else b2 = (2.0*sigproportion)/l2;
  if (sigvalue<=0.0) b1 = b2;
  else b1 = (2.0*(((1.0-sigproportion)/sigvalue)-b2)) + b2;
  // determiny parameters in the domain of the cdf
  y_sig = sigvalue*((0.5*(b2-b1)) + (b1+b2));
  a1 = (b2-b1)/sigvalue; // *** BEWARE
  a2 = -b2/l2; // *** BEWARE
  C = (1.0-sigproportion) + (0.5*a2*sigvalue*sigvalue) - (b2*sigvalue);
  max_y = C + (0.5*a2*max_x*max_x) + ((b2-sigvalue)*max_x);
  // cached values: max_y, y_sig, -(b1+b2), (b1+b2)^2, 2a1, a1,
  //                -(b2-sigvalue), (b2-sigvalue)^2, 2a2, C, a2
  negsumb1b2 = -(b1+b2);
  sqsumb1b2 = negsumb1b2*negsumb1b2;
  twoa1 = 2.0*a1;
  negdifb2sigvalue = -(b2-sigvalue);
  sqdifb2sigvalue = negdifb2sigvalue*negdifb2sigvalue;
  twoa2 = 2.0*a2;
  twomax_y = 2.0*max_y;
}

spline_normal_pdf::spline_normal_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni) {
  double max_x = 1.0, sigvalue = 0.25, sigproportion = 0.1;
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"max_x"))>=0) max_x = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"significance_threshold"))>=0) sigvalue = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"proportion_significant"))>=0) {
    double sigprop = atof(clp.ParValue(n));
    if ((sigprop<0.0) || (sigprop>1.0)) warning("Warning: In spline_normal_pdf the proportion_significant must be between 0 and 1, setting to 0.1 instead of "+clp.ParValue(n)+".\n");
    else sigproportion = sigprop;
  }
  setup(max_x,sigvalue,sigproportion);
}

String spline_normal_pdf::report_local_parameters() {
  String res(_max_x,"max_x=%.3f ");
  res += String(_sigvalue,"sigvalue=%.3f ");
  res += String(_sigproportion,"sigproportion=%.3f");
  return res;
}

String spline_normal_pdf::report_parameters() {
  return String("spline_normal_pdf: "+report_local_parameters());
}

int splinetest = 0;
double spline_normal_pdf::random_positive() {
  // implement the solution of step 4 in the paper notes
  double y = max_y*X_uni->get_rand_real1();
  switch (splinetest) {
  case 0:
    if (y<=y_sig) return (negsumb1b2+sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return (negdifb2sigvalue+sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
    break;;
  case 1:
    if (y<=y_sig) return (negsumb1b2-sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return (negdifb2sigvalue-sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
    break;;
  case 2:
    if (y<=y_sig) return (negsumb1b2+sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return (negdifb2sigvalue-sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
    break;;
  default:
    if (y<=y_sig) return (negsumb1b2-sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return (negdifb2sigvalue+sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
  }
  return 0.0;
}

double spline_normal_pdf::random_selection() {
  double y = twomax_y*X_uni->get_rand_real1();
  if (y>max_y) {
    return random_positive();
    //    y -= max_y;
    //    if (y<=y_sig) return (negsumb1b2+sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    //    else return (negdifb2sigvalue+sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
    
  } else {
    return -random_positive();
    //    if (y<=y_sig) return -(negsumb1b2+sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    //    else return -(negdifb2sigvalue+sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
  }
}

/*
double spline_normal_pdf::random_positive() {
  // implement the solution of step 4 in the paper notes
  double y = max_y*X_uni->get_rand_real1();
  if (y<=y_sig) return (negsumb1b2+-sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
  else return (negdifb2sigvalue+-sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
}

double spline_normal_pdf::random_selection() {
  double y = twomax_y*X_uni->get_rand_real1();
  if (y>max_y) {
    y -= max_y;
    if (y<=y_sig) return (negsumb1b2+-sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return (negdifb2sigvalue+-sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
    
  } else {
    if (y<=y_sig) return -(negsumb1b2+-sqrt(sqsumb1b2 + twoa1*y))/a1; // *** BEWARE
    else return -(negdifb2sigvalue+-sqrt(sqdifb2sigvalue + twoa2*(y-C)))/a2; // *** BEWARE
  }
}
*/

spline_normal_pdf_with_min::spline_normal_pdf_with_min(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): spline_normal_pdf(thislabel,clp,_X_uni), min_x(0.0) {
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  int n;
  if ((n=clp.Specifies_Parameter(thislabel+"min_x"))>=0) min_x = atof(clp.ParValue(n));
}

String spline_normal_pdf_with_min::report_local_parameters() {
  return spline_normal_pdf::report_local_parameters() + String(min_x," min_x=%.3f");
}

String spline_normal_pdf_with_min::report_parameters() {
  return String("spline_normal_pdf_with_min: "+report_local_parameters());
}

normal_pdf::normal_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni), mean(0.0), std(1.0), trunc(DBL_MAX) {
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  int n;
  //if (thislabel==String("all_axons.eri.PDF.")) cout << "LOOKING FOR SPECIFIC NORMAL PARAMETERS ";
  if ((n=clp.Specifies_Parameter(thislabel+"mean"))>=0) mean = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"std"))>=0) std = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(thislabel+"trunc"))>=0) trunc = atof(clp.ParValue(n));
  //if (thislabel==String("all_axons.eri.PDF.")) cout << "mean = " << mean << " std = " << std << " trunc = " << trunc << '\n';
  zigset();
  set_mean_std(mean,std,trunc);
}

// [***BEWARE] I might have to make the definition below explicitly "(signed)".
#define SHR3 (X_uni->get_rand())
#define UNI (X_uni->get_rand_real3())

float normal_pdf::nfix(void) { /*provides RNOR if #define cannot */
  const float r = 3.442620f; static float x, y;
  for(;;){
    x=hz*wn[iz];
    if(iz==0) {
      do {
	x=-log(UNI)*0.2904764;
	y=-log(UNI);
      } while(y+y<x*x);
      return (hz>0)? r+x : -r-x;
    }
    if ( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
    hz=SHR3;
    iz=hz&127;
    if (((unsigned) abs(hz))<kn[iz]) return (hz*wn[iz]);
  }
}

// Obligatory global initialization of static member variables
uint32_t normal_pdf::kn[128];
float normal_pdf::wn[128];
float normal_pdf::fn[128];
bool normal_pdf::tablesinitialized = false;

void normal_pdf::zigset() {
  // This function prepares the static tables kn[], wn[], and fn[] used by all instances of normal_pdf.
  if (tablesinitialized) return;
  const double m1 = 2147483648.0;
  double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
  int i;
  /* Tables for RNOR: */
  q=vn/exp(-.5*dn*dn);
  kn[0]=(unsigned) ((dn/q)*m1);
  kn[1]=0;
  wn[0]=q/m1;
  wn[127]=dn/m1;
  fn[0]=1.;
  fn[127]=exp(-.5*dn*dn);
  for(i=126;i>=1;i--) {
    dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
    kn[i+1]=(unsigned) ((dn/tn)*m1);
    tn=dn;
    fn[i]=exp(-.5*dn*dn);
    wn[i]=dn/m1;
  }
  tablesinitialized = true; // do this initialization only once
}

double normal_pdf::random_positive() {
#ifdef INCLUDE_PDF_SAMPLING
  double X = abs(random_selection());
  positivesamples=sample(positivesamples,X);
  return X;
#else
  return abs(random_selection());
#endif
}

double normal_pdf::random_selection() {
#ifdef INCLUDE_PDF_SAMPLING
  double X = random_selection_calculation();
  selectionsamples=sample(selectionsamples,X);
  return X;
}
double normal_pdf::random_selection_calculation() {
#endif
  hz=SHR3;
  iz=hz&127;
  if (standard) return ( (((unsigned) abs(hz))<kn[iz])? hz*wn[iz] : nfix() );
  else for (;;) {
    register double r = ( (((unsigned) abs(hz))<kn[iz])? hz*wn[iz] : nfix() )*std + mean;
    if (abs(r)<=trunc) return r;
    hz=SHR3;
    iz=hz&127;
  }
}

String normal_pdf::report_local_parameters() {
  String res(mean,"mean=%.5f ");
  res += String(std,"std=%.5f ");
  res += String(trunc,"trunc=%.5f");
  return res;
}

String normal_pdf::report_parameters() {
  return String("normal_pdf: "+report_local_parameters());
}

exponential_pdf::exponential_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni) {
  zigset();
}

float exponential_pdf::efix(void) { /*provides REXP if #define cannot */
  float x;
  for(;;) {
    if(iz==0) return (7.69711-log(UNI));
    x=jz*we[iz];
    if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);
    jz=SHR3;
    iz=(jz&255);
    if(jz<ke[iz]) return (jz*we[iz]);
  }
}

/*--------This procedure sets the seed and creates the tables------*/
void exponential_pdf::zigset() {
  const double m2 = 4294967296.;
  double q;
  double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
  int i;
  /* Tables for REXP */
  q = ve/exp(-de);
  ke[0]=(unsigned) ((de/q)*m2); ke[1]=0;
  we[0]=q/m2; we[255]=de/m2;
  fe[0]=1.; fe[255]=exp(-de);
  for(i=254;i>=1;i--) {
    de=-log(ve/de+exp(-de));
    ke[i+1]= (unsigned) ((de/te)*m2); te=de;
    fe[i]=exp(-de); we[i]=de/m2;
  }
}

double exponential_pdf::random_positive() {
  return abs(random_selection());
}

double exponential_pdf::random_selection() {
  jz = SHR3;
  iz=jz&255;
  return ( (jz <ke[iz])? jz*we[iz] : efix() );
}

String exponential_pdf::report_local_parameters() {
  return String();
}

String exponential_pdf::report_parameters() {
  return String("exponential_pdf: "+report_local_parameters());
}

discinv_pdf::discinv_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): probability_distribution_function(_X_uni), p(0), q(0), r(0), n(0), r_min(0.0), r_max(1.0), r_mean(0.0), r_std(0.0) {
  if (!thislabel.empty()) if (thislabel.lastchar()!='.') thislabel += '.';
  int n_cmd;
  int N = 0;
  if ((n_cmd=clp.Specifies_Parameter(thislabel+"N"))>=0) N = clp.get_int(n_cmd,1,STRVECBUFSIZE,"Warning: The Discrete PDF requires 1<N<=256.\n");
  double * P = 0; int P_len = 0;
  if ((n_cmd=clp.Specifies_Parameter(thislabel+"P"))>=0) P = str2vec(clp.ParValue(n),&P_len);
  if (P_len<N) {
    warning("Warning: List of P values shorter than N="+String((long) N)+" in discinv_pdf::discinv_pdf(), setting N="+String((long) P_len)+".\n");
    N = P_len;
  }
  if (P_len>N) warning("Warning: Only using N="+String((long) N)+" values from list of P values with length "+String((long) P_len)+"in discinv_pdf::discinv_pdf().\n");
  double R_min = 0.0, R_max = 1.0;
  if ((n_cmd=clp.Specifies_Parameter(thislabel+"R_min"))>=0) R_min = atof(clp.ParValue(n));
  if ((n_cmd=clp.Specifies_Parameter(thislabel+"R_max"))>=0) R_max = atof(clp.ParValue(n));
  if (N>0) set_discrete_pdf(P,N,R_min,R_max);
  if (P) delete[] P;
}

void discinv_pdf::set_discrete_pdf(double * P, int N, double R_min, double R_max) {
  if ((!P) || (N<=0)) return;
  if (p) delete[] p;
  if (q) delete[] q;
  if (r) delete[] r;
  p = new double[N];
  q = new double[N];
  r = new double[N];
  n = N;
  r_min = R_min;
  r_max = R_max;
  double P_sum = 0.0;
  for (int i = 0; i<n; i++) P_sum += P[i];
  double p_sum = 0.0;
  for (int i = 0; i<n; i++) { p[i] = P[i]/P_sum; p_sum += p[i]; q[i] = p_sum; }
  r[0] = r_min;
  if (n>1) {
    double r_step = (r_max-r_min)/((double) (n-1));
    for (int i = 1; i<n; i++) r[i] = r[i-1] + r_step;
    r[n-1] = r_max; // just to make sure it is exact
  }
  r_mean = 0.0;
  double E_rsq = 0.0;
  for (int i = 0; i<n; i++) {
    double rp = r[i]*p[i];
    r_mean += rp; // calculate the discrete mean
    E_rsq += r[i]*rp; // calculate the expectation of r squared
  }
  double var_r = E_rsq - (r_mean*r_mean); // calculate the variance
  if (var_r<0) error("Error: Negative variance computed in discinv_pdf().\n");
  r_std = sqrt(var_r);
}

double discinv_pdf::random_positive() {
  return abs(random_selection());
}

double discinv_pdf::random_selection() {
  if (n<1) return r_min;
  double u = X_uni->get_rand_real1();
  for (int i=0; i<n; i++) if (q[i]>=u) return r[i];
  return r_max;
}

String discinv_pdf::report_local_parameters() {
  String res("N="+String((long) n));
  res += ", R_min="+String(r_min,"%.2f");
  res += ", R_max="+String(r_max,"%.2f");
  res += ", P value list with discrete mean "+String(r_mean,"%.2f");
  res += " and standard deviation "+String(r_std,"%.2f");
  return res;
}

String discinv_pdf::report_parameters() {
  return String("discinv_pdf: "+report_local_parameters());
}

probability_distribution_function * pdfselection(String basecommand, Command_Line_Parameters & clp, rng & _X_uni) {
  // This function aids in the selection of a specific PDF class.
  // basecommand specifies that command that is used to specify a PDF class, e.g. a growth model.
  // The label "[.]PDF" is appended to basecommand to find the corresponding command line parameter.
  // Further parameters of individual PDF classes are then specified with the prefix "[basecommand.]PDF.".
  // E.g. synapse_formation.PDF=delta; synapse_formation.PDF.value=0.0;
  if (basecommand.empty()) basecommand = "PDF";
  else if (basecommand.lastchar()=='.') basecommand += "PDF";
  else basecommand += ".PDF";
  //if (basecommand==String("all_axons.eri.PDF")) cout << "REGISTERED PDF COMMAND\n";
  int n;
  if ((n=clp.Specifies_Parameter(basecommand))>=0) {
    String PDFstr(downcase(clp.ParValue(n)));
    if (PDFstr.empty()) {
      warning("Warning: Empty value specified for "+basecommand+'\n');
      return NULL;
    } else {
      if (PDFstr.matches("delta")) return new delta_pdf(basecommand,clp,_X_uni);
      else if (PDFstr.matches("uniform")) return new uniform_pdf(_X_uni);
      else if (PDFstr.matches("linear")) return new linear_pdf(basecommand,clp,_X_uni);
      else if (PDFstr.matches("spline_normal")) return new spline_normal_pdf(basecommand,clp,_X_uni);
      else if (PDFstr.matches("spline_normal_with_min")) return new spline_normal_pdf_with_min(basecommand,clp,_X_uni);
      else if (PDFstr.matches("normal")) return new normal_pdf(basecommand,clp,_X_uni);
      else if (PDFstr.matches("exponential")) return new exponential_pdf(basecommand,clp,_X_uni);
      else if (PDFstr.matches("discrete")) return new discinv_pdf(basecommand,clp,_X_uni);
    }
    warning("Warning: Unrecognized "+basecommand+" class "+clp.ParValue(n)+'\n');
    return NULL;
  }
  return NULL; // No warning is given if no PDF specification is found
}

double * str2vec(const char * vecstr, int * veclen) {
  // This function takes a string and looks for a list of numbers to
  // return in a newly allocated array of doubles. The numbers may
  // be separated by any non-numerical characters. Consequently,
  // this function does not detect "INF" or "NAN".
  double tmpbuf[STRVECBUFSIZE];
  int numfound = 0;
  const char * vecptr = vecstr;
  char * endptr = NULL;
  if (veclen) *veclen = 0;
  for (int i = 0; i<STRVECBUFSIZE; i++) {
    tmpbuf[i] = strtod(vecptr,&endptr);
    if (endptr==vecptr) break;
    numfound++;
    while ((*endptr!='\0') && (*endptr!='+') && (*endptr!='-') && ((*endptr<'0') || (*endptr>'9'))) endptr++;
    if (*endptr=='\0') break;
    vecptr=endptr;
  }
  if (numfound<1) return NULL;
  double * vec = new double[numfound];
  for (int i = 0; i<numfound; i++) vec[i] = tmpbuf[i];
  *veclen = numfound;
  return vec;
}
