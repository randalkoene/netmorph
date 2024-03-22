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
// global.hh
// Randal A. Koene, 20040830
//
// Miscellaneous global definitions and declarations.

#ifndef __GLOBAL_HH
#define __GLOBAL_HH

#include <cmath>
#include "BigString.hh"
#include "Command_Line_Parameters.hh"
#include "mtprng.hh"
#include <stdint.h>

// global defines

//typedef unsigned long uint32_t;
//typedef long int32_t;

//#define TESTING_AEM_RANDOM_DISTRIBUTION
#ifdef TESTING_AEM_RANDOM_DISTRIBUTION
extern unsigned long int aemrandomdist[];
#endif

// Testing new code
#define STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION

#define FIGSCALED(figcoord) ((long) (figscale*(figcoord)))

#define INV_FIGSCALED(scaledcoord) (((double) scaledcoord)/figscale)

// bit pattens for neuron::figattr flags
// cleared bit = visible
#define FIGATTR_SHOWALL     0x00
#define FIGMASK_SHOWALL     0x0F
#define FIGATTR_SHOWPRESYN  0x0E
#define FIGMASK_SHOWPRESYN  0x01
#define FIGATTR_SHOWPOSTSYN 0x0D
#define FIGMASK_SHOWPOSTSYN 0x02
#define FIGATTR_SHOWCELL    0x0B
#define FIGMASK_SHOWCELL    0x04
#define FIGATTR_SHOWNOCONN  0x0F

// set bit = visible
#define FIGATTR_REPORTDEPTH 0x10
#define FIGMASK_REPORTDEPTH 0x10

// default physiological parameters
#define CONNECTION_STATISTICS_MAX_FIBRE_DISTANCE_DEFAULT 0.1

// global type enumerations

enum growth_cone_competition_type { gcc_within_tree, gcc_within_all_AorD, gcc_whole_neuron, gcc_NUM };

enum fibre_structure_class_id { axon_fs, dendrite_fs, NUM_fs };

enum branching_model { van_Pelt_bm, polynomial_O1_bm, polynomial_O3_bm, NUM_bm };
enum TSBM { van_Pelt_tsbm, van_Pelt_specBM_tsbm, NUM_tsbm };
enum TSTM { tm_none, tm_each_dt, tm_linear_rate, tm_branch_coupled, tm_NUM };
enum branch_angle_model { balanced_forces_bam, fork_main_daughter_bam, NUM_bam };
enum turning_model { independent_rate_tm, elongation_locked_tm, bifurcation_coupled_tm, NUM_tm };
enum direction_model { radial_dm, segment_history_tension_dm, cell_attraction_dm, vector_dm, NUM_dm };
enum arbor_elongation_model { van_Pelt_aem, polynomial_O1_aem, polynomial_O3_aem, NUM_aem };
enum terminal_segment_elongation_model { simple_tsem, inertia_tsem, second_order_tsem, constrained_second_order_tsem, initialized_CSO_tsem, decaying_second_order_tsem, BESTL_tsem, nonnorm_BESTL_tsem, pyrAD_BESTLNN_tsem, NUM_tsem };
enum elongation_rate_initialization_model { length_distribution_eri, nonnorm_BESTL_length_distribution_eri, pure_stochastic_eri, zero_eri, unitary_eri, continue_defaults_eri, NUM_eri };

// PROTOCOL: [See TL#200709121534.] As for all defined sets, in the case of natural sets also, parsing specifications must occur for more general
//           sets before proceeding to more specific (narrow) sets. This insures that the more general model selections are available for
//           inheritance when such a selection is not made explicitly for a more specific set. Note that this is similar to object inheritance
//           in languages such as C++, and the order of constructor calls.
enum natural_schema_parent_set {
  universal_sps,
  all_axons_sps,
  all_dendrites_sps,
  all_multipolar_axons_sps,
  all_multipolar_dendrites_sps,
  all_bipolar_axons_sps,
  all_bipolar_dendrites_sps,
  all_pyramidal_axons_sps,
  all_pyramidal_dendrites_sps,
  all_interneuron_axons_sps,
  all_interneuron_dendrites_sps,
  all_apical_pyramidal_dendrites_sps,
  NUM_NATURAL_SPS };
extern natural_schema_parent_set superior_natural_set[NUM_NATURAL_SPS];

// global variables

extern unsigned int random_seed;
extern bool outattr_show_progress;
extern bool outattr_show_figure;
extern bool outattr_show_stats;
extern bool outattr_track_synaptogenesis;
extern bool outattr_track_nodegenesis;
extern bool outattr_make_full_Txt;
extern bool outattr_make_full_X3D;
extern bool outattr_make_full_Catacomb;
extern fibre_structure_class_id statsattr_fan_in_analysis;
extern bool outattr_Txt_sequence;
extern bool outattr_Txt_separate_files;
extern bool outattr_synapse_distance_frequency;
extern bool outattr_connection_distance_frequency;
extern double outattr_distance_frequency_distbinsize;
extern bool figattr_show_electrodes;
extern bool figattr_show_neurons;
extern bool figattr_show_abstract_connections;
extern bool figattr_show_presynaptic_structure;
extern bool figattr_show_postsynaptic_structure;
extern bool figattr_show_synapse_structure;
extern bool figattr_show_spatial_segment_subsets;
extern bool figattr_show_connection_evaluation;
extern bool figattr_use_color;
extern bool figattr_fill_somas;
extern bool figattr_update_terminal_segments_visibly;
extern bool figattr_zoom; // set by set_focus_box in nibr.cc
extern bool figattr_fibres_nobox;
extern bool figattr_box_fibre_independently;
extern bool figattr_make_full_Fig;
extern bool figattr_make_zoom_Fig;
extern bool figattr_make_connections_Fig;
extern bool figattr_make_abstract_Fig;
extern bool figattr_make_neurons_Figs;
extern int figattr_sample_show_progress;
extern bool figattr_show_scale;
extern bool figattr_show_axis_arrows;

extern int _figvar_connection_depth;
extern long _figvar_connection_axon_x;
extern long _figvar_connection_axon_y;
extern long _figvar_connection_color;

class Fig_Group;
extern Fig_Group * _figgroup;

class Sampled_Growth_Output;
extern Sampled_Growth_Output * sampled_output; // originally in Sampled_Output.hh

extern double max_growth_time;
extern bool reduce_rand_calls;

// === begin parameter section ===
// The following parameter section should be moved into a Multipolar_Statistics
// object so that each neuron can link to its own specific instance.
extern bool random_orientation;
extern double random_orientation_amplitude;
extern bool apical_specified_overrides_centroid;
extern bool use_specified_basal_direction;
extern double specified_basal_direction_theta;
extern double specified_basal_direction_phi;
extern double multipolar_initangledevratio; // *** this should be in Multipolar_Statistics
extern double multipolar_axon_initangledevratio; // *** this should be in Multipolar_Statistics

extern double pyramidaltreesmaxdeviation; // *** this should be in Pyramidal_statistics or in the pspreadmin/pspreadmax initialization
extern double bipolartreesmaxdeviation;
extern double multipolartreesmaxdeviation;

extern double branchanglemin; // *** this should be in Multipolar_Statistics
extern double branchanglemax; // *** this should be in Multipolar_Statistics
extern double branchangleothermax; // *** this should be in Multipolar_Statistics
extern double branchIIanglemin; // *** this should be in Multipolar_Statistics
extern double branchIIanglemax; // *** this should be in Multipolar_Statistics
extern double branchIIangleothermax; // *** this should be in Multipolar_Statistics
extern double branchIIprobability; // *** this should be in Multipolar_Statistics
extern bool dendrite_two_branch_types; // *** this should be in Multipolar_Statistics
extern bool axon_two_branch_types; // *** this should be in Multipolar_Statistics
extern double turnanglemin; // *** this should be in Multipolar_Statistics
extern double turnanglemax; // *** this should be in Multipolar_Statistics
// === end parameter section ===

extern bool fibreswithturns;
extern bool initiallengthatactualfirstbranch;
extern bool branchsomewhereinsegment;

extern int spatialsegmentsubsetsizelimit;
extern unsigned int maxnumberofspatialsegmentsubsets;
extern double minimum_span_spatialsegmentsubset; // (default 2 micron)
extern double specified_max_spatial_segment_coordinate; 

extern bool evaluate_connectivity_immediately; // ** this could be in dendritic_growth_model, accessed via net->agm (if net is stored in Spatial_Segment_Subset and agm in network)

extern unsigned long warning_fibre_segment_net_Fig_zero_length;
extern unsigned long warning_Spatial_Segment_Subset_intersects_zero_length;
#define WARN_OFF        0
#define WARN_STDOUT     1
#define WARN_STDOUTFILE 2
#define WARN_FILE       3
extern int warnings_on; // OFF, STDOUT, STDOUT & FILE, FILE
extern int reports_on; // OFF, STDOUT, STDOUT & FILE, FILE
extern int progress_on; // OFF, STDOUT, STDOUT & FILE, FILE

extern String outputdirectory;
extern String absoluteoutputURL;

extern String warningfile;
extern String reportfile;
extern String progressfile;

//#define TESTQUOTAS
#ifdef TESTQUOTAS
extern String quotas;
#endif

// function declarations

void error(String msg);
void warning(String msg);
void report(String msg);
void progress(String msg);

void global_parse_CLP(Command_Line_Parameters & clp);
String global_report_parameters();

// PRECISE_RNG or FAST_RNG
#ifdef PRECISE_RNG

#define rng libcoyote::mtprng

#else
#ifdef FAST_RNG

#define rng shiftrng
class shiftrng { // see mtprng for an explanation of the member functions
public:
  //! Defines the 32-bit integer type.
  typedef uint32_t int_type;  // needs to be a 32-bit type
  //typedef unsigned long int int_type;  // needs to be a 32-bit type
  //! Defines the 64-bit IEEE-754 type.
  typedef double        real_type; // needs to be a 64-bit or larger type
private:
  // Working storage
  int_type m_seed;
  int_type jsr;
  //int_type m_mt[N];
  size_t   m_mti;
  int_type m_multiplier;
public:
  shiftrng(int_type seed = 19650218UL) { init(seed); }
  void init(int_type seed) { m_seed = seed; jsr = seed; }
  int_type get_seed() { return m_seed; }
  int_type get_max() { return (int_type) 4294967295UL; }
  int_type get_rand() {
    register int_type jz = jsr;
    jsr^=(jsr<<13);
    jsr^=(jsr>>17);
    jsr^=(jsr<<5);
    return jz+jsr;
  }
  int_type get_rand_range(int_type lo, int_type hi) {
    real_type range = hi - lo + 1.0;
    return lo + int_type(floor(range * get_rand_real2()));
  }
  real_type get_rand_real1() {
    // privides a granularity of approx. 2.3E-10
    return real_type(get_rand()) * (1.0 / 4294967295.0);
  }
  real_type get_rand_range_real1(real_type lo, real_type hi) {
    return range_scale(get_rand_real1(),lo,hi);
  }
  real_type get_rand_real2() {
    // privides a granularity of approx. 2.3E-10
    return real_type(get_rand()) * (1.0 / 4294967296.0);
  }
  real_type get_rand_real3() {
    // privides a granularity of approx. 2.3E-10
    return real_type((double(get_rand()) + 0.5) * (1.0 / 4294967296.0));
  }
  real_type get_rand_real53() { return get_rand_real2(); }
};
#else

#define rng gnucrng
class gnucrng { // see mtprng for an explanation of the member functions
public:
  //! Defines the 32-bit integer type.
  //typedef uint32_t int_type;  // needs to be a 32-bit type
  //typedef unsigned long int int_type;  // needs to be a 32-bit type
  //! Defines the 64-bit IEEE-754 type.
  typedef double        real_type; // needs to be a 64-bit or larger type
private:
  // Working storage
  unsigned int m_seed;
public:
  shiftrng(unsigned int seed = 19650218U) { init(seed); }
  void init(unsigned int seed) { m_seed = seed; srand(seed); }
  unsigned int get_seed() { return m_seed; }
  unsigned int get_max() { return (unsigned int) RAND_MAX; }
  unsigned int get_rand() { return rand(); }
  unsigned int get_rand_range(int_type lo, int_type hi) {
    real_type range = hi - lo + 1.0;
    return lo + ((unsigned int) (floor(range * get_rand_real2())));
  }
  real_type get_rand_real1() {
    return real_type(get_rand())/real_type(RAND_MAX);
  }
  real_type get_rand_range_real1(real_type lo, real_type hi) {
    return range_scale(get_rand_real1(),lo,hi);
  }
  real_type get_rand_real2() {
    return real_type(get_rand())/(real_type(RAND_MAX)+1.0);
  }
  real_type get_rand_real3() {
    return real_type(double(get_rand()) + 0.5)/(real_type(RAND_MAX)+1.0);
  }
  real_type get_rand_real53() { return get_rand_real2(); }
};

#endif
#endif

void set_random_seed(unsigned int seed);
void set_time_random_seed();

//double random_double(double rmin, double rmax);
//double random_double();
//int random_int(int rmin, int rmax);
//int random_int(int rmax);

double square_distance(double x1, double y1, double x2, double y2);

// class declarations

class probability_distribution_function {
  // This is the base class for a variety of probability distribution
  // functions.
  // There are certain functions that each PDF must support:
  //   mean_X()
  //   (std_X())
  //   random_positive()
  //   random_selection()
  //   amplitude()
private:
  probability_distribution_function() { error("Error: Default constructor probability_distribution_function() should never be called.\n"); }
protected:
  rng * X_uni;
#ifdef INCLUDE_PDF_SAMPLING
  unsigned long * selectionsamples;
  unsigned long * positivesamples;
  double samplemax;
  unsigned long * sample(unsigned long * samplearray, double s);
#endif
public:
#ifdef INCLUDE_PDF_SAMPLING
  probability_distribution_function(rng & _X_uni): X_uni(&_X_uni), selectionsamples(NULL), positivesamples(NULL), samplemax(1.0) {}
  probability_distribution_function(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): X_uni(&_X_uni), selectionsamples(NULL), positivesamples(NULL), samplemax(1.0) {}
  virtual ~probability_distribution_function() { if (selectionsamples) delete[] selectionsamples; if (positivesamples) delete[] positivesamples; }
#else
  probability_distribution_function(rng & _X_uni): X_uni(&_X_uni) {}
  probability_distribution_function(String thislabel, Command_Line_Parameters & clp, rng & _X_uni): X_uni(&_X_uni) {}
  virtual ~probability_distribution_function() {}
#endif
  virtual double mean_X() = 0;
  virtual double random_positive() = 0;
  virtual double random_selection() = 0;
  virtual double amplitude() = 0;
  virtual String report_parameters() = 0;
  virtual probability_distribution_function * clone() = 0;
#ifdef INCLUDE_PDF_SAMPLING
  void set_sample_max(double x) { samplemax = x; }
  String sample_histogram_data();
#endif
};

class delta_pdf: public probability_distribution_function {
  // This PDF always returns the same value when
  // random_positive() or random_selection() is called.
protected:
  double v;
public:
  delta_pdf(rng & _X_uni, double _v = 0.0): probability_distribution_function(_X_uni), v(_v) {}
  delta_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  virtual double mean_X() { return v; }
  virtual double random_positive() { if (v>=0.0) return v; else return -v; }
  virtual double random_selection() { return v; }
  virtual double amplitude() { return 0.0; } // I corrected this from 1.0.
  virtual String report_parameters() { return String("delta PDF: "+String(v,"value=%.3f")); }
  virtual probability_distribution_function * clone() { return new delta_pdf(*X_uni,v); }
};

class uniform_pdf: public probability_distribution_function {
public:
  uniform_pdf(rng & _X_uni): probability_distribution_function(_X_uni) {}
  virtual double mean_X() { return 0.0; }
  virtual double random_positive() { return X_uni->get_rand_real1(); }
  virtual double random_selection() { return X_uni->get_rand_range_real1(-1.0,1.0); }
  virtual double amplitude() { return 1.0; }
  virtual String report_parameters() { return String("uniform PDF"); }
  virtual probability_distribution_function * clone() { return new uniform_pdf(*X_uni); }
};

class linear_pdf: public probability_distribution_function {
  // A linear or triangular probability distribution function with a mean
  // at 0.0. For random selection, the maximum is set to max_x, while the
  // height of the probability distribution at the mean is given by b.
protected:
  double a, max_y, negb, SQb, twoa, twomax_y;
public:
  linear_pdf(rng & _X_uni, double max_x = 1.0, double b = 1.0): probability_distribution_function(_X_uni) { setup(max_x,b); }
  linear_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  virtual void setup(double max_x = 1.0, double b = 1.0) {
    if (max_x<=0.0) error("Error: max_x must be >= 0.0 in linear_pdf\n");
    a = -b/max_x;
    max_y = 0.5*b*max_x;
    negb = -b;
    SQb = b*b;
    twoa = 2.0*a;
    twomax_y = 2.0*max_y;    
  }
  virtual double mean_X() { return 0.0; }
  virtual double random_positive() {
    double y = max_y*X_uni->get_rand_real1();
    return (negb + sqrt(SQb + (twoa*y)))/a;
  }
  virtual double random_selection() {
    double y = twomax_y*X_uni->get_rand_real1();
    if (y>max_y) {
      y -= max_y;
      return (negb + sqrt(SQb + (twoa*y)))/a;
    } else return -((negb + sqrt(SQb + (twoa*y)))/a);
  }
  virtual double amplitude() { return negb/a; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new linear_pdf(*X_uni,amplitude(),-negb); }
};

class spline_normal_pdf: public probability_distribution_function {
  // A spline-based approximation of a normal distribution function by linear
  // sections with a mean at 0.0.
  // Parameters:
  // sigvalue is the value beyond which the slowly sloping linear approximation
  // is used (e.g. for significant turns).
  // sigproportion is the proportion of randomly selected values that should
  // fall in the range beyond sigvalue.
  // max_x is the maximum value that random selection can return.
  // Note: If you wish to have a range from 0.0 to a value min_x for which the
  // probability is zero, then (a) supply sigvalue = actualsigvalue - min_x,
  // (b) supply max_x = actualmax_x - min_x, and (c) after obtaining a value
  // with random_selection(), either add or subtract min_x, depending on the
  // sign of the random selection.
protected:
  double a1, a2, y_sig, max_y, C, negsumb1b2, sqsumb1b2, twoa1, negdifb2sigvalue, sqdifb2sigvalue, twoa2, twomax_y;
  double _max_x, _sigvalue, _sigproportion;
public:
  spline_normal_pdf(rng & _X_uni, double max_x, double sigvalue, double sigproportion): probability_distribution_function(_X_uni) { setup(max_x,sigvalue,sigproportion); }
  spline_normal_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  virtual void setup(double max_x, double sigvalue, double sigproportion);
  virtual double mean_X() { return 0.0; }
  virtual double random_positive();
  virtual double random_selection();
  virtual double amplitude() { return _max_x; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new spline_normal_pdf(*X_uni,_max_x,_sigvalue,_sigproportion); }
};

class spline_normal_pdf_with_min: public spline_normal_pdf {
  // This is the same as spline_normal_pdf, except it takes care of all the
  // necessary conversions noted above when a min_x that defines a range
  // from 0 to min_x in which the probabilities are 0 during random selection.
protected:
  double min_x;
public:
  spline_normal_pdf_with_min(rng & _X_uni, double _min_x, double max_x, double sigvalue, double sigproportion): spline_normal_pdf(_X_uni,max_x-_min_x,sigvalue-_min_x,sigproportion), min_x(_min_x) {}
  spline_normal_pdf_with_min(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  virtual double mean_X() { return 0.0; }
  virtual double random_positive() { return min_x + spline_normal_pdf::random_positive(); }
  virtual double random_selection() { double y = spline_normal_pdf::random_selection(); if (y<0) return y-min_x; return y+min_x; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new spline_normal_pdf_with_min(*X_uni,min_x,_max_x,_sigvalue,_sigproportion); }
};

class normal_pdf: public probability_distribution_function {
  // A random generator with normal distribution and default mean 0.0 and std 1.0.
  // Parameters:
  //   mean = the mean of the normal distribution, default 0.0
  //   std  = the standard deviation of the normal distribution, default 1.0
  //   trunc = the truncation constraint of output random variables, default DBL_MAX
/* This is a slightly modified version of the Ziggurat Normal and Exponential
   Random Number generator published by Marsaglia and Tsang [MARSAGLIA;ZIGGURAT].

   The short period problem that was noted in published comments on the
   Ziggurat method is solved here by using the Mersenne Twister random
   number generator. */
  // [***NOTE] This Ziggurat method of generating random variables with a normal
  // distribution can be further improved by switching to doubles and using 1.25
  // uniform random variables to decide an x,y point, as described in 
  // Jurgen A. Doornik's "An Improved Ziggurat Method to Generate Normal Random Samples",
  // presently (20070704) available at http://www.doornik.com/research/ziggurat.pdf.
  // The improvement with a better uniform random number generator mentioned there
  // and by P.H.W. Leong et al. (2005) is already done with the mtprng code when
  // PRECISE_RNG is set.
  // Usage example:
  // 1. Somewhere in the preceding code declare a random variables generator:
  //    rng X;
  // 2. Declare a normal distributed random variables generator that uses X:
  //    normal_PDF X_norm(X);
  // All instances of PDFs that use X will sample random variables from the
  // same series of generated random variables.
protected: // static tables (Note: These do not add to the memory footprint of each instance of this class.)
  static uint32_t kn[128];
  static float wn[128];
  static float fn[128];
  static bool tablesinitialized;
protected:
  uint32_t iz;
  int32_t hz;
  double mean, std, trunc;
  bool standard;
  float nfix(void);
  void zigset();
#ifdef INCLUDE_PDF_SAMPLING
  double random_selection_calculation();
#endif
public:
  normal_pdf(rng & _X_uni): probability_distribution_function(_X_uni), mean(0.0), std(1.0), trunc(DBL_MAX), standard(true) { zigset(); }
  normal_pdf(rng & _X_uni, double _mean, double _std, double _trunc = DBL_MAX): probability_distribution_function(_X_uni) { zigset(); set_mean_std(_mean,_std,_trunc); }
  normal_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  double mean_X() { return mean; }
  double std_X() { return std; }
  void set_mean_std(double _mean, double _std, double _trunc = DBL_MAX) { mean=_mean; std=_std; trunc=_trunc; standard=((mean==0.0) && (std==1.0) && (trunc==DBL_MAX)); }
  virtual double random_positive();
  virtual double random_selection();
  virtual double amplitude() { return trunc-mean; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new normal_pdf(*X_uni,mean,std,trunc); }
};

class exponential_pdf: public probability_distribution_function {
  // A random generator with exponential distribution of default mean 0.0.
  // Parameters:
  //   ...to be defined...
/* This is a slightly modified version of the Ziggurat Normal and Exponential
   Random Number generator published by Marsaglia and Tsang [MARSAGLIA:ZIGGURAT].

   The short period problem that was noted in published comments on the
   Ziggurat method is solved here by using the Mersenne Twister random
   number generator. */
protected:
  rng * X_uni;
  uint32_t iz,ke[256];
  int32_t hz;
  uint32_t jz;
  //unsigned long int jz;
  //unsigned long int iz,ke[256];
  //long int hz;
  float we[256],fe[256];
  float efix(void);
  void zigset();
public:
  exponential_pdf(rng & _X_uni): probability_distribution_function(_X_uni) { zigset(); }
  exponential_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  virtual double mean_X() { return 0.0; }
  virtual double random_positive();
  virtual double random_selection();
  virtual double amplitude() { return 1.0; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new exponential_pdf(*X_uni); }
};

class discinv_pdf: public probability_distribution_function {
  // A random generator with a distribution that is obtained by using the inverse
  // of a discrete cdf.
  // Parameters:
  //   P = a vector of values
  //   N = the length of p
  //   R_min = lowest output value in the random value output range
  //   R_max = highest output value in the random value output range
  // This is done as explained in ~/octave/arbrnd.m, i.e.:
  // 1. provide discrete P
  // 2. compute normalized discrete p
  // 3. compute and store discrete cumulative q=cumsum(P)
  // 4. compute the output conversion vector r
  // 5. get a random variable u between 0 and 1
  // 6. walk through q to find the first value >= u
  // 7. at that index, look into r, return that value
protected:
  double * p;
  double * q;
  double * r;
  int n;
  double r_min, r_max;
  double r_mean, r_std;
#ifdef INCLUDE_PDF_SAMPLING
  double random_selection_calculation();
#endif
public:
  discinv_pdf(rng & _X_uni): probability_distribution_function(_X_uni), p(0), q(0), r(0), n(0), r_min(0.0), r_max(1.0), r_mean(0.0), r_std(0.0) {}
  discinv_pdf(rng & _X_uni, double * P, int N, double R_min = 0.0, double R_max = 1.0): probability_distribution_function(_X_uni) { set_discrete_pdf(P, N, R_min, R_max); }
  discinv_pdf(String thislabel, Command_Line_Parameters & clp, rng & _X_uni);
  ~discinv_pdf() { if (p) delete[] p; if (q) delete[] q; if (r) delete[] r; } 
  double mean_X() { return r_mean; }
  double std_X() { return r_std; }
  void set_discrete_pdf(double * P, int N, double R_min = 0.0, double R_max = 1.0);
  virtual double random_positive();
  virtual double random_selection();
  virtual double amplitude() { return 1.0; }
  virtual String report_local_parameters();
  virtual String report_parameters();
  virtual probability_distribution_function * clone() { return new discinv_pdf(*X_uni,p,n,r_min,r_max); }
};

probability_distribution_function * pdfselection(String basecommand, Command_Line_Parameters & clp, rng & _X_uni);

double * str2vec(const char * vecstr, int * veclen);

// inline function definitions

//[*** REMOVE]inline double square_distance(double x1, double y1, double x2, double y2) {
//[*** REMOVE]  double h = x1-x2; double v = y1-y2;
//[*** REMOVE]  return (h*h)+(v*v);
//[*** REMOVE]}

/* The C++ FAQ Lite explains that non-member inline functions should be
   declared normally in the header file and that their definition should
   be prepended with the keyword "inline" and contained within the header
   file. Otherwise a linker error occurs.
*/

extern int splinetest;

extern rng X_misc;
extern rng X_elongate;
extern rng X_branch;
extern rng X_turn;
extern rng X_synapses;

#endif
