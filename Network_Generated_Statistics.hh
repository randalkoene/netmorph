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
// Network_Generated_Statistics.hh
// Randal A. Koene, 20041118
//
// Classes that provide network statistics after simulation.

/* Documentation of network_statistics_sample/set/base:
   (This is the new and improved method.)

   To add a new type of statistics for all data:
   1) Add a label to stat_type_labels.
   2) In network_statistics_sample, add a calculation and caching function,
      and add a call to that function in fill_cache().
   3) Add necessary code to Octave_Output().

   To add a new type of data for which to obtain statistics:
   1) Add a label to stat_labels. When compiling, that will automatically
      require a number of other definitions (e.g. a string label).
   2) Make sure the corresponding data is collected by net->Collect_Data()
      when that call propagates throughout the network.
      OR make sure the data is collected during network development
      functions and that the data is declared static.
 */

/* Documentation (older method):
   1) An object of class Network_Generated_Statistics contains a buffer
      of type stat_array (a list of double arrays, sizes STAT_ARRAY_LABELS_NUM
      by arraysize). The elements are indexed by stat_array_labels and number.
   2) The functions get_value() and add_value() are used to read and write
      the buffer. Specialized functions also exist as a legacy interface for
      specific statistics.
      The add_value() function treats the buffer like a stack, automatically
      incrementing the array index if a value has already been stored in the
      current indexed element of an array.
   3) The function Octave_Output() writes a file with an Octave script that
      can be used to examine the statistical data in the buffer.
      This function is generally called by the Output() function of an
      object derived from Sampled_Growth_Output.
   4) The object network_statistics_data is used to collect standard
      statistics through compatible functions (such as those of the network
      object).

   This is used in Sampled_Output.hh. Objects with base class
   Sampled_Growth_Output inherit Network_Generated_Statistics, so that the
   statistics buffer is generally initialized when a Sampled_Output object
   is allocated (e.g. in nibr.cc). Initialize_Terminals_Development_
   Statistics() is called in the constructor of Sampled_Growth_Output().
   The initialized values are all -1.0.

   (!) In order to add statistics:
   A - Labels must be added to stat_array_labels (standard labels include
       mean values, standard deviation, minimum and maximum values).
   B - Functions must exist (e.g. in network objects) that collect the data.
   C - Those functions are called when a Sampled_.*Output::Sample() is taken.
   D - Sampled_.*Output::Output() should be produced, at least in terms of
       Network_Generated_Statistics::Octave_Output().
 */

#ifndef __NETWORK_GENERATED_STATISTICS_HH
#define __NETWORK_GENERATED_STATISTICS_HH

//#define VERBOSE_EMPTY_DATA_SAMPLE_WARNINGS

#include <values.h>
#include "global.hh"
#include "BigString.hh"
#include "Command_Line_Parameters.hh"
#include "templates.hh"

typedef double * stat_array;

#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
enum stat_array_labels {arrayoftime,arrayofdendritesmean_nt,arrayofdendritesstd_nt,arrayofaxonsmean_nt,arrayofaxonsstd_nt,arrayofdendritearbors,arrayofaxonarbors,arrayofdendritelength,arrayofaxonlength,arrayofdendriteterminallength,arrayofaxonterminallength,arrayofdendritelengthstd,arrayofaxonlengthstd,arrayofdendriteterminallengthstd,arrayofaxonterminallengthstd,arrayofdendritesmin_nt,arrayofdendritesmax_nt,arrayofaxonsmin_nt,arrayofaxonsmax_nt,arrayofdendritelengthmin,arrayofdendritelengthmax,arrayofaxonlengthmin,arrayofaxonlengthmax,arrayofdendriteterminallengthmin,arrayofdendriteterminallengthmax,arrayofaxonterminallengthmin,arrayofaxonterminallengthmax,arrayofdendriteintermediatelength,arrayofdendriteintermediatelengthstd,arrayofdendriteintermediatelengthmin,arrayofdendriteintermediatelengthmax,arrayofaxonintermediatelength,arrayofaxonintermediatelengthstd,arrayofaxonintermediatelengthmin,arrayofaxonintermediatelengthmax,STAT_ARRAY_LABELS_NUM};
#endif

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
enum stat_labels {
  dendrites_arbor_len, // per arbor
  axons_arbor_len,
  dendrites_termsegs_n,
  axons_termsegs_n,
  // dendrites_intersegs_n,
  // axons_intersegs_n,
  dendrites_termsegs_len, // per fibre segment
  axons_termsegs_len,
  dendrites_intersegs_len,
  axons_intersegs_len,
  dendrites_len_between_bifurcations,
  axons_len_between_bifurcations,
  dendrites_term_len_since_bifurcation,
  axons_term_len_since_bifurcation,
  dendrites_term_len_since_soma,
  axons_term_len_since_soma,
  dendrites_turns_between_bifurcations,
  axons_turns_between_bifurcations,
  dendrites_branch_angles,
  axons_branch_angles,
  dendrites_turn_angles,
  axons_turn_angles,
  dendrites_cartratio_between_bifurcations, // cartesian to fiber ratio
  axons_cartratio_between_bifurcations,
  dendrites_cartratio_bifurcationtoterm,
  axons_cartratio_bifurcationtoterm,
  dendrites_cartratio_somatoterm,
  axons_cartratio_somatoterm,
  dendrites_seven_micron_branch_angles, // measured with the 7 micron method
  axons_seven_micron_branch_angles,
  STAT_LABELS_NUM};

extern const char stat_label_str[STAT_LABELS_NUM][30];

extern const bool stat_is_static[STAT_LABELS_NUM];

enum stat_type_labels {stat_n,stat_mean,stat_std,stat_min,stat_max,STAT_TYPE_LABELS_NUM};
#endif

class network_statistics_data {
  // This class collects typical information needed to compute a full set
  // of statistics about sample of specific parameters: mean, standard
  // deviation, minimum and maximum, and the number of samples.
public:
  network_statistics_data(): n(0.0), mean(0.0), std(0.0), min(MAXDOUBLE), max(MINDOUBLE) {}
  double n;
  double mean;
  double std;
  double min;
  double max;
};

/* Documentation:
   The new network generated statistics construct involves arrays of statistics
   labels and statistics objects. Each object contains variables and functions
   that enable calculation of the statistics from raw data. Each object also
   contains the root of a list of raw data collector objects. The list is
   automatically expanded as more array space is needed.
   The statistics objects are known entities (either globally declared,
   or assigned to globally declared pointers).
   In network, neuron, connection, synapse and other objects, collect_data()
   functions are called to traverse the network, and each sends its data to
   a function of the relevant statistics object. There, the data is sent to
   the root of the raw data collector objects, which sends it to the head of
   the list.
   The statistics and distributions are calculated once all raw data has
   been collected. Optionally, raw data may be stored as well. Raw data
   can then be discarded.
   At the beginning of data collection at a specific age t, the network
   statistics base adds new network statistics samples, and those network
   statistics set heads are then used for data collection.

   Alternative considerations can lead to the use of modified derived
   versions of these classes. One such derived version could rely on two
   data collection events in order to calculate statistics immediately,
   so that raw data need not be stored during data collection.
 */

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
class raw_data_collector: public PLLHandle<raw_data_collector> {
#define RAW_DATA_BLOCK_SIZE 500
protected:
  double data[RAW_DATA_BLOCK_SIZE];
  long numdata;
public:
  raw_data_collector(): numdata(0) {}
  void Collect_Data(double d) {
    if (numdata<RAW_DATA_BLOCK_SIZE) {
      data[numdata] = d;
      numdata++;
    } else {
      raw_data_collector * rdc(new raw_data_collector());
      Root()->link_before(rdc);
      rdc->Collect_Data(d);
    }
  }
  long N() { return numdata; }
  double Sum() {
    double sum = 0.0;
    for (int i = 0; i<numdata; i++) sum += data[i];
    return sum;
  }
  double Min() {
    double min = MAXDOUBLE;
    for (int i = 0; i<numdata; i++) if (data[i]<min) min = data[i];
    return min;
  }
  double Max() {
    double max = MINDOUBLE;
    for (int i = 0; i<numdata; i++) if (data[i]>max) max = data[i];
    return max;
  }
  double SumofSquares(double m) {
    double ssq = 0.0, sq;
    for (int i = 0; i<numdata; i++) {
      sq = data[i] - m;
      ssq += (sq*sq);
    }
    return ssq;
  }
  void str(String & s, char sep = '\n') {
    for (int i = 0; i<numdata; i++) s += String(data[i],"%e")+sep;
  }
};

class raw_data_collector_root: public PLLRoot<raw_data_collector> {
public:
  raw_data_collector_root() {}
  void Collect_Data(double d) {
    if (!tail()) link_before(new raw_data_collector());
    tail()->Collect_Data(d);
  }
  long N() {
    long n = 0;
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) n += e->N();
    return n;
  }
  double Sum() {
    double sum = 0.0;
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) sum += e->Sum();
    return sum;
  }
  double Min() {
    double min = MAXDOUBLE, emin;
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) {
      emin = e->Min();
      if (emin<min) min = emin;
    }
    return min;
  }
  double Max() {
    double max = MINDOUBLE, emax;
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) {
      emax = e->Max();
      if (emax>max) max = emax;
    }
    return max;
  }
  double SumofSquares(double m) {
    double ssq = 0.0;
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) ssq += e->SumofSquares(m);
    return ssq;
  }
  void str(String & s, char sep = '\n') {
    PLL_LOOP_FORWARD(raw_data_collector,head(),1) e->str(s,sep);
  }
};

extern unsigned long empty_data_samples;

class network_statistics_sample: public PLLHandle<network_statistics_sample> {
  // This is a linked list approach to the collection of many types of
  // statistical information, including probability distributions.
  // For each data label, there is one of these sample objects per sample
  // time (age).
protected:
  double t;
  raw_data_collector_root rdcr;
  // flags indicate which values have already been calculated and cached
  bool cached[STAT_TYPE_LABELS_NUM];
  bool allcached;
  // caches (once statistics have been cached, raw data may be deleted)
  double statscache[STAT_TYPE_LABELS_NUM];
public:
  network_statistics_sample(): t(-1.0) { reset_cache(); }
  //network_statistics_data nsd; // sample statistics
  void set_age(double sample_t) { t = sample_t; }
  double age() { return t; }
  raw_data_collector_root & Raw_Data() { return rdcr; }
  double N() {
    if (!cached[stat_n]) {
      statscache[stat_n] = (double) rdcr.N();
      cached[stat_n] = true;
    }
    return statscache[stat_n];
  }
  double Mean() { // mean of data at this sample time (age)
    if (!cached[stat_mean]) {
      statscache[stat_mean] = rdcr.Sum() / N();
      cached[stat_mean] = true;
    }
    return statscache[stat_mean];
  }
  double Min() {
    if (!cached[stat_min]) {
      statscache[stat_min] = rdcr.Min();
      cached[stat_min] = true;
    }
    return statscache[stat_min];
  }
  double Max() {
    if (!cached[stat_max]) {
      statscache[stat_max] = rdcr.Max();
      cached[stat_max] = true;
    }
    return statscache[stat_max];
  }
  double Std() {
    if (!cached[stat_std]) {
      if (N()<=1.0) statscache[stat_std] = 0.0; // *** was -1.0 tag, but that caused trouble in the evaluation of stats.m (see TL#200606020621.1)
      else statscache[stat_std] = sqrt(rdcr.SumofSquares(Mean())/(N()-1.0));
      cached[stat_std] = true;
    }
    return statscache[stat_std];
  }
  void reset_cache() {
    for (int i = 0; i<STAT_TYPE_LABELS_NUM; i++) cached[i] = false;
    allcached = false;
  }
  void fill_cache() {
    if (N()<=0.0) {
      empty_data_samples++;
#ifdef VERBOSE_EMPTY_DATA_SAMPLE_WARNINGS
       warning("Empty data sample at t="+String(t,"%.3f"));
#endif
      for (int i = 0; i<STAT_TYPE_LABELS_NUM; i++) statscache[i] = 0.0;
      allcached = true;
      return;
    }
    Mean();
    Min();
    Max();
    Std();
    allcached = true;
  }
  void str(String & s, char sep = ' ') {
    if (!allcached) fill_cache();
    for (int i = 0; i<STAT_TYPE_LABELS_NUM; i++) s += String(statscache[i],"%.7f")+sep;
  }
  void Octave_Output(String & s, char sep = ' ') {
    /* Documentation:
       The current version prepends the output with N so that the actual
       number of samples used to calculate the statistics is explicitly
       specified. This can be dealt with in the nibr.m script to enable
       backward compatibility with the previous version that omitted N.
     */
    if (!allcached) fill_cache();
    //for (int i=stat_n; i<stat_max; i++) { cout << statscache[i] << '/'; cout.flush();}
    for (int i = stat_n; i<=stat_max; i++) {
      //cout << statscache[i] << ','; cout.flush();
      s += String(statscache[i],"%.7f")+sep; // Note: This can crash if the value is MAXDOUBLE
    }
  }
};

class network_statistics_set: public PLLRoot<network_statistics_sample>, public CLP_Modifiable {
  friend class network_statistics_base;
  // Each set contains the statistics for all sample times (ages) for one
  // data label.
protected:
  String label;
  bool collect;
public:
  network_statistics_set(): collect(true) {} //(String l): label(l) {}
  virtual ~network_statistics_set() {}
  void Collect_Data(double d) { if (collect) tail()->Raw_Data().Collect_Data(d); }
  String Label() { return label; }
  void set_collect(bool c) { collect = c; }
  bool Collect() { return collect; }
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
};

class network_statistics_base {
  // This provides access to all network generated statistics and is called
  // during and after data collection to obtain raw data and to calculate
  // statistics at each sample time (age).
protected:
  long agesamples; // initialized by derived class (Sampled_Growth_Output)
  bool collect_statistics;
public:
  network_statistics_base();
  /* Documentation:
     Static verus Variable data sets: A data set is variable if the data of
     sampled components that exist at two different sample times may be
     different. A data set is static if existing components have unchanged
     sampled data, changing the set only by adding the data for new
     components. The data of each component in a variable data sets must be
     sampled during each sample interval. Previously collected raw data of
     a static data set is reused and appended to as new components are
     created.
     Notes: (1) If a static data set is not collected during object creation
     then later parsing must either identify previously existing objects
     (i.e. store an object pointer) or resample the data anyway. (2) Existing
     static data is currently simply relinked, so that the full war data set
     is not available for all ages if raw data is not deleted after each
     sample interval. To obtain that, either copies must be made or a time
     value must be stored with each raw datum.
   */
  network_statistics_set netstats[STAT_LABELS_NUM];
  bool Collecting_Statistics() { return collect_statistics; }
  void Allocate_New_Data();
  void Set_Tail_Age(double sample_t) {
    for (int i=0; i<STAT_LABELS_NUM; i++) if (netstats[i].tail()) netstats[i].tail()->set_age(sample_t);
  }
  // Sample time (age) increases from head() to tail() so that an increasing
  // ordered retrieval is possible by using the PLL_LOOP_FORWARD macro.
  void Store_Tail_Raw_Data(String datafile);
  void Cache_Tail_Statistics();
  void Remove_Tail_Raw_Data(bool evenstatic = false);
  void Octave_Output(String filename, bool autoplot = false);
};
// [***INCOMPLETE] I may wish to create a virtual class
// network_statistics_data_source, which forces those objects that inherit
// the class to define Collect_Data() functions. (See TL#200602120629.2.)
/* [***INCOMPLETE --- THE FOLLOWING MAY BE OUTDATED INSTRUCTIONS (20060212)]
  Steps toward a new approach toward the gathering of statistical information:
  1. Create a new Network_Generated_Statistics class that uses
     networks_statistics_set: (a) make Network_Generated_Statistics a base
     class that can't be initialized directly, (b) make the current version a
     derivative of that with virtual functions, but make some functions only
     available in that specific class so that references to them are found
     and can be replaced during compilation, (c) make a new derivative that
     uses the new approach.
  2. Create a new Octave output function for the new approach. Much of the
     same output should be created for backward compatibility.
  3. Replace the old approach with the new approach in the nibr simulation
     functions, compile and test.
  4. Add new parts (and even previously provided parts), such as probability
     distributions in a manner where the Octave output demands their
     availability, forcing me to provide the necesssary data gathering
     functions.
 */
#endif

#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
class Network_Generated_Statistics {
#define NGS_ADDSAMPLE(arrayname,valuename) { \
  if (!statistics[arrayname]) return; \
  if ((statistics[arrayname])[arrayindex]>=0.0) { \
    if ((arrayindex+1)>=arraysize) return; \
    arrayindex++; \
  } \
  (statistics[arrayname])[arrayindex] = valuename; \
}
protected:
  stat_array statistics[STAT_ARRAY_LABELS_NUM];
  int arraysize, arrayindex;
  void clear() { for (int i = 0; i < STAT_ARRAY_LABELS_NUM; i++) delete[] (statistics[i]); arraysize = 0; arrayindex = -1; }
public:
  Network_Generated_Statistics(): arraysize(0), arrayindex(-1) { for (int i = 0; i < STAT_ARRAY_LABELS_NUM; i++) statistics[i] = NULL; }
  ~Network_Generated_Statistics() { clear(); }
  void Initialize_Terminals_Development_Statistics(int numdatapoints); // , bool std = false);
  int array_size() { return arraysize; }
  int number_of_data_points() { return arrayindex; }
  double t(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayoftime])[n]; else return -1.0; }
  double recent_t() { if (statistics[arrayoftime]) return (statistics[arrayoftime])[arrayindex]; else return -1.0; }
  double get_value(stat_array_labels parlabel, int n) { if ((parlabel<STAT_ARRAY_LABELS_NUM) && (n<=arrayindex) && (arrayindex>=0)) return (statistics[parlabel])[n]; else return -1.0; }
  double Dendrites_Mean_nt(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendritesmean_nt])[n]; else return -1.0; }
  double Dendrites_Std_nt(int n) { if ((arrayofdendritesstd_nt) && (n<=arrayindex)) return (statistics[arrayofdendritesstd_nt])[n]; else return -1.0; }
  double Axons_Mean_nt(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonsmean_nt])[n]; else return -1.0; }
  double Axons_Std_nt(int n) { if ((arrayofaxonsstd_nt) && (n<=arrayindex)) return (statistics[arrayofaxonsstd_nt])[n]; else return -1.0; }
  double Dendrites_Arbors(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendritearbors])[n]; else return -1.0; }
  double Axons_Arbors(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonarbors])[n]; else return -1.0; }
  double Dendrites_Mean_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendritelength])[n]; else return -1.0; }
  double Axons_Mean_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonlength])[n]; else return -1.0; }
  double Dendrites_Mean_Terminal_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendriteterminallength])[n]; else return -1.0; }
  double Axons_Mean_Terminal_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonterminallength])[n]; else return -1.0; }
  double Dendrites_Std_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendritelengthstd])[n]; else return -1.0; }
  double Axons_Std_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonlengthstd])[n]; else return -1.0; }
  double Dendrites_Std_Terminal_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofdendriteterminallengthstd])[n]; else return -1.0; }
  double Axons_Std_Terminal_Length(int n) { if ((n<=arrayindex) && (arrayindex>=0)) return (statistics[arrayofaxonterminallengthstd])[n]; else return -1.0; }
  void add_t(double T) { NGS_ADDSAMPLE(arrayoftime,T) }
  void add_value(stat_array_labels parlabel, double v) { if (parlabel<STAT_ARRAY_LABELS_NUM) NGS_ADDSAMPLE(parlabel,v) }
  void add_dendrites_mean_nt(double mnt) { NGS_ADDSAMPLE(arrayofdendritesmean_nt,mnt) }
  void add_dendrites_std_nt(double std) { NGS_ADDSAMPLE(arrayofdendritesstd_nt,std) }
  void add_axons_mean_nt(double mnt) { NGS_ADDSAMPLE(arrayofaxonsmean_nt,mnt) }
  void add_axons_std_nt(double std) { NGS_ADDSAMPLE(arrayofaxonsstd_nt,std) }
  void add_dendrites_arbors(double arborsamples) { NGS_ADDSAMPLE(arrayofdendritearbors,arborsamples) }
  void add_axons_arbors(double arborsamples) { NGS_ADDSAMPLE(arrayofaxonarbors,arborsamples) }
  void add_dendrites_mean_length(double dendritelen) { NGS_ADDSAMPLE(arrayofdendritelength,dendritelen) }
  void add_axons_mean_length(double axonlen) { NGS_ADDSAMPLE(arrayofaxonlength,axonlen) }
  void add_dendrites_mean_terminal_length(double dendritetermlen) { NGS_ADDSAMPLE(arrayofdendriteterminallength,dendritetermlen) }
  void add_axons_mean_terminal_length(double axontermlen) { NGS_ADDSAMPLE(arrayofaxonterminallength,axontermlen) }
  void add_dendrites_std_length(double dendritelen) { NGS_ADDSAMPLE(arrayofdendritelengthstd,dendritelen) }
  void add_axons_std_length(double axonlen) { NGS_ADDSAMPLE(arrayofaxonlengthstd,axonlen) }
  void add_dendrites_std_terminal_length(double dendritetermlen) { NGS_ADDSAMPLE(arrayofdendriteterminallengthstd,dendritetermlen) }
  void add_axons_std_terminal_length(double axontermlen) { NGS_ADDSAMPLE(arrayofaxonterminallengthstd,axontermlen) }
  void add_dendrites_axons_nt(double t, network_statistics_data & dnumnsd, network_statistics_data & anumnsd);
  void add_dendrites_axons_lengths(double t, network_statistics_data & dlennsd, network_statistics_data & alennsd, network_statistics_data & dtermnsd, network_statistics_data & atermnsd);
  void add_dendrites_axons_intermediate_lengths(double t, network_statistics_data & dinternsd, network_statistics_data & ainternsd);
  void Octave_Output(String filename, bool autoplot = false);
  // public cache variables
};
#endif

#endif
