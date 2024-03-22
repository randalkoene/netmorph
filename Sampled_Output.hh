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
// Sampled_Output.hh
// Randal A. Koene, 20041118
//
// Classes that provide simulation output.

/* Documentation:
   1) An external Sampled_Output object sampled_output can point to any
      derived type of Sampled_Output (e.g. Sampled_Growth_Fig_Output) in
      order to provide all the desired types of sampled output.
   2) The Output() function is called to produce the sampled output.
   3) A command line interface is provided by inheriting CLP_Modifiable
      and through the functions parse_CLP() and report_parameters().
   4) The Sampled_Growth_Output (and derived objects) inherit
      Network_Generated_Statistics. The number of samples (numsamples) is
      initialized as (max_time/sample_interval)+1.
   5) Samples are taken during calls to a Sample() function. This can call
      statistics gathering functions, such as those defined in the network
      object (e.g. network::mean_number_of_input_terminal_segments()).
 */

#ifndef __SAMPLED_OUTPUT_HH
#define __SAMPLED_OUTPUT_HH

#include <limits.h>
#include "Network_Generated_Statistics.hh"
#include "Command_Line_Parameters.hh"
#include "Spatial_Presentation.hh"
#include "spatial.hh"
//#include "network.hh" // include in .cc if needed to avoid problems

// classes

#ifndef __NEURON_HH
class neuron;
#endif

class network; // forward declaration resolved in network.hh

#define FANIN_ANGLE_SLICE (M_PI/50.0)
struct fanin_histogram {
  // one bar per 3.6 degrees, i.e. 50 in a half circle
  unsigned long n[50];
  double len[50];
  bool collect;
  fanin_histogram(): collect(false) {
    for (int i = 0; i<50; i++) {
      n[i] = 0;
      len[i] = 0.0;
    }
  }
  void means() {
    for (int i = 0; i<50; i++) if (n[i]>1) len[i] /= (double) n[i];
  }
  String octave_output();
  Fig_Group * fig_output();
};

class Sampled_Output: public CLP_Modifiable {
protected:
  network * net;
  double max_time;
  double sample_interval;
  int numsamples;
  double t_nextmark;
  bool maximum_resolution;
public:
  Sampled_Output(network * netptr, double mt, double si, bool std = false);
  virtual ~Sampled_Output() {}
  void Null_Net() { net = NULL; } // *** call to insure segfault if net is used
  void max_res_sampling_on() { maximum_resolution = true; }
  void max_res_sampling_off() { maximum_resolution = false; }
  bool max_res() { return maximum_resolution; }
  network * sampled_network() { return net; }
  double sampled_time() { return max_time; }
  double interval() { return sample_interval; }
  int number_of_samples() { return numsamples; }
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual bool Sample(double t) { return false; }
  virtual bool Max_Resolution_Sample(double t, const void * id) { return false; }
  virtual void Output(bool show_stats = true) {}
};

class Sampled_Growth_Output: public Sampled_Output
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
			   , public Network_Generated_Statistics
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
			     , public network_statistics_base
#endif
{
protected:
  // *** (inherited instead) Network_Generated_Statistics collect_statistics;
  String Txtname;
  String statsname;
  bool autoplot;
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  bool store_raw_data;
#endif
  bool segmentcount; // 0 = no, 1 = at end of simulation, 2 = every sample
  bool synapsesearchprofile; // 0 = no, 1 = at end of simulation, 2 = every sample
  bool synapsegenesisandloss; // 0 = no, 1 = at end of simulation, 2 = every sample
public:
  Sampled_Growth_Output(network * netptr, double mt, double si, String tn, String sn, bool ap = false, bool std = false): Sampled_Output(netptr,mt,si,std), Txtname(tn), statsname(sn), autoplot(ap)
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
													       , store_raw_data(false)
#endif
													       , segmentcount(1), synapsesearchprofile(1), synapsegenesisandloss(1)
 {
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
    Initialize_Terminals_Development_Statistics(numsamples);
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
    agesamples = numsamples;
#endif
  }
  virtual ~Sampled_Growth_Output() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual bool Sample(double t);
  virtual bool Max_Resolution_Sample(double t, const void * id) { return true; }
  virtual void Output(bool show_stats = true);
};

class Sampled_Growth_Fig_Output: public Sampled_Growth_Output {
protected:
  String figname;
  double figwidth;
  double combinemagnification;
  long sequencedelay;
  bool combine;
  bool autorotatesequence;
  int netfigidx;
  int topleftx, toplefty, bottomrightx, bottomrighty;
  spatial zoomcenter;
  double zoomdisttoedge; // -1.0 = show full figure (default)
  double zoomwidth; // -1.0 = use zoomdisttoedge
  double zoomheight;
  double zoomdepth;
  String combinetype;
public:
  Sampled_Growth_Fig_Output(network * netptr, double mt, double si, String tn, String sn, bool ap, String fn, double fw, bool c = false, bool ar = false, bool std = false): Sampled_Growth_Output(netptr,mt,si,tn,sn,ap,std), figname(fn), figwidth(fw), combinemagnification(1.0), combine(c), autorotatesequence(ar), netfigidx(0), topleftx(INT_MAX), toplefty(INT_MAX), bottomrightx(INT_MIN), bottomrighty(INT_MIN), zoomdisttoedge(-1.0), zoomwidth(-1.0), zoomheight(-1.0), zoomdepth(-1.0), combinetype("gif") {
    sequencedelay = 2500/numsamples; // default total simulation time 25 seconds
#ifdef VECTOR3D
    if (autorotatesequence) spat_pres->set_rotate_speed(numsamples);
#endif
  }
  bool Combine() { return combine; }
  String CombineType() { return combinetype; }
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  virtual bool Sample(double t);
  virtual bool Max_Resolution_Sample(double t, const void * id);
  virtual void Output(bool show_stats = true);
};

class Activity_Results_Output {
  bool postop_done;
public:
  Activity_Results_Output(): postop_done(false) {}
  virtual ~Activity_Results_Output() { if (!postop_done) postop(); }
  virtual void spike(neuron * n, double t);
  virtual void psp(neuron * pre, neuron * post, double t);
  virtual void postop() { progress("Activity simulation completed.\n"); postop_done = true; }
};

class ARO_to_File: public Activity_Results_Output {
protected:
  String fname;
  String spikedata;
  String pspdata;
  long numspikes;
  long numpsps;
public:
  ARO_to_File(String _fname): fname(_fname), numspikes(0), numpsps(0) {}
  virtual ~ARO_to_File() { if (!postop_done) postop(); }
  virtual void spike(neuron * n, double t);
  virtual void psp(neuron * pre, neuron * post, double t);
  virtual void postop();
};

extern Activity_Results_Output * actres;

// global functions

void set_focus_box(spatial & center, double sizer);
void set_focus_volume(spatial & center, double width, double height, double depth);
void parse_focus_box_CLP(String zoomlabel, Command_Line_Parameters & clp, spatial & center, double & sizer);
void parse_focus_volume_CLP(String zoomlabel, Command_Line_Parameters & clp, double & width, double & height, double & depth);
#ifdef VECTOR2D
bool fig_in_zoom(double x, double y);
#endif
#ifdef VECTOR3D
bool fig_in_zoom(double x, double y, double z);
#endif
bool fig_in_zoom(const spatial & p);

// global variables

extern fanin_histogram net_fanin_histogram;

extern double figattr_focus_box_x1;
extern double figattr_focus_box_y1;
extern double figattr_focus_box_x2;
extern double figattr_focus_box_y2;
#ifdef VECTOR3D
extern double figattr_focus_box_z1;
extern double figattr_focus_box_z2;
#endif
extern bool camera;

//extern Sampled_Output * sampled_output; // *** moved to global.hh

#endif
