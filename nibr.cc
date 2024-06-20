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
// nibr.cc
// Randal A. Koene, 20040830

/*
   The ways to determine network connectivity are described in
   ~/doc/tex/generation-framework/generation-framework.tex.
 */

/* Strategy and future to-dos:
- Instead of focusing on any speed/efficiency related items during the
  initial implementation, focus on producing output that is useful for
  research results. Then, when speed-ups become necessary, identify the
  main offenders (from an overall perspective) and tackle those.
- I won't focus on this now, but in a future rewrite there will certainly
  be files and classes that can best be split or combined (e.g.
  presynaptic_structure and postsynaptic_structure may be represented by
  a single class or one may be derived from the other to reduce copied
  code).
*/

#include <sys/types.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <cstddef>
using namespace std;
#include "diagnostic.hh"
#include "Sampled_Output.hh"
#include "Txt_Object.hh"
#include "file.hh"
#include "nibr.hh"
#include "prepost_structure.hh"
#include "synapse.hh"
#include "synapse_formation_model.hh"
#include "axon_direction_model.hh"
#include "environment_physics.hh"
#include "slice.hh"
#include "event.hh"
//#include "generic_synapse.hh"
#ifdef TESTING_SPIKING
#include "IFactivity.hh"
#endif
#include "global.hh"
#include "mtprng.hh"

void reliability_checklist() {
  // This function is used to check some essential variables
  // before attempting a simulation. This is necessary, because
  // some issues are not caught by the compiler. For example,
  // if a new neuron type was added, but the synapse formation
  // maxtrices were not initialized with the new number of
  // neurons, then remaining array entries are automatically
  // initialized to 0, which can cause unexpected results.
  if (superior_natural_set[NUM_NATURAL_SPS-1]<=universal_sps) error("Compile-Time Error: Not all elements of global:superior_natural_set[] have been correctly initialized.\n");
  if (natural_subset_idstr[NUM_NATURAL_SPS-1][0]=='\0') error("Compile-Time Error: Not all elements of neuron:natural_subset_idstr[] have been correctly initialized.\n");
  if (min_basal[UNTYPED_NEURON]<1) error("Compile-Time Error: Not all elements of neuron:min_basal[] have been correctly initialized.\n");
  if (max_basal[UNTYPED_NEURON]<1) error("Compile-Time Error: Not all elements of neuron:max_basal[] have been correctly initialized.\n");
  if (max_angle[UNTYPED_NEURON]<=0.0) error("Compile-Time Error: Not all elements of neuron:max_angle[] have been correctly initialized.\n");
  if (max_axons[UNTYPED_NEURON]<1) error("Compile-Time Error: Not all elements of neuron:max_axons[] have been correctly initialized.\n");
  if (possible_syntypes[UNTYPED_NEURON][UNTYPED_NEURON][syntype_candidate-1]>-9.0) error("Compile-Time Error: Not all elements of synapse_formation_model:possible_syntypes[][][] have been correctly initialized.\n");
  if (axons_most_specific_natural_set[UNTYPED_NEURON]<=universal_sps) error("Compile-Time Error: Not all elements of prepost_structure:axons_most_specific_natural_set[] have been correctly initialized.\n");
  if (dendrites_most_specific_natural_set[UNTYPED_NEURON]<=universal_sps) error("Compile-Time Error: Not all elements of prepost_structure:dendrites_most_specific_natural_set[] have been correctly initialized.\n");
}

electrodearray * make_electrodes_hexagon(const spatial & center) {
  // produces a specific electrode design
  // electrode design in which every electrod is 70 micrometers from the
  // nearest electrode (i.e. sets of three electrodes form equilateral
  // triangles)
  // vertical distance = sqrt(70^2 - (70/2)^2) = 60.622
#define NUMHEXAGONELECTRODES 61
#define ELEC_h 60.622
  const double xpos[NUMHEXAGONELECTRODES] = {0.0,70.0,140.0,210.0,280.0,-70.0,-140.0,-210.0,-280.0,
			 35.0,105.0,175.0,245.0,-35.0,-105.0,-175.0,-245.0,
			 35.0,105.0,175.0,245.0,-35.0,-105.0,-175.0,-245.0,
			 0.0,70.0,140.0,210.0,-70.0,-140.0,-210.0,
			 0.0,70.0,140.0,210.0,-70.0,-140.0,-210.0,
			 35.0,105.0,175.0,-35.0,-105.0,-175.0,
			 35.0,105.0,175.0,-35.0,-105.0,-175.0,
			 0.0,70.0,140.0,-70.0,-140.0,
			 0.0,70.0,140.0,-70.0,-140.0};
  const double ypos[NUMHEXAGONELECTRODES] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			 ELEC_h,ELEC_h,ELEC_h,ELEC_h,ELEC_h,ELEC_h,ELEC_h,ELEC_h,
			 -ELEC_h,-ELEC_h,-ELEC_h,-ELEC_h,-ELEC_h,-ELEC_h,-ELEC_h,-ELEC_h,
			 2.0*ELEC_h,2.0*ELEC_h,2.0*ELEC_h,2.0*ELEC_h,2.0*ELEC_h,2.0*ELEC_h,2.0*ELEC_h,
			 -2.0*ELEC_h,-2.0*ELEC_h,-2.0*ELEC_h,-2.0*ELEC_h,-2.0*ELEC_h,-2.0*ELEC_h,-2.0*ELEC_h,
			 3.0*ELEC_h,3.0*ELEC_h,3.0*ELEC_h,3.0*ELEC_h,3.0*ELEC_h,3.0*ELEC_h,
			 -3.0*ELEC_h,-3.0*ELEC_h,-3.0*ELEC_h,-3.0*ELEC_h,-3.0*ELEC_h,-3.0*ELEC_h,
			 4.0*ELEC_h,4.0*ELEC_h,4.0*ELEC_h,4.0*ELEC_h,4.0*ELEC_h,
			 -4.0*ELEC_h,-4.0*ELEC_h,-4.0*ELEC_h,-4.0*ELEC_h,-4.0*ELEC_h};
  electrodearray * els = new electrodearray(NUMHEXAGONELECTRODES);
  for (int i=0; i<NUMHEXAGONELECTRODES; i++) {
    els->PLLRoot<neuron>::el(i)->set_position_in_Z_plane(xpos[i]+center.X(),ypos[i]+center.Y());
  }
  els->set_center(center);
  return els;
}

void running_instance_preop() {
  // set PID file
  pid_t thispid = getpid();
  String thispidstr((long) thispid);
  if (!write_file_from_String(outputdirectory+"nibr.pid",thispidstr)) {
    cout << "Content-Type: text/html\n\n<HTML>\n<BODY>\n<H1>Network Generation: Simulation</H1>\n\n<FONT COLOR=\"#FF0000\"><B>Error: Unable to store process ID in .pid file.</B></FONT>\n\n</BODY>\n</HTML>\n";
    exit(1);
  }
}

void running_instance_postop() {
  // remove PID file
  unlink(outputdirectory+"nibr.pid");
}

void single_instance_restriction() {
  String pidstr;
  if (read_file_into_String(outputdirectory+"nibr.pid",pidstr,false)) { // another instance may be running
#ifdef LINUX
    pid_t instancepid = (pid_t) atol(pidstr);
    sigval sv; sv.sival_int = 0;
    int res = sigqueue(instancepid,0,sv);
    if (res!=0) if (errno!=ESRCH) res = 0; // test if nibr.pid contains stale PID
    if (res==0) {
#endif
      cout << "Content-Type: text/html\n\n<HTML>\n<BODY>\n<H1>Network Generation: Simulation</H1>\n\n<FONT COLOR=\"#FF0000\"><B>Another simulation is currently running. To conserve resources on the server, only one instance is permitted through the internet form interface. To run multiple instances concurrently, please contact Randal A. Koene as indicated at <A HREF=\"http://rak.minduploading.org\">http://rak.minduploading.org</A>.</B></FONT>\n\n</BODY>\n</HTML>\n";
      exit(1);
#ifdef LINUX
    }
#endif
  }
  running_instance_preop();
}

void report_form_aware(Command_Line_Parameters & clp, CLP_Modifiable * clpm) {
  // Handles parameter reporting for different output targets such as
  // a command line terminal or HTML output following form input.
  // The reporting text is properly encapsulated and where necessary,
  // control characters are replaced.
  if (clp.IsFormInput()) report("</PRE>\n<!-- ");
  report(clpm->report_parameters());
  if (clp.IsFormInput()) report("-->\n<PRE>");
}

void report_form_aware(Command_Line_Parameters & clp, String s) {
  // The same as above, but can present any string.
  if (clp.IsFormInput()) report("</PRE>\n<!-- ");
  report(s);
  if (clp.IsFormInput()) report("-->\n<PRE>");
}

void devsim_present_HTML_graphical_output(String & statsdatafile, String & netfigfile, String & zoomfigfile, String & abstractfigfile, String & statsfile, bool figsequenceflag, bool combineflag, String & combinetype) {
  String statsdatabase(statsdatafile.before('.',-1));
  String statsdataURL(statsdatabase.after('/',-1));
  String statsdatadirectory(statsdatabase.before('/',-1));
  if (statsdataURL.empty()) statsdataURL=statsdatabase;
  statsdataURL.prepend(absoluteoutputURL);

#ifndef COMPILE_WITHOUT_FIG2DEV
  // prepare network structure PDF files
  if (outattr_show_figure && figattr_make_full_Fig) system(String(FIG2DEV)+" -L pdf "+netfigfile+' '+outputdirectory+"net.pdf");
  if (outattr_show_figure && figattr_make_zoom_Fig) system(String(FIG2DEV)+" -L pdf "+zoomfigfile+' '+outputdirectory+"zoom.pdf");
  if (outattr_show_figure && figattr_make_abstract_Fig) system(String(FIG2DEV)+" -L pdf "+abstractfigfile+' '+outputdirectory+"abstract.pdf");
  //system(String(FIG2DEV)+" -L pdf "+connectionexamplefigfile+' '+outputdirectory+"connection-example.pdf");

  // prepare processed network statistics output
  if (outattr_show_stats)
    system(
	   "octave -q "+statsfile+' '+statsdatadirectory+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendritent.fig "+statsdatabase+"-dendritent.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axonnt.fig "+statsdatabase+"-axonnt.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendritelength.fig "+statsdatabase+"-dendritelength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axonlength.fig "+statsdatabase+"-axonlength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendritetermlength.fig "+statsdatabase+"-dendritetermlength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axontermlength.fig "+statsdatabase+"-axontermlength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendriteinterlength.fig "+statsdatabase+"-dendriteinterlength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axoninterlength.fig "+statsdatabase+"-axoninterlength.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendritebranchangles.fig "+statsdatabase+"-dendritebranchangles.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axonbranchangles.fig "+statsdatabase+"-axonbranchangles.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-dendriteturnangles.fig "+statsdatabase+"-dendriteturnangles.png"+" && "+
	   String(FIG2DEV)+" -L png -m 0.5 "+statsdatabase+"-axonturnangles.fig "+statsdatabase+"-axonturnangles.png"
	   );

  // network structure development animated sequence
  if (figsequenceflag && combineflag) progress("<BR><IMG SRC=\""+absoluteoutputURL+"sample."+combinetype+"\">\n");

  // network statistics figures
  if (outattr_show_stats)
    progress("<TABLE><TR><TD><A HREF=\""+statsdataURL+"-dendritent.fig\">Mean dendrite terminal segment numbers</A>:<BR><IMG SRC=\""+statsdataURL+"-dendritent.png\">\n<TD><A HREF=\""+statsdataURL+"-axonnt.fig\">Mean axon terminal segment numbers</A>:<BR><IMG SRC=\""+statsdataURL+"-axonnt.png\">\n<TR><TD><A HREF=\""+statsdataURL+"-dendritelength.fig\">Mean dendrite length</A>:<BR><IMG SRC=\""+statsdataURL+"-dendritelength.png\">\n<TD><A HREF=\""+statsdataURL+"-axonlength.fig\">Mean axon length</A>:<BR><IMG SRC=\""+statsdataURL+"-axonlength.png\">\n<TR><TD><A HREF=\""+statsdataURL+"-dendritetermlength.fig\">Mean dendrite terminal segment length</A>:<BR><IMG SRC=\""+statsdataURL+"-dendritetermlength.png\">\n<TD><A HREF=\""+statsdataURL+"-axontermlength.fig\">Mean axon terminal segment length</A>:<BR><IMG SRC=\""+statsdataURL+"-axontermlength.png\">\n<TR><TD><A HREF=\""+statsdataURL+"-dendriteinterlength.fig\">Mean dendrite intermediate segment length</A>:<BR><IMG SRC=\""+statsdataURL+"-dendriteinterlength.png\">\n<TD><A HREF=\""+statsdataURL+"-axoninterlength.fig\">Mean axon intermediate segment length</A>:<BR><IMG SRC=\""+statsdataURL+"-axoninterlength.png\">\n<TR><TD><A HREF=\""+statsdataURL+"-dendritebranchangles.fig\">Mean dendrite branch angles</A>:<BR><IMG SRC=\""+statsdataURL+"-dendritebranchangles.png\">\n<TD><A HREF=\""+statsdataURL+"-axonbranchangles.fig\">Mean axon branch angles</A>:<BR><IMG SRC=\""+statsdataURL+"-axonbranchangles.png\">\n<TR><TD><A HREF=\""+statsdataURL+"-dendriteturnangles.fig\">Mean dendrite turn angles</A>:<BR><IMG SRC=\""+statsdataURL+"-dendriteturnangles.png\">\n<TD><A HREF=\""+statsdataURL+"-axonturnangles.fig\">Mean axon turn angles</A>:<BR><IMG SRC=\""+statsdataURL+"-axonturnangles.png\">\n</TABLE>");
#else
  progress("<P><B>COMPILED WITHOUT FIG2DEV - FIGURES NOT EMBEDDED IN HTML PAGE.</B><P>");
#endif
}

// The following function is a built-in example of command sequences as they
// would be used in command scripts.

void developmental_simulation(Command_Line_Parameters & clp) {
  // A temporary function that mimics a command script for a simulation of
  // dendritic and axonal growth with equations by van Pelt and Uylings.

  // setting up parameters
  // Defaults (that differ from those compiled into the program) are given
  // in the .nibrrc file.
  int n;
  // Kludge: Eventually replace this with object
  // initialization and an object::parse_CLP() function call.
  // OUTPUT FORMAT PARAMETERS
  if ((n=clp.Specifies_Parameter("outattr_show_progress"))>=0) outattr_show_progress = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_show_figure"))>=0) outattr_show_figure = (downcase(clp.ParValue(n))==String("true"));
  // [***NOTE] The two lines below may be outdated. This is now done through
  // stats_attr_collect_statistics.
  //bool collect_statistics = true;
  //if ((n=clp.Specifies_Parameter("collect_statistics"))>=0) collect_statistics = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_show_stats"))>=0) outattr_show_stats = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_track_synaptogenesis"))>=0) outattr_track_synaptogenesis = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_track_nodegenesis"))>=0) outattr_track_nodegenesis = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_make_full_Txt"))>=0) outattr_make_full_Txt = (downcase(clp.ParValue(n))==String("true"));
  if ((outattr_track_nodegenesis) || (outattr_make_full_Txt)) outattr_track_synaptogenesis = true;
  if (outattr_track_synaptogenesis) SynaptoGenesis_Data = new synaptogenesis_data;
  if (outattr_track_nodegenesis) NodeGenesis_Data = new nodegenesis_data;
  if (outattr_make_full_Txt) {
    if ((n=clp.Specifies_Parameter("outattr_Txt_sequence"))>=0) outattr_Txt_sequence = (downcase(clp.ParValue(n))==String("true"));
  } else outattr_Txt_sequence = false;
  if ((n=clp.Specifies_Parameter("outattr_Txt_separate_files"))>=0) outattr_Txt_separate_files = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_make_full_X3D"))>=0) outattr_make_full_X3D = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("outattr_make_full_Catacomb"))>=0) outattr_make_full_Catacomb = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("statsattr_fan_in_analysis"))>=0) {
    if (downcase(clp.ParValue(n))==String("axons")) statsattr_fan_in_analysis = axon_fs;
    else if (downcase(clp.ParValue(n))==String("dendrites")) statsattr_fan_in_analysis = dendrite_fs;
    else warning("Warning: Unknown fiber structure reference '"+clp.ParValue(n)+"', assuming statsattr_fan_in_analysis=none.\n");
  }
  if ((n=clp.Specifies_Parameter("figattr_use_color"))>=0) figattr_use_color = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_tsupd_visibly"))>=0) figattr_update_terminal_segments_visibly = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_neurons"))>=0) figattr_show_neurons = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_connections"))>=0) figattr_show_abstract_connections = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_presynaptic"))>=0) figattr_show_presynaptic_structure = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_postsynaptic"))>=0) figattr_show_postsynaptic_structure = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_synapses"))>=0) figattr_show_synapse_structure = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_AMPAR"))>=0) figattr_receptor[syntype_AMPAR] = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_NMDAR"))>=0) figattr_receptor[syntype_NMDAR] = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_GABAR"))>=0) figattr_receptor[syntype_GABAR] = (downcase(clp.ParValue(n))==String("true"));
  if (figattr_show_synapse_structure) { // report the display parameters (possibly in HTML comments)
    report_form_aware(clp,"figattr_AMPAR="+String(figattr_receptor[syntype_AMPAR])+"\nfigattr_NMDAR="+String(figattr_receptor[syntype_NMDAR])+"\nfigattr_GABAR="+String(figattr_receptor[syntype_GABAR]));
  }
  if ((n=clp.Specifies_Parameter("figattr_partitions"))>=0) figattr_show_spatial_segment_subsets = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_connection_eval"))>=0) figattr_show_connection_evaluation = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_progress"))>=0) {
    String figattrstr(downcase(clp.ParValue(n)));
    if (figattrstr==String("none")) figattr_sample_show_progress = 0;
    else if (figattrstr==String("text")) figattr_sample_show_progress = 1;
    else if (figattrstr==String("bar")) figattr_sample_show_progress = 2;
  }
  switch (figattr_sample_show_progress) {
  case 0: report("\nFigure Attribute: No progress indicator.\n"); break;;
  case 1: report("\nFigure Attribute: Textual progress indicator.\n"); break;;
  case 2: report("\nFigure Attribute: Bar progress indicator.\n"); break;;
  }
  if ((n=clp.Specifies_Parameter("figattr_make_full_Fig"))>=0) figattr_make_full_Fig = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_make_zoom_Fig"))>=0) figattr_make_zoom_Fig = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_make_connections_Fig"))>=0) figattr_make_connections_Fig = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_make_abstract_Fig"))>=0) figattr_make_abstract_Fig = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figattr_make_neurons_Figs"))>=0) figattr_make_neurons_Figs = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figureTeXwidth"))>=0) tex_textwidth = atof(clp.ParValue(n));
  bool figsequenceflag=false;
  if ((n=clp.Specifies_Parameter("figuresequence"))>=0) figsequenceflag = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("figureTeXwidth"))>=0) tex_textwidth = atof(clp.ParValue(n));
  double sampleinterval = 24.0*60.0*60.0;
  if ((n=clp.Specifies_Parameter("sample_dt"))>=0) {
    sampleinterval = atof(clp.ParValue(n));
    if (sampleinterval==0.0) {
      sampleinterval = 24.0*60.0*60.0;
      warning("Warning: Parameter sample_dt was set to 0.0, defaulting to "+String(sampleinterval,"%.1f")+" seconds.\n");
    }
  }
  colortable->parse_CLP(clp);
  report_form_aware(clp,colortable);
  report('\n');
  // NETWORK PARAMETERS
  if ((n=clp.Specifies_Parameter("days"))>=0) max_growth_time = atof(clp.ParValue(n))*24.0*60.0*60.0;
  if ((n=clp.Specifies_Parameter("seconds"))>=0) max_growth_time = atof(clp.ParValue(n));
  unsigned int rseed = 0;
  if ((n=clp.Specifies_Parameter("randomseed"))>=0) rseed = atol(clp.ParValue(n));

  // *** (1) Put in all the neuron number setting code
  // *** (2) Move it to Network_Statistics or some other appropriate place
  // *** --> Can move to Network_Statistics::initialize();
  // Neuron type proportions set by approximate proportion parameters
  String parname;
  double approxproportion[UNTYPED_NEURON+1];
  for (int i = 0; i<=UNTYPED_NEURON; i++) approxproportion[i] = 0.0;
  approxproportion[PYRAMIDAL] = 0.7; // default value
  approxproportion[INTERNEURON] = 1.0 - approxproportion[PYRAMIDAL]; // default value
  for (int i = 0; i<=UNTYPED_NEURON; i++) {
    parname = "approxproportion";
    parname += neuron_type_name[i];
    parname.gsub(" ","_");
    parname.downcase();
    if ((n=clp.Specifies_Parameter(parname))>=0) {
      approxproportion[i] = atof(clp.ParValue(n));
      //      cout << "approxproportion[" << i << "] = " << approxproportion[i] << '\n';
      if (approxproportion[i] < 0.0) {
	approxproportion[i] = 0.0;
	warning("Warning: The approximate proportion of neurons that are "+String(neuron_type_name[i])+" was set to < 0, defaulting to 0.0.\n");
      }
    }
  }
  for (int j = UNTYPED_NEURON; j>=0; j--) {
    double approxproportiontotal = approxproportion[0];
    for (int i = 1; i<=UNTYPED_NEURON; i++) approxproportiontotal += approxproportion[i];
    if (approxproportiontotal == 1.0) break;
    double newtotal = approxproportiontotal - approxproportion[j];
    approxproportion[j] = 1.0 - newtotal;
    if (approxproportion[j]<0.0) approxproportion[j] = 0.0;
    warning("Warning: The sum of the proportions of specific types of neurons is "+String(approxproportiontotal,"%.3f")+". adjusting the proportion of "+neuron_type_name[j]+" to "+String(approxproportion[j],"%.3f")+'\n');
  }
  // Neuron type numbers set directly
  int numneurons = 0;
  int populationsize[UNTYPED_NEURON+1];
  for (int i = 0; i<=UNTYPED_NEURON; i++) populationsize[i] = 0;
  for (int i = 0; i<=UNTYPED_NEURON; i++) {
    parname = "populationsize";
    parname += neuron_type_name[i];
    parname.gsub(" ","_");
    parname.downcase();
    //cout << "Seeking " << parname << '\n';
    if ((n=clp.Specifies_Parameter(parname))>=0) {
      populationsize[i] = atoi(clp.ParValue(n));
      if (populationsize[i] < 0) {
	populationsize[i] = 0;
	warning("Warning: The population size for neurons that are "+String(neuron_type_name[i])+" was set to < 0, defaulting to 0.\n");
      }
      numneurons += populationsize[i];
      report("  Absolute size of "+String(neuron_type_name[i])+" population requested = "+String((long) populationsize[i])+" (total="+String((long) numneurons)+")\n");
    }
  }
  bool exactpopulationsizes = (numneurons > 0);
  // Number of neurons of all types
  numneurons_is_default = false;
  if (numneurons == 0) {
    numneurons = 9; // default value
    if ((n=clp.Specifies_Parameter("neurons"))>=0) {
      numneurons = atoi(clp.ParValue(n));
      if (numneurons<1) {
	numneurons = 9;
	numneurons_is_default = true;
	warning("Warning: The number of neurons was set to < 1, defaulting to "+String((long) numneurons)+" neurons.\n");
      }
    } else {
      numneurons_is_default = true;
      report("Note: Using default network size (9 neurons).\n");
    }
  }
  if ((n=clp.Specifies_Parameter("fibreswithturns"))>=0) fibreswithturns = (downcase(clp.ParValue(n))==String("true"));
  double timesteps = 0.0;
  if ((n=clp.Specifies_Parameter("dt"))>=0) timesteps = atof(clp.ParValue(n));
  // GROWTH FUNCTION PARAMETERS
  if ((n=clp.Specifies_Parameter("partitionsubsetsize"))>=0) {
    spatialsegmentsubsetsizelimit = atoi(clp.ParValue(n));
    if (spatialsegmentsubsetsizelimit<1) {
      spatialsegmentsubsetsizelimit = 65;
      warning("Warning: The preferred maximum number of segments per spatial partition was\nset to < 1, defaulting to 65 fibre segments.\n");
    }
  }
  if ((n=clp.Specifies_Parameter("branchinsegment"))>=0) branchsomewhereinsegment = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("branchatinitlength"))>=0) initiallengthatactualfirstbranch = (downcase(clp.ParValue(n))==String("true"));
  global_parse_CLP(clp);
  if (maxnumberofspatialsegmentsubsets<INT_MAX) {
    unsigned int sssmaxmem = (sizeof(Spatial_Segment_Subset)+(sizeof(fibre_segment_ptr)*spatialsegmentsubsetsizelimit))*maxnumberofspatialsegmentsubsets;
    report("\nMemory consumption by Spatial Segment Subsets limited to: "+String(((double) sssmaxmem)/(1024.0*1024.0),"%.3f")+" Mbytes + memory allocation overhead\n");
  }

  // PARAMETER RECORDING BY PRINTING TO STDOUT
  report("Generating a network with "+String((long) numneurons)+" general population neurons.\n(See below for additional region specific numbers of neurons of specific type added.)\n");
  report("Creating presynaptic and postsynaptic structure over "+String(max_growth_time/(24.0*60.0*60.0),"%.2f")+" days.\n");
  if (fibreswithturns) report("Dendrite and axon fibre segments include turns at continuation points.\n");
  else report("Dendrite and axon fibre segments are straight between bifurcation points.\n");
  if (branchsomewhereinsegment) report("The branch points at bifurcating segments occur anywhere between steps of the growth function.\n");
  else report("The branch points at bifurcating segments occur only at the end of a terminal segment during a growth function step.\n");
  if (initiallengthatactualfirstbranch) report("A branch point must occur at the initial segment length.\n");
  else report("Probabilistic branch point calculations apply at the initial segment length.\n");
  report(global_report_parameters());
  if (figsequenceflag) report("A sequence of figures is produced during simulation.\n");
  //if (figsequenceflag && combineflag) cout << "The sequence is subsequently combined into an animation.\n";
  // [***INCOMPLETE] add clp test for evaluate_connectivity_immediately

  if (rseed<1) set_time_random_seed();
  else set_random_seed(rseed);

  spatial smin(-3.0,-3.0,-3.0), smax(3.0,3.0,3.0);
  principal pspreadmin(smin); pspreadmin.set_radius(4.0); // sample sets minimums of parameter spread
  principal pspreadmax(smax); pspreadmax.set_radius(6.0); // sample sets maximums of parameter spread
  // *** --> Can move to Network_Statistics::initialize();
  Network_Statistics_Root nsr(exactpopulationsizes);
  if (nsr.UseExactPopulationSizes()) for (int i = 0; i<=UNTYPED_NEURON; i++) nsr.link_before(new Network_Statistics(static_cast<neuron_type>(i),populationsize[i]));
  else for (int i = 0; i<=UNTYPED_NEURON; i++) nsr.link_before(new Network_Statistics(static_cast<neuron_type>(i),approxproportion[i]));
  report(nsr.head()->all_str()+'\n');
  // [***NOTE] The Connection Statistics, i.e. statistics that indicate the
  // desired connectivity, are currently not in use.
  Connection_Statistics_Root cstats;
  cstats.link_before(new Connection_Statistics(PRINCIPAL_NEURON,PRINCIPAL_NEURON,20,28,37.0)); 
  cstats.link_before(new Connection_Statistics(PRINCIPAL_NEURON,INTERNEURON,17,23,232.0)); 
  cstats.link_before(new Connection_Statistics(INTERNEURON,PRINCIPAL_NEURON,60,88,82.0)); 
  cstats.link_before(new Connection_Statistics(INTERNEURON,INTERNEURON,5,7,82.0));
  // selecting growth functions
  van_Pelt_dendritic_growth_model * vPdgm = new van_Pelt_dendritic_growth_model(true,clp);
  van_Pelt_dendritic_growth_model * vPagm = new van_Pelt_dendritic_growth_model(false,clp);
  if (timesteps>0.0) {
    vPdgm->set_fixed_time_step_size(timesteps);
    vPagm->set_fixed_time_step_size(timesteps);
  }
  report_form_aware(clp,vPdgm);
  report_form_aware(clp,vPagm);
  general_fibre_structure_parameters.parse_CLP(clp);
  report_form_aware(clp,&general_fibre_structure_parameters);
  general_direction_model_parameters.parse_CLP(clp); // detect some direction model stuff
  report_form_aware(clp,&general_direction_model_parameters);
  general_synapse_formation_model_parameters.parse_CLP(clp);
  report_form_aware(clp,&general_synapse_formation_model_parameters);
#ifdef VECTOR3D
  general_slice_parameters.parse_CLP(clp);
  report_form_aware(clp,&general_slice_parameters);
#endif
  report("The simulation time step size is "+String(vPdgm->get_fixed_time_step_size(),"%.1")+" seconds.\n");
  report("The output sample interval size is "+String(sampleinterval,"%.1")+" seconds.\n");
  report("The preferred maximum number of segments per spatial partition is "+String((long) spatialsegmentsubsetsizelimit)+".\n");
  //reduce_rand_calls = true;
  bool autoplot = false;
  if (clp.IsFormInput()) autoplot = true;
  String samplefigfiles(outputdirectory+"sample"), netfigfile(outputdirectory+"net.fig"), zoomfigfile(outputdirectory+"zoom.fig"), connectionexamplefigfile(outputdirectory+"connections-example.fig"), abstractfigfile(outputdirectory+"abstract.fig"), neuronfigfilesbase(outputdirectory+"neuron-XXXXX.fig"), statsdatafile(outputdirectory+"stats.m"), statsfile(outputdirectory+"nibr.m"), nettxtfile(outputdirectory+"net.txt"), netvrmlfile(outputdirectory+"net.x3d"), netccmfile(outputdirectory+"net.ccm"), slicefile(outputdirectory+"slice"), netdatasyndistfile(outputdirectory+"netdata-synapsedistance.m"), netdataconndistfile(outputdirectory+"netdata-connectiondistance.m"), faninfigfile(outputdirectory+"fanin.fig"), faninhistfile(outputdirectory+"fanin.m");

  // LOAD TEST FOR OPEN FORM ACCESS
  if (clp.IsFormInput()) {
    double load = ((double) numneurons) * max_growth_time / vPdgm->get_fixed_time_step_size();
    if (load > 2000000) {
      cout << "</PRE><P><B><FONT COLOR=\"FF0000\">The computational requirements for the specified combination of the number of neurons, simulated duration of development and step size exceed the load permitted through the open internet form interface.<P>To run simulations of this or greater load, please contact Randal A. Koene as indicated at <A HREF=\"http://rak.minduploading.org\">rak.minduploading.org</A>.</FONT></B><PRE>\n";
      return;
    }
  }
  // create network and electrode array
  network * net = new network(numneurons,&pspreadmin,&pspreadmax,nsr); // ,false); // 1st call that can create neurons
  eq = net;

  // PARAMETERS THAT REQUIRE KNOWLEDGE OF THE NETWORK
  net->parse_CLP(clp);
  report_form_aware(clp,net);
  Shape_Result res = net->shape_network(clp,&pspreadmin,&pspreadmax); // 2nd call that can create neurons
  if (net->PLLRoot<neuron>::length()<1) error("Error: The network does not appear to contain any neurons!\n");
  general_neuron_parameters.parse_CLP(clp,*net); // detect main direction/elongation/neuron model stuff
  report_form_aware(clp,&general_neuron_parameters);
  /* *** uniform random connectivity or growth choice */
  report(res.message);
  progress("Network creation completed.\n");
  spatial net_zoom_center(net->Center());
  double net_zoom_disttoedge = 80.0, net_zoomwidth=0.0, net_zoomheight=1.0, net_zoomdepth=1.0;
  parse_focus_box_CLP("net",clp,net_zoom_center,net_zoom_disttoedge);
  parse_focus_volume_CLP("net",clp,net_zoomwidth,net_zoomheight,net_zoomdepth);
  if (net_zoomwidth>0.0) net_zoom_disttoedge = net_zoomwidth;
#ifdef VECTOR3D
  spat_map = spat_pres = new Spatial_Presentation(net->Center());
  spat_pres->parse_CLP(clp);
  report_form_aware(clp,spat_pres);
#endif
  bool electrodes = false;
  if ((n=clp.Specifies_Parameter("electrodes"))>=0) electrodes = (downcase(clp.ParValue(n))==String("true"));
  electrodearray * els = NULL;
  if (electrodes) {
    els = make_electrodes_hexagon(net->Center());
    progress("Electrode placement completed.\n");
  }
  if (figsequenceflag) sampled_output = new Sampled_Growth_Fig_Output(net,max_growth_time,sampleinterval,samplefigfiles,statsdatafile,autoplot,samplefigfiles,res.displaywidth,false,false,true);
  else sampled_output = new Sampled_Growth_Output(net,max_growth_time,sampleinterval,samplefigfiles,statsdatafile,autoplot,true);
  sampled_output->parse_CLP(clp);
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  if (!sampled_output->Collecting_Statistics()) { // for sanity
    outattr_show_stats = false;
  }
#endif
  if ((outattr_track_nodegenesis) && (outattr_make_full_Txt)) {
    int nlsize = net->PLLRoot<neuron>::length()*64;
    if (Txt_tuftrootbranchnodelist) delete Txt_tuftrootbranchnodelist;
    if (Txt_obliquerootbranchnodelist) delete Txt_obliquerootbranchnodelist;
    if (outattr_Txt_separate_files) {
      Txt_tuftrootbranchnodelist = new String(COLUMN_LABELS_TUFT_NODES);
      Txt_obliquerootbranchnodelist = new String(COLUMN_LABELS_OBLIQUE_NODES);
    } else {
      Txt_tuftrootbranchnodelist = new String();
      Txt_obliquerootbranchnodelist = new String();
    }
    Txt_tuftrootbranchnodelist->alloc(nlsize);
    Txt_obliquerootbranchnodelist->alloc(nlsize*8);
  }
  bool combineflag = false;
  if (figsequenceflag) combineflag = static_cast<Sampled_Growth_Fig_Output *>(sampled_output)->Combine();
  report_form_aware(clp,sampled_output);
  system("rm -f "+samplefigfiles+"*fig "+samplefigfiles+"*png");
  if (pb.select_physical_boundary("",clp)>0) report_form_aware(clp,&pb);
  
  // select synapse formation model (*** currently only one available)
  sfm = new synapse_formation_model(net,&cstats);

  // generate presynaptic and postsynaptic structure
  net->develop_connection_structure(cstats,*vPdgm,*vPagm);
  if (warning_fibre_segment_net_Fig_zero_length>0) warning("nibr Warning: "+String((long) warning_fibre_segment_net_Fig_zero_length)+" occurrences of fibre_segment with length == 0.0 in fibre_segment::net_Fig()\n");
  if (warning_Spatial_Segment_Subset_intersects_zero_length>0) warning("nibr Warning: "+String((long) warning_Spatial_Segment_Subset_intersects_zero_length)+" occurrences of zero length segment in Spatial_Segment_Subset::intersects()\n");
  progress("Connections generated (total="+String((long) net->number_of_connections())+").\n");

#ifdef TESTING_SPIKING
  // [***INCOMPLETE] The following temporary, as the proper implementation
  // should allow specific neuron types or even specific neurons to accept
  // appropriate options for activity models.
  //actres = new Activity_Results_Output();
  if (eq->T_End()>max_growth_time) {
    report("Simulating "+String((eq->T_End()-max_growth_time),"%.1f")+" additional seconds of ACTIVITY in the generated network.\n");
    actres = new ARO_to_File(outputdirectory+"activity.ascii");
    net->abstract_connections();
    PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->set_activity(new IFactivity(e));
    
    eq->add(new Random_Spikes()); // first random spike is immediate
    eq->variablesteps();
    actres->postop();
    delete actres;
  } else report("NOT simulating additional ACTIVITY events following network generation (parameters: event_seconds <= seconds).\n");
#endif

  // output results data
  progress("Counting segments in final network structure...\n");
  int totaldendriticsegments = 0, totalaxonalsegments = 0;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    totaldendriticsegments += e->total_input_segments();
    totalaxonalsegments += e->total_output_segments();
  }
  progress("Total number of dendritic segments = "+String((long) totaldendriticsegments)+" (average per neuron = "+String(((double) totaldendriticsegments)/((double) net->PLLRoot<neuron>::length()),"%.2f")+")\n");
  progress("Total number of axonal segments    = "+String((long) totalaxonalsegments)+" (average per neuron = "+String(((double) totalaxonalsegments)/((double) net->PLLRoot<neuron>::length()),"%.2f")+")\n");
  if (net->Seek_Candidate_Synapses()) report("Proximity threshold used to establish candidate synapses = "+String(sfm->Proximity_Threshold(),"%.3f")+" microns\n");
  else report("NO candidate synaptic sites were sought in this simulation, due to the chosen parameter value 'candidate_synapses=false'.\n");
#ifdef EVALUATE_POSSIBLE_CONNECTION_PROFILING
  if (net->Seek_Candidate_Synapses()) {
    double brute_force = ((double) totaldendriticsegments)*((double) totalaxonalsegments);
    progress("Evaluated possible candidate synapses = "+String((long) evaluate_possible_connection_calls)+" (compared to "+String(brute_force,"%.3f")+" by brute force, i.e. "+String(((((double) evaluate_possible_connection_calls)/brute_force)*100.0),"%.3f")+" percent)\n");
    progress("Created candidate synapses = "+String((long) candidate_synapses_created)+"\n");
    progress("Minimum distance found during candidate synapses tests = "+String(candidate_synapses_minimum_distance,"%.5f")+"\n");
    progress("Number of candidate synapses disstance checks = "+String((long) candidate_synapses_num_distance_checks)+"\n");
    String hist_str;
    for (size_t i = 0; i<10; i++) hist_str += String((long) candidate_synapse_distances_histogram[i])+" ";
    progress("Histogram of candidate synapse distances (steps of 0.1) = "+hist_str+"\n");
    hist_str = "";
    for (size_t i = 0; i<10; i++) hist_str += String((long) threshold_distances_histogram[i])+" ";
    progress("Histogram of threshold distances (steps of 0.1) = "+hist_str+"\n");
  }
#endif
#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
  if (net->Seek_Candidate_Synapses()) {
    progress("Inventory of synapse genesis and loss:\n");
    for (int i=0; i<syntype_IDs; i++) progress("  "+String(synapse_type_name[i])+"\tgenesis="+String((long) synapse_inventory[i][SYNGENESIS])+"\tloss="+String((long) synapse_inventory[i][SYNLOSS])+'\n');
    progress("Attempted conversions of candidate synapses = "+String((long) candidates_conversion_attempted)+'\n');
  }
#endif
  if (empty_data_samples>0) warning("During the collection of simulation statistics, "+String((long) empty_data_samples)+" samples were empty.\n");
  int showprogresscache = 0;
#ifdef TEST_FOR_NAN
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps)
      PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
      if (isnan(ts->TerminalSegment()->P0.X())) {cout << "NAN before Figure calls!\n"; cout.flush();}
  }
#endif
  if (outattr_show_figure) {
    if (figattr_show_abstract_connections) net->abstract_connections();
    if (figattr_make_full_Fig) {
      set_focus_box(net->Center(),-1.0);
      net->Fig_Output(netfigfile,res.displaywidth,true);
      if (electrodes) {
	showprogresscache=figattr_sample_show_progress;
	figattr_sample_show_progress = 0;
	els->Fig_Output(netfigfile);
	figattr_sample_show_progress=showprogresscache;
      }
#ifdef VECTOR3D
      if (general_slice_parameters.show_slice_outlines()) net->Slice_Outlines_Output(netfigfile,general_slice_parameters);
#endif
      if (figattr_connections_thresholdlistlength>1) { // make a set of connection figures
	bool fsncache = figattr_show_neurons, fspscache = figattr_show_presynaptic_structure;
	figattr_show_neurons = true; figattr_show_presynaptic_structure = false;
	net->set_figattr(FIGATTR_SHOWCELL & FIGATTR_SHOWPRESYN);
	bool fsaccache = figattr_show_abstract_connections;
	double fctcache = figattr_connections_threshold;
	figattr_show_abstract_connections = true;
	String netfigfilebase(netfigfile.before(".fig"));
	for (long i = 0; i<figattr_connections_thresholdlistlength; i++) {
	  figattr_connections_threshold = figattr_connections_thresholdlist[i];
	  net->Fig_Output(netfigfilebase+"-threshold-"+String(figattr_connections_threshold,"%08.2f.fig"),res.displaywidth,true);
	}
	figattr_show_abstract_connections = fsaccache;
	figattr_connections_threshold = fctcache;
	figattr_show_neurons = fsncache; figattr_show_presynaptic_structure = fspscache;
	net->set_figattr(FIGATTR_SHOWALL);
      }
    }
#ifdef TEST_FOR_NAN
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1)
    PLL_LOOP_FORWARD_NESTED(fibre_structure,e->OutputStructure()->head(),1,ps)
      PLL_LOOP_FORWARD_NESTED(terminal_segment,ps->TerminalSegments()->head(),1,ts) {
      if (isnan(ts->TerminalSegment()->P0.X())) {cout << "NAN after net.fig call!\n"; cout.flush();}
  }
#endif
    if (figattr_make_connections_Fig) {
      neuron * zoomneuron = net->nearest(net->Center());
      if (zoomneuron) {
	net->set_figattr(FIGATTR_SHOWNOCONN);
	net->visible_pre_and_target_post(*zoomneuron);
	net->Fig_Output(connectionexamplefigfile,res.displaywidth,true);
	if (electrodes) {
	  showprogresscache=figattr_sample_show_progress;
	  figattr_sample_show_progress = 0;
	  els->Fig_Output(connectionexamplefigfile);
	  figattr_sample_show_progress=showprogresscache;
	}
	net->set_figattr(FIGATTR_SHOWALL);
      }
    }
    if (figattr_make_zoom_Fig) {
      if (net_zoomwidth>0.0) set_focus_volume(net_zoom_center,net_zoomwidth,net_zoomheight,net_zoomdepth);
      else set_focus_box(net_zoom_center,net_zoom_disttoedge);
      //usewidth = figattr_focus_box_x2-figattr_focus_box_x1;
      net->Fig_Output(zoomfigfile,figattr_focus_box_x2-figattr_focus_box_x1,true);
      if (electrodes) {
	showprogresscache=figattr_sample_show_progress;
	figattr_sample_show_progress = 0;
	els->Fig_Output(zoomfigfile);
	figattr_sample_show_progress=showprogresscache;
      }
      if (figattr_connections_thresholdlistlength>1) { // make a set of connection figures
	bool fsncache = figattr_show_neurons, fspscache = figattr_show_presynaptic_structure;
	figattr_show_neurons = true; figattr_show_presynaptic_structure = false;
	net->set_figattr(FIGATTR_SHOWCELL & FIGATTR_SHOWPRESYN);
	bool fsaccache = figattr_show_abstract_connections;
	double fctcache = figattr_connections_threshold;
	figattr_show_abstract_connections = true;
	String netfigfilebase(zoomfigfile.before(".fig"));
	for (long i = 0; i<figattr_connections_thresholdlistlength; i++) {
	  figattr_connections_threshold = figattr_connections_thresholdlist[i];
	  net->Fig_Output(netfigfilebase+"-threshold-"+String(figattr_connections_threshold,"%08.2f.fig"),res.displaywidth,true);
	}
	figattr_show_abstract_connections = fsaccache;
	figattr_connections_threshold = fctcache;
	figattr_show_neurons = fsncache; figattr_show_presynaptic_structure = fspscache;
	net->set_figattr(FIGATTR_SHOWALL);
      }
    }
    if (figattr_make_abstract_Fig) {
      if (!figattr_show_abstract_connections) net->abstract_connections();
      set_focus_box(net->Center(),-1.0);
      net->Fig_Abstract_Connections(abstractfigfile,res.displaywidth,true);
      if (figattr_connections_thresholdlistlength>1) { // make a set of connection figures
	double fctcache = figattr_connections_threshold;
	String abstractfigfilebase(abstractfigfile.before(".fig"));
	for (long i = 0; i<figattr_connections_thresholdlistlength; i++) {
	  figattr_connections_threshold = figattr_connections_thresholdlist[i];
	  net->Fig_Abstract_Connections(abstractfigfilebase+"-threshold-"+String(figattr_connections_threshold,"%08.2f.fig"),res.displaywidth,true);
	}
	figattr_connections_threshold = fctcache;
      }
    }
    if (figattr_make_neurons_Figs) {
      set_focus_box(net->Center(),-1.0);
      net->Fig_Neurons(neuronfigfilesbase,res.displaywidth);
    }
  }
  if (outattr_make_full_Txt) net->Txt_Output(nettxtfile);
  if (outattr_make_full_X3D) {
    if (net_zoomwidth>0.0) set_focus_volume(net_zoom_center,net_zoomwidth,net_zoomheight,net_zoomdepth);
    else set_focus_box(net_zoom_center,net_zoom_disttoedge);
    net->VRML_Output(netvrmlfile);
  }
  if (outattr_make_full_Catacomb) net->Catacomb_Output(netccmfile);
  if (outattr_synapse_distance_frequency) net->Data_Output_Synapse_Distance(netdatasyndistfile);
  if (outattr_connection_distance_frequency) net->Data_Output_Connection_Distance(netdataconndistfile);
  if (max_abstract_strength>0.0) progress("Maximum abstract connection strength computed = "+String(max_abstract_strength,"%.3f")+'\n');
#ifdef VECTOR3D
  if (general_slice_parameters.get_slice()) net->Slice_Output(slicefile,general_slice_parameters);
#endif
  // [***BEWARE] Unless a network copy is made, the following optional analysis changes the network.
  if (statsattr_fan_in_analysis!=NUM_fs) { // [***NOTE] This should be moved into a separate function.
    spatial * origpos = new spatial[net->PLLRoot<neuron>::length()];
    int n_pos = 0;
    PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) origpos[n_pos++] = e->Pos();
    spatial cpos;
    net->move(cpos);
    net->fanin_rot();
    n_pos = 0;
    PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) e->move(origpos[n_pos++]);
    delete[] origpos;
    net_fanin_histogram.means();
    write_file_from_String(faninhistfile,net_fanin_histogram.octave_output());
    //set_focus_box(cpos,-1.0);
    set_focus_box(net->Center(),-1.0);
    net->Fig_Output(faninfigfile,res.displaywidth,true);
    append_file_from_String(faninfigfile,net_fanin_histogram.fig_output()->str());
  }

  //#ifdef VECTOR3D
  //  PLL_LOOP_FORWARD(neuron,net->head(),1) cout << e->Pos().X() << ',' << e->Pos().Y() << ',' << e->Pos().Z() << '\n'; 
  //#endif

  /*  int nn = 0;
  PLL_LOOP_FORWARD(neuron,net->PLLRoot<neuron>::head(),1) {
    if (e->activity()==&NullActivity) {
      cout << "Neuron " << nn << " uses NullActivity\n";
      cout.flush();
    }
    nn++;
    }*/

  if (zerocoordinateconversions>0) progress("Conversion from Euler to Spherical coordinates was attempted for "+String((long) zerocoordinateconversions)+" ZERO vectors.\n");

  // clean up to make more memory available
  delete els;
  delete net;
  delete vPdgm;
  delete vPagm;
  sampled_output->Null_Net();

  // Place HTML output in advance
  String combinetype;
  if (clp.IsFormInput()) {
    if (outattr_show_figure && figattr_make_full_Fig) progress("The resulting network: <A HREF=\""+absoluteoutputURL+"net.pdf\">net.pdf</A>\n");
    if (outattr_show_figure && figattr_make_zoom_Fig)  progress("An enlargement of the center of the network: <A HREF=\""+absoluteoutputURL+"zoom.pdf\">zoom.pdf</A>\n");
    if (outattr_show_figure && figattr_make_abstract_Fig) progress("The abstracted connections: <A HREF=\""+absoluteoutputURL+"abstract.pdf\">abstract.pdf</A>\n");
    if (figsequenceflag && combineflag) {
      combinetype = static_cast<Sampled_Growth_Fig_Output *>(sampled_output)->CombineType();
      progress("The combined and animated sequence: <A HREF=\""+absoluteoutputURL+"sample."+combinetype+"\">sample."+combinetype+"</A>\n");
    }
    if (outattr_show_figure && figattr_make_neurons_Figs)  progress("A set of morphology .fig plots of individual neurons has also been prepared with base name "+neuronfigfilesbase+".\n");
    //if (outattr_show_stats) cout << "Resulting network statistics: <A HREF=\""+absoluteoutputURL+"stats.png\">stats.png</A>\n";
    progress("(Please wait for network output to complete before accessing files at the links above.)\n");
  }

  sampled_output->Output(outattr_show_stats);

  // more clean up
  if (Txt_tuftrootbranchnodelist) delete Txt_tuftrootbranchnodelist;
  if (Txt_obliquerootbranchnodelist) delete Txt_obliquerootbranchnodelist;
  delete sampled_output;
#ifdef VECTOR3D
  delete spat_pres;
#endif

#ifdef TESTING_RANDOM_DISTRIBUTION
  cout << "TESTING_RANDOM_DISTRIBUTION:\n";
  for (int ndi = 0; ndi<11; ndi++) cout << randomdist[ndi] << ' ';
  cout << '\n';
#endif
#ifdef TESTING_AEM_RANDOM_DISTRIBUTION
  cout << "TESTING_AEM_RANDOM_DISTRIBUTION:\n";
  for (int ndi = 0; ndi<21; ndi++) cout << aemrandomdist[ndi] << ' ';
  cout << '\n';
#endif

  if (clp.IsFormInput()) devsim_present_HTML_graphical_output(statsdatafile,netfigfile,zoomfigfile,abstractfigfile,statsfile,figsequenceflag,combineflag,combinetype);
  // some additional output

  // *** To create a truly realistic implementation of the van Pelt
  // algorithm, I would have to move away from the approach in which
  // all neurons and terminal segments are treated at the same time points.
  // Afterall, if a terminal branches at a different time than other
  // terminals, the new branches continue to elongate, while others are
  // still approaching their branch points. It would be necessary to
  // set up a stack of terminal segments in the order of their time step
  // treatment, along with necessary information about the tree they are
  // in. In this case, it may be best to focus immediately on a method that
  // does away with time-steps and instead tries to predict branching
  // events in a manner simialr to the event prediction method.
  // (Perhaps for now, the current implementation produces output that is
  // not significantly less realistic than that produced by the van Pelt
  // algorithms.)
  // *** Ask van Pelt if it is the meaning of the algorithm that it be
  // applied given a certain time for the initial segments, or that it
  // only start after those segments, i.e. that the time at the onset of
  // growth past the initial segments is 0.
  // *** There is probably a standard deviation of about 50% on nu0 that
  // needs to be built into the growth functions by multiplying their
  // return values with random_double(-std_nu0,std_nu0), but ask van Pelt
  // if that variance applies for each segment, or for a tree, or for
  // a neuron (i.e. is std_nu0 determined at a specific level, how variable
  // is it over the network, and what is the function that determines
  // the growth speed?). Note that there may be a resource supply function
  // that could be affected by input such as nutrition, blood flow, etc.
  // *** If growth includes periodic motility that leads to changes in
  // the angle or length of segments within the tree (other than terminal
  // segments) then coordinates or normal angles of the segments need
  // to be refreshed before drawing them or before computing a new
  // normal angle at a continuation point.
  // *** normal_random():
  // mean + sqrt(variance)*randn() uses F77 normal_dist

  progress("DEVELOPMENTAL (MORPHOGENESIS) SIMULATION: Successful, clean exit.\n");
}

const char helpstr[] = "\
Usage: netmorph [parameter=value]\n\
       netmorph2D [parameter=value]\n\
\n\
  parameters:\n\
\n\
    Simulation Context (global):\n\
\n\
    seconds=<decimal#>                   duration of simulated growth time\n\
    days=<decimal#>                      duration of simulated growth time\n\
    events_seconds=<decimal#>            duration of simulated events\n\
    events_append_seconds=<decimal#>     append event time to growth time\n\
    maxrandomspikeinterval=<decimal#>    max. interval between random spikes\n\
    randomseed=<integer#/0>              random seed (0 = use time value)\n\
    dt=<decimal#>                        time step size (seconds)\n\
\n\
    Complete simulated network:\n\
\n\
    neurons=<integer#>                   number of neurons in network\n\
    approxproportion<type>=<decimal#>    approximate proportion of neurons\n\
                                         that are of type:\n\
                                         principal_neuron, interneuron,\n\
                                         pyramidal, multipolar_nonpyramidal,\n\
                                         bipolar, untyped_neuron\n\
    populationsize<type>=<integer#>      number of neurons of type:\n\
                                         principal_neuron, interneuron,\n\
                                         pyramidal, multipolar_nonpyramidal,\n\
                                         bipolar, untyped_neuron\n\
    random_orientation=<true/decimal#>   orientation of axon and dendrites\n\
    use_specified_basal_direction=<true/false> align as specified\n\
    specified_basal_direction_theta=<decimal#>\n\
    specified_basal_direction_phi=<decimal#>\n\
    specified_basal_direction_relx=<decimal#>\n\
    specified_basal_direction_rely=<decimal#>\n\
    specified_basal_direction_relz=<decimal#>\n\
    apical_specified_overrides_centroid=<true/false> use negated specified\n\
    shape=<shape name>                   a shape selection: rectangle, box,\n\
                                         regions, circle, hexagon\n\
\n\
    Network environments:\n\
\n\
    shape_horizontal=<decimal#>          x dimension length of network shape\n\
    shape_vertical=<decimal#>            y dimension length of network shape\n\
    shape_depth=<decimal#>               z dimension length of network shape\n\
    shape_sidelength=<decimal#>          length of a network (hexagon) side\n\
    shape_radius=<decimal#>              radius of a (circular) network\n\
    shape_randompolar=<true/false>       randomize allocation in polar crds.\n\
\n\
    regions=<space separated list>       list of region labels\n\
    <region>.<neuron-type>=<integer#>    specify additional neurons of type:\n\
                                         principal_neuron, interneuron,\n\
                                         pyramidal, multipolar_nonpyramidal,\n\
                                         bipolar\n\
    <region>.shape=<disc/box/sphere>     shape of a specific region\n\
    <region>.neurons=<integer#>          number of neurons in region\n\
    <region>.minneuronseparation=<decimal#> minimum separation > cell sizes\n\
    <region>.center<X/Y/Z>=<decimal#>    region center location\n\
    [region.]shape.radius=<decimal#>     radius of disc/sphere region shape\n\
    [region.]shape.thickness=<decimal#>  thickness of disc region shape\n\
    [region.]shape.width=<decimal#>      horizontal size of box region shape\n\
    [region.]shape.height=<decimal#>     vertical size of box region shape\n\
    [region.]shape.depth=<decimal#>      z dimension size of box region shape\n\
\n\
    electrodes=<true/false>              include electrodes\n\
    pia_attraction_repulsion_hypothesis=<true/false> axon/apical locations\n\
    <neuron-type>.min_basal=<integer#>   minimum number of basal dendrites\n\
    <neuron-type>.max_basal=<integer#>   maximum number of basal dendrites\n\
    <neuron-type>.basal.minangle=<decimal#> minimum initial angular deviation\n\
    <neuron-type>.basal.maxangle=<decimal#> maximum initial angular deviation\n\
    <neuron-type>.basal.force_model=<label> placement: unrestricted,\n\
                                         surface_division\n\
    multipolar.max_axons=<integer#>      maximum number of axons\n\
\n\
    [label.]environment_physics=<label[ ... ]>   physical constraints\n\
      <label>.physical_boundary=<label>  spherical, point_attractor\n\
      spherical boundary:\n\
        <label>.spherical_center[X/Y/Z]=<neuron id#/decimal#>\n\
        <label>.spherical_radius=<decimal#>\n\
        <label>.spherical_maxrange=<decimal#>\n\
        <label>.spherical_c_repulse=<decimal#>\n\
      point attractor:\n\
        <label>.attractor_point[X/Y/Z]=<neuron id#/decimal#>\n\
        <label>.attractor_maxrange=<decimal#>\n\
        <label>.attractor_c_attract=<decimal#>\n\
\n\
    Morphological Development -  Growth Cone Elongation:\n\
\n\
    <schema.>arbor_elongation_model=     van_Pelt, polynomial_O1,\n\
      <model-id>                         polynomial_O3\n\
      van_Pelt parameters:\n\
        [set.][contributing.]growth_F,\n\
        [set.][contributing.]growth_nu0,\n\
        [set.][contributing.]F_competes_with  whole_neuron, all_axons,\n\
                                         all_dendrites, same_arbor\n\
      polynomial_O1 parameters:\n\
        [set.][contributing.]growth_F,\n\
        [set.][contributing.]growth_nu0\n\
      polynomial_O3 parameters:\n\
        [set.][contributing.]growth_F,\n\
        [set.][contributing.]growth_nu0,\n\
        [set.][contributing.]growth_nu1,\n\
        [set.][contributing.]growth_nu2\n\
      .aem.PDF=<PDF-id>                  delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
      .aem_weight\n\
      .arbor_elongation_model_label\n\
      .aem_label\n\
    <schema.>terminal_segment_elongation_model=  simple, inertia, second_order,\n\
      <model-id>                         constrained_second_order,\n\
                                         decaying_second_order,\n\
                                         initialized_cso, BESTL, nonnorm_BESTL,\n\
                                         pyrAD_BESTLNN\n\
      second_order parameters:\n\
        .tsem_inertia=<decimal#>         first order (speed) inertia\n\
      decaying_second_order parameters:\n\
        .tsem_decay=<decimal#>           decay time constant (seconds)\n\
      initialized_cso parameters:\n\
        .tsem_initial_acceleration=<decimal#> initial quota acceleration\n\
      BESTL parameters:\n\
        .tsem.branch.PDF=<PDF-id>        initial elongation at branches\n\
      pyrAD_BESTLNN parameters:\n\
        .tsem.trunklength.PDF=<PDF-id>   length of apical dendrite trunk\n\
        .tsem.obliques.PDF=<PDF-id>      number of oblique branches\n\
        .tsem.prefix                     label for oblique and tuft parameters\n\
          <prefix.>oblique.<model-parameters> various model parameters for\n\
                                         oblique branches, e.g.\n\
                                         terminal_segment_elongation_model\n\
          <prefix.>tuft.<model-parameters> various model parameters for\n\
                                         tuft branches, e.g.\n\
                                         terminal_segment_elongation_model\n\
      .tsem.PDF=<PDF-id>                 delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
      .tsem_weight\n\
      .terminal_segment_elongation_model_label\n\
      .tsem_label\n\
    <schema.>elongation_rate_initialization_model =\n\
      <model-id>                         length_distribution,\n\
                                         nonnorm_BESTL_length_distribution,\n\
                                         pure_stochastic, zero, unitary,\n\
                                         continue_defaults\n\
      unitary parameters:\n\
        .eri_initialquota                elongation quota to initialize each to\n\
      .eri.PDF=<PDF-id>                  delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
      .eri_weight\n\
      .elongation_rate_initialization_model_label\n\
      .eri_label\n\
\n\
    Morphological Development - Growth Cone Bifurcation:\n\
\n\
    <schema.>L0=<#decimal[,#decimal]>    range of initial arbor lengths\n\
\n\
    branchinsegment=<true/false>         branches occur between steps\n\
    branchatinitlength=<true/false>      first branch at initial length\n\
    Abranchesatturns=<true/false>        branches may occur at existing turns\n\
    Dbranchesatturns=<true/false>        branches may occur at existing turns\n\
\n\
    <schema.>branching_model=            van_Pelt, polynomial_O1, polynomial_O3\n\
      <model-id>                         \n\
      general parameters:\n\
        .min_node_interval=<decimal#>    min.length (microns) of segment pieces\n\
      van_Pelt parameters:\n\
        .B_inf=<decimal#>                asymptotic value\n\
        .tau=<decimal#>                  time coefficient (seconds)\n\
        .E=<decimal#>                    competition coefficient\n\
        .E_competes_with=                whole_neuron, all_axons,\n\
                                         all_dendrites, same_arbor\n\
      .bm_weight\n\
      .branching_model_label\n\
      .bm_label\n\
    <schema.>TSBM=<model-id>             van_Pelt, van_Pelt_specBM\n\
      van_Pelt parameters:\n\
        .S=<decimal#>                    Centrifugal order dependency\n\
      .tsbm_weight\n\
      .TSBM_label\n\
      .tsbm_label\n\
    <schema.>branch_angle_model=         Balanced_Forces, Fork_Main_Daughter\n\
      <model-id>                         \n\
      Balanced_Forces parameters:\n\
        .bam.bfbam.PDF=<PDF-id>          angle in the branching parallelogram\n\
      .bam.PDF=<PDF-id>                  delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
      .bam_weight\n\
      .branch_angle_model_label\n\
      .bam_label\n\
\n\
    Morphological Development - Growth Cone Change of Direction:\n\
\n\
    fibreswithturns=<true/false>         fibre segments with turns\n\
    <schema.>TSTM=<model-id>             terminal segment turn event selection\n\
                                         model: none, each_dt, linear_rate,\n\
                                         branch_coupled\n\
      linear_rate parameters:\n\
        .turn_rate                       turns per micrometer (or second)\n\
        .turn_separation                 microns (or seconds) between turns\n\
      branch_coupled parameters:\n\
        .tf=<decimal#>                   turnfactor likelihood multiplier\n\
      .tstm_weight\n\
      .TSTM_label\n\
      .tstm_label\n\
\n\
    <schema.>direction_model=<model-id>  segment_history_tension,\n\
                                         cell_attraction, radial, vector\n\
      .veeranglemin=<decimal#>           minimum and maximum pitch veer angles\n\
      .veeranglemax=<decimal#>           that perturb expected turn direction\n\
      .dm_weight\n\
      .direction_model_label\n\
      .dm_label\n\
    tension parameters:\n\
      .history_power=<decimal#>          history has increasing (<0),equal (0)\n\
                                         or decreasing (>0) influence\n\
    radial parameters:\n\
      .Samsonovich_hypothesis=<true/false> do turns return to a normal vector\n\
    vector parameters:\n\
      .direction=<vector>                expected direction vector\n\
    dirhistory_selection=                dynamic changes of the segment history\n\
      <none/random_truncation>           used in the direction model\n\
    general_rtdhs_probability=<decimal#> probability of history truncation\n\
    general_rtdhs_minfraction=<decimal#> min. truncated history fraction\n\
    general_rtdhs_maxfraction=<decimal#> max. truncated history fraction\n\
    turnanglemin=<decimal#>              minimum angle at continuation node\n\
    turnanglemax=<decimal#>              maximum angle at continuation node\n\
\n\
    Morphological Development - Neurite Fiber Diameter:\n\
\n\
    fibrediameter=<true/false>           use specific diameters\n\
    neurite_diameter_model=<label>       none, asymptotic, rall\n\
      Asymptotic Soma Proportion:\n\
      ndm.asymptotic_proportion=<decimal#>  proportion of soma diameter\n\
      ndm.asymptotic_lambda=<decimal#>      exponential lambda decay factor\n\
      Rall Power Law:\n\
      ndm.e_power.PDF=<PDF-id>           delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
      ndm.d_term.PDF=<PDF-id>            delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
\n\
    Connectivity Development - Synapse formation:\n\
\n\
    candidate_synapses=<true/false>      seek candidate synapses\n\
    no_autapses=<true/false>             (dis)allow synapses onto same neuron\n\
    partition_to_synapses=<true/false>   seek candidate syn. once partitioned\n\
    synapses_during_development=<true/false>     synapses during or after\n\
    partitionsubsetsize=<integer#>       max.segments preferred per partition\n\
    max_spatial_segment_subsets=<integer#> max. number of spat. seg. subsets\n\
    min_spatial_segment_span=<decimal#>  min. span of spatial segment subset\n\
    max_spatial_segment_coordinate=<decimal#> explicit max. coordinates\n\
    synapse_formation.PDF=<delta/uniform/  PDF used to decide if an actual\n\
       linear/spline_normal/               synapse is formed at a candidate\n\
       spline_normal_with_min/normal/      site\n\
       exponential>\n\
    D_synmax.<pre>.<post>=<decimal#>     Max. dist. between fibers at synapse\n\
                                         pre/post: principal, interneuron,\n\
                                         multipolar, pyramidal, bipolar,\n\
                                         untyped\n\
    P_receptor.<pre>.<post>.<receptor type>=<decimal#> receptor probability\n\
                                         pre/post: principal, interneuron,\n\
                                         multipolar, pyramidal, bipolar,\n\
                                         untyped\n\
                                         recentor type: AMPAR, NMDAR, GABAR\n\
\n\
    Simulation Output - Histological Slice Generation:\n\
\n\
    slice=<label>                        generate histological slices\n\
      <label>: multi-level labels (1.2.3.4), <id>: deepest label level\n\
      <label>.batchsize=<integer#>       create batch if greater than 1\n\
      <label.<id>>.slice=<new-id>        rename slice/batch <id> to <new-id>\n\
      <label>.<vertex>[x/y/z]=           vertex coordinates for vertices:\n\
        <decimal#>[,decimal#,decimal#]   <S/Z/I/z><A/Y/P/y><D/X/S/x>\n\
      <label>.relativecoords=<subindex>  coordinates relative to element\n\
                                         in this-<subindex> batch position\n\
      <label>.batchrelativecoords=<subindex>\n\
      <label><.nalpha/.ntransverse>=<decimal#>\n\
      <label><.nbeta/.ncoronal>=<decimal#>\n\
      <label><.ngamma/.nsagittal>=<decimal#>\n\
      <label>.origin=<vertex>            use <vertex> for relative offset\n\
      <label>.width=<decimal#>           width of slice (uses origin)\n\
      <label>.height=<decimal#>          height of slice (uses origin)\n\
      <label>.depth=<decimal#>           depth of slice (uses origin)\n\
\n\
    Simulation Output - Network Data:\n\
\n\
    outattr_make_full_Txt=<true/false>   create .txt of network structure\n\
    outattr_Txt_sequence=<true/false>    produce sequence of .txt data\n\
    outattr_Txt_separate_files=<true/false> store full_Txt in separate files\n\
    outattr_track_synaptogenesis=<true/false> track times of synapse creation\n\
    outattr_track_nodegenesis=<true/false> track times of node creation\n\
    outattr_make_full_X3D=<true/false>   create .x3d of network structure\n\
    outattr_make_full_Catacomb=<true/false> create .ccm of network structure\n\
    outattr_synapse_distance_frequency=<true/false> synapses by distance\n\
    outattr_connection_distance_frequency=<true/false> connection by distance\n\
    outattr_distance_frequency_distbinsize=<decimal#> distance bin size\n\
\n\
    Simulation Output - Statistical Data:\n\
\n\
    statsattr_collect_statistics=<true/false> collect statistics data\n\
    statsattr_<specific>=<true/false>    collect a specific statistic:\n\
                                         Dlength, Alength, Dtermsegsperarbor,\n\
                                         Atermsegsperarbor,Dterminallength,\n\
                                         Aterminallength,Dintermediatelength,\n\
                                         Aintermediatelength,\n\
                                         Dlengthbetweenbifurcations,\n\
                                         Alengthbetweenbifurcations,\n\
                                         Dtermlensincebifurcation,\n\
                                         Atermlensincebifurcation,\n\
                                         Dturnsbetweenbifurcations,\n\
                                         Aturnsbetweenbifurcations,\n\
                                         Dbranchangles, Abranchangles,\n\
                                         Dturnangles, Aturnangles,\n\
                                         Dtermlensincesoma, Atermlensincesoma,\n\
                                         Dcartratiobetweenbifurcations,\n\
                                         Acartratiobetweenbifurcations,\n\
                                         Dcartratiobifurcationtoterm,\n\
                                         Acartratiobifurcationtoterm,\n\
                                         Dcartratiosomatoterm,\n\
                                         Acartratiosomatoterm,\n\
                                         Dsevenmicronbranchangles\n\
                                         Asevenmicronbranchangles\n\
    statsattr_autoplot=<true/false>      preset output format of statistics\n\
    statsattr_store_raw_data=<true/false> store raw data for each sample age\n\
    statsattr_fan_in_analysis=<true/false> perform fan-in angle analysis\n\
    outattr_show_stats=<true/false>      show statistics of resulting network\n\
    outattr_count_segments=<true/false/sample> \n\
    outattr_count_synapse_search=<true/false/sample>\n\
    outattr_synapse_genesis_and_loss=<true/false/sample>\n\
\n\
    Simulation Output - Runtime Options:\n\
\n\
    warnings_on=<off/stdout/stdoutfile/file> where to send warning messages\n\
    reports_on=<off/stdout/stdoutfile/file> where to send report messages\n\
    progress_on=<off/stdout/stdoutfile/file> where to send progress messages\n\
    outattr_show_progress=<true/false>   verbose progress indication\n\
    outattr_directory=<directory path>   path to output directory\n\
    outattr_URL=<output URL>             absolute output URL\n\
    figattr_tsupd_visibly=<true/false>   update terminal segments visibly\n\
\n\
    Simulation Output - Graphical Visualization:\n\
\n\
    outattr_show_figure=<true/false>     create figures of resulting network\n\
    figattr_use_color=<true/false>       presentation or publication figures\n\
    figattr_fill_somas=<true/false>      draw filled circles for somata\n\
    figattr_neurons=<true/false>         show neurons\n\
    figattr_connections=<true/false>     show abstract connections\n\
    figattr_connections_threshold=<decimal#> (list of) threshold values,\n\
                                         smallest is used in full figure,\n\
                                         a list creates -threshold- figures\n\
    figattr_presynaptic=<true/false>     show presynaptic structure\n\
    figattr_postsynaptic=<true/false>    show postsynaptic structure\n\
    figattr_synapses=<true/false>        show synapse structure\n\
    figattr_<receptor>=<true/false>      specify types: AMPAR, NMDAR, GABAR\n\
    figattr_partitions=<true/false>      show spatial partitions\n\
    figattr_connection_eval=<true/false> show partitions being evaluated\n\
    figattr_progress=<none/text/bar>     progress indicator\n\
    figattr_show_scale=<true/false>      draw scale bar\n\
    figattr_make_full_Fig=<true/false>   create .fig(.pdf) of full network\n\
    figattr_make_zoom_Fig=<true/false>   create .fig in specified focus box\n\
    figattr_make_connections_Fig=<true/false> plot connections of specific cell\n\
    figattr_make_abstract_Fig=<true/false> create .fig(.pdf) of abstract conn.\n\
    figattr_make_neurons_Figs=<true/false> create .fig for each neuron\n\
    net_zoom_disttoedge                  net zoom box radius (default=-1)\n\
    net_zoom_<width/height/depth>        net zoom volume extents\n\
    net_zoom_center<X/Y/Z>               net zoom box center coordinates\n\
    figattr_fibres_nobox=<true/false>    draw fibres regardless of focus box\n\
    figattr_box_fibre_independently=<true/false> show fibres of unshown n.\n\
    figattr_show_axis_arrows=<true/fasle> show axes in 3D\n\
    figureTeXwidth=<decimal#>            figure .fig width in inches\n\
    figattr_Fig_rescale=<true/false>     rescale .fig output for figureTeXwidth\n\
    color_table=<default/culturestain>   select a standard color table\n\
    <CT_id>=<[0x]integer#>               color table RGB entry with CT_id:\n\
                                         CT_background, CT_neuron_untyped,\n\
                                         CT_neuron_principal,\n\
                                         CT_neuron_interneuron,\n\
                                         CT_connection_excitatory,\n\
                                         CT_connection_inhibitory,\n\
                                         CT_synapses, CT_dendrites,\n\
                                         CT_axon_excitatory,\n\
                                         CT_axon_inhibitory,\n\
                                         CT_partition_border,\n\
                                         CT_partition_overfull,\n\
                                         CT_partition_evaluated,\n\
                                         CT_progress_text,\n\
                                         CT_AMPAR, CT_NMDAR, CT_GABAR\n\
    <schema.>color=<[0x]integer#>        A display color RGB specification\n\
\n\
    Simulation Output - Sequences and Animation:\n\
\n\
    figuresequence=<true/false>          produce a sequence of output figures\n\
    combinesequence=<true/false>         combine a sequence of figures\n\
    combinetype=<type extension>         default: gif\n\
    sequence_total_time=<decimal#>       seconds (default=25)\n\
    sequence_zoom_disttoedge             sequence zoom box radius (default=-1)\n\
    sequence_zoom_<width/height/depth>   sequence zoom volume extents\n\
    sequence_zoom_center<X/Y/Z>          sequence zoom box center coordinates\n\
    autorotatesequence=<true/false>      rotate the 3D view in a sequence\n\
    ROT_origin[_x/_y/_z]=<v/d#>          origin of rotation by samples\n\
    ROT[_x/_y/_z]=<vector/decimal#>      initial rotation angles around axes\n\
    ROT_interval[_x/_y/_z]=<v/d#>        total interval of rotation by samples\n\
    camera[_x/_y/_z]=<vector/decimal#>   projection camera location\n\
    camera_ROT[_x/_y/_z]=<v/d#>          projection camera rotation\n\
    viewer[_x/_y/_z]=<vector/decimal#>   projection viewer position\n\
    sample_dt=<decimal#>                 output sample interval (seconds)\n\
    max_res_samples=<true/false>         maximum resolution sample output\n\
    combinemagnification=<decimal#>      magnificant applied when combining\n\
    camera=<true/false>                  use camera specifications\n\
\n\
    Simulation Output - Output for NES:\n\
\n\
    NES_output=<true/false>              Produce output for NES\n\
\n\
    General:\n\
\n\
    substitute=<regular-expression>      RegEx substitution applied to commands\n\
    include=<file-name>                  include commands from file\n\
\n\
    type IDs of <neuron-type>:\n\
      principal, interneuron, multipolar, bipolar, pyramidal, untyped\n\
\n\
    natural sets of <schema.>:\n\
      (empty universal set), all_axons., all_dendrites., all_pyramidal_axons.,\n\
      all_pyramidal_dendrites., all_interneuron_axons.,\n\
      all_interneuron_dendrites., all_apical_pyramidal_dendrites.\n\
\n\
    PDF specific parameters:\n\
      delta: <label>.PDF.value\n\
      uniform: (none)\n\
      linear: <label>.PDF.max_x, <label>.PDF.height_b\n\
      spline_normal: <label>.PDF.max_x, <label>.PDF.significance_threshold,\n\
        <label>.PDF.proportion_significant\n\
      spline_normal_with_min: <label>.PDF.max_x,\n\
        <label>.PDF.significance_threshold, <label>.PDF.proportion_significant,\n\
        <label>.PDF.min_x\n\
      normal: <label>.PDF.mean, <label>.PDF.std, <label>.PDF.trunc\n\
      exponential: (none)\n\
      discrete: <label>.PDF.N, <label>.PDF.P, <label>.PDF.R_min,\n\
        <lable>.PDF.R_max\n\
\n\
  files:\n\
\n\
    .nibr(2D)rc contains commands included at initialization if the file exists.\n\
";
/* Former commands:
    sidelength=<decimal#>                length of side of network shape\n\
    Agrowth_n_inf
    Dgrowth_n_inf
    <axon/dendrite>_direction_model=     specific direction of growth model\n\
      <none/segment_history_tension/     to use at continuation nodes\n\
      cell_attraction>\n\
    <axon/dendrite>_direction_model_label identifier at root when chained\n\
    dendrite_two_branch_types=<true/false> use main and daughter branching\n\
    axon_two_branch_types=<true/false>     use main and daughter branching\n\

    <Agrowth_id>=<decimal#>              Axon growth parameters:\n\
                                         Agrowth_E, Agrowth_tau,\n\
                                         Agrowth_B_inf,\n\
                                         Agrowth_c,\n\
                                         Agrowth_c1, Agrowth_c2,\n\
    <Dgrowth_id>=<decimal#>              Dendrite growth parameters:\n\
                                         Dgrowth_E, Dgrowth_tau,\n\
                                         Dgrowth_B_inf,\n\
                                         Dgrowth_c,\n\
                                         Dgrowth_c1, Dgrowth_c2,\n\
    <Dapical_growth_id>=<decimal#>       Apical dendrite growth parameters:\n\
                                         Dgrowth_E_apicaltrunc,\n\
                                         Dgrowth_tau_apicaltrunc,\n\
                                         Dgrowth_B_inf_apicaltrunc,\n\
                                         Dgrowth_c_apicaltrunc,\n\
                                         Dgrowth_E_apicaltuft,\n\
                                         Dgrowth_tau_apicaltuft,\n\
                                         Dgrowth_B_inf_apicaltuft,\n\
                                         Dgrowth_c_apicaltuft\n\

    The following does not yet appear to exist in the model implementation:
      .bm.PDF=<PDF-id>                   delta, uniform, linear, spline_normal,\n\
                                         spline_normal_with_min, normal,\n\
                                         exponential\n\
*/

/* An example of tests to perform with INCLUDE_PDF_SAMPLING:
  rng X;
  normal_pdf X_norm(X);
  for (int i=0; i<1000000; i++) X_norm.random_selection();
  for (int i=0; i<1000000; i++) X_norm.random_positive();
  cout << X_norm.sample_histogram_data();
  normal_pdf X2_norm(X,0.0,0.3,1.0);
  X2_norm.set_sample_max(2.0);
  for (int i=0; i<1000000; i++) X2_norm.random_selection();
  cout << X2_norm.sample_histogram_data();
  exit(0);
 */

void copyright() {
  cout << "NETMORPH Copyright (C) 2008 Randal A. Koene <randalk@netmorph.org>\nThis program comes with ABSOLUTELY NO WARRANTY; for details see the\nincluded LICENSE file or <http://www.gnu.org/licenses/>.\nThis is free software, and you are welcome to redistribute it\nunder certain conditions; see the included LICENSE file or\n<http://www.gnu.org/licenses/> for details.\n\n";
}


/*!
  \mainpage NETMORPH Class Documentation and UML

  \section Entry function main()

  \image html /home/randalk/src/nnmodels/nibr/doc/main.png "[UML] Actions of the NETMORPH program entry function main()."
  \image latex /home/randalk/src/nnmodels/nibr/doc/main.png "[UML] Actions of the NETMORPH program entry function main()."

  \section The simulation loop
 */
int main(int argc, char * argv[]) {
  copyright();
#ifdef USE_RCFILE
  Command_Line_Parameters clp(argc,argv,helpstr,RCFILE);
#else
  Command_Line_Parameters clp(argc,argv,helpstr);
#endif
  main_clp = &clp;

  if (clp.IsFormInput()) outputdirectory = "../nibr/output/"; // the default when called by a form
  int n;
  if ((n=clp.Specifies_Parameter("outattr_directory"))>=0) outputdirectory = clp.URI_unescape_ParValue(n);
  if ((n=clp.Specifies_Parameter("outattr_URL"))>=0) absoluteoutputURL = clp.URI_unescape_ParValue(n);
  if (!CLP_recentinclude.empty()) {
    String scriptfilename;
    if (outputdirectory.empty()) scriptfilename = CLP_recentinclude;
    else {
      scriptfilename = CLP_recentinclude.after('/',-1);
      if (scriptfilename.empty()) scriptfilename = CLP_recentinclude;
    }
    String outputnameprepend(scriptfilename.before('.',-1));
    if (outputnameprepend.empty()) outputnameprepend = scriptfilename;
    outputdirectory += outputnameprepend + '_';
  }
  char dstr[80];
  time_t t = time(NULL);
  strftime(dstr,80,"%Y%m%d%H%M_",localtime(&t));
  outputdirectory += dstr;
  // At this point "outputdirectory" is fully available

  if ((n=clp.Specifies_Parameter("warnings_on"))>=0) {
    if (downcase(clp.ParValue(n))==String("off")) {
      warnings_on = WARN_OFF;
    } else if (downcase(clp.ParValue(n))==String("stdout")) {
      warnings_on = WARN_STDOUT;
    } else if (downcase(clp.ParValue(n))==String("stdoutfile")) {
      warnings_on = WARN_STDOUTFILE;
    } else if (downcase(clp.ParValue(n))==String("file")) {
      warnings_on = WARN_FILE;
    }
  }
  if (warnings_on>WARN_STDOUT) {
    warningfile = outputdirectory + "warnings";
  }
  if ((n=clp.Specifies_Parameter("reports_on"))>=0) {
    if (downcase(clp.ParValue(n))==String("off")) {
      reports_on = WARN_OFF;
    } else if (downcase(clp.ParValue(n))==String("stdout")) {
      reports_on = WARN_STDOUT;
    } else if (downcase(clp.ParValue(n))==String("stdoutfile")) {
      reports_on = WARN_STDOUTFILE;
    } else if (downcase(clp.ParValue(n))==String("file")) {
      reports_on = WARN_FILE;
    }
  }
  if (reports_on>WARN_STDOUT) {
    reportfile = outputdirectory + "report";
  }
  if ((n=clp.Specifies_Parameter("progress_on"))>=0) {
    if (downcase(clp.ParValue(n))==String("off")) {
      progress_on = WARN_OFF;
    } else if (downcase(clp.ParValue(n))==String("stdout")) {
      progress_on = WARN_STDOUT;
    } else if (downcase(clp.ParValue(n))==String("stdoutfile")) {
      progress_on = WARN_STDOUTFILE;
    } else if (downcase(clp.ParValue(n))==String("file")) {
      progress_on = WARN_FILE;
    }
  }
  if (progress_on>WARN_STDOUT) {
    progressfile = outputdirectory + "progress";
  }
  report_compiler_directives();
  reliability_checklist();

  write_file_from_String(outputdirectory+LASTCLPFILE,clp.str(";\n")); // log command line
  if (clp.IsFormInput()) {
    single_instance_restriction();
    progress("Content-Type: text/html\n\n<HTML>\n<BODY>\n<H1>Network Generation: Simulation</H1>\n\n<PRE>\n");
  }
  if (clp.IncludeFileErrors()>0) {
    if (clp.IsFormInput()) progress("\n<B>");
    warning("Warning: The simulation was unable to read "+String((long) clp.IncludeFileErrors())+" include file(s).\n");
    if (clp.IsFormInput()) progress("\n</B>");
  }

  //read_commands();
  developmental_simulation(clp);

  global_debug_report();
  if (clp.IsFormInput()) {
    progress("</PRE>\n\n</BODY>\n</HTML>\n");
    running_instance_postop();
  }

#ifdef TRACK_RECOGNIZED_COMMANDS
  warning(clp.unrecognized_str());
#endif

#ifdef TESTQUOTAS
  warning("BEWARE: TESTQUOTAS is still on and included!\n");
  write_file_from_String("quotas.txt",quotas);
#endif

  warnings_on = WARN_OFF; // avoids the unnecessary warning about attempting to destruct the Null_Activity object
  exit(0);
}
