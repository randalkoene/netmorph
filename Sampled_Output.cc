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
// Sampled_Output.cc
// Randal A. Koene, 20041118

#include "global.hh"
#include "Color_Table.hh"
#include "file.hh"
#include "network.hh"
#include "Sampled_Output.hh"

//#define CONVERT_READS_FIG_FILES

// global variables

fanin_histogram net_fanin_histogram;

Sampled_Growth_Output * sampled_output = NULL;
Activity_Results_Output * actres = NULL;

bool camera = false;

// classes

String fanin_histogram::octave_output() {
  String res("fanin_histogram_means = [");
  for (int i=0; i<50; i++) res += String(len[i],"%.3f\n");
  res += "];\n\nn_fanin = [";
  for (int i=0; i<50; i++) res += String((long) n[i]) + '\n';
  res += "];\n";
  return res;
}

Fig_Group * fanin_histogram::fig_output() {
  Fig_Group * faninfig = new Fig_Group();
  double angle = 0.0;
  long x[3];
  long y[3];
  x[0] = 0;
  y[0] = 0;
  double cosa = 2.0;
  double sina = 0.0;
  for (int i=0; i<50; i++) {
    x[1] = (long) (len[i]*cosa);
    y[1] = (long) (len[i]*sina);
    angle += FANIN_ANGLE_SLICE;
    cosa = -2.0*cos(angle);
    sina = 2.0*sin(angle);
    x[2] = (long) (len[i]*cosa);
    y[2] = (long) (len[i]*sina);
    faninfig->Add_Fig_Object(new Fig_Polygon(0,1,0,7,50,20,2.0,0,0,y,x,3));
  }
  for (int i=0; i<10; i++) faninfig->Add_Fig_Object(new Fig_Circle(2,1,0,0,49,-1,0.0,0.0,0,0,400*i,0,0,400,200));
  return faninfig;
}

Sampled_Output::Sampled_Output(network * netptr, double mt, double si, bool std): net(netptr), max_time(mt), sample_interval(si), numsamples(((int) (mt/si))+1), t_nextmark(si), maximum_resolution(false) {
  if (!netptr) error("nibr Error: No network specified in Sampled_Output::Sampled_Output()\n");
}

void Sampled_Output::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("max_res_samples"))>=0) maximum_resolution = (downcase(clp.ParValue(n))==String("true"));
}

String Sampled_Output::report_parameters() {
  String res("Sample Output: max_res_samples=");
  if (maximum_resolution) res += "true\n";
  else res += "false\n";
  return res;
}

int parse_CLP_samplesetting(Command_Line_Parameters & clp, String label, int defaultvalue) {
  // detects "false", "true" and "sample" settings and translates those to 0, 1 and 2 values
  int n;
  if ((n=clp.Specifies_Parameter(label))>=0) {
    String vstr(downcase(clp.ParValue(n)));
    if (vstr.matches("sample")) return 2;
    if (vstr.matches("true")) return 1;
    return 0;
  }
  return defaultvalue;
}

bool collect_statistics_for_fig_sequence = false;

void Sampled_Growth_Output::parse_CLP(Command_Line_Parameters & clp) {
  Sampled_Output::parse_CLP(clp);
  int n;
  if ((n=clp.Specifies_Parameter("statsattr_autoplot"))>=0) autoplot = (downcase(clp.ParValue(n))==String("true"));
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  if ((n=clp.Specifies_Parameter("statsattr_collect_statistics"))>=0) collect_statistics = (downcase(clp.ParValue(n))==String("true"));
  if (collect_statistics_for_fig_sequence) collect_statistics = true;
  if (collect_statistics) {
    for (int i = 0; i<STAT_LABELS_NUM; i++) netstats[i].parse_CLP(clp);
  } else {
    for (int i = 0; i<STAT_LABELS_NUM; i++) netstats[i].set_collect(false);
  }
  if ((n=clp.Specifies_Parameter("statsattr_store_raw_data"))>=0) store_raw_data = (downcase(clp.ParValue(n))==String("true"));
#endif
  segmentcount = parse_CLP_samplesetting(clp,"outattr_count_segments",segmentcount);
  synapsesearchprofile = parse_CLP_samplesetting(clp,"outattr_count_synpase_search",synapsesearchprofile);
  synapsegenesisandloss = parse_CLP_samplesetting(clp,"outattr_synapse_genesis_and_loss",synapsegenesisandloss);
}

String Sampled_Growth_Output::report_parameters() {
  const char boolstr[3][8] = {"false\n","true\n","sample\n"};
  String res(Sampled_Output::report_parameters());
  res += "Sample Output: statsattr_autoplot="; res += boolstr[autoplot];
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  res += "Sample Output: statsattr_collect_statistics="; res += boolstr[collect_statistics];
  if (collect_statistics) {
    res += "  Selected statistics: ";
    for (int i = 0; i<STAT_LABELS_NUM; i++) res += netstats[i].report_parameters();
    res += "\n";
  }
  res += "Sample Output: statsattr_store_raw_data="; res += boolstr[store_raw_data];
  if (collect_statistics) res += "(Note: Collecting statistics causes periodic segment vector updates.)\n";
#endif
  res += "Sample Output: outattr_count_segments="; res += boolstr[segmentcount];
  res += "Sample Output: outattr_count_synapse_search="; res += boolstr[synapsesearchprofile];
  res += "Sample Output: outattr_synapse_genesis_and_loss="; res += boolstr[synapsegenesisandloss];
  return res;
}

bool Sampled_Growth_Output::Sample(double t) {
  if (t<t_nextmark) return false;
  t_nextmark += sample_interval;
  if (t_nextmark>max_time) t_nextmark = max_time;
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  net->Collect_Data(*this);
  if (store_raw_data) Store_Tail_Raw_Data(statsname+"-raw-data");
  //cout << 'S' << t << '\n';
  Cache_Tail_Statistics();
  Remove_Tail_Raw_Data((t>=max_time));
  if (t<max_time) Allocate_New_Data();
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
  network_statistics_data dnumnsd, anumnsd, dlennsd, alennsd, dtermnsd, atermnsd, dinternsd, ainternsd;
  net->mean_number_of_input_terminal_segments(&dnumnsd);
  net->mean_number_of_output_terminal_segments(&anumnsd);
  net->mean_presynaptic_structure_length(atermnsd,alennsd,&ainternsd);
  net->mean_postsynaptic_structure_length(dtermnsd,dlennsd,&dinternsd);
  add_dendrites_axons_nt(t,dnumnsd,anumnsd);
  add_dendrites_axons_lengths(t,dlennsd,alennsd,dtermnsd,atermnsd);
  add_dendrites_axons_intermediate_lengths(t,dinternsd,ainternsd);
#endif
  //if ((segmentcount==2) || ((segmentcount==1) && (t>=max_time))) Count_Segments();
  //if ((synapsesearchprofile==2) || ((synapsesearchprofile==1) && (t>=max_time))) Count_Synapse_Search();
  //if ((synapsegenesisandloss==2) || ((synapsegenesisandloss==1) && (t>=max_time))) Synapse_Genesis_and_Loss();
  if (outattr_Txt_sequence) net->Txt_Output(Txtname+'-'+String(t,"%012.3f.txt"));
  return true;
}

void Sampled_Growth_Output::Output(bool show_stats) {
  progress("Sampled growth output.\n");
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  if (show_stats) network_statistics_base::Octave_Output(statsname,autoplot);
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
  if (show_stats) Network_Generated_Statistics::Octave_Output("ngs"+statsname,autoplot);
#endif
}

void Sampled_Growth_Fig_Output::parse_CLP(Command_Line_Parameters & clp) {
  collect_statistics_for_fig_sequence=true;
  Sampled_Growth_Output::parse_CLP(clp);
  int n;
  if ((n=clp.Specifies_Parameter("combinesequence"))>=0) combine = (downcase(clp.ParValue(n))==String("true"));
  if ((n=clp.Specifies_Parameter("combinemagnification"))>=0) combinemagnification = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter("combinetype"))>=0) combinetype = clp.ParValue(n);
  if ((n=clp.Specifies_Parameter("sequence_total_time"))>=0) {
    double seqtottime = atof(clp.ParValue(n));
    sequencedelay = lround((seqtottime*100.0)/((double) numsamples));
  }
  parse_focus_box_CLP("sequence",clp,zoomcenter,zoomdisttoedge);
  parse_focus_volume_CLP("sequence",clp,zoomwidth,zoomheight,zoomdepth);
  if (zoomwidth>0.0) zoomdisttoedge = zoomwidth; // This is necessary in order to pass a test in Sampled_Growth_Fig_Output::Sample()
#ifdef VECTOR3D
  if ((n=clp.Specifies_Parameter("autorotatesequence"))>=0) autorotatesequence = (downcase(clp.ParValue(n))==String("true"));
  if (autorotatesequence) spat_pres->set_rotate_speed(numsamples);
  if ((n=clp.Specifies_Parameter("camera"))>=0) camera = (downcase(clp.ParValue(n))==String("true"));
#endif
}

String Sampled_Growth_Fig_Output::report_parameters() {
  String res(Sampled_Growth_Output::report_parameters());
  res += "(Collecting statistics is auto-enabled by the option to create a sequence of sample figures.)\n";
  res += "Sample Output: combinesequence=";
  if (combine) res += "true (format: "+combinetype+")\n";
  else res += "false\n";
  if (combine) res += "Sample Output: combinemagnification="+String(combinemagnification,"%.2f\n");
  if (combine) res += "Sample Output: sequence_total_time="+String((((double) sequencedelay)*((double) numsamples))/100.0,"%.2f")+" seconds (sequencedelay="+String(sequencedelay)+")\n";
#ifdef VECTOR3D
  res += "Sample Output: autorotatesequence=";
  if (autorotatesequence) res += "true\n";
  else res += "false\n";
#endif
  if (zoomwidth>0.0) {
    res += "Sample Output: zoom width="+String(zoomwidth,"%.2f, ")+"zoom height="+String(zoomheight,"%.2f, ")+"zoom depth="+String(zoomdepth,"%.2f\n");
  } else res += "Sample Output: zoom distance to edge="+String(zoomdisttoedge,"%.2f\n");
  
  res += "Sample Output: zoom center=("+String(zoomcenter.X(),"%.2f")+','+String(zoomcenter.Y(),"%.2f")
#ifdef VECTOR3D
    +','+ String(zoomcenter.Z(),"%.2f")+")\n";
#endif
#ifdef VECTOR2D
  +")\n";
#endif
  return res;
}

bool Sampled_Growth_Fig_Output::Sample(double t) {
  if (!Sampled_Growth_Output::Sample(t)) return false;
  // *** to do this faster, I can capture all the Fig_Groups in one and
  // use its box size for scaling, then produce output from the component
  // groups, to save memory do as follows
  double usewidth = figwidth;
  if (zoomdisttoedge>0.0) {
    if (zoomwidth>0.0) set_focus_volume(zoomcenter,zoomwidth,zoomheight,zoomdepth);
    else set_focus_box(zoomcenter,zoomdisttoedge);
    usewidth = figattr_focus_box_x2-figattr_focus_box_x1;
  }
  sampling = true;
  //cout << 'F' << t << '\n';
  Fig_Group * figgroup = net->Fig_Output(figname+'-'+String(t,"%012.3f.fig"),usewidth,true);
  sampling = false;
  if (figgroup->TopLeftX()<topleftx) topleftx = figgroup->TopLeftX();
  if (figgroup->TopLeftY()<toplefty) toplefty = figgroup->TopLeftY();
  if (figgroup->BottomRightX()>bottomrightx) bottomrightx = figgroup->BottomRightX();
  if (figgroup->BottomRightY()>bottomrighty) bottomrighty = figgroup->BottomRightY();
  delete figgroup;
  netfigidx++;
#ifdef VECTOR3D
  if (autorotatesequence) spat_pres->update_rotation();
#endif
  return true;
}

bool Sampled_Growth_Fig_Output::Max_Resolution_Sample(double t, const void * id) {
  // *** if (!Sampled_Growth_Output::Max_Resolution_Sample(t,id)) return false;
  //     I may reinstate the above if we need statistics at maximum
  //     resolution in addition to network plots. But if so, then
  //     arrayoftime must be dealt with properly to avoid overruns due to
  //     many calls of Max_Resolution_Sample().
  // *** to do this faster, I can capture all the Fig_Groups in one and
  // use its box size for scaling, then produce output from the component
  // groups, to save memory do as follows
  Fig_Group * figgroup = net->Fig_Output(figname+'-'+String(t,"%012.3f-")+String((long) id)+".fig",figwidth,true);
  /*if (figgroup->TopLeftX()<topleftx) topleftx = figgroup->TopLeftX();
  if (figgroup->TopLeftY()<toplefty) toplefty = figgroup->TopLeftY();
  if (figgroup->BottomRightX()>bottomrightx) bottomrightx = figgroup->BottomRightX();
  if (figgroup->BottomRightY()>bottomrighty) bottomrighty = figgroup->BottomRightY();*/
  delete figgroup;
  // *** (see first comment line in this function) netfigidx++;
  return true;
}

void Sampled_Growth_Fig_Output::Output(bool show_stats) {
  Sampled_Growth_Output::Output(show_stats);
  String fnamebase;
  //  if (!((figattr_show_spatial_segment_subsets) && (net->spatial_segment_subset()))) {
  if (netfigidx>0) progress("Aligning sequence of sample figures.\n");
  //cout << topleftx << ',' << toplefty << ',' << bottomrightx << ',' << bottomrighty << '\n';
  Fig_Rectangle figrectangle(0,1,colortable->colnum(CT_background),colortable->colnum(CT_background),999,-1,0.0,0,0,topleftx-2,toplefty-2,bottomrightx+2,bottomrighty+2);
  String figrectanglestr(figrectangle.str());
  if (figattr_show_scale) {
    Fig_Group * scalebarfig = network::scale_bar_Fig(topleftx,bottomrighty);
    figrectanglestr += scalebarfig->str();
  }
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  if (netfigidx>0) PLL_LOOP_FORWARD(network_statistics_sample,netstats[0].head(),1) append_file_from_String(figname+'-'+String(e->age(),"%012.3f.fig"),figrectanglestr);
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
  for (int i=0; i<netfigidx; i++) append_file_from_String(figname+'-'+String(t(i),"%012.3f.fig"),figrectanglestr);
#endif
  //  }
  if (combine) {
    progress("Converting sample figures for combined representation:\n");
    String magnificationstr(" ");
    if ((combine) && (combinemagnification!=1.0)) magnificationstr = " -m "+String(combinemagnification,"%.2f ");
#ifndef CONVERT_READS_FIG_FILES
#ifndef COMPILE_WITHOUT_FIG2DEV
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
    if (netfigidx>0) PLL_LOOP_FORWARD(network_statistics_sample,netstats[0].head(),1) {
      fnamebase = figname+'-'+String(e->age(),"%012.3f");
#else
    for (int i=0; i<netfigidx; i++) {
      fnamebase = figname+'-'+String(t(i),"%012.3f");
#endif
      system(String(FIG2DEV)+String(" -L png -g#"+RGBstr(colortable->coldef(CT_background))+magnificationstr+fnamebase+".fig > "+fnamebase+".png").chars());
      progress('.');
    }
#else
    progress("NETMORPH WAS COMPILED WITHOUT FIG2DEV - NO SEQUENCE CONVERSION FROM .fig TO .png.\n");
#endif
#endif
    // duplicate the last figure in order to linger a bit longer on the end result before looping back to the beginning
    // [***NOTE] It now appears that "-pause" can be used for the delay between
    // loops.
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
    if (netfigidx>0) fnamebase = figname+'-'+String(netstats[0].tail()->age(),"%012.3f");
#else
    fnamebase = figname+'-'+String(t(netfigidx-1),"%012.3f");
#endif
#ifdef CONVERT_READS_FIG_FILES
    system(String("cp -f "+fnamebase+".fig "+fnamebase+"last.fig"));
#else
    system(String("cp -f "+fnamebase+".png "+fnamebase+"last.png"));
#endif
    progress("\nCombining samples into "+figname+'.'+combinetype+'\n');
#ifdef CONVERT_READS_FIG_FILES
    system(String("convert -delay "+String(sequencedelay)+' '+figname+"-*fig "+figname+'.'+combinetype).chars());
#else
    system(String("convert -delay "+String(sequencedelay)+' '+figname+"-*png "+figname+'.'+combinetype).chars());
#endif
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
    if (netfigidx>0) PLL_LOOP_FORWARD(network_statistics_sample,netstats[0].head(),1) {
      fnamebase = figname+'-'+String(e->age(),"%012.3f");
#else
    for (int i=0; i<netfigidx; i++) {
      fnamebase = figname+'-'+String(t(i),"%012.3f");
#endif
      remove(String(fnamebase+".fig").chars());
#ifndef CONVERT_READS_FIG_FILES
      remove(String(fnamebase+".png").chars());
#endif
    }
#ifdef CONVERT_READS_FIG_FILES
    remove(String(fnamebase+"last.fig").chars());
#else
    remove(String(fnamebase+"last.png").chars());
#endif
  }
}

void Activity_Results_Output::spike(neuron * n, double t) {
  progress(String(t,"%.3f")+':'+String((long) n)+'\n');
}

void Activity_Results_Output::psp(neuron * pre, neuron * post, double t) {
  progress(String(t,"%.3f")+':'+String((long) pre)+"->"+String((long) post)+'\n');
}

void ARO_to_File::spike(neuron * n, double t) {
  spikedata += String(t,"%f ");
  spikedata += String(n->numerical_ID());
  spikedata += '\n';
  numspikes++;
}

void ARO_to_File::psp(neuron * pre, neuron * post, double t) {
  pspdata += String(t,"%f ");
  pspdata += String(pre->numerical_ID());
  pspdata += ' ';
  pspdata += String(post->numerical_ID());
  pspdata += '\n';
  numpsps++;
}

void ARO_to_File::postop() {
  String asciidata("# Created by NETMORPH: network IFactivity\n#name: spikes\n# type: matrix\n# rows: ");
  asciidata += String(numspikes);
  asciidata += "\n# columns: 2\n";
  asciidata += spikedata;
  asciidata += "\n# name: psps\n# type: matrix\n# rows: ";
  asciidata += String(numpsps);
  asciidata += "\n# columns: 3\n";
  asciidata += pspdata;
  write_file_from_String(fname,asciidata);
  postop_done = true;
}

void set_focus_box(spatial & center, double sizer) {
  if (sizer<=0.0) {
    figattr_zoom = false;
    return;
  }
  figattr_zoom = true;
  figattr_focus_box_x1 = center.X()-sizer;
  figattr_focus_box_x2 = center.X()+sizer;
  figattr_focus_box_y1 = center.Y()-sizer;
  figattr_focus_box_y2 = center.Y()+sizer;
#ifdef VECTOR3D
  figattr_focus_box_z1 = center.Z()-sizer;
  figattr_focus_box_z2 = center.Z()+sizer;
#endif
}

void set_focus_volume(spatial & center, double width, double height, double depth) {
  if ((width<=0.0) || (height<=0.0)
#ifdef VECTOR3D
      || (depth<=0.0)
#endif
      ) {
    figattr_zoom = false;
    return;
  }
  figattr_zoom = true;
  figattr_focus_box_x1 = center.X()-(0.5*width);
  figattr_focus_box_x2 = center.X()+(0.5*width);
  figattr_focus_box_y1 = center.Y()-(0.5*height);
  figattr_focus_box_y2 = center.Y()+(0.5*height);
#ifdef VECTOR3D
  figattr_focus_box_z1 = center.Z()-(0.5*depth);
  figattr_focus_box_z2 = center.Z()+(0.5*depth);
#endif
}

void parse_focus_box_CLP(String zoomlabel, Command_Line_Parameters & clp, spatial & center, double & sizer) {
  // This is used both by Sampled_Growth_Fig_Output::parse_CLP() and
  // in nibr.cc.
  int n;
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_disttoedge"))>=0) sizer = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_center"))>=0) {
    neuron * nptr = static_cast<PLLRoot<neuron> *>(eq->Net())->el(atoi(clp.ParValue(n)));
    if (nptr) center = nptr->Pos();
    else warning("Parameter ERROR, "+zoomlabel+"_zoom_center: Neuron index "+clp.ParValue(n)+" does not exist in network.\n");
  }
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_centerX"))>=0) center.set_X(atof(clp.ParValue(n)));
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_centerY"))>=0) center.set_Y(atof(clp.ParValue(n)));
#ifdef VECTOR3D
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_centerZ"))>=0) center.set_Z(atof(clp.ParValue(n)));
#endif
  double x,y,z;
  center.get_all(x,y,z);
  report("Zoom center of "+zoomlabel+" = "+String(x,"%.2f,")+String(y,"%.2f,")+String(z,"%.2f\n"));
}

void parse_focus_volume_CLP(String zoomlabel, Command_Line_Parameters & clp, double & width, double & height, double & depth) {
  // This should be used both by Sampled_Growth_Fig_Output::parse_CLP() and
  // in nibr.cc.
  int n;
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_width"))>=0) width = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_height"))>=0) height = atof(clp.ParValue(n));
  if ((n=clp.Specifies_Parameter(zoomlabel+"_zoom_depth"))>=0) depth = atof(clp.ParValue(n));
}

#ifdef VECTOR2D
bool fig_in_zoom(double x, double y) {
#endif
#ifdef VECTOR3D
bool fig_in_zoom(double x, double y, double z) {
#endif
  if (!figattr_zoom) return true;
  if (x<figattr_focus_box_x1) return false;
  if (x>figattr_focus_box_x2) return false;
  if (y<figattr_focus_box_y1) return false;
  if (y>figattr_focus_box_y2) return false;
#ifdef VECTOR3D
  if (z<figattr_focus_box_z1) return false;
  if (z>figattr_focus_box_z2) return false;
#endif
  return true;
}

bool fig_in_zoom(const spatial & p) {
  if (!figattr_zoom) return true;
  if (p.X()<figattr_focus_box_x1) return false;
  if (p.X()>figattr_focus_box_x2) return false;
  if (p.Y()<figattr_focus_box_y1) return false;
  if (p.Y()>figattr_focus_box_y2) return false;
#ifdef VECTOR3D
  if (p.Z()<figattr_focus_box_z1) return false;
  if (p.Z()>figattr_focus_box_z2) return false;
#endif
  return true;  
}
