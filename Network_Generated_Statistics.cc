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
// Network_Generated_Statistics.cc
// Randal A. Koene, 20041230

#include "file.hh"
#include "Network_Generated_Statistics.hh"


const char stat_labels_str[STAT_LABELS_NUM][30] = {
  "Dlength",
  "Alength",
  "Dtermsegsperarbor",
  "Atermsegsperarbor",
  //  "Dintermediatesegsperarbor",
  //  "Aintermediatesegsperarbor",
  "Dterminallength",
  "Aterminallength",
  "Dintermediatelength",
  "Aintermediatelength",
  "Dlengthbetweenbifurcations",
  "Alengthbetweenbifurcations",
  "Dtermlensincebifurcation",
  "Atermlensincebifurcation",
  "Dtermlensincesoma",
  "Atermlensincesoma",
  "Dturnsbetweenbifurcations",
  "Aturnsbetweenbifurcations",
  "Dbranchangles",
  "Abranchangles",
  "Dturnangles",
  "Aturnangles",
  "Dcartratiobetweenbifurcations",
  "Acartratiobetweenbifurcations",
  "Dcartratiobifurcationtoterm",
  "Acartratiobifurcationtoterm",
  "Dcartratiosomatoterm",
  "Acartratiosomatoterm",
  "Dsevenmicronbranchangles",
  "Asevenmicronbranchangles"
  // "Dtermsomadistance",
  // "Atermsomadistance",
  // "Dtermbranchdistance",
  // "Atermbranchdistance",
  // "Dbifurcationpointseparation",
  // "Abifurcationpointseparation"
};

const bool stat_is_static[STAT_LABELS_NUM] = {
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  true,
  true,
  true,
  true,
  false,
  false,
  false,
  false,
  false,
  false,
  false,
  false
  //false,
  //false,
  //false,
  //false,
  //false,
  //false
};

unsigned long empty_data_samples = 0; // This counts empty data samples taken.

void Octave_Autoplot(String & s, String filebase) {
  s += "function res = standard_auto_plots()\n\
global terminal_segments_per_dendrite; global Eageidx; global Emeanidx; global terminal_segments_per_axon; global total_dendrite_length; global total_axon_length; global dendrites_mean_terminal_segment_length; global dendrites_mean_intermediate_segment_length; global axons_mean_terminal_segment_length; global axons_mean_intermediate_segment_length; global DNT; global ANT; global DL; global AL; global DTL; global ATL; global DIL; global AIL; global DLBB; global ALBB; global DTLB; global ATLB; global DTBB; global ATBB; global ABA; global DBA; global ATA; global DTA; global ASTR; global DSTR; global dendrites_seven_micron_branch_angles; global axons_seven_micron_branch_angles; global E7ageidx; global E7meanidx; global D7BA; global A7BA;\n\
gset term fig big\n\
gset output \""+filebase+"-dendritent.fig\"\ngplot terminal_segments_per_dendrite using Eageidx:Emeanidx, DNT using 1:2\ncloseplot\n\
gset output \""+filebase+"-axonnt.fig\"\ngplot terminal_segments_per_axon using Eageidx:Emeanidx, ANT using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendritelength.fig\"\ngplot total_dendrite_length using Eageidx:Emeanidx, DL using 1:2\ncloseplot\n\
gset output \""+filebase+"-axonlength.fig\"\ngplot total_axon_length using Eageidx:Emeanidx, AL using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendritetermlength.fig\"\ngplot dendrites_mean_terminal_segment_length using Eageidx:Emeanidx, DTLB using 1:2\ncloseplot\n\
gset output \""+filebase+"-axontermlength.fig\"\ngplot axons_mean_terminal_segment_length using Eageidx:Emeanidx, ATLB using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendriteinterlength.fig\"\ngplot dendrites_mean_intermediate_segment_length using Eageidx:Emeanidx, DLBB using 1:2\ncloseplot\n\
gset output \""+filebase+"-axoninterlength.fig\"\ngplot axons_mean_intermediate_segment_length using Eageidx:Emeanidx, ALBB using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendritebranchangles.fig\"\ngplot dendrites_seven_micron_branch_angles using E7ageidx:E7meanidx, D7BA using 1:2\ncloseplot\n\
gset output \""+filebase+"-axonbranchangles.fig\"\ngplot axons_seven_micron_branch_angles using E7ageidx:E7meanidx, A7BA using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendriteturnangles.fig\"\ngplot DTA using 1:2\ncloseplot\n\
gset output \""+filebase+"-axonturnangles.fig\"\ngplot ATA using 1:2\ncloseplot\n\
gset output \""+filebase+"-dendritesomatoterminalratio.fig\"\ngplot DSTR using 1:2\ncloseplot\n\
gset output \""+filebase+"-axonsomatoterminalratio.fig\"\ngplot ASTR using 1:2\ncloseplot\n\
endfunction\n\
\n";
  /*gset output \""+filebase+"-dendritebranchangles.fig\"\ngplot DBA using 1:2\ncloseplot\n\
    gset output \""+filebase+"-axonbranchangles.fig\"\ngplot ABA using 1:2\ncloseplot\n\ */
}

#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
void network_statistics_set::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("statsattr_"+label))>=0) collect = (downcase(clp.ParValue(n))==String("true"));
}

String network_statistics_set::report_parameters() {
  String res(label+'(');
  if (collect) res += "true) ";
  else res += "false) ";
  return res;
}

network_statistics_base::network_statistics_base(): collect_statistics(true) {
  for (int i = 0; i < STAT_LABELS_NUM; i++) netstats[i].label = stat_labels_str[i];
  // [***INCOMPLETE] THIS IS TEMPORARY! Since seven micron angles are not yet
  // measured, we turn off the collection of those statistics by default:
  netstats[dendrites_seven_micron_branch_angles].collect = false;
  netstats[axons_seven_micron_branch_angles].collect = false;
  // -------------------
  Allocate_New_Data();
}

void network_statistics_base::Allocate_New_Data() {
  // Note: Even if we are not collecting statistics, this will allocated
  // one sample to each type of statistics, so that there is no danger of
  // runtime errors caused by accessing NULL pointers.
  if ((!collect_statistics) && (netstats[0].head())) return;
  for (int i = 0; i<STAT_LABELS_NUM; i++) if (netstats[i].collect) {
    network_statistics_sample * nss = new network_statistics_sample();
    netstats[i].link_before(nss);
    if ((stat_is_static[i]) && (nss->Prev()) ) {
      // take raw data from data collected previously
      nss->Raw_Data().link_after(&(nss->Prev()->Raw_Data()));
    }
  }
}

void network_statistics_base::Store_Tail_Raw_Data(String datafile) {
  if ((!collect_statistics) || (!netstats[0].tail())) return;
  double t = netstats[0].tail()->age();
  String datastr("# Created by nibr(2D) network_statistics_base::Store_Tail_Raw_Data()\n\
# Raw sample data obtained with network::Collect_Data()\n\
# name: sample_t\n\
# type: matrix\n\
# rows: 1\n\
# columns: 1\n");  
  datastr += String(t,"%.3f\n");
  network_statistics_sample * nss;
  for (int i = 0; i<STAT_LABELS_NUM; i++) if (netstats[i].collect) {
    nss = netstats[i].tail();
    if (nss) {
      datastr += "# name: raw_"; datastr += stat_labels_str[i];
      datastr += "\n# type: matrix\n# rows: "+String(nss->Raw_Data().N());
      datastr += "\n# columns: 1\n";
      nss->Raw_Data().str(datastr);
    }
  }
  write_file_from_String(datafile+'-'+String(t,"%012.3f.ascii"),datastr);
}

void network_statistics_base::Cache_Tail_Statistics() {
  for (int i = 0; i<STAT_LABELS_NUM; i++) if (netstats[i].collect) netstats[i].tail()->fill_cache();
}

void network_statistics_base::Remove_Tail_Raw_Data(bool evenstatic) {
  for (int i = 0; i<STAT_LABELS_NUM; i++)
    if (evenstatic || (!stat_is_static[i])) if (netstats[i].tail()) netstats[i].tail()->Raw_Data().clear();
}

void network_statistics_base::Octave_Output(String filename, bool autoplot) {
  String filebase(filename.before('.',-1));
  if (filebase.empty()) filebase = filename;
  progress("Creating Octave file for statistics: "+filename+'\n');
  String s("#! /usr/bin/octave\n# Randal A. Koene\n# nibr statistics output\n\n# LOADPATH=[LOADPATH,':',system('printf $HOME',1),'/octave/common-m//'];\n\n# Statistical data of simulated network development\n\nT = [\n");
  PLL_LOOP_FORWARD(network_statistics_sample,netstats[0].head(),1) s += String(e->age()/(24.0*60.0*60.0),"%.3f\n");
  s += "];\nDarborssampled = [\n";
  PLL_LOOP_FORWARD(network_statistics_sample,netstats[dendrites_arbor_len].head(),1) s += String(e->N(),"%.7f\n");
  s += "];\nAarborssampled = [\n";
  PLL_LOOP_FORWARD(network_statistics_sample,netstats[axons_arbor_len].head(),1) s += String(e->N(),"%.7f\n");
  s += "];\n";
  //cout << "A\n" << s << '\n'; cout.flush();
  for (int i = 0; i<STAT_LABELS_NUM; i++) if (netstats[i].collect) {
    //cout << '/'; cout.flush();
    s += stat_labels_str[i]; s += " = [\n";
    PLL_LOOP_FORWARD(network_statistics_sample,netstats[i].head(),1) {
      //cout << e << '\n'; cout.flush();
      e->Octave_Output(s);
      //cout << "r\n"; cout.flush();
      s[s.length()-1] = '\n';
    }
    s += "];\n";
    //cout << '+'; cout.flush();
  }
  //cout << 'B'; cout.flush();
  s += "filebase = \""+filebase+"\";\n";
  if (autoplot) s += "autoplot = 1;\n";
  else s += "autoplot = 0;\n";
  if (autoplot) Octave_Autoplot(s,filebase);
  write_file_from_String(filename,s);
  //cout << 'C'; cout.flush();
}
#endif

#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
void Network_Generated_Statistics::Initialize_Terminals_Development_Statistics(int numdatapoints) { // , bool std) { // [*** NOTE] we could deallocate the
  // std buffers if that is desired when std values are not collected
  clear();
  arraysize = numdatapoints;
  arrayindex = 0;
  for (int i = 0; i < STAT_ARRAY_LABELS_NUM; i++) {
    statistics[i] = new double[numdatapoints];
    for (int j = 0; j < arraysize; j++) (statistics[i])[j] = -1.0;
  }
}

void Network_Generated_Statistics::add_dendrites_axons_nt(double t, network_statistics_data & dnumnsd, network_statistics_data & anumnsd) {
  if (recent_t()<t) add_t(t);
  add_dendrites_mean_nt(dnumnsd.mean);
  add_dendrites_std_nt(dnumnsd.std);
  add_axons_mean_nt(anumnsd.mean);
  add_axons_std_nt(anumnsd.std);
  add_dendrites_arbors(dnumnsd.n);
  add_axons_arbors(anumnsd.n);
  add_value(arrayofdendritesmin_nt,dnumnsd.min);
  add_value(arrayofdendritesmax_nt,dnumnsd.max);
  add_value(arrayofaxonsmin_nt,anumnsd.min);
  add_value(arrayofaxonsmax_nt,anumnsd.max);
}

void Network_Generated_Statistics::add_dendrites_axons_lengths(double t, network_statistics_data & dlennsd, network_statistics_data & alennsd, network_statistics_data & dtermnsd, network_statistics_data & atermnsd) {
   if (recent_t()<t) add_t(t);
   add_dendrites_mean_length(dlennsd.mean);
   add_axons_mean_length(alennsd.mean);
   add_dendrites_mean_terminal_length(dtermnsd.mean);
   add_axons_mean_terminal_length(atermnsd.mean);
   add_dendrites_std_length(dlennsd.std);
   add_axons_std_length(alennsd.std);
   add_dendrites_std_terminal_length(dtermnsd.std);
   add_axons_std_terminal_length(atermnsd.std);
   add_value(arrayofdendritelengthmin,dlennsd.min);
   add_value(arrayofdendritelengthmax,dlennsd.max);
   add_value(arrayofaxonlengthmin,alennsd.min);
   add_value(arrayofaxonlengthmax,alennsd.max);
   add_value(arrayofdendriteterminallengthmin,dtermnsd.min);
   add_value(arrayofdendriteterminallengthmax,dtermnsd.max);
   add_value(arrayofaxonterminallengthmin,atermnsd.min);
   add_value(arrayofaxonterminallengthmax,atermnsd.max);
}

void Network_Generated_Statistics::add_dendrites_axons_intermediate_lengths(double t, network_statistics_data & dinternsd, network_statistics_data & ainternsd) {
   if (recent_t()<t) add_t(t);
   add_value(arrayofdendriteintermediatelength,dinternsd.mean);
   add_value(arrayofdendriteintermediatelengthstd,dinternsd.std);
   add_value(arrayofdendriteintermediatelengthmin,dinternsd.min);
   add_value(arrayofdendriteintermediatelengthmax,dinternsd.max);
   add_value(arrayofaxonintermediatelength,ainternsd.mean);
   add_value(arrayofaxonintermediatelengthstd,ainternsd.std);
   add_value(arrayofaxonintermediatelengthmin,ainternsd.min);
   add_value(arrayofaxonintermediatelengthmax,ainternsd.max);
}

void Network_Generated_Statistics::Octave_Output(String filename, bool autoplot) {
  String filebase(filename.before('.',-1));
  if (filebase.empty()) filebase = filename;
  progress("Creating Octave file for statistics: "+filename+'\n');
  String s("#! /usr/bin/octave\n# Randal A. Koene\n# nibr statistics output\n\n# LOADPATH=[LOADPATH,':',system('printf $HOME',1),'/octave/common-m//'];\n\n# Statistical data of simulated network development\n\nT=[");
  for (int i=0; i<=arrayindex; i++) s += String(t(i)/(24.0*60.0*60.0),"%.3f\n");
  //cout << "6"; cout.flush();
  s[s.length()-1] = ']'; s += ";\nDarborssampled = [";
  for (int i=0; i<=arrayindex; i++) s += String(Dendrites_Arbors(i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nAarborssampled = [";
  for (int i=0; i<=arrayindex; i++) s += String(Axons_Arbors(i),"%.7f\n");
  // [*** NOTE] I could simplify the following by (a) having an array of
  // octave output names such as Dlength and (b) by sorting the
  // stat_array_labels so that I can simply use an incrementing index to
  // produce all the proper output.
  s[s.length()-1] = ']'; s += ";\nDlength = [";
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofdendritelength,i),"%.7f ") + String(get_value(arrayofdendritelengthstd,i),"%.7f ") + String(get_value(arrayofdendritelengthmin,i),"%.7f ") + String(get_value(arrayofdendritelengthmax,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nAlength = [";
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofaxonlength,i),"%.7f ") + String(get_value(arrayofaxonlengthstd,i),"%.7f ") + String(get_value(arrayofaxonlengthmin,i),"%.7f ") + String(get_value(arrayofaxonlengthmax,i),"%.7f\n");
  //cout << "7"; cout.flush();
  s[s.length()-1] = ']'; s += ";\nDterminallength = [";
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofdendriteterminallength,i),"%.7f ") + String(get_value(arrayofdendriteterminallengthstd,i),"%.7f ") + String(get_value(arrayofdendriteterminallengthmin,i),"%.7f ") + String(get_value(arrayofdendriteterminallengthmax,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nAterminallength = [";
  //cout << "A"; cout.flush();
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofaxonterminallength,i),"%.7f ") + String(get_value(arrayofaxonterminallengthstd,i),"%.7f ") + String(get_value(arrayofaxonterminallengthmin,i),"%.7f ") + String(get_value(arrayofaxonterminallengthmax,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nDtermsegsperarbor = [";
  //cout << "B"; cout.flush();
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofdendritesmean_nt,i),"%.7f ") + String(get_value(arrayofdendritesstd_nt,i),"%.7f ") + String(get_value(arrayofdendritesmin_nt,i),"%.7f ") + String(get_value(arrayofdendritesmax_nt,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nAtermsegsperarbor = [";
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofaxonsmean_nt,i),"%.7f ") + String(get_value(arrayofaxonsstd_nt,i),"%.7f ") + String(get_value(arrayofaxonsmin_nt,i),"%.7f ") + String(get_value(arrayofaxonsmax_nt,i),"%.7f\n");
  //cout << "C"; cout.flush();
  s[s.length()-1] = ']'; s += ";\nDintermediatelength = [";
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofdendriteintermediatelength,i),"%.7f ") + String(get_value(arrayofdendriteintermediatelengthstd,i),"%.7f ") + String(get_value(arrayofdendriteintermediatelengthmin,i),"%.7f ") + String(get_value(arrayofdendriteintermediatelengthmax,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\nAintermediatelength = [";
  //cout << "D"; cout.flush();
  for (int i=0; i<=arrayindex; i++) s += String(get_value(arrayofaxonintermediatelength,i),"%.7f ") + String(get_value(arrayofaxonintermediatelengthstd,i),"%.7f ") + String(get_value(arrayofaxonintermediatelengthmin,i),"%.7f ") + String(get_value(arrayofaxonintermediatelengthmax,i),"%.7f\n");
  s[s.length()-1] = ']'; s += ";\n\
\n";
  s += "filebase = \""+filebase+"\";\n";
  //cout << "8"; cout.flush();
  if (autoplot) s += "autoplot = 1;\n";
  else s += "autoplot = 0;\n";
  if (autoplot) Octave_Autoplot(s,filebase);
  //cout << "9"; cout.flush();
  write_file_from_String(filename,s);
}

#endif
