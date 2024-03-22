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
// Txt_Object.cc
// Randal A. Koene, 20070201

#include "diagnostic.hh"
#include "global.hh"
#include "Txt_Object.hh"

String * Txt_neuronlist = NULL;
String * Txt_synapselist = NULL;
String * Txt_fiberrootlist = NULL;
String * Txt_continuationnodelist = NULL;
String * Txt_bifurcationnodelist = NULL;
String * Txt_terminalgrowthconelist = NULL;
String * Txt_tuftrootbranchnodelist = NULL;
String * Txt_obliquerootbranchnodelist = NULL;
long Txt_neuronindex = 0;
long Txt_synapseindex = 0;
long Txt_nodeindex = -1;

Txt_Header::Txt_Header(): text("#NETMORPH Txt; Generated Structure\n#Format:\n#1. List of neurons\n#2. List of synapses\n") {
  if (outattr_track_nodegenesis) text += "#3. List of fiber structure (axon/dendrite) roots\n#4. List of continuation nodes (turn points)\n#5. List of bifurcation nodes (branch points)\n#6. List of terminal growth cones (end points)\n#7. List of apical dendrite tuft root nodes (branch points)\n#8. List of apical dendrite oblique root nodes (branch points)\n";
  text += "#\n#List of neurons contains:\n";
  text += COLUMN_LABELS_NEURONS;
  text += "#List of synapses contains:\n";
  text += COLUMN_LABELS_SYNAPSES;
  text += "#(If multiple synapses have the same location then there are multiple different types of receptor populations at that synaptic location.)\n";
  if (outattr_track_nodegenesis) {
    text += "#List of root nodes contains:\n";
    text += COLUMN_LABELS_ROOT_NODES;
    text += "#List of continuation nodes contains:\n";
    text += COLUMN_LABELS_FIBER_NODES;
    text += "#List of bifurcation nodes contains:\n";
    text += COLUMN_LABELS_FIBER_NODES;
    text += "#List of terminal growth cones contains:\n";
    text += COLUMN_LABELS_FIBER_NODES;
    text += "#List of apical dendrite tuft root nodes contains:\n";
    text += COLUMN_LABELS_TUFT_NODES;
    text += "#List of apical dendrite oblique root nodes contains:\n";
    text += COLUMN_LABELS_OBLIQUE_NODES;
    text += "#(Fiber structure types are: axon, dendrite, apical dendrite. The parent label refers to the first centripetal point, i.e. towards the soma.)\n";
  }
  text += '\n';
}

String Txt_Group::str() {
  String s;
  if (!comment.empty()) {
    if (comment[comment.length()-1]!='\n') comment += '\n';
    comment[comment.length()-1] = '.';
    comment.gsub("\n","\n# ");
    comment[comment.length()-1] = '\n';
    s = "# " + comment;
  }
  // *** Do something here to print the desired lists
  PLL_LOOP_FORWARD(Txt_Object,txtobjects.head(),1) s+=e->str();
  return s;
}

void Txt_Group::Add_Txt_Object(Txt_Object * txtobject) {
  if (!txtobject) return;
  txtobjects.link_before(txtobject);
}
