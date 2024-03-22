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
// Color_Table.cc
// Randal A. Koene, 20050101

#include "Color_Table.hh"
#include "global.hh"
#include "Fig_Object.hh"

Color_Table defaultcolortable;
Color_Table * colortable = &defaultcolortable;
long definedcolors = 32; // the base colors pre-defined in XFig, colors numbered 32 and above are user-defined colors

const char ctname[CT_ARRAY_SIZE][25] = {
  "CT_background",
  "CT_neuron_untyped",
  "CT_neuron_principal",
  "CT_neuron_interneuron",
  "CT_connection_excitatory",
  "CT_connection_inhibitory",
  "CT_synapses",
  "CT_dendrites",
  "CT_axon_excitatory",
  "CT_axon_inhibitory",
  "CT_partition_border",
  "CT_partition_overfull",
  "CT_partition_evaluated",
  "CT_progress_text",
  "CT_AMPAR",
  "CT_NMDAR",
  "CT_GABAR"
};

// classes

void Color_Table::init() {
  for (int i = 0; i<CT_ARRAY_SIZE; i++) { ct[i].colnum = definedcolors; definedcolors++; }
  ct[CT_background].coldef = 0xffffff;
  ct[CT_neuron_untyped].coldef = 0x000000;
  ct[CT_neuron_principal].coldef = 0x0000ae;
  ct[CT_neuron_interneuron].coldef = 0xcf0000;
  ct[CT_connection_excitatory].coldef = 0x0000ff;
  ct[CT_connection_inhibitory].coldef = 0xff0000;
  ct[CT_synapses].coldef = 0xff00ff;
  ct[CT_dendrites].coldef = 0x0000ff;
  ct[CT_axon_excitatory].coldef = 0x009200;
  ct[CT_axon_inhibitory].coldef = 0xff0000;
  ct[CT_partition_border].coldef = 0x000000;
  ct[CT_partition_overfull].coldef = 0xffff00;
  ct[CT_partition_evaluated].coldef = 0x00ffff;
  ct[CT_progress_text].coldef = 0x000000;
  ct[CT_AMPAR].coldef = 0x007F00;
  ct[CT_NMDAR].coldef = 0x00007F;
  ct[CT_GABAR].coldef = 0x7F0000;
}

void Color_Table::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  if ((n=clp.Specifies_Parameter("color_table"))>=0) {
    if (downcase(clp.ParValue(n))==String("default")) init();
    else if (downcase(clp.ParValue(n))==String("culturestain")) set_culture_stain();
    else warning("Warning: Ignoring unknown color table specification '"+clp.ParValue(n)+"'.\n");
  }
  for (int i = 0; i<CT_ARRAY_SIZE; i++)
    if ((n=clp.Specifies_Parameter(ctname[i]))>=0) ct[i].coldef = strtoul(clp.ParValue(n),NULL,0);
}

String Color_Table::report_parameters() {
  String res("Color table: ");
  report(String(ctname[7])+String(ctname[8])+String(ctname[9])+'\n');
  for (int i = 0; i<CT_ARRAY_SIZE; i++) {
    res += ctname[i];
    res += "=0x";
    res += RGBstr(ct[i].coldef);
    if (i==7) report(res+'\n');
  }
  return res;
}

void Color_Table::set_culture_stain() {
  ct[CT_background].coldef = 0x000000;
  ct[CT_neuron_untyped].coldef = 0xffffff;
  ct[CT_neuron_principal].coldef = 0x007f00;
  ct[CT_neuron_interneuron].coldef = 0x7f0000;
  ct[CT_connection_excitatory].coldef = 0xffff00;
  ct[CT_connection_inhibitory].coldef = 0xff0000;
  ct[CT_synapses].coldef = 0xffffff;
  ct[CT_dendrites].coldef = 0x007f00;
  ct[CT_axon_excitatory].coldef = 0xffff00;
  ct[CT_axon_inhibitory].coldef = 0xff0000;
  ct[CT_partition_border].coldef = 0xffffff;
  ct[CT_partition_overfull].coldef = 0x00003f;
  ct[CT_partition_evaluated].coldef = 0x00005f;
  ct[CT_progress_text].coldef = 0xffffff;
  ct[CT_AMPAR].coldef = 0x00FF00;
  ct[CT_NMDAR].coldef = 0x0000FF;
  ct[CT_GABAR].coldef = 0xFF0000;
}

unsigned int Color_Table::get_color(long cn) {
  // get a color definition from the additional list, where
  // cn is the assigned numerical color ID
  if ((cn<(CT_ARRAY_SIZE+32)) || (cn>=definedcolors) || (cn>extracolorsmax)) {
    warning("Warning: Requested color outside of 'extra colors' range with color ID number "+String(cn)+", setting to black.\n");
    return 0;
  }
  return extracolor[cn-(CT_ARRAY_SIZE+32)];
}

void Color_Table::set_color(long cn, unsigned int cd) {
  // set a color or add a color definition according to
  // the assigned numerical ID in the additional colors list
  if (cn>255) error("Error: Attempted to define more than 256 colors.\n");
  if (cn<(CT_ARRAY_SIZE+32)) {
    warning("Warning: Ignoring attempt to set an 'extra color' at a color ID number below the 'extra colors' range.\n");
    return;
  }
  if (cn<=extracolorsmax) {
    extracolor[cn-(CT_ARRAY_SIZE+32)] = cd;
    return;
  }
  unsigned int * newextracolor = new unsigned int[(cn-(CT_ARRAY_SIZE+32))+1];
  if (extracolor) {
    for (long i = (extracolorsmax-(CT_ARRAY_SIZE+32)); i>=0; i--) newextracolor[i] = extracolor[i];
    delete[] extracolor;
  }
  extracolor = newextracolor;
  extracolorsmax = cn;
  definedcolors = extracolorsmax+1;
  extracolor[cn-(CT_ARRAY_SIZE+32)] = cd;
}

String Color_Table::Fig_str() {
  String figstr;
  for (int i=0; i<CT_ARRAY_SIZE; i++) figstr += Fig_Color(ct[i].colnum,ct[i].coldef).str();
  for (long i=(CT_ARRAY_SIZE+32); i<=extracolorsmax; i++) figstr += Fig_Color(i,extracolor[i-(CT_ARRAY_SIZE+32)]).str();
  return figstr;
}

