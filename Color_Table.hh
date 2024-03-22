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
// Color_Table.hh
// Randal A. Koene, 20050101
//
// Classes that provide and manipulate color tables.

#ifndef __COLOR_TABLE_HH
#define __COLOR_TABLE_HH

#include "BigString.hh"
#include "Command_Line_Parameters.hh"

enum color_table_name { CT_background, CT_neuron_untyped, CT_neuron_principal, CT_neuron_interneuron, CT_connection_excitatory, CT_connection_inhibitory, CT_synapses, CT_dendrites, CT_axon_excitatory, CT_axon_inhibitory, CT_partition_border, CT_partition_overfull, CT_partition_evaluated, CT_progress_text, CT_AMPAR, CT_NMDAR, CT_GABAR, CT_ARRAY_SIZE};

extern const char ctname[][25];

// classes

class colorpair {
  // This is used instead of Fig_Color so that the Color_Table.hh header
  // file does not depend on the Fig_Object.hh header file and so that the
  // color table can be used with multiple output formats. This is also why
  // the Color_Table::Fig_str() function is not simply called
  // Color_Table::str().
public:
  colorpair(): colnum(0), coldef(0) {}
  colorpair(long cn, unsigned int cd): colnum(cn), coldef(cd) {}
  long colnum;
  unsigned int coldef;
};

class Color_Table: public CLP_Modifiable {
protected:
  colorpair ct[CT_ARRAY_SIZE];
  unsigned int * extracolor;
  long extracolorsmax; // Note that once non-zero, this does not count from 1
  void init();
public:
  Color_Table(): extracolor(NULL), extracolorsmax(0) { init(); }
  ~Color_Table() { if (extracolor) delete[] extracolor; }
  long colnum(color_table_name n) { return ct[n].colnum; }
  unsigned int coldef(color_table_name n) { return ct[n].coldef; } // from list of named colors
  unsigned int get_color(long cn); // from list of additional colors
  void set_color(color_table_name n, long cn, unsigned int cd) { ct[n].colnum = cn; ct[n].coldef = cd; } // list of named colors
  void set_color(long cn, unsigned int cd); // list of additional colors
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  void set_defaults() { init(); }
  void set_culture_stain();
  String Fig_str();
};

extern Color_Table defaultcolortable;
extern Color_Table * colortable;
extern long definedcolors;

#endif
