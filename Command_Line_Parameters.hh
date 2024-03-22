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
// Command_Line_Parameters.hh
// Randal A. Koene, 20041215
//
// These classes parse the command line or HTML form input and create
// an object with command line parameter data.

#ifndef __COMMAND_LINE_PARAMETERS_HH
#define __COMMAND_LINE_PARAMETERS_HH

#include <limits>
#include <cfloat>
#include <climits>
//STL#include <vector>
#include "StringList.hh"

#define TRACK_RECOGNIZED_COMMANDS

class Command_Line_Parameters {
protected:
  bool calledbyforminput;
  int clpargc;
  char ** clpargv;
  void Add_to_Substitutions(String subs);
  void Apply_Substitutions(String & pname);
#ifdef TRACK_RECOGNIZED_COMMANDS
  void Add_Parameter(String pname, String pvalue, String co);
#else
  void Add_Parameter(String pname, String pvalue);
#endif
  void Parse_File(const char filename[]);
  bool Detect_Form_Input();
  bool Parse_Command_Line();
  StringList parname, parvalue;
  int numparameters;
  int includefileerrors;
#ifdef TRACK_RECOGNIZED_COMMANDS
  //STLstd::vector<int> * recognized;
  int * recognized;
  StringList origin;
  int recognizedsize;
#endif
  StringList substitutionpat, substitutionstr;
  int numsubstitutions;
#ifdef TRACK_RECOGNIZED_COMMANDS
protected:
  void Expand_Recognized(unsigned int oldlen, unsigned int newlen);
#endif
public:
  Command_Line_Parameters(int _clpargc, char * _clpargv[], const char helpstr[], const char rcfile[] = "");
  int NumParameters() { return numparameters; }
  String ParName(int n) { if (n<numparameters) return parname[n]; else return String(""); }
  String ParValue(int n) { if (n<numparameters) return parvalue[n]; else return String(""); }
  String URI_unescape_ParValue(int n);
#ifdef TRACK_RECOGNIZED_COMMANDS
  //STLlong Specifies_Parameter(String pname) { long n = parname.iselement(pname); if ((n>=0) && ((unsigned long) n < recognized->size())) (*recognized)[n]++; return n; }
  long Specifies_Parameter(String pname) { long n = parname.iselement(pname); if ((n>=0) && (n<recognizedsize) && (recognized!=NULL)) recognized[n]++; return n; }
  int Recognized(unsigned int n) { if ((recognized) && (n<(unsigned int) numparameters)) return recognized[n]; return -1; }
  String Origin(int n) { if (n<numparameters) return origin[n]; else return String("unknown"); }
  int NumUnrecognized();
  String unrecognized_str(String sep = "\n");
#else
  long Specifies_Parameter(String pname) { long n = parname.iselement(pname); if (n<numparameters) return n; return -1; }
#endif
  String Par(int n) { if ((n>=0) && (n<numparameters)) return ParName(n)+'='+ParValue(n); else return String(""); }
  String str(String sep = " ");
  bool IsFormInput() { return calledbyforminput; }
  int IncludeFileErrors() { return includefileerrors; }
  int get_int(int n, int min = INT_MIN, int max = INT_MAX, const char * warn = NULL);
  int init_int(int & parnum, String parlbl, int defaultval, int min = INT_MIN, int max = INT_MAX, const char * warn = NULL);
  double get_double(int n, double min = -DBL_MAX, double max = DBL_MAX, const char * warn = NULL);
  double init_double(int & parnum, String parlbl, double defaultval, double min = -DBL_MAX, double max = DBL_MAX, const char * warn = NULL);
};

extern Command_Line_Parameters * main_clp;

extern String CLP_recentinclude;

class CLP_Modifiable {
  // This is a base class inherited by classes that contain parameters
  // that are modifiable through the CLP interface. This way, all such
  // objects may be automatically registered and called to parse the CLP.
  // (Note 20050101: Registration is not yet implemented.)
public:
  CLP_Modifiable() {}
  virtual ~CLP_Modifiable() {}
  virtual void parse_CLP(Command_Line_Parameters & clp) = 0; // read parameter settings from a CLP object
  virtual String report_parameters() = 0; // human readable explanatory report of parameter settings
  //***  virtual String parameter_instructions() = 0; // list of instructions for parse_CLP() that recreate the parameter settings
};

void report_compiler_directives();

#endif

/*
The concept:

Command line or FORM parameters are a set of individual parameter-value
pairs. The parameters can belong to any valid component within the tree
of functions called as the program is executed. For instance, if running
the nibr program results in execution of the temporary function that
represents a van Pelt developmental simulation then components involved
are: output format and output generation, the network object, parameter
spread sample neuron objects, Network Statistics objects, Connection
Statistics objects, Dendritic and Axonal Growth Model objects.

Parameters can be:
  Default
  Specified

Parameters should be:
  Initialized
  Recorded

Default parameter are included with each object class.
Parameter specification is achieved by:
  Specific initialization of the object, part of the constructor arguments
  Function call of a _set function in the object
  Calling a parse_parameter function in the object for textual specification
    (that is used when parsing with read_commands())
  Retrieval from the Command_Line_Parameters during object construction
    (the CLP is an argument to the constructor and its parameter-value pairs
    can be parsed by the same parse_parameter function as above)

Parameter initialization is taken care of as noted above.

Parameters used during a simulation should be recorded so that the
simulation can be understood and repeated. This can be done by:
  Adding comments to output figures
  Adding text to output figures
  Printing output to standard output
  Storing parameters in a file
The output for all of these can be provided during:
  Initialization and specification of parameters
    (e.g. this allows a command line to be reconstructed later)

To ensure an up to date help text, every object that can parse parameters
should supply a help function.
And those objects should be derived from a base class that guarantees
certain functions through pure virtual functions.

Development progress:

  (1) Status: Most objects contain reasonable default parameter values,
      but some demand that values be specified (non-optimal). Only the
      constructor arguments and _set function call methods are currently
      supported. Objects are not yet able to obtain their own parameters
      by parsing.
      Focus: Complete the Command_Line_Parameters class basics so that
      parameter-value pairs are stored and available. Complete the FORM
      entry interface. Obtain parameters from the CLP (even if the
      method is still a kludge).
      Date: 20041220
  (2) Status: Basic Command_Line_Parameters parsing is in use. A start
      has been made on a CLP_Modifiable base class to implement parsing,
      reporting and storage for reproduction by derived objects. Only the
      parsing member function is currently active. A FORM interface exists.
      Focus: Complete the parse_CLP() function, turn some objects into
      CLP_Modifiable derived objects, add rc-file and other file input
      forms of parameter specification, make multiple specification of the
      same parameter possible (so that the last specification is the active
      one).
      Date: 20050101

 */
