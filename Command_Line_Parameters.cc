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
// Command_Line_Parameters.cc
// Randal A. Koene, 20041215

#include "Command_Line_Parameters.hh"
#include "file.hh"
#include "global.hh"

Command_Line_Parameters * main_clp;

String CLP_recentinclude;

#ifdef TRACK_RECOGNIZED_COMMANDS
//#define FIXED_TRACKING_ARRAY
#ifdef FIXED_TRACKING_ARRAY
int fixedrecognizedarray[10240];
#endif
//STLstd::vector<int> nullrecognized;
bool initbatchaddparams = false; // This flag is used to prevent inefficient array expansion.

Command_Line_Parameters::Command_Line_Parameters(int _clpargc, char * _clpargv[], const char helpstr[], const char rcfile[]): calledbyforminput(false), clpargc(_clpargc), clpargv(_clpargv), numparameters(0), includefileerrors(0), recognized(NULL), recognizedsize(0), numsubstitutions(0) {
  //STLrecognized(&nullrecognized) {
  initbatchaddparams = true;
#ifdef FIXED_TRACKING_ARRAY
  recognized = fixedrecognizedarray;
#endif
#else
Command_Line_Parameters::Command_Line_Parameters(int _clpargc, char * _clpargv[], const char helpstr[], const char rcfile[]): calledbyforminput(false), clpargc(_clpargc), clpargv(_clpargv), numparameters(0), includefileerrors(0), numsubstitutions(0) {
#endif
  Parse_File(rcfile);
  if ((!Parse_Command_Line()) || (!Detect_Form_Input())) {
    cout << helpstr;
    cout << '\n';
    report_compiler_directives();
    exit(0);
  }
#ifdef TRACK_RECOGNIZED_COMMANDS
  if (numparameters>0) {
#ifndef FIXED_TRACKING_ARRAY
    recognized = new int[numparameters];
    recognizedsize = numparameters;
    //STLcout << "Vector size: " << recognized->size() << '\n'; cout.flush();
    //STLrecognized = new std::vector<int>(numparameters,(int) 0);
    //STLcout << "Vector size: " << recognized->size() << '\n'; cout.flush();
#else
    if (numparameters<10240) recognizedsize=numparameters;
    else recognizedsize=10240;
#endif
    for (int i=0; i<recognizedsize; i++) recognized[i] = 0;
  }
  initbatchaddparams = false;
#endif
}

void remove_preceding_whitespace(String & s) {
  int i = 0, len = s.length();
  for (; i<len; i++) if ((s[i]=='_') || ((s[i]>='A') && (s[i]<='Z')) || ((s[i]>='a') && (s[i]<='z'))) break;
  if ((i<len) && (i>0)) s = s.from(i);
}

#ifdef TRACK_RECOGNIZED_COMMANDS
void Command_Line_Parameters::Expand_Recognized(unsigned int oldlen, unsigned int newlen) {
  if (newlen<=oldlen) return;
  int * r = new int[newlen];
  if (recognized) {
    for (unsigned int i=0; (i<oldlen) && (i<(unsigned int) recognizedsize); i++) r[i] = recognized[i];
    delete[] recognized;
  }
  for (unsigned int i=oldlen; i<newlen; i++) r[i] = 0;
  recognized = r;
  recognizedsize = newlen;
}
#endif

void Command_Line_Parameters::Add_to_Substitutions(String subs) {
  // Add another line to the list of substitutions to be applied before parsing commands.
  if (subs.index(':')<0) error("Error: A substitution string ("+subs+") must consist of two parts separated by a colon ':'.\n");
  substitutionpat[numsubstitutions] = subs.before(':');
  substitutionstr[numsubstitutions] = subs.after(':');
  numsubstitutions++;
}

void Command_Line_Parameters::Apply_Substitutions(String & pname) {
  // A list of substitutions is applied to command names before parsing them for addition
  // to the list of commands.
  for (int i=0; i<numsubstitutions; i++) pname.gsub(substitutionpat[i],substitutionstr[i]);
}

#ifdef TRACK_RECOGNIZED_COMMANDS
void Command_Line_Parameters::Add_Parameter(String pname, String pvalue, String co) {
#else
void Command_Line_Parameters::Add_Parameter(String pname, String pvalue) {
#endif
  long n;
  if (pvalue.empty()) return;
  if (pname=="include") {
    Parse_File(pvalue);
    return;
  }
  if (pname=="substitute") {
    Add_to_Substitutions(pvalue);
    return;
  }
  Apply_Substitutions(pname);
  if ((n=Specifies_Parameter(pname))>=0) {
    parvalue[n] = pvalue; // new value takes precedence over old value
#ifdef TRACK_RECOGNIZED_COMMANDS
    origin[n] = co; // new origin of new value
#endif
  } else {
    parname[numparameters] = pname;
    parvalue[numparameters] = pvalue;
    origin[numparameters] = co;
    numparameters++;
#ifdef TRACK_RECOGNIZED_COMMANDS
#ifdef FIXED_TRACKING_ARRAY
    if (!initbatchaddparams) {
      // [***INCOMPLETE] This breaks if numparameters>10240.
      recognized[numparameters-1]=0;
      if (numparameters<10240) recognizedsize=numparameters;
      else recognizedsize=10240;
    }
#else
    //STLif (!initbatchaddparams) recognized->push_back(0);
    if (!initbatchaddparams) Expand_Recognized(numparameters-1,numparameters);
#endif
#endif
  }
}

void Command_Line_Parameters::Parse_File(const char filename[]) {
  String fname(filename);
  if (fname.empty()) return;
  String parstr;
  if (!read_file_into_String(fname,parstr)) {
    includefileerrors++;
    return;
  }
  CLP_recentinclude = fname;
  // remove comment lines, which may not terminate with ';'
  BigRegex inithashcomments("^[ \t]*#[^\n]*\n"); // (init)hashcomment must be separated, since gsub treats the whole string as one line
  BigRegex hashcomments("\n[ \t]*#[^\n]*\n");
  while (parstr.gsub(hashcomments,"\n")>0) ;
  parstr.gsub(inithashcomments,"");
  BigRegex initCPPcomments("^[ \t]*[/][/][^\n]*\n"); // (init)CPPcomment must be separated, since gsub treats the whole string as one line
  BigRegex CPPcomments("\n[ \t]*[/][/][^\n]*\n");
  while (parstr.gsub(CPPcomments,"\n")>0) ;
  parstr.gsub(initCPPcomments,"");
  BigRegex Ccomments("[/][*].*[*][/]"); // [***UNTESTED] This may not work as expected, as it may find the last instance of a comment end.
  parstr.gsub(Ccomments,"");
  // make a list of parameter=value pairs
  StringList plist(parstr,';'); int pnum = plist.length()-1; // *** empty element at end of StringList
  String pname;
  for (int i = 0; i<pnum; i++) {
    pname = plist[i].before('=');
    remove_preceding_whitespace(pname); // this also removes empty lines
#ifdef TRACK_RECOGNIZED_COMMANDS
    Add_Parameter(pname,plist[i].after('='),fname);
#else
    Add_Parameter(pname,plist[i].after('='));
#endif
  }
}

bool Command_Line_Parameters::Detect_Form_Input() {
  const int LLEN = 10240;
  char lbuf[LLEN];
  char * qstr = getenv("QUERY_STRING"); // GET method form-data
  char * clenstr = getenv("CONTENT_LENGTH"); // POST method form-data
  if ((!qstr) && (!clenstr)) return true;
  calledbyforminput = true; // set flag
  String querystr;
  if (qstr) querystr = qstr; // get GET form-data
  if (clenstr) { // get POST form-data
    int clen = atoi(clenstr);
    if ((clen>0) && (clen<(LLEN-1))) {
      cin.get(lbuf,clen+1);
	 if (!querystr.empty()) querystr += '&';
	 querystr += lbuf;
    }
  }
  if (querystr.length()==0) return true;
  StringList qlist(querystr,'&'); int qnum = qlist.length()-1; // *** empty element at end of StringList
#ifdef TRACK_RECOGNIZED_COMMANDS
  for (int i=0; i<qnum; i++) Add_Parameter(qlist[i].before('='),qlist[i].after('='),"form input");
#else
  for (int i=0; i<qnum; i++) Add_Parameter(qlist[i].before('='),qlist[i].after('='));
#endif
  return true;
}

bool Command_Line_Parameters::Parse_Command_Line() {
  if (clpargc>1) {
    String pname(clpargv[1]);
    if ((pname == "-h") || (pname == "-?")) return false;
#ifdef TRACK_RECOGNIZED_COMMANDS
    Add_Parameter(pname.before('='),pname.after('='),"command line");
#else
    Add_Parameter(pname.before('='),pname.after('='));
#endif
    for (int i = 2; i<clpargc; i++) {
      pname = clpargv[i];
#ifdef TRACK_RECOGNIZED_COMMANDS
      Add_Parameter(pname.before('='),pname.after('='),"command_line");
#else
      Add_Parameter(pname.before('='),pname.after('='));
#endif
    }
  }
  return true;
}

String Command_Line_Parameters::URI_unescape_ParValue(int n) {
  if (n>=numparameters) return String("");
  String escapedURI(parvalue[n]);
  String s; char c; int d1,d2;
  for (unsigned int i = 0; i<escapedURI.length(); i++) {
    c = escapedURI[i];
    if (c=='%') {
      i++; c = escapedURI[i];
      if ((c>='0') && (c<='9')) d1 = (int) (c - '0');
      else if ((c>='A') && (c<='F')) d1 = ((int) (c - 'A')) + 10;
      else d1 = ((int) (c - 'a')) + 10;
      i++; c = escapedURI[i];
      if ((c>='0') && (c<='9')) d2 = (int) (c - '0');
      else if ((c>='A') && (c<='F')) d2 = ((int) (c - 'A')) + 10;
      else d2 = ((int) (c - 'a')) + 10;
      d1 *=16; d1 += d2;
      c = (char) d1;
    }
    s += c;
  }
  return s;
}

String Command_Line_Parameters::str(String sep) {
  String s;
  for (int i=0; i<numparameters; i++) s += parname[i]+'='+parvalue[i]+sep;
  return s;
}

int Command_Line_Parameters::get_int(int n, int min, int max, const char * warn) {
  int res = atoi(ParValue(n));
  if (res<min) {
    res = min;
    if (warn) warning(String(warn)+", setting to minimum value "+String((long) res)+".\n");
  } else if (res>max) {
    res = max;
    if (warn) warning(String(warn)+", setting to maximum value "+String((long) res)+".\n");
  }
  return res;
}

double Command_Line_Parameters::get_double(int n, double min, double max, const char * warn) {
  char * endptr = NULL;
  const char * dblstr = ParValue(n).chars();
  double res = 0.0;
  if (!ParValue(n).empty()) res = strtod(dblstr,&endptr);
  if ((res<min) || (endptr!=(dblstr+ParValue(n).length()))) {
    res = min;
    if (warn) warning(String(warn)+", setting to minimum value "+String(res,"%f")+".\n");
  } else if (res>max) {
    res = max;
    if (warn) warning(String(warn)+", setting to maximum value "+String(res,"%f")+".\n");
  }
  return res;
}

int Command_Line_Parameters::init_int(int & parnum, String parlbl, int defaultval, int min, int max, const char * warn) {
  int res = defaultval;
  if ((parnum=Specifies_Parameter(parlbl))>=0) {
    res = get_int(parnum,min,max,warn);
  }
  return res;
}

double Command_Line_Parameters::init_double(int & parnum, String parlbl, double defaultval, double min, double max, const char * warn) {
  double res = defaultval;
  if ((parnum=Specifies_Parameter(parlbl))>=0) {
    res = get_double(parnum,min,max,warn);
  }
  return res;
}

#ifdef TRACK_RECOGNIZED_COMMANDS
int Command_Line_Parameters::NumUnrecognized() {
  int res = 0;
  //STLfor (unsigned int i=0; i<recognized->size(); i++) if ((*recognized)[i]<1) res++;
  if (recognized) for (int i=0; i<recognizedsize; i++) if (recognized[i]<1) res++;
  return res;
}

String Command_Line_Parameters::unrecognized_str(String sep) {
  String s;
  if (!recognized) return s;
  //STLcout << "GOT INTO UNRECOGNIZED_STR - vector size is " << recognized->size() << " and numparameters=" << numparameters << '\n'; cout.flush();
  //STLfor (unsigned int i=0; ((i<recognized->size()) && (i<(unsigned int) numparameters)); i++) {
  for (int i=0; (i<recognizedsize); i++) {
     //STLcout << (*recognized)[i] << '\n'; cout.flush();
    //STLif ((*recognized)[i]<1) {
    if (recognized[i]<1) {
      s += "Warning: Unused command ";
      s += parname[i];
      s += '=';
      s += parvalue[i];
      s += sep;
      s += "  [origin: ";
      s += origin[i];
      s += ']';
      s += sep;
    }
    //cout << s; cout.flush();
  }
  return s;
}

void report_compiler_directives() {
  String res("Compiler directives:\n");
  res += "  USE_RCFILE: ";
#ifdef USE_RCFILE
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  MEMTESTS: ";
#ifdef MEMTESTS
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  PRECISE_RNG / FAST_RNG / GNUC_RNG:\n";
#ifdef PRECISE_RNG
  res += "     yes           no         no\n";
#else
#ifdef FAST_RNG
  res += "     no            yes        no\n";
#else
  res += "     no            no         yes\n";
#endif
#endif
  res += "  EVALUATE_POSSIBLE_CONNECTION_PROFILING: ";
#ifdef EVALUATE_POSSIBLE_CONNECTION_PROFILING
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  SYNAPTOGENESIS_AND_LOSS_INVENTORY: ";
#ifdef SYNAPTOGENESIS_AND_LOSS_INVENTORY
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS: ";
#ifdef SAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE: ";
#ifdef SAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  ADM_PRIORITIZE_SPEED_OVER_SPACE: ";
#ifdef ADM_PRIORITIZE_SPEED_OVER_SPACE
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  TESTING_SIMPLE_ATTRACTION: ";
#ifdef TESTING_SIMPLE_ATTRACTION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  TESTING_SPIKING: ";
#ifdef TESTING_SPIKING
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  ENABLE_FIXED_STEP_SIMULATION: ";
#ifdef ENABLE_FIXED_STEP_SIMULATION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  FIXED_STEP_SIMULATION_START_AT_DT: ";
#ifdef FIXED_STEP_SIMULATION_START_AT_DT
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  LEGACY_ELONGATION: ";
#ifdef LEGACY_ELONGATION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  TESTING_ELONGATION_TOTAL: ";
#ifdef TESTING_ELONGATION_TOTAL
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  DEBUGGING_ELONGATION: ";
#ifdef DEBUGGING_ELONGATION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  DEBUGGING_DIRECTION: ";
#ifdef DEBUGGING_DIRECTION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  DEBUG_SPHERICAL_BOUNDARY: ";
#ifdef DEBUG_SPHERICAL_BOUNDARY
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  INCLUDE_PDF_SAMPLING: ";
#ifdef INCLUDE_PDF_SAMPLING
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  TRACK_RECOGNIZED_COMMANDS: ";
#ifdef TRACK_RECOGNIZED_COMMANDS
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  CHECK_XFIG_RANGES: ";
#ifdef CHECK_XFIG_RANGES
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION: ";
#ifdef STRICTLY_ADHERE_TO_TAU_C_BINF_RELATION
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  DIVIDE_INITLENGTHS_OVER_ARBORS: ";
#ifdef DIVIDE_INITLENGTHS_OVER_ARBORS
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += "  INCLUDE_SIMPLE_FIBER_DIAMETER:";
#ifdef INCLUDE_SIMPLE_FIBER_DIAMETER
  res += "yes\n";
#else
  res += "no\n";
#endif
  res += '\n';
  report(res);
}

#endif
