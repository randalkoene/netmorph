/*
  © Copyright 2008-2024 Randal A. Koene <randalk@netmorph.org>
  
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
// Netmorph.cpp
// Randal A. Koene, 20240628

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
#include "Netmorph.h"

/*!
  \mainpage NETMORPH Class Documentation and UML

  \section Entry function main()

  \image html /home/randalk/src/nnmodels/nibr/doc/main.png "[UML] Actions of the NETMORPH program entry function main()."
  \image latex /home/randalk/src/nnmodels/nibr/doc/main.png "[UML] Actions of the NETMORPH program entry function main()."

  \section The simulation loop
 */
int main2(std::string & clpcontent) {
  copyright();
  Command_Line_Parameters clp(clpcontent, helpstr);
// #ifdef USE_RCFILE
//   Command_Line_Parameters clp(argc,argv,helpstr,RCFILE);
// #else
//   Command_Line_Parameters clp(argc,argv,helpstr);
// #endif
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
  return 0;
}

/**
 * @brief Run netmorph on the provided modelfile. Will write the current status to the provided int pointer (0-100%).
 * 
 * Note that _Modelfile is not a path to a file but rather contains the entire content of a model file.
 * 
 * @param _StatusPercent 
 * @return NetmorphResult Struct containing all produced data. You must check Status is equal to true, otherwise all other data may not be initialized.
 */
NetmorphResult Netmorph(int* _StatusPercent, std::string _Modelfile) {

  if (_StatusPercent != nullptr) {
    progress_percentage = _StatusPercent;
  }

  NetmorphResult res;
  int res_int = main2(_Modelfile);
  res.Status = (res_int == 0);

  return res;
}
