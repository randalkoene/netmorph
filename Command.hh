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
// Command.hh
// Randal A. Koene, 20041118
//
// Classes that perform command parsing.

#ifndef __COMMAND_HH
#define __COMMAND_HH

//#include <unistd.h>
//#include "global.hh"
#include "BigString.hh"
//#include "BigRegex.hh"
//#include "templates.hh"
//#include "network.hh"

/* more about command parsing:
grammer_parse(str) {
-> divide str into LH, RH, operation
-> operations & functions (e.g. +, .something()) call parsers WITHIN objects
-> variable assignments set up arrays with object pointers
}
automatically give each function a command method - that way being able
to create a connection automatically means specific connections can be
set
setting connection distribution example:
  # functions:
    distfuncab = probdist_numbydistance(20,1,0.6)
  # set selection:
    a = all[interneuron]
    b = all[principal]
  # connections:
    a.make_connections_with_distribution_function(b,distfuncab)
 */

enum cmd_type { CMDUNDEFINED, CMDEXIT, CMDEMPTY, CMDCOMMENT, CMDFIGOUTPUT, CMDMAKENETWORK, CMDSHAPE };

// classes

class Command {
protected:
  String cmd_init_str;
  cmd_type cmd_number;
  void parse_command(String cmdstr);
  void remove_whitespace(String & cmdstr);
  String get_command_left_hand(String & cmdstr);
  String get_command_right_hand(String & cmdstr);
  String get_command_name(String & cmdstr);
  void split_object_command(String & cmdname, String & objname);
  bool cmd_Undefined();
  bool cmd_Exit();
public:
  Command(String cmdstr): cmd_init_str(cmdstr), cmd_number(CMDUNDEFINED) { parse_command(cmdstr); }
  bool execute();
};

// function declarations

void read_commands();

#endif
