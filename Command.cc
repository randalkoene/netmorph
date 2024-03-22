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
// Command.cc
// Randal A. Koene, 20041230

#include <iostream>
using namespace std;
#include "Command.hh"
#include "global.hh"

void Command::parse_command(String cmdstr) {
  remove_whitespace(cmdstr);
  if (cmdstr.empty()) { cmd_number = CMDEMPTY; return; }
  if (cmdstr[0]=='#') { cmd_number = CMDCOMMENT; return; }
  String cmdlh(get_command_left_hand(cmdstr));
  String cmdrh(get_command_right_hand(cmdstr));
  String cmdname(get_command_name(cmdrh));
  String objname; split_object_command(cmdname,objname);
  if (objname.empty()) { // global commands
    if (cmdname=="exit") { cmd_number = CMDEXIT; return; }
    if (cmdname=="figoutput") {
      cmd_number = CMDFIGOUTPUT;
      // process content of brackets here or in execute?
      return;
    }
    if (cmdname=="make_network") {
      cmd_number = CMDMAKENETWORK;
      // process content of brackets here or in execute?
      return;
    }
  } else { // object specific commands
    // get type of object
    // process commands according to type
    if (cmdname=="shape") {
      cmd_number = CMDSHAPE;
      // proces content of brackets here or in execute?
      return;
    }
  }
}

void Command::remove_whitespace(String & cmdstr) {
  cmdstr.gsub(BRXwhite,"");
}

String Command::get_command_left_hand(String & cmdstr) {
  return String(cmdstr.before('='));
}

String Command::get_command_right_hand(String & cmdstr) {
  int rhstart = cmdstr.index('=');
  if (rhstart<0) return cmdstr;
  return String(cmdstr.after('='));
}

String Command::get_command_name(String & cmdstr) {
  int cmdend = cmdstr.index('(');
  if (cmdend<0) return cmdstr;
  return String(cmdstr.before('('));
}

void Command::split_object_command(String & cmdname, String & objname) {
  int objend = cmdname.index('.');
  if (objend<0) return;
  objname = cmdname.before(objend);
  cmdname = cmdname.after(objend);
}

bool Command::cmd_Undefined() {
  warning("Warning: Undefined command ignored ("+cmd_init_str+")\n");
  return true;
}

bool Command::cmd_Exit() {
  progress("Command Script Completed\n");
  return false;
}

bool Command::execute() {
  switch (cmd_number) {
  case CMDUNDEFINED: return cmd_Undefined();
  case CMDEXIT: return cmd_Exit();
  case CMDEMPTY: case CMDCOMMENT: return true;
  case CMDFIGOUTPUT: return true;
  case CMDMAKENETWORK: return true;
  case CMDSHAPE: return true;
  }
  return false;
}

void read_commands() {
  // this reads a command script from standard input
  // Protocol: A script consists of individual commands, one on each line
  // of input delimited with the '\n' character. A multiline command is
  // created by making '\' the last caracter before the line delimiter.
  // Note that a script may be generated dynamically, so that new commands
  // are produced in response to simulation output.
  const int CMDBUFFERSIZE = 1024;
  char cmdbuffer[CMDBUFFERSIZE];
  String nextcommand("");
  while (!cin.eof()) {
    cin.getline(cmdbuffer,CMDBUFFERSIZE);
    if (cin.gcount()==0) break;
    nextcommand += cmdbuffer;
    if (!nextcommand.empty())
      if (nextcommand.lastchar()!='\\') { // execute command
	Command cmd(nextcommand);
	if (!cmd.execute()) break;
	nextcommand = ""; // clear for following command
      }
  }
}

