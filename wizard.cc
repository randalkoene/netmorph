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
// wizard.cc
// Randal A. Koene, 20080204

// See TL#200802040900.1.

#include "wizard.hh"
#include "StringList.hh"

StringList script;

char binary_choice(char _default, char alt) {
  String bc("?");
  while ((bc[0] != _default) && (bc[0] != alt) && (!bc.empty())) {
    readline(cin,bc);
    bc = downcase(bc);
  }
  if (bc.empty()) bc = _default;
  return bc[0];
}

String integer_value() {
  String iv;
  do {
    readline(cin,iv);
  } while(!iv.matches(BRXint));
  return iv;
}

void terminal() {
  // Terminal version of the wizard, which simply runs through a series of
  // questions and answers to create a valid NETMORPH script.
  cout << "Welcome to the NETMORPH Script Wizard (version " << version << ")\n\n\
The Terminal I/O function of this wizard will guide you through a\n\
series of questions and answers to create a valid NETMORPH script\n\
customized to your simulation aims.\n\n\
Please note that options in capital letters are default choices made if you\n\
simply press ENTER.\n\n\
SECTION 1: GENERAL SIMULATION PARAMETERS\n\n\
Would you like to [S]pecify a random seed to make simulation results\n\
replicable, or would you like to use a [d]ifferent random seed on each\n\
run? (S/d) => ";
  String randomseedstr("randomseed=-1;");
  if (binary_choice('s','d')=='s') {
    cout << "\nPlease enter an integer number as random seed: ";
    randomseedstr = "randomseed="+integer_value()+';';
    script.append(script.length(),randomseedstr);
  } else script.append(script.length(),randomseedstr);

  cout << "Wizard data entry has been completed.\n\n";
  String scripttext(script.concatenate("\n"));
  cout << scripttext;
}

int main(int argc, char * argv[]) {

  terminal();

  exit(0);
}
