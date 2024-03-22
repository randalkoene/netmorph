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
// file.cc
// Randal A. Koene, 20041230

#include "file.hh"
#include "global.hh"

bool read_file_into_String(ifstream & sfl, String & filestr, bool warnnoopen) {
// read from an open file stream into the string filestr,
// optionally warning about the inability to open the file
// returns true if file content was read into filestr
  if (sfl) {
    filestr = sfl;
    sfl.close();
    return true;
  } else {
    if (warnnoopen) warning("Warning: Unable to open file in read_file_into_String()\n");
    return false;
  }
}

bool read_file_into_String(String filename, String & filestr, bool warnnoopen) {
// read the file called filename into the string filestr,
// optionally warning about the inability to open the file
// returns true if file content was read into filestr
  ifstream sfl(filename);
  if (read_file_into_String(sfl,filestr,false)) return true;
  else {
    if (warnnoopen) warning("Warning: Unable to open "+filename+" in read_file_into_String()\n");
    return false;
  }
}

bool write_file_from_String(String filename, const String & filestr) {
// Write a file from string content, with automatic backup of
// a previous version (if one exists).
  String newfilename = filename;
  ofstream ntl(newfilename);
  if (!ntl) {
    warning("Warning: Unable to create "+newfilename+" in write_file_from_String()\n");
    return false;
  }
  ntl << filestr;
  ntl.close();
  return true;
}

bool append_file_from_String(String filename, const String & filestr) {
// Append to a file from string content, with automatic backup of
// a previous version (if one exists).
  String newfilename = filename;
  ofstream ntl(newfilename,ios_base::out | ios_base::ate | ios_base::app);
  if (!ntl) {
    warning("Warning: Unable to open "+newfilename+" in append_file_from_String()\n");
    return false;
  }
  ntl << filestr;
  ntl.close();
  return true;
}
