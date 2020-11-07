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
// file.hh
// Randal A. Koene, 20040830
//
// File read and write functions.

#ifndef __FILE_HH
#define __FILE_HH

#include <fstream>
#include "BigString.hh"

bool read_file_into_String(ifstream & sfl, String & filestr, bool warnnoopen = true);
bool read_file_into_String(String filename, String & filestr, bool warnnoopen = true);
bool write_file_from_String(String filename, const String & filestr);
bool append_file_from_String(String filename, const String & filestr);

#endif
