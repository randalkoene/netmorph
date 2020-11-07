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
// StringList.hh
// Randal A. Koene, 20041119

#ifndef __STRINGLIST_HH
#define __STRINGLIST_HH

#include "BigString.hh"
#include "BigRegex.hh"

class StringList {
protected:
  StringList * Next;
  String s;
public:
  ~StringList() { if (Next) delete Next; }
  StringList() { Next=NULL; }
  StringList(unsigned long n) { if (n) Next = new StringList(n-1); else Next = NULL; }
  StringList(const String cs, StringList * cn = NULL) { s=cs; Next=cn; }
  StringList(String cs, const char sep);
  String & operator[](unsigned long n); // (StringList.cc)
  StringList * get_Next() { return Next; }
  StringList * List_Index(long n); // (StringList.cc)
  // append to n (without breaking link to n+1)
  void append(unsigned long n, String & ins); // (StringList.cc)
  // insert at n (n becomes n+1)
  void insert(unsigned long n, String & ins); // (StringList.cc)
  // count number of items following this one
  unsigned long length() { if (Next) return 1+Next->length(); else return 1; }
  // find in list
  long iselement(String el, unsigned long n = 0); // (StringList.cc)
  // find substring in list
  long contains(String substr, unsigned long n = 0); // (StringList.cc)
  // find first in list that is substring of source string
  long first_contained_in(String & sourcestr, unsigned long n = 0); // (StringList.cc)
  // split a string into list elements
  unsigned long split(unsigned long n, String cs, const char sep = '\n');
  // concatenate list
  String concatenate(String sepstr = ""); // (StringList.cc)
};

#endif
