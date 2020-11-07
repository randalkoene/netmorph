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
// StringList.cc
// Randal A. Koene, 20041119

#include "StringList.hh"

StringList::StringList(String cs, const char sep) {
// initialize list from a string with a specified separator
// (Note: In the current implementation this leaves an unused
// list element at the end of the list.
// Beware! If you fix this to leave no unused list element at
// the end of the list you must check all source code where
// this initialization function is used to insure that the
// code is modified so as not to depend on this "feature"!)
	Next=NULL;
	split(0,cs,sep);
}

String & StringList::operator[](unsigned long n) {
	if (!n) return s;
	if (!Next) Next = new StringList();
	return (*Next)[n-1];
}

StringList * StringList::List_Index(long n) {
  if (n<0) return NULL;
  if (!n) return this;
  if (!Next) return NULL;
  return Next->List_Index(n-1);
}

void StringList::append(unsigned long n, String & ins) {
	if (!n) Next = new StringList(ins,Next);
	else
		if (Next) Next->append(n-1,ins);
		else {
			Next = new StringList();
			Next->append(n-1,ins);
		}
}

void StringList::insert(unsigned long n, String & ins) {
	if (!n) { Next = new StringList(s,Next); s = ins; }
	else
		if (Next) Next->insert(n-1,ins);
		else 
			if (n>1) {
				Next = new StringList();
				Next->insert(n-1,ins);
			} else Next = new StringList(ins);
}

long StringList::iselement(String el, unsigned long n) {
	if (el == s) return n;
	else {
		if (Next) return Next->iselement(el,n+1);
		else return -1;
	}
}

long StringList::contains(String substr, unsigned long n) {
	if (s.contains(substr)) return n;
	else {
		if (Next) return Next->contains(substr,n+1);
		else return -1;
	}
}

long StringList::first_contained_in(String & sourcestr, unsigned long n) {
  if (sourcestr.contains(s)) return n;
  else {
    if (Next) return Next->first_contained_in(sourcestr,n+1);
    else return -1;
  }
}

unsigned long StringList::split(unsigned long n, String cs, const char sep) {
// insert list elements from a string with a specified separator
	if (cs.length()<=0) return length();
	String el;
	if (cs.index(sep)>=0) el = cs.before(sep); else el = cs;
	// *** this is not a speed optimized implementation
	insert(n,el);
	split(n+1,cs.after(sep),sep);
	return length();
}

String StringList::concatenate(String sepstr) {
	if (Next) return s+sepstr+Next->concatenate(sepstr);
	else return s;
}
