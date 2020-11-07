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
// Results.hh
// Randal A. Koene, 20041118
//
// Classes that provide information about process results.

#ifndef __RESULTS_HH
#define __RESULTS_HH

#include <iostream>
#include "BigString.hh"

// classes

// <A NAME="class-result">class Result: flexible processing of function results</A>
class Result {
public:
  Result() {}
  String message;
  inline friend ostream & operator<<(ostream & s, const Result & x);
};

class Shape_Result: public Result {
public:
  Shape_Result(): displaywidth(1.0) {};
  double displaywidth;
};

class Shape_Hexagon_Result: public Shape_Result {
public:
  Shape_Hexagon_Result() {};
  double rectangulararea;
  int neuronsperline;
};

class Shape_Rectangle_Result: public Shape_Result {
public:
  Shape_Rectangle_Result() {};
  double rectangulararea;
  int neuronsperline;
};

class Shape_Box_Result: public Shape_Result {
public:
  Shape_Box_Result() {};
  double boxvolume;
  //[***REMOVE]double rectangulararea;
  int neuronsperline;
  int neuronsdeep;
  int neuronsperxzplane;
};

class Shape_Regions_Result: public Shape_Result {
  // Required parameters are already included in the base class.
  // Parameters defined here are largely used for calculations
  // during shaping.
public:
  Shape_Regions_Result() {};
};

class Shape_Circle_Result: public Shape_Result {
public:
  Shape_Circle_Result() {};
  double circlearea;
  int maxneuronsperline;
};

// function declarations

ostream & operator<<(ostream & s, const Result & x);

// inline function definitions

inline ostream & operator<<(ostream & s, const Result & x) {
  s << x.message.chars();
  return s;
}

#endif
