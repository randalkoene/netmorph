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
// slice.cc
// Randal A. Koene, 20070223

#include <limits>
#include <cfloat>
#include "slice.hh"
#include "StringList.hh"
#include "BigRegex.hh"
#include "Fig_Object.hh"
#include "global.hh"

#ifdef VECTOR3D

// The minimum width, height or depth of a slice volume is 0.001 micron (1 nanometer).
#define MINSLICESPAN 0.001

general_slice_parameters_interface general_slice_parameters;

const char slicevertexID[NUM_slice_vertices][4] = {
  "SAD", 
  "SAS", 
  "SPD", 
  "SPS", 
  "IAD", 
  "IAS", 
  "IPD", 
  "IPS"
};
StringList slicevertexIDstr("SAD SAS SPD SPS IAD IAS IPD IPS",' ');

String spatial_str(const spatial & p, char sep = ',') {
#ifdef VECTOR3D
  return String(p.X(),"%.3f")+sep+String(p.Y(),"%.3f")+sep+String(p.Z(),"%.3f");
#endif
#ifdef VECTOR2D
  return String(p.X(),"%.3f")+sep+String(p.Y(),"%.3f");
#endif
}

void range_check(double & res, double min, double max, char warnchar, const char * warn = NULL) {
  if (res<min) {
    res = min;
    if (warn) warning(String(warn)+", setting "+String(warnchar)+" to minimum value "+String(res,"%.3f")+".\n");
  } else if (res>max) {
    res = max;
    if (warn) warning(String(warn)+", setting "+String(warnchar)+" to maximum value "+String(res,"%.3f")+".\n");
  }
}

void get_spatial(Command_Line_Parameters & clp, int n, spatial & res, double min = -DBL_MAX, double max = DBL_MAX, const char * warn = NULL) {
  const BigRegex BRXspatialxy("\\(-?\\(\\([0-9]+\\.[0-9]*\\)\\|\\([0-9]+\\)\\|\\(\\.[0-9]+\\)\\)\\([eE][---+]?[0-9]+\\)?\\)[ \t,;]+\\(-?\\(\\([0-9]+\\.[0-9]*\\)\\|\\([0-9]+\\)\\|\\(\\.[0-9]+\\)\\)\\([eE][---+]?[0-9]+\\)?\\)",1,400); // [***NOTE] We cannot check directly for x, y and z, since regex-gnu or BigRegex does not seem to be able to handle more than 16 match registers.
  String v(clp.ParValue(n));
  if (v.index(BRXspatialxy)<0) { // check x and y
    if (warn) warning(String(warn)+", spatial coordinates ("+spatial_str(res)+") unmodified\n");
    return;
#ifdef VECTOR3D
  } else if (v.index(BRXdouble,BRXspatialxy.subpos(0)+BRXspatialxy.sublen(0))<0) { // check z
    if (warn) warning(String(warn)+", spatial coordinates ("+spatial_str(res)+") unmodified\n");
    return;
#endif
  }
  double x = atof(v.sub(BRXspatialxy,1));
  range_check(x,min,max,'x',warn);
  double y = atof(v.sub(BRXspatialxy,7));
  range_check(y,min,max,'y',warn);
#ifdef VECTOR3D
  double z = atof(v.sub(BRXdouble,0));
  range_check(z,min,max,'z',warn);
  res.set_all(x,y,z);
#endif
#ifdef VECTOR2D
  res.set_all(x,y);
#endif
}

void init_spatial(Command_Line_Parameters & clp, int & parnum, String parlbl, spatial & res, double min = -DBL_MAX, double max = DBL_MAX, const char * warn = NULL) {
  // It is assumed that res already contains its default values.
  if ((parnum=clp.Specifies_Parameter(parlbl))>=0) {
    get_spatial(clp,parnum,res,min,max,warn);
  }
}

Slice::Slice(Slice * _parentbatch, long _id, Command_Line_Parameters * clp): id(_id), parentbatch(_parentbatch), batchsize(1) {
  if (clp) parse_CLP(*clp);
}

String Slice::Prefix() {
  if (!label.empty()) return label+'.';
  else if (parentbatch) return parentbatch->Prefix()+String(id)+'.';
  else return label; // no prefix or prefix separating period at unlabeled root
}

void Slice::set_relative_to(double width, double height, double depth) {
  // This is always relative to vertex SAD.
  vertex[SAS] = vertex[SAD]; vertex[SAS].set_X(vertex[SAS].X()-width);
  vertex[IAD] = vertex[SAD];
#ifdef VECTOR3D
  vertex[IAD].set_Z(vertex[IAD].X()-depth);
#endif
  vertex[SPD] = vertex[SAD]; vertex[SPD].set_Y(vertex[SPD].Y()-height);
  vertex[IAS] = vertex[SAS];
#ifdef VECTOR3D
  vertex[IAS].set_Z(vertex[IAS].Z()-depth);
#endif
  vertex[IPD] = vertex[IAD]; vertex[IPD].set_Y(vertex[IPD].Y()-height);
  vertex[SPS] = vertex[SAS]; vertex[SPS].set_Y(vertex[SPS].Y()-height);
  vertex[IPS] = vertex[SPS];
#ifdef VECTOR3D
  vertex[IPS].set_Z(vertex[IPS].Z()-depth);
#endif
}

void Slice::parse_CLP(Command_Line_Parameters & clp) {
  int n;
  String thisslice(Prefix());
  if ((n=clp.Specifies_Parameter(thisslice+"slice"))>=0) {
    if (parentbatch) label=parentbatch->Prefix()+clp.ParValue(n);
    else label=clp.ParValue(n);
    thisslice = Prefix();
  }
  batchsize = clp.init_int(n,thisslice+"batchsize",batchsize,1,INT_MAX,String("Warning: Specified "+thisslice+"batchsize invalid"));
  if (batchsize>1) {
    batch.clear();
    for (unsigned int i=0; i<batchsize; i++) batch.link_before(new Slice(this,i,&clp));
  }
  for (slice_vertices i = SAD; i<NUM_slice_vertices; i = slice_vertices(i+1)) {
    String pname(thisslice+slicevertexID[i]);
    init_spatial(clp,n,pname,vertex[i],-DBL_MAX,DBL_MAX,String("Warning: Specified "+pname+" invalid"));
    double x = clp.init_double(n,pname+"x",vertex[i].X(),-DBL_MAX,DBL_MAX,String("Warning: Specified "+pname+"x invalid"));
    double y = clp.init_double(n,pname+"y",vertex[i].Y(),-DBL_MAX,DBL_MAX,String("Warning: Specified "+pname+"y invalid"));
#ifdef VECTOR3D
    double z = clp.init_double(n,pname+"z",vertex[i].Z(),-DBL_MAX,DBL_MAX,String("Warning: Specified "+pname+"z invalid"));
    vertex[i].set_all(x,y,z);
#endif
#ifdef VECTOR2D
    vertex[i].set_all(x,y);
#endif
  }
  slice_vertices relto = SAD; // default for relative slice volume extension
  if ((n=clp.Specifies_Parameter(thisslice+"origin"))>=0) {
    String originstr(upcase(clp.ParValue(n)));
    for (slice_vertices i = SAD; i<NUM_slice_vertices; i = slice_vertices(i+1)) if (originstr==slicevertexIDstr[(long unsigned int) i]) { relto = i; break; }
  }
  double swidth = clp.init_double(n,thisslice+"width",0.0,MINSLICESPAN,DBL_MAX,String("Warning: Specified "+thisslice+"width invalid"));
  double sheight = clp.init_double(n,thisslice+"height",0.0,MINSLICESPAN,DBL_MAX,String("Warning: Specified "+thisslice+"height invalid"));
#ifdef VECTOR3D
  double sdepth = clp.init_double(n,thisslice+"depth",0.0,MINSLICESPAN,DBL_MAX,String("Warning: Specified "+thisslice+"depth invalid"));
#endif
#ifdef VECTOR2D
  double sdepth = 0.0; // placeholder
  if ((swidth+sheight)>0.0) sdepth = 1.0; // so that sdepth is not the cause of a warning
#endif
  if ((swidth+sheight+sdepth)>0.0) { // attempt relative slice volume extension
    if ((swidth<=0.0) || (sheight<=0.0) || (sdepth<=0.0)) warning("Warning: Relative slice volume extension for "+thisslice+" did not specify all extents\n");
    else {
      // Adjust to origin
      switch (relto) {
      case SAS: swidth = -swidth; break;;
      case IAD: sdepth = -sdepth; break;;
      case SPD: sheight = -sheight; break;;
      case IAS: swidth = -swidth; sdepth=-sdepth; break;;
      case IPD: sheight = -sheight; sdepth=-sdepth; break;;
      case SPS: sheight = -sheight; swidth=-swidth; break;;
      case IPS: sheight = -sheight; swidth=-swidth; sdepth=-sdepth; break;;
      default: break;;
      }
      set_relative_to(swidth,sheight,sdepth);
    }
  }
  inner_normals();
}

void Slice::inner_normals() {
  // SASSPSSPDSAD
  spatial e1(vertex[SPS]); e1 -= vertex[SAS];
  spatial e2(vertex[SPD]); e2 -= vertex[SPS];
  innernormal[SASSPSSPDSAD] = e1; innernormal[SASSPSSPDSAD] ^= e2;
  //innernormal[SASSPSSPDSAD] = e1 ^ e2;
  // Testing to see if the innernormal is properly pointing inward (see figure 4.4 of http://www.doc.ic.ac.uk/~dfg/graphics/GraphicsLecture04.pdf)
  spatial tst(vertex[IPS]); tst -= vertex[SAS]; // This is "(B-A)"
  if ((innernormal[SASSPSSPDSAD]*tst) <= 0) innernormal[SASSPSSPDSAD].Negate();
  // IASIADIPDIPS
  e1 = vertex[IAD]; e1 -= vertex[IAS];
  e2 = vertex[IPD]; e2 -= vertex[IAD];
  innernormal[IASIADIPDIPS] = e1; innernormal[IASIADIPDIPS] ^= e2;
  tst = vertex[SPS]; tst -= vertex[IAS];
  if ((innernormal[IASIADIPDIPS]*tst) <= 0) innernormal[IASIADIPDIPS].Negate();
  // SASIASIPSSPS
  e1 = vertex[IAS]; e1 -= vertex[SAS];
  e2 = vertex[IPS]; e2 -= vertex[IAS];
  innernormal[SASIASIPSSPS] = e1; innernormal[SASIASIPSSPS] ^= e2;
  tst = vertex[SAD]; tst -= vertex[IAS];
  if ((innernormal[SASIASIPSSPS]*tst) <= 0) innernormal[SASIASIPSSPS].Negate();
  // SASSADIADIAS
  e1 = vertex[SAD]; e1 -= vertex[SAS];
  e2 = vertex[IAD]; e2 -= vertex[SAD];
  innernormal[SASSADIADIAS] = e1; innernormal[SASSADIADIAS] ^= e2;
  tst = vertex[IPS]; tst -= vertex[SAS];
  if ((innernormal[SASSADIADIAS]*tst) <= 0) innernormal[SASSADIADIAS].Negate();
  // SADSPDIPDIAD
  e1 = vertex[SPD]; e1 -= vertex[SAD];
  e2 = vertex[IPD]; e2 -= vertex[SPD];
  innernormal[SADSPDIPDIAD] = e1; innernormal[SADSPDIPDIAD] ^= e2;
  tst = vertex[IAS]; tst -= vertex[SAD];
  if ((innernormal[SADSPDIPDIAD]*tst) <= 0) innernormal[SADSPDIPDIAD].Negate();
  // SPSIPSIPDSPD
  e1 = vertex[IPS]; e1 -= vertex[SPS];
  e2 = vertex[IPD]; e2 -= vertex[IPS];
  innernormal[SPSIPSIPDSPD] = e1; innernormal[SPSIPSIPDSPD] ^= e2;
  tst = vertex[IAS]; tst -= vertex[SPS];
  if ((innernormal[SPSIPSIPDSPD]*tst) <= 0) innernormal[SPSIPSIPDSPD].Negate();
}

Slice * Slice::iterator_first() {
  // This seeks into depth to find the first end-point along this set of branches.
  // See the example: ~/src/nnmodels/nibr/slice-iteration.fig.
  if (batch.head()) return batch.head()->iterator_first();
  else return this;
}

Slice * Slice::iterator_next() {
  // This seeks the next at this level, then finds the first end-point along there.
  // See the example: ~/src/nnmodels/nibr/slice-iteration.fig.
  if (Next()) return Next()->iterator_first();
  else if (parentbatch) return parentbatch->iterator_next(); // try it one level up
  else return NULL; // this was the last slice in the tree
}

bool Slice::contains(const spatial & centerpoint, double radius) {
  // True if a sphere (or circle) with the given radius intersects the slice.
#ifdef INCLUDE_INCOMPLETE_SLICE_CODE
  Segment centersegment(centerpoint,centerpoint); // [***INCOMPLETE] See TL#200705011444.
#endif
  // Compare the radius and distance to all slice boundaries.

  return false;
}

Fig_Group * Slice::Outlines_Fig() { // [***INCOMPLETE]
  Fig_Group * fg = NULL;
  return fg;
}

void general_slice_parameters_interface::parse_CLP(Command_Line_Parameters & clp) {
  Slice::parse_CLP(clp);
  slice = (!label.empty());
  int n;
  if ((n=clp.Specifies_Parameter("outattr_slice_outlines"))>=0) showsliceoutlines = (downcase(clp.ParValue(n))==String("true"));
}

String general_slice_parameters_interface::report_parameters() {
  String res("Slice specification: ");
  return res;
}

// VECTOR3D
#endif
