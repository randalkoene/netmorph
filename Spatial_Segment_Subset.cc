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
// Spatial_Segment_Subset.cc
// Randal A. Koene, 20041206

//#include "nibr.hh"
#include "diagnostic.hh"
#include "global.hh"
#include "Color_Table.hh"
#include "Fig_Object.hh"
#include "synapse_formation_model.hh"
#include "Spatial_Segment_Subset.hh"

// normally pass the following define with make
//#define __SPATIAL_SEGMENT_SUBSET_TEST
#ifdef __SPATIAL_SEGMENT_SUBSET_TEST
  int __recursiondepth = 0;
  #define __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_INIT \
    if (!part0) { cout << "ERROR unable to allocate segments\n"; cout.flush(); return true; } \
    int __localrecursiondepth = __recursiondepth; \
    __recursiondepth++;
  #define __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_RETURN \
    __recursiondepth = __localrecursiondepth;
  #define __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_DENDRITIC \
    cout << "Accessing P from " << this << " (dendritic)\n"; cout.flush();
  #define __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_AXONAL \
    cout << "Accessing P from " << this << " (axonal)\n"; cout.flush();
  #define __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_PARTITION \
    cout << "P at " << __recursiondepth << " (" << this << ") LRTB=[" << minvertex.X() << ',' << maxvertex.X() << ',' << minvertex.Y() << ',' << ,maxvertex.Y() << "]\n"; cout.flush(); \
    cout << "  subset_array_size = " << subset_array_size << '\n'; cout.flush(); \
    if ((numdendritesegments<subset_array_size) && (numaxonsegments<subset_array_size)) { cout << "  ERROR partitioning not needed!\n"; cout.flush(); return; }
  #define __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_SUBDIVISIONS \
    cout << "  subdivided by (" << part0 << ',' << part1 << ")\n"; cout.flush(); \
    cout << "  subdividing: numdendritesegments = " << numdendritesegments << " and numaxonsegments = " << numaxonsegments << '\n'; cout.flush();
#else
  #define __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_INIT
  #define __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_RETURN
  #define __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_DENDRITIC
  #define __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_AXONAL
  #define __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_PARTITION
  #define __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_SUBDIVISIONS
#endif

// variables, caches and flags
bool partitioningflag = false;

#ifdef TEST_SSS_DEPTH
unsigned int sss_greatest_depth = 0;
#endif

unsigned int sss_number = 0;

Spatial_Segment_Subset::Spatial_Segment_Subset(Spatial_Segment_Subset * rt, unsigned char divd, int sizelim, double minspan, spatial & minv, spatial maxv, unsigned int d): depth(d), minvertex(minv), maxvertex(maxv), subset_size_limit(sizelim), subset_array_size(sizelim), minimum_span(minspan), numdendritesegments(0), numaxonsegments(0), root(rt), part0(NULL), part1(NULL), dividedimension(divd), segment_state_flag(0), state_flag(0) {
  sss_number++;
#ifdef TEST_SSS_DEPTH
  if (depth>sss_greatest_depth) {
    sss_greatest_depth = depth;
    //    cout << "D = " << depth << '\n'; cout.flush();
  }
#endif
    if (!rt) root = this;
    if (divd>=PDIMSIZE) dividedimension=0; // rotating dimensions
    // [*] Note: This implementation does not include a test to insure that
    // all elements of minvertex are indeed smaller than the corresponding
    // elements in maxvertex.
    DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
    dendritesegments = new fibre_segment_ptr[sizelim];
    axonsegments = new fibre_segment_ptr[sizelim];
    DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,2*sizelim,sizeof(dendritesegments[0]));
  }

Spatial_Segment_Subset::~Spatial_Segment_Subset() {
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
  delete[] dendritesegments;
  delete[] axonsegments;
  DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,2*subset_array_size,-sizeof(fibre_segment_ptr));
}

bool Spatial_Segment_Subset::intersects(Segment & s) {
  // This function uses segment with 3D polyhedron intersection algorithms
  // as described by Dan Sunday in accordance with published algorithms
  // [CYRUS:CLIPPING,FOLEY:CLIPPING,HAINES:RAYCONVEXINTERSECTION,
  // HILL:PERPDOT,LIANG:CLIPPING]
  // Extended line: P(t) = s.P0 + t*(s.P1-s.P0) = s.P0 + t*v
  // with v = s.P1-s.P0
  // The segment s contains those points P(t) with 0<=t<=1.
  // Note that for less obvious normal vectors (here they are unit
  // vectors of a box aligned with the coordinate axes) it is possible
  // to derive plane normal vectors as the cross product:
  // n = (V1-V0) x (V2-V0), whre V0, V1 and V2 are points in the plane,
  // in counter-clockwise order to define an outward n.
  // Here, we merely use the two points minvertex and maxvertex, and
  // the standard normal vectors of the coordinate axes aligned box.
  // (See ./partition-box-normal-vectors.fig.)
  
  // [*] Note: I omit testing for point inclusion if s.P0==s.P1
  if (s.P0==s.P1) {
    warning_Spatial_Segment_Subset_intersects_zero_length++;
    return false;
  }
  double tE = 0; // the maximum entering segment parameter
  double tL = 1; // the minimum leaving segment parameter
  double t, N, D; // intersect parameter t = N / D
  Vector dS = s.P1 - s.P0; // the segment direction vector
  
  //Vector e; // edge vector
  // Vector ne; // edge outward normal (not explicit in code)
  
  // [*] Note: This is specialized for partition boxes with the standard
  // normal vectors, thus N = -dot(ni, s.P0-V[i]) and D = dot(ni,dS) are
  // converted to a rapid switch statement, since most elements of ni
  // are zero.
  // [***INCOMPLETE] I must check what happens with this code when a
  // segment is entirely within (not intersecting) a partition (as in the
  // case of the root partition). If this function returns true, then
  // there is no need for the outside test.
#ifdef VECTOR3D
  for (int i=0; i < 6; i++) { // n=6, process polyhedron faces
#endif
#ifdef VECTOR2D
  for (int i=0; i < 4; i++) { // n=4, process polygon edges
#endif
    switch (i) {
    case 0: // plane; minvertex, (-1,0,0)
      N = s.P0.X()-minvertex.X();
      D = -dS.X();
      break;
    case 1: // plane: maxvertex, (1,0,0)
      N = minvertex.X()-s.P0.X();
      D = dS.X();
      break;
    case 2: // plane: minvertex, (0,-1,0)
      N = s.P0.Y()-minvertex.Y();
      D = -dS.Y();
      break;
    case 3: // plane: maxvertex, (0,1,0)
      N = minvertex.Y()-s.P0.Y();
      D = dS.Y();
      break;
#ifdef VECTOR3D
    case 4: // plane: minvertex, (0,0,-1)
      N = s.P0.Z()-minvertex.Z();
      D = -dS.Z();
      break;
    case 5: // plane: maxvertex, (0,0,1)
      N = minvertex.Z()-s.P0.Z();
      D = dS.Z();
      break;
#endif
    default:
      N = 0.0; D = 1.0;
    }
    if (fabs(D) < SMALL_NUM) { // s is nearly parallel to this face/edge
      if (N < 0) return false; // P0 is outside this face/edge, so s is outside the polyhedron/polygon                               
      else continue; // s cannot cross this face/edge, so ignore this face/edge
    }
    t = N / D;
    if (D < 0) { // segment s is entering through this face/edge
      if (t > tE) { // new max tE
	tE = t;
	if (tE > tL) return false; // s enters after leaving polyhedron/polygon
      }
    } else { // segment s is leaving through this face/edge
      if (t < tL) { // new min tL
	tL = t;
	if (tL < tE) return false; // s leaves before entering polyhedron/polygon
      }
    }
  }
  // tE <= tL implies that there is a valid intersection subsegment
  //IS->P0 = S.P0 + tE * dS;   // = P(tE) = point where S enters polygon
  //IS->P1 = S.P0 + tL * dS;   // = P(tL) = point where S leaves polygon
  return true;
}

/* The following is an alternative implementation of the intersects()
   function for 2D:
    // This is a separate and possibly faster implementation specifically
    // for 2D.
    double x,y, xdiff = x2-x1, ydiff = y2-y1;
    if (xdiff!=0.0) {
      double m = ydiff/xdiff;
      y = (m*(left-x1)) + y1;
      if ((y>=miny) && (y<=maxy) && (y>=top) && (y<=bottom)) return true;
      y = (m*(right-x1)) + y1;
      if ((y>=miny) && (y<=maxy) && (y>=top) && (y<=bottom)) return true;
    }
    if (ydiff!=0.0) {
      double minv = xdiff/ydiff;
      x = (minv*(top-y1)) + x1;
      if ((x>=minx) && (x<=maxx) && (x>=left) && (x<=right)) return true;
      x = (minv*(bottom-y1)) + x1;
      if ((x>=minx) && (x<=maxx) && (x>=left) && (x<=right)) return true;
    }
    return false;
  }
*/

void Spatial_Segment_Subset::partition() {
  // [*** INCOMPLETE/IMPROVE] Partitioning can be made more intelligent if
  // it takes into account the mean orientation of segments and/or the
  // variance of their center points.
  partitioningflag = true; // This changes some of the behavior of add_dendritic_segment() and add_axonal_segment().
  __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_PARTITION
  double middim = (minvertex[dividedimension] + maxvertex[dividedimension])/2.0;
  spatial minv(minvertex);
  spatial maxv(maxvertex);
  maxv[dividedimension] = middim;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
  part0 = new Spatial_Segment_Subset(root,dividedimension+1,subset_size_limit,minimum_span,minv,maxv,depth+1);
  maxv[dividedimension] = maxvertex[dividedimension];
  minv[dividedimension] = middim;
  part1 = new Spatial_Segment_Subset(root,dividedimension+1,subset_size_limit,minimum_span,minv,maxv,depth+1);
  DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,2,sizeof(*part0));
  __SPATIAL_SEGMENT_SUBSET_TEST_REPORT_SUBDIVISIONS
  for (int i=0; i<numdendritesegments; i++) {
    // [***NOTE] I could add a test here, if ((!part0->...) && (!part1->...)) then warning("Lost fiber during partitioning...");
    part0->add_dendritic_segment(dendritesegments[i]);
    part1->add_dendritic_segment(dendritesegments[i]);
  }
  for (int i=0; i<numaxonsegments; i++) {
    // [***NOTE] I could add a test here, if ((!part0->...) && (!part1->...)) then warning("Lost fiber during partitioning...");
    part0->add_axonal_segment(axonsegments[i]);
    part1->add_axonal_segment(axonsegments[i]);
  }
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
  delete[] dendritesegments;
  delete[] axonsegments;
  DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,-2*subset_array_size,sizeof(fibre_segment_ptr));
  dendritesegments = NULL;
  axonsegments = NULL;
  partitioningflag = false;
}

void Spatial_Segment_Subset::expand() {
  DIAGNOSTIC_TEMP_VARIABLES(int,diagoldarraysize,subset_array_size);
  subset_array_size += subset_size_limit;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
  fibre_segment_ptr * tmpdendritesegments = new fibre_segment_ptr[subset_array_size];
  fibre_segment_ptr * tmpaxonsegments = new fibre_segment_ptr[subset_array_size];
  DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,2*subset_array_size,sizeof(fibre_segment_ptr));
  for (int i = 0; i<numdendritesegments; i++) tmpdendritesegments[i] = dendritesegments[i];
  for (int i = 0; i<numaxonsegments; i++) tmpaxonsegments[i] = axonsegments[i];
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Spatial_Segment_Subset);
  delete[] dendritesegments;
  delete[] axonsegments;
  DIAGNOSTIC_AFTER_ALLOCATION(new_Spatial_Segment_Subset,-2*diagoldarraysize,sizeof(fibre_segment_ptr));
  dendritesegments = tmpdendritesegments;
  axonsegments = tmpaxonsegments;
}

void Spatial_Segment_Subset::axon_this_subset_connectivity(fibre_segment * fs) {
  // Determine possible synapses from this axon to dendrites in the subset
  // for this spatial partition.
  state_flag++;
  for (int i=0; i<numdendritesegments; i++) if (dendritesegments[i]->connprocid!=fs) sfm->evaluate_possible_connection(fs,dendritesegments[i]);
}

void Spatial_Segment_Subset::dendrite_this_subset_connectivity(fibre_segment * fs) {
  // Determine possible synapses to this dendrite from axons in the subset
  // for this spatial partition.
  state_flag++;
  for (int i=0; i<numaxonsegments; i++) if (axonsegments[i]->connprocid!=fs) sfm->evaluate_possible_connection(axonsegments[i],fs);
}

void Spatial_Segment_Subset::axon_neighbor_subset_connectivity(fibre_segment * fs, const Spatial_Segment_Subset * fp) {
  // Determine if this subset could contain further subsets that could be
  // neighbors or if it might itself be an immediate neighbor subset
  // (without further spatial divisions) of the partition fp in which the
  // fibre segment fs resides.
  // (This is part of an immediate neighbor search that is usually started
  // at the root of the spatial partition hierarchy.)
  if (this==fp) return;
  if (this!=root) { // if this is not the root test if neighboring is possible
    if (maxvertex.X()<fp->minvertex.X()) return;
    if (minvertex.X()>fp->maxvertex.X()) return;
    if (maxvertex.Y()<fp->minvertex.Y()) return;
    if (minvertex.Y()>fp->maxvertex.Y()) return;
#ifdef VECTOR3D
    if (maxvertex.Z()<fp->minvertex.Z()) return;
    if (minvertex.Z()>fp->maxvertex.Z()) return;
#endif
  }
  if (part0) { // if this has further spatial divisions
    part0->axon_neighbor_subset_connectivity(fs,fp);
    part1->axon_neighbor_subset_connectivity(fs,fp);
    return;
  }
  // [*** INCOMPLETE] Check: The test below is probably unnecessary, except if
  // this is the root, since the tests above assure that it is not unconnected
  // outside fp, and it cannot be within fp, since fp is a leaf of the
  // partition tree. Since it was already determined that this!=fp, then the
  // root partition must have subdivisions, so that if this is the root then
  // the function already took a different path.
  if ((maxvertex.X()==fp->minvertex.X()) ||
      (minvertex.X()==fp->maxvertex.X()) ||
      (maxvertex.Y()==fp->minvertex.Y()) ||
      (minvertex.Y()==fp->maxvertex.Y())
#ifdef VECTOR3D
      || (maxvertex.Z()==fp->minvertex.Z()) ||
      (minvertex.Z()==fp->maxvertex.Z())
#endif
      ) axon_this_subset_connectivity(fs);
}

void Spatial_Segment_Subset::dendrite_neighbor_subset_connectivity(fibre_segment * fs, const Spatial_Segment_Subset * fp) {
  // Determine if this subset could contain further subsets that could be
  // neighbors or if it might itself be an immediate neighbor subset
  // (without further spatial divisions) of the partition fp in which the
  // fibre segment fs resides.
  // (This is part of an immediate neighbor search that is usually started
  // at the root of the spatial partition hierarchy.)
  if (this==fp) return;
  if (this!=root) { // if this is not the root test if neighboring is possible
    if (maxvertex.X()<fp->minvertex.X()) return;
    if (minvertex.X()>fp->maxvertex.X()) return;
    if (maxvertex.Y()<fp->minvertex.Y()) return;
    if (minvertex.Y()>fp->maxvertex.Y()) return;
#ifdef VECTOR3D
    if (maxvertex.Z()<fp->minvertex.Z()) return;
    if (minvertex.Z()>fp->maxvertex.Z()) return;
#endif
  }
  if (part0) { // if this has further spatial divisions
    part0->dendrite_neighbor_subset_connectivity(fs,fp);
    part1->dendrite_neighbor_subset_connectivity(fs,fp);
    return;
  }
  // [*** INCOMPLETE] Check: The test below is probably unnecessary, except if
  // this is the root, since the test above assure that it is not unconnected
  // outside fp, and it cannot be within fp, since fp is a leaf of the
  // partition tree. Since it was already determined that this!=fp, then the
  // root partition must have subdivisions, so that if this is the root then
  // the function already took a different path.
  if ((maxvertex.X()==fp->minvertex.X()) ||
      (minvertex.X()==fp->maxvertex.X()) ||
      (maxvertex.Y()==fp->minvertex.Y()) ||
      (minvertex.Y()==fp->maxvertex.Y())
#ifdef VECTOR3D
      || (maxvertex.Z()==fp->minvertex.Z()) ||
      (minvertex.Z()==fp->maxvertex.Z())
#endif
  ) dendrite_this_subset_connectivity(fs);
}

void Spatial_Segment_Subset::axon_connectivity(fibre_segment * fs) {
  // Finds dendrites in this subset or immediately neighboring subsets that
  // may form a synapse with a specific axon.
  // This is called when an axon is added to a subset array or when all
  // axon segments are evaluated for connectivity after the growth
  // simulation.
  ///if ((eq->T()>=101000.0) && (eq->T()<102000.0)) cout << "\nSUBSET: subset=" << (unsigned long) this << "axonsegment=" << (unsigned long) fs << '\n';
  segment_state_flag |= 0x02;
  // evaluate dendrites in the same partition as the axon
  axon_this_subset_connectivity(fs);
  // find immediately neighboring partitions
  if (this!=root) root->axon_neighbor_subset_connectivity(fs,this);
}

void Spatial_Segment_Subset::dendrite_connectivity(fibre_segment * fs) {
  // Finds axons in this subset or immediately neighboring subsets that
  // may form a synapse with a specific dendrite.
  // This is called when a dendrite is added to a subset array and
  // evaluate_connectivity_immediately is true.
  segment_state_flag |= 0x01;
  // evaluate axons in the same partition as the dendrite
  dendrite_this_subset_connectivity(fs);
  // find immediately neighboring partitions
  if (this!=root) root->dendrite_neighbor_subset_connectivity(fs,this);
}

bool Spatial_Segment_Subset::add_dendritic_segment(fibre_segment * fs) {
  // Note: Assumes that fs->x and fs->y are valid!
  // if (!fs) return false;
  if ((outside(*fs)) && (!intersects(*fs))) return false; // no part of the segment is within this partition
  if (!part0) { // is not yet subdivided
    if (numdendritesegments<subset_array_size) {
      dendritesegments[numdendritesegments++] = fs;
      if (evaluate_connectivity_immediately && (!partitioningflag)) dendrite_connectivity(fs);
      return true;
    }
    if ((sss_number>=maxnumberofspatialsegmentsubsets) || ((maxvertex.X()-minvertex.X())<=minimum_span)
	|| ((maxvertex.Y()-minvertex.Y())<=minimum_span)
#ifdef VECTOR3D
	|| ((maxvertex.Z()-minvertex.Z())<=minimum_span)
#endif
	) {
      expand();
      dendritesegments[numdendritesegments++] = fs;
      if (evaluate_connectivity_immediately && (!partitioningflag)) dendrite_connectivity(fs);
      return true;
    }
    __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_DENDRITIC
    partition(); // subdivide
  }
  // is subdivided
  // Note: You don't need to evaluate the return values here, since we know
  // the segment will be in or pass through at least one partition.
  __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_INIT
  part0->add_dendritic_segment(fs);
  part1->add_dendritic_segment(fs);
  __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_RETURN
  numdendritesegments++; // still gives a valid number of individual segments when subdivided
  return true;
}

bool Spatial_Segment_Subset::add_axonal_segment(fibre_segment * fs) {
  // Note: Assumes that fs->x and fs->y are valid!
  // if (!fs) return false;
  if ((outside(*fs)) && (!intersects(*fs))) return false; // no part of the segment is within this partition
  if (!part0) { // is not yet subdivided
    if (numaxonsegments<subset_array_size) {
      axonsegments[numaxonsegments++] = fs;
      if (evaluate_connectivity_immediately && (!partitioningflag)) axon_connectivity(fs);
      return true;
    }
    if ((sss_number>=maxnumberofspatialsegmentsubsets) || ((maxvertex.X()-minvertex.X())<=minimum_span)
	|| ((maxvertex.Y()-minvertex.Y())<=minimum_span)
#ifdef VECTOR3D
	|| ((maxvertex.Z()-minvertex.Z())<=minimum_span)
#endif
	) {
      expand();
      axonsegments[numaxonsegments++] = fs;
      if (evaluate_connectivity_immediately && (!partitioningflag)) axon_connectivity(fs);
      return true;
    }
    __SPATIAL_SEGMENT_SUBSET_TEST_PARTITION_ACCESS_AXONAL
    partition(); // subdivide
  }
  // is subdivided
  // Note: You don't need to evaluate the return values here, since we know
  // the segment will be in or pass through at least one partition.
  __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_INIT
  part0->add_axonal_segment(fs);
  part1->add_axonal_segment(fs);
  __SPATIAL_SEGMENT_SUBSET_TEST_RECURSE_RETURN
  numaxonsegments++; // still gives a valid number of individual segments when subdivided
  return true;
}

Fig_Object * Spatial_Segment_Subset::net_Fig() {
#ifdef VECTOR3D
  Fig_Group * sssfig = NULL;
  if (part0) {
    long fillstyle = -1;
    long fillcolor = colortable->colnum(CT_background);
    long linecolor = colortable->colnum(CT_partition_border);
    if (subset_array_size>subset_size_limit) fillcolor = colortable->colnum(CT_partition_overfull);
    if ((figattr_show_connection_evaluation) && (state_flag>0)) { 
      //if (segment_state_flag>0) fillstyle = 20;
      //else fillstyle = 19;
      linecolor = colortable->colnum(CT_partition_evaluated);
    }
    DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
    sssfig = new Fig_Group();
    DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(Fig_Group));
    spatial p[4];
    double x_double[4], y_double[4];
    long x[4], y[4];
    for (int j=0; j<4; j++) { // front, back, top and bottom polygons
      switch (j) {
      case 0: // back
	p[0] = minvertex;
	p[1].set_all(maxvertex.X(),minvertex.Y(),minvertex.Z());
	p[2].set_all(maxvertex.X(),maxvertex.Y(),minvertex.Z());
	p[3].set_all(minvertex.X(),maxvertex.Y(),minvertex.Z());
	break;
      case 1: // bottom
	p[0].set_all(minvertex.X(),maxvertex.Y(),minvertex.Z());
	p[1].set_all(maxvertex.X(),maxvertex.Y(),minvertex.Z());
	p[2] = maxvertex;
	p[3].set_all(minvertex.X(),maxvertex.Y(),maxvertex.Z());	
	break;
      case 2: // top
	p[0] = minvertex;
	p[1].set_all(maxvertex.X(),minvertex.Y(),minvertex.Z());
	p[2].set_all(maxvertex.X(),minvertex.Y(),maxvertex.Z());
	p[3].set_all(minvertex.X(),minvertex.Y(),maxvertex.Z());	
	break;
      default: // front
	p[0].set_all(minvertex.X(),minvertex.Y(),maxvertex.Z());
	p[1].set_all(maxvertex.X(),minvertex.Y(),maxvertex.Z());
	p[2] = maxvertex;
	p[3].set_all(minvertex.X(),maxvertex.Y(),maxvertex.Z());	
      }
      for (int i=0; i<4; i++) {
	p[i].plane_mapped(x_double[i],y_double[i]);
	x[i] = FIGSCALED(x_double[i]); y[i] = FIGSCALED(y_double[i]);
      }
      DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
      sssfig->Add_Fig_Object(new Fig_Polygon(2,1,linecolor,fillcolor,998,fillstyle,2.0,0,0,x,y,4));
      DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(Fig_Polygon));
    }
  }
#endif
#ifdef VECTOR2D
  long fillstyle = -1;
  long fillcolor = colortable->colnum(CT_partition_overfull);
  if (subset_array_size>subset_size_limit) fillstyle = 20;
  if ((figattr_show_connection_evaluation) && (state_flag>0)) { 
    if (segment_state_flag>0) fillstyle = 20;
    else fillstyle = 19;
    fillcolor = colortable->colnum(CT_partition_evaluated);
  }
  double left, right, top, bottom;
  minvertex.plane_mapped(left,top);
  maxvertex.plane_mapped(right,bottom);
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
  Fig_Rectangle * sssfig = new Fig_Rectangle(2,1,colortable->colnum(CT_partition_border),fillcolor,998,fillstyle,2.0,0,0,FIGSCALED(left),FIGSCALED(top),FIGSCALED(right),FIGSCALED(bottom));
  DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(Fig_Rectangle));
#endif
  // [***INCOMPLETE] total state_flag values can be collected here
  segment_state_flag = 0;
  state_flag = 0;
  if (!part0) return sssfig;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
  Fig_Group * sssgroup = new Fig_Group();
  DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,1,sizeof(Fig_Group));
  sssgroup->Add_Fig_Object(sssfig);
  sssgroup->Add_Fig_Object(part0->net_Fig());
  sssgroup->Add_Fig_Object(part1->net_Fig());
  return sssgroup;
}

/* 

[***INCOMPLETE]
Problem: When I look at maximum resolution there are sometimes more than
one partition that are the source of an evaluation, which is only possible
if a subdivision occurs during the addition of a segment - but in the
problematic cases there is no segment that actually extends into both of
the source partitions!

The command used for testing is:

./nibr neurons=10 days=1 figuresequence=true figattr_partitions=true figattr_connection_eval=true max_res_samples=true partitionsubsetsize=4

*/
