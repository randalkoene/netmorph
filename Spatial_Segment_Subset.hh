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
// Spatial_Segment_Subset.hh
// Randal A. Koene, 20041206
//
// Classes that enable the grouping of subsets of dendritic or axonal
// fibre segments according to their spatial location. This is used to
// generate connectivity.

#ifndef __SPATIAL_SEGMENT_SUBSET_HH
#define __SPATIAL_SEGMENT_SUBSET_HH

#include <unistd.h>
#include "templates.hh"
#include "spatial.hh"
#include "fibre_structure.hh"

//#define TEST_SSS_DEPTH
#ifdef TEST_SSS_DEPTH
extern unsigned int sss_greatest_depth;
#endif

extern unsigned int sss_number;

class Fig_Object;

class Spatial_Segment_Subset {
  //#define SSS_STATE_CLEAR '\0'
  //#define SSS_STATE_
protected:
  unsigned int depth;
  spatial minvertex, maxvertex; // these describe the spatial region of the partition
  int subset_size_limit; // preferred limit of array size
  int subset_array_size; // actual array size (may be greater)
  double minimum_span; // don't make subsets below regions smaller than this (set this so that you need to test only within a segment and its immediate neighbors)
  fibre_segment_ptr * dendritesegments;
  fibre_segment_ptr * axonsegments;
  int numdendritesegments, numaxonsegments; // these continue to track when subdivided
  Spatial_Segment_Subset * root; // NULL initialization sets root to this
  Spatial_Segment_Subset * part0; // partition pointers leading to binary
  Spatial_Segment_Subset * part1; // subtrees or leafs
  unsigned char dividedimension; // which dimension to subdivide next
  // Note that an alternative implementation would keep segments that
  // do not fit entirely within one of the spatial subdivisions at this
  // level, with connectivity equations to be applied to that segment
  // and those in subdivisions. The drawback of such a method is that the
  // number of such segments is unpredictable, so that a variable size
  // list would need to be used (with its consequent increase in memory
  // usage and allocation time delays). The advantage would be less
  // spatial segmentation caused by a sparse number of long fibres.
  // The same segmentation should be applied to dendrites and axons!
  inline bool outside(Segment & s) {
    return ((Pmax(s.P0.X(),s.P1.X())<minvertex.X()) ||
	    (Pmin(s.P0.X(),s.P1.X())>maxvertex.X()) ||
	    (Pmax(s.P0.Y(),s.P1.Y())<minvertex.Y()) ||
#ifdef VECTOR2D
	    (Pmin(s.P0.Y(),s.P1.Y())>maxvertex.Y()));
#endif
#ifdef VECTOR3D
	    (Pmin(s.P0.Y(),s.P1.Y())>maxvertex.Y()) ||
	    (Pmax(s.P0.Z(),s.P1.Z())<minvertex.Z()) ||
	    (Pmin(s.P0.Z(),s.P1.Z())>maxvertex.Z()));
#endif
  }
  bool intersects(Segment & s);
  void partition();
  void expand();
  void axon_connectivity(fibre_segment * fs);
  void dendrite_connectivity(fibre_segment * fs);
  void axon_this_subset_connectivity(fibre_segment * fs);
  void dendrite_this_subset_connectivity(fibre_segment * fs);
  void axon_neighbor_subset_connectivity(fibre_segment * fs, const Spatial_Segment_Subset * fp);
  void dendrite_neighbor_subset_connectivity(fibre_segment * fs, const Spatial_Segment_Subset * fp);
public:
  Spatial_Segment_Subset(Spatial_Segment_Subset * rt, unsigned char divd, int sizelim, double minspan, spatial & minv, spatial maxv, unsigned int d = 0);
  ~Spatial_Segment_Subset();
  bool add_dendritic_segment(fibre_segment * fs);
  bool add_axonal_segment(fibre_segment * fs);
  Fig_Object * net_Fig();
  // public cache variables
  unsigned char segment_state_flag; // bit0=dendrite/bit1=axon segment added
  int state_flag; // >0 means partition evaluated since last reset of flag
};

// *** If I want to trade space for speed, I can do the following:
// cache x2, y2 and the line gradient m for each fibre segment in the
// fibre_segment class. The caches x, y, x2, y2 and m would be protected
// and automatically updated when the length or angle of the segment are
// changed. Also, changes should then propagate to branches in a rubber
// band manner. All instances of refresh_branch_xy() then become
// unnecessary.

#endif
