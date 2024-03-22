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
// slice.hh
// Randal A. Koene, 20070223
//
// Classes that contain sample slice generation functions.

/* Slices and Batches:

   Each slice provides a minimum number of parameters needed to be able to compute the data that it should show.
   Additional parameters are either specified in the slice or in a batch of slices that it belongs to.

   At the root level, the special general_slice_parameters_interface class is used to set up a specific static
   slice that can be used to initialize histological slice genereation. Slice generation is enabled if the
   root is given a slice label. E.g.
     "slice=SL"

   The batchsize parameter is used to create batches instead of individual slices. If the batchsize is greater
   than 1 then a batch is created. E.g.
     "SL.batchsize=3"

   The slice label changes how specific slices or batches are addressed. E.g.
     "SL.0.slice=dorsal" means that the slice or batch "SL.0" MUST now be addressed as "SL.dorsal".
     "SL.dorsal.batchsize=100"
     "SL.1.slice=medial SL.medial.batchsize=1"
     "SL.2.slice=ventral"

   Within each batch, elements of a list may themselves be individual slices or batches, e.g.:
     "SL.dorsal.99.batchsize=20"

   The items (single or in a list) contain basic slice information, e.g:
     "SL.medial.SASx=-100.0" "SL.dorsal.0.SPD=100.0,-50.0,25.0"

   Slice vertices are defined as shown in ~/doc/tex/generation-framework/slice-vertices.fig.

   There are many different ways to specify slice vertices. The coordinates of the vertices may be
   identified by anatomical terminology, which translates as follows when viewing in the direction
   of an animal's gaze when the brain is located within the animal's skull, for the identifying
   characters <.123>:
     .S--,.Z--: Superior (upper), z_max
     .I--,.z--: Inferior (lower), z_min
     .-A-,.-Y-; Anterior (front), y_max
     .-P-,.-y-; Posterior (rear), y_min
     .--D,.--X; Dexter (right),   x_max
     .--S,.--x; Sinister (left),  x_min

   The eight vertices that define a slice are SAD, SAS, SPD, SPS, IAD, IAS, IPD and IPS.

   For these vertices, appended axis identifiers x, y or z may be used to specify individual number
   values for the corresponding element of the vertex coordinates.

   Alternatively, all coordinate values in 3D may be specified in a comma-separated list.

   Coordinates may also be obtained in relation to other slices within a batch, if the parameter
   "relativecoords" is set, which may be specified explicitly or by the batch-root:
     SL.dorsal.2.relativecoords=1 SL.dorsal.2.SADz=-1.0
       This example means: Generate coordinates relative to the list element 1 index prior to this one,
       copying all coordinates, maintaining width, height and length, but shifting the Superior
       location (and therefore also the Inferior location) by one thickness of the previous slice.
     SL.dorsal.batchrelativecoords=1 SL.dorsal.nalpha=-1.0
       This example means: Use the batch coordinates as the location for the first slice in the batch,
       then compute following slice coordinates relative to that one, by shifting one slice thickness
       along the alpha normal vector. Note that the coordinates of the batch itself may also be
       specified in relation to another slice/batch in the superior list to which the dorsal batch
       belongs.

   The normal vectors are defined as follows:
     .nalpha: The normal vector through the Superior plane, directed from Inferior to Superior.
     .nbeta : The normal vector through the Anterior plane, directed from Posterior to Anterior.
     .ngamma: The normal vector through the Dexter plane, directed from Sinister to Dexter.

   The normal vectors are also recognized by equivalent anatomical descriptors of the normal vectors:
     .nalpha is .ntransverse
     .nbeta is .ncoronal
     .ngamma is .nsagittal

   The slice definitions also store calculated "inner normals" to each slice surface. These are useful
   vectors to store, since they are used in many volumetric calculations, such as wether a point lies
   within a slice volume. Normals are retained for all six surfaces, in case the slice is not a
   parallelogram. The surface identifiers are derived from vertex identifiers concatenated clockwise
   from an inner perspective as shown in ~/doc/tex/generation-framework/slice-vertices.fig:
     SASSPSSPDSAD
     IASIADIPDIPS
     SASIASIPSSPS
     SASSADIADIAS
     SADSPDIPDIAD
     SPSIPSIPDSPD

   When slice output is generated from by the call in the nibr.cc file, the net->Slice_Output() call
   in turn calls Slice() functions in network, neuron, fibre_structure and other components. To be
   more efficient with fiber and synapses, it does so only for those fibers that lie within Spatial
   Segments that either contain or intersect with defined Slices.

   The first implementation of slice output is a proof of concept and relies on functions already
   implemented for the XFig file format. The resulting files are scalable, can be converted to a
   bitmap format, and can be post-processed as desired.
   In future work output in a 3D file format may be added, as may be output processing within
   NETMORPH.

 */

#ifndef __SLICE_HH
#define __SLICE_HH

#ifdef VECTOR3D

#include "templates.hh"
#include "Command_Line_Parameters.hh"
#include "spatial.hh"
#include "BigString.hh"

#ifndef __FIG_OBJECT_HH
class Fig_Group;
#endif

enum slice_vertices { SAD, SAS, SPD, SPS, IAD, IAS, IPD, IPS, NUM_slice_vertices };

enum slice_innernormals { SASSPSSPDSAD, IASIADIPDIPS, SASIASIPSSPS, SASSADIADIAS, SADSPDIPDIAD, SPSIPSIPDSPD, NUM_slice_innernormals };

extern const char slicevertexID[NUM_slice_vertices][4];

class Slice: public PLLHandle<Slice> {
protected:
  long id;
  Slice * parentbatch;
  unsigned int batchsize;
  PLLRoot<Slice> batch; // if there are no connected Slices, then this is a single slice
  String label;
  spatial vertex[NUM_slice_vertices]; // *** Does this need to be dynamic in order to destruct correctly?
  spatial innernormal[NUM_slice_innernormals]; // *** Does this need to be dynamic in order to destruct correctly?
protected:
  void set_relative_to(double width, double height, double depth);
  void inner_normals();
public:
  Slice(Slice * _parentbatch, long _id, Command_Line_Parameters * clp);
  //Slice(const char * _label, unsigned int _batchsize): batchsize(_batchsize), label(_label) {}
  ~Slice() { }
  void parse_CLP(Command_Line_Parameters & clp);
  String Prefix();
  Slice * iterator_first();
  Slice * iterator_next();
  bool contains(const spatial & centerpoint, double radius);
  Fig_Group * Outlines_Fig();
};

class general_slice_parameters_interface: public CLP_Modifiable, public Slice {
  // This class serves to create the general_slice_parameters
  // object for parameters that should apply to all sample slice generation.
  // A global reference object is initialized in nibr.cc.
protected:
  bool slice; // true if slice output is requested
  bool showsliceoutlines; // true if slice outlines and slice output are requested
public:
  general_slice_parameters_interface(): Slice(NULL,0,NULL), slice(false), showsliceoutlines(false) {}
  virtual ~general_slice_parameters_interface() {}
  virtual void parse_CLP(Command_Line_Parameters & clp);
  virtual String report_parameters();
  bool get_slice() { return slice; }
  bool show_slice_outlines() { return showsliceoutlines; }
};

extern general_slice_parameters_interface general_slice_parameters;

// VECTOR3D
#endif

#endif
