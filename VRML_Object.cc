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
// VRML_Object.cc
// Randal A. Koene, 20080226

#include "diagnostic.hh"
#include "global.hh"
#include "spatial.hh"
#include "VRML_Object.hh"

String * VRML_neuronlist = NULL;
String * VRML_synapselist = NULL;
String * VRML_fiberlist = NULL;
String * VRML_continuationnodelist = NULL;
String * VRML_bifurcationnodelist = NULL;
String * VRML_terminalgrowthconelist = NULL;
long VRML_neuronindex = 0;
long VRML_synapseindex = 0;
long VRML_nodeindex = -1;

VRML_Header::VRML_Header(): text("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.2//EN\"   \"http://www.web3d.org/specifications/x3d-3.2.dtd\">\n\
<X3D version='3.0' profile='Interchange'>\n\
  <head>\n\
  </head>\n\
  <Scene>\n\
    <ProtoDeclare name='Soma'>\n\
      <ProtoInterface>\n\
        <field name='cellcolor' type='SFColor' value='.0 .0 .8' accessType='initializeOnly'/>\n\
      </ProtoInterface>\n\
      <ProtoBody>\n\
        <Shape>\n\
          <Appearance>\n\
            <Material DEF='SomaMaterial'>\n\
              <IS>\n\
                <connect nodeField='diffuseColor' protoField='cellcolor'/>\n\
              </IS>\n\
            </Material>\n\
          </Appearance>\n\
          <Sphere radius='5.0'/>\n\
        </Shape>\n\
      </ProtoBody>\n\
    </ProtoDeclare>\n\
    <ProtoDeclare name='NeuritePiece'>\n\
      <ProtoInterface>\n\
        <field name='neuritecolor' type='SFColor' value='.0 .8 .0' accessType='initializeOnly'/>\n\
        <field name='neuriteheight' type='SFFloat' value='10.0' accessType='initializeOnly'/>\n\
        <field name='neuriteradius' type='SFFloat' value='1.0' accessType='initializeOnly'/>\n\
      </ProtoInterface>\n\
      <ProtoBody>\n\
        <Shape>\n\
          <Appearance>\n\
            <Material DEF='NeuriteMaterial'>\n\
              <IS>\n\
                <connect nodeField='diffuseColor' protoField='neuritecolor'/>\n\
              </IS>\n\
            </Material>\n\
          </Appearance>\n\
          <Cylinder>\n\
            <IS>\n\
              <connect nodeField='height' protoField='neuriteheight'/>\n\
              <connect nodeField='radius' protoField='neuriteradius'/>\n\
            </IS>\n\
          </Cylinder>\n\
        </Shape>\n\
      </ProtoBody>\n\
    </ProtoDeclare>\n\
") {
}

String VRML_Group::str() {
  String s;
  // [***NOTE] A comment is probably not allowed in an VRML file.
  /*  if (!comment.empty()) {
    if (comment[comment.length()-1]!='\n') comment += '\n';
    comment[comment.length()-1] = '.';
    comment.gsub("\n","\n# ");
    comment[comment.length()-1] = '\n';
    s = "# " + comment;
    }*/
  // *** Do something here to print the desired lists
  PLL_LOOP_FORWARD(VRML_Object,vrmlobjects.head(),1) s+=e->str();
  return s;
}

void VRML_Group::Add_VRML_Object(VRML_Object * vrmlobject) {
  if (!vrmlobject) return;
  vrmlobjects.link_before(vrmlobject);
}

void X3D_cylinder_center(spatial & p0, spatial & p1, spatial & center) {
  // Computes the center of a cylinder object in X3D from segment
  // start and end points p0 and p1. The vector is:
  // center = 0.5*(p0+p1)
  center = p0;
  center += p1;
  center *= 0.5;
}

bool X3D_SFrotation_to_vector(spatial & p0, spatial & p1, spatial & rotaxis, double & theta, double & len) {
  // This code was inspired by the txt2xml.c code at
  // http://people.ccmr.cornell.edu/~veit/semiproteins/main.html
  // Unlike that code, here we explicitly normalize before determining
  // the vector rotation, and false is returned for zero-length
  // vectors.
  spatial comp(p1);
  comp -= p0;
  if (comp.isZero()) return false;
  len = comp.normalize(); // comp now has length 1.0
#ifdef VECTOR3D
  // [***NOTE] The following suffices for the cross product, because an
  // X3D cylinder is oriented along (0,1,0) by default.
  rotaxis.set_all(comp.Z(),0.0,-comp.X());
#endif
#ifdef VECTOR2D
  // [***Note] This is not from the original code and is presently untested.
  rotaxis.set_all(0.0,0.0,0.0);
#endif
  double mag = rotaxis.len();
  theta = asin(mag);
  if (comp.Y()<0.0) theta = -theta;
  return true;
}
