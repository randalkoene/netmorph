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
// Fig_Object.cc
// Randal A. Koene, 20041118

#include "diagnostic.hh"
#include "Fig_Object.hh"

unsigned long int Fig_output_modified = 0;
long Fig_output_min_excess = 0;
long Fig_output_max_excess = 0;

double tex_textwidth = 6.5;
double figscale = 1.0;
bool figattr_Fig_rescale = true;

void XFIG_RANGE_CHECK_INIT() {
  Fig_output_modified = 0;
  Fig_output_min_excess = 0;
  Fig_output_max_excess = 0;
}

String XFIG_RANGE_CHECK_WARNING(String fname, double outwidth) {
  String warnstr("WARNING: Output in "+fname+" has been MODIFIED "+String((long) Fig_output_modified)+" times to fit maximum XFig coordinates.\n");
  warnstr += "  Greatest negative excess value encountered = " + String(Fig_output_min_excess) + '\n';
  warnstr += "  Greatest positive excess value encountered = " + String(Fig_output_max_excess) + '\n';
  warnstr += "  Current network cell population output width requested (figureTeXwidth): " + String(outwidth,"%f\n");
  double minfactor = Fig_output_min_excess/((double) XFIG_RANGE_MIN);
  double maxfactor = Fig_output_max_excess/((double) XFIG_RANGE_MAX);
  if (minfactor>maxfactor) maxfactor = minfactor;
  warnstr += "  Suggested setting: figureTeXwidth=" + String(outwidth/maxfactor,"%f;\n");
  return warnstr;
}

#ifdef CHECK_XFIG_RANGES
long XFIG_RANGE_CHECK(long v) {
  if (v<XFIG_RANGE_MIN) {
    if (v<Fig_output_min_excess) Fig_output_min_excess = v;
    Fig_output_modified++;
    return XFIG_RANGE_MIN;
  } else if (v>XFIG_RANGE_MAX) {
    if (v>Fig_output_max_excess) Fig_output_max_excess = v;
    Fig_output_modified++;
    return XFIG_RANGE_MAX;
  }
  return v;
}
#endif

bool istoplevelfiggroup = true;

String Fig_Group::str() {
  String s;
  bool tlfgcache = istoplevelfiggroup;
  if (istoplevelfiggroup) {
    istoplevelfiggroup = false; // prevent enclosed Fig_Group objects from acting as the top level Fig_Group
    if (figattr_Fig_rescale) {
      double w = (double) (bottomrightx - topleftx);
      double w_target = tex_textwidth*1200.0;
      double _rescale = w_target/w;
      if ((_rescale>1.33) || (_rescale<0.67)) {
	rescale(_rescale);
	figscale *= _rescale;
	s += "# rescaled figscale = "+String(figscale,"%.3f\n");
      }
    }
  }
  if (!comment.empty()) {
    if (comment[comment.length()-1]!='\n') comment += '\n';
    comment[comment.length()-1] = '.';
    comment.gsub("\n","\n# ");
    comment[comment.length()-1] = '\n';
    s += "# " + comment;
  }
  s += "6 "+String(topleftx)+' '+String(toplefty)+' '+String(bottomrightx)+' '+String(bottomrighty)+'\n';
  PLL_LOOP_FORWARD(Fig_Object,figobjects.head(),1) s+=e->str();
  s+="-6\n";
  istoplevelfiggroup = tlfgcache;
  return s;
}

void Fig_Group::rescale(double r) {
  topleftx = (long) (((double) topleftx)*r);
  toplefty = (long) (((double) toplefty)*r);
  bottomrightx = (long) (((double) bottomrightx)*r);
  bottomrighty = (long) (((double) bottomrighty)*r);
  PLL_LOOP_FORWARD(Fig_Object,figobjects.head(),1) e->rescale(r);
}

void Fig_Group::Add_Fig_Object(Fig_Object * figobject) {
  if (!figobject) return;
  figobjects.link_before(figobject);
  if (figobject->TopLeftX()<topleftx) topleftx = figobject->TopLeftX();
  if (figobject->TopLeftY()<toplefty) toplefty = figobject->TopLeftY();
  if (figobject->BottomRightX()>bottomrightx) bottomrightx = figobject->BottomRightX();
  if (figobject->BottomRightY()>bottomrighty) bottomrighty = figobject->BottomRightY();
}

String Fig_Line::str() {
  String s("2 1 "+String(linestyle)+' '+String(linewidth)+' '+String(linecolor)+' '+String(fillcolor)+' '+String(depth)+" 0 "+String(fillstyle)+' '+String(dashlengthdotgap,"%.3f")+" 0 0 -1 "+String(leftarrow)+' '+String(rightarrow)+' '+String(pointxysfollow)+"\n\t");
  //if (leftarrow) "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  //if (rightarrow) ; "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  s += String(x1)+' '+String(y1)+' '+String(x2)+' '+String(y2)+'\n';
  return s;
}

void Fig_Line::rescale(double r) {
  x1 = (long) (((double) x1)*r);
  y1 = (long) (((double) y1)*r);
  x2 = (long) (((double) x2)*r);
  y2 = (long) (((double) y2)*r);
}

String Fig_Rectangle::str() {
  String s("2 2 "+String(linestyle)+' '+String(linewidth)+' '+String(linecolor)+' '+String(fillcolor)+' '+String(depth)+" 0 "+String(fillstyle)+' '+String(dashlengthdotgap,"%.3f")+" 0 0 -1 "+String(leftarrow)+' '+String(rightarrow)+' '+String(pointxysfollow)+"\n\t");
  //if (leftarrow) "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  //if (rightarrow) ; "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  s += String(x1)+' '+String(y1)+' '+String(x2)+' '+String(y1)+' '+String(x2)+' '+String(y2)+' '+String(x1)+' '+String(y2)+' '+String(x1)+' '+String(y1)+'\n';
  return s;
}

Fig_Polygon::Fig_Polygon(long ls, long lw, long lc, long fc, long d, long fs, double dldg, long larr, long rarr, const long * x, const long * y, int n): Fig_Line(ls,lw,lc,fc,d,fs,dldg,larr,rarr,LONG_MAX,LONG_MAX,LONG_MIN,LONG_MIN), numpoints(n) {
  pointxysfollow = n+1;
  DIAGNOSTIC_BEFORE_ALLOCATION(new_Fig_element);
  X = new long[n];
  Y = new long[n];
  DIAGNOSTIC_AFTER_ALLOCATION(new_Fig_element,2*n,sizeof(long));
  for (int i=0; i<n; i++) {
    X[i] = XFIG_RANGE_CHECK(x[i]);
    Y[i] = XFIG_RANGE_CHECK(y[i]);
    if (X[i]<x1) x1 = X[i];
    if (X[i]>x2) x2 = X[i];
    if (Y[i]<y1) y1 = Y[i];
    if (Y[i]>y2) y2 = Y[i];
    }
}

String Fig_Polygon::str() {
  String s("2 3 "+String(linestyle)+' '+String(linewidth)+' '+String(linecolor)+' '+String(fillcolor)+' '+String(depth)+" 0 "+String(fillstyle)+' '+String(dashlengthdotgap,"%.3f")+" 0 0 -1 "+String(leftarrow)+' '+String(rightarrow)+' '+String(pointxysfollow)+"\n\t");
  //if (leftarrow) "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  //if (rightarrow) ; "0 0 thickness.00 width.00 height.00\n"; // in steps of 15.00
  for (int i=0; i<numpoints; i++) s += String(X[i])+' '+String(Y[i])+' ';
  s += String(X[0])+' '+String(Y[0])+'\n';
  return s;
}

void Fig_Polygon::rescale(double r) {
  for (int i=0; i<numpoints; i++) {
    X[i] = (long) (((double) X[i])*r);
    Y[i] = (long) (((double) Y[i])*r);
  }
}

void Fig_Circle::rescale(double r) {
  x = (long) (((double) x)*r);
  y = (long) (((double) y)*r);
  radius = (long) (((double) radius)*r);
  centerptx = (long) (((double) centerptx)*r);
  centerpty = (long) (((double) centerpty)*r);
  edgetipx = (long) (((double) edgetipx)*r);
  edgetipy = (long) (((double) edgetipy)*r);
}

long Fig_Text::TopLeftX() {
  switch (justification) {
  case 1: return x - (textlength/2) - 30;
  case 2: return x - textlength - 60;
  default: return x;
  }
}

long Fig_Text::BottomRightX() {
  switch (justification) {
  case 1: return x + (textlength/2) + 30;
  case 2: return x;
  default: return x + textlength + 60;;
  }
}

void Fig_Text::rescale(double r) {
  x = (long) (((double) x)*r);
  y = (long) (((double) y)*r);
  textheight = (long) (((double) textheight)*r);
  textlength = (long) (((double) textlength)*r);
}

String RGBstr(unsigned int coldef) {
  char rgbstr[7];
  for (int i=5; i>=0; i--) {
    rgbstr[i] = '0'+(coldef & 0x00000f);
    if (rgbstr[i]>'9') rgbstr[i] += 'a' - '9' - 1;
    coldef >>= 4;
  }
  rgbstr[6] = '\0';
  return rgbstr;
}
