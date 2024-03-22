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
// Fig_Object.hh
// Randal A. Koene, 20041118
//
// Classes the define Fig_Object and objects that inherit it.

#ifndef __FIG_OBJECT_HH
#define __FIG_OBJECT_HH

#include <limits.h>
#include "BigString.hh"
#include "templates.hh"

extern unsigned long int Fig_output_modified;
extern long Fig_output_min_excess;
extern long Fig_output_max_excess;

extern double tex_textwidth;
extern double figscale;
extern bool figattr_Fig_rescale;

// definitions
void XFIG_RANGE_CHECK_INIT();
String XFIG_RANGE_CHECK_WARNING(String fname, double outwidth);
#ifdef CHECK_XFIG_RANGES
#define XFIG_RANGE_MIN -40000
#define XFIG_RANGE_MAX 40000
long XFIG_RANGE_CHECK(long v);
#else
  #define XFIG_RANGE_CHECK(v) v
#endif

// functions

String RGBstr(unsigned int coldef);

// classes

class Fig_Object: public PLLHandle<Fig_Object> {
public:
  Fig_Object() {}
  virtual String str() { return String(""); }
  virtual long TopLeftX() { return LONG_MAX; }
  virtual long TopLeftY() { return LONG_MAX; }
  virtual long BottomRightX() { return LONG_MIN; }
  virtual long BottomRightY() { return LONG_MIN; }
  virtual void rescale(double r) = 0;
};

class Fig_Header: public Fig_Object {
protected:
  String text;
public:
  Fig_Header(): text("#FIG 3.2\nPortrait\nCenter\nInches\nLetter\n100.00\nSingle\n-2\n1200 2\n") {} // Note: used to be Landscape
  virtual String str() { return text; }
  virtual void rescale(double r) {}
};

class Fig_Color: public Fig_Object {
protected:
  long colnum; // color number (32 and above)
  unsigned int coldef; // RGB color definition, e.g. #FF0000 for red
public:
  Fig_Color(long cn, unsigned int cd): colnum(cn), coldef(cd) {}
  virtual String str() { return String("0 "+String(colnum)+" #"+RGBstr(coldef)+'\n'); }
  virtual void rescale(double r) {}
};

class Fig_Group: public Fig_Object {
protected:
  long topleftx, toplefty, bottomrightx, bottomrighty;
  String comment;
  PLLRoot<Fig_Object> figobjects;
public:
  Fig_Group(long tlx, long tly, long brx, long bry, String c): topleftx(XFIG_RANGE_CHECK(tlx)), toplefty(XFIG_RANGE_CHECK(tly)), bottomrightx(XFIG_RANGE_CHECK(brx)), bottomrighty(XFIG_RANGE_CHECK(bry)), comment(c) {}
  Fig_Group(String c): topleftx(LONG_MAX), toplefty(LONG_MAX), bottomrightx(LONG_MIN), bottomrighty(LONG_MIN), comment(c) {}
  Fig_Group(): topleftx(LONG_MAX), toplefty(LONG_MAX), bottomrightx(LONG_MIN), bottomrighty(LONG_MIN) {}
  virtual String str();
  void Add_Fig_Object(Fig_Object * figobject);
  virtual long TopLeftX() { return topleftx; }
  virtual long TopLeftY() { return toplefty; }
  virtual long BottomRightX() { return bottomrightx; }
  virtual long BottomRightY() { return bottomrighty; }
  virtual void rescale(double r);
};

class Fig_Line: public Fig_Object {
protected:
  long linestyle;
  long linewidth;
  long linecolor;
  long fillcolor;
  long depth;
  long fillstyle;
  double dashlengthdotgap;
  //double angle;
  long leftarrow, rightarrow; // 0 = no, 1 = yes
  //long radius;
  long pointxysfollow;
  long x1, y1;
  long x2, y2;
public:
  Fig_Line(long ls, long lw, long lc, long fc, long d, long fs, double dldg, long larr, long rarr, long xa, long ya, long xb, long yb): linestyle(ls), linewidth(lw), linecolor(lc), fillcolor(fc), depth(d), fillstyle(fs),dashlengthdotgap(dldg), leftarrow(larr), rightarrow(rarr), pointxysfollow(2), x1(XFIG_RANGE_CHECK(xa)), y1(XFIG_RANGE_CHECK(ya)), x2(XFIG_RANGE_CHECK(xb)), y2(XFIG_RANGE_CHECK(yb)) {}
  virtual String str();
  virtual long TopLeftX() { if (x1<x2) return x1; return x2; }
  virtual long TopLeftY() { if (y1<y2) return y1; return y2; }
  virtual long BottomRightX() { if (x1>x2) return x1; return x2; }
  virtual long BottomRightY() { if (y1>y2) return y1; return y2; }
  //void set_leftarrow(double thickness, double width, double height) {}
  //void set_rightarrow() {}
  virtual void rescale(double r);
};

class Fig_Rectangle: public Fig_Line {
public:
  Fig_Rectangle(long ls, long lw, long lc, long fc, long d, long fs, double dldg, long larr, long rarr, long xa, long ya, long xb, long yb): Fig_Line(ls,lw,lc,fc,d,fs,dldg,larr,rarr,xa,ya,xb,yb) { pointxysfollow = 5; }
  virtual String str();
  virtual long TopLeftX() { if (x1<x2) return x1; return x2; }
  virtual long TopLeftY() { if (y1<y2) return y1; return y2; }
  virtual long BottomRightX() { if (x1>x2) return x1; return x2; }
  virtual long BottomRightY() { if (y1>y2) return y1; return y2; }
  //void set_leftarrow(double thickness, double width, double height) {}
  //void set_rightarrow() {}
};

class Fig_Polygon: public Fig_Line {
protected:
  long * X, * Y;
  int numpoints;
public:
  Fig_Polygon(long ls, long lw, long lc, long fc, long d, long fs, double dldg, long larr, long rarr, const long * x, const long * y, int n);
  ~Fig_Polygon() { delete[] X; delete[] Y; }
  virtual String str();
  virtual long TopLeftX() { return x1; }
  virtual long TopLeftY() { return y1; }
  virtual long BottomRightX() { return x2; }
  virtual long BottomRightY() { return y2; }
  //void set_leftarrow(double thickness, double width, double height) {}
  //void set_rightarrow() {}
  virtual void rescale(double r);
};

class Fig_Circle: public Fig_Object {
protected:
  long linestyle;
  long linewidth;
  long linecolor;
  long fillcolor;
  long depth;
  long fillstyle;
  double dashlengthdotgap;
  double angle;
  long x, y;
  long radius;
  long centerptx, centerpty;
  long edgetipx, edgetipy;
public:
  Fig_Circle(long ls, long lw, long lc, long fc, long d, long fs, double dldg, double a, long X, long Y, long r, long cx, long cy, long ex, long ey): linestyle(ls), linewidth(lw), linecolor(lc), fillcolor(fc), depth(d), fillstyle(fs),dashlengthdotgap(dldg), angle(a), x(X), y(Y), radius(r), centerptx(XFIG_RANGE_CHECK(cx)), centerpty(XFIG_RANGE_CHECK(cy)), edgetipx(XFIG_RANGE_CHECK(ex)), edgetipy(XFIG_RANGE_CHECK(ey)) {}
  virtual String str() { return "1 3 "+String(linestyle)+' '+String(linewidth)+' '+String(linecolor)+' '+String(fillcolor)+' '+String(depth)+" 0 "+String(fillstyle)+' '+String(dashlengthdotgap,"%.3f")+" 1 "+String(angle,"%.4f")+' '+String(x)+' '+String(y)+' '+String(radius)+' '+String(radius)+' '+String(centerptx)+' '+String(centerpty)+' '+String(edgetipx)+' '+String(edgetipy)+'\n'; }
  virtual long TopLeftX() { return x-radius; }
  virtual long TopLeftY() { return y-radius; }
  virtual long BottomRightX() { return x+radius; }
  virtual long BottomRightY() { return y+radius; }
  virtual void rescale(double r);
};

class Fig_Text: public Fig_Object {
protected:
  long justification;
  long textcolor;
  long depth;
  long textfont;
  long textsize; // in points
  double textangle;
  long textheight;
  long textlength; // approximated for variable size fonts
  long x, y;
  String text;
public:
  Fig_Text(long tj, long tc, long d, long tf, long ts, double a, long X, long Y, String s): justification(tj), textcolor(tc), depth(d), textfont(tf), textsize(ts), textangle(a), textheight(15*ts), textlength(70*s.length()), x(XFIG_RANGE_CHECK(X)), y(XFIG_RANGE_CHECK(Y)), text(s) {}
  virtual String str() { return "4 "+String(justification)+' '+String(textcolor)+' '+String(depth)+" 0 "+String(textfont)+' '+String(textsize)+' '+String(textangle,"%.4f")+" 4 "+String(textheight)+' '+String(textlength)+' '+String(x)+' '+String(y)+' '+text+"\\001\n"; }
  virtual long TopLeftX();
  virtual long TopLeftY() { return y-((textheight*2)/3); }
  virtual long BottomRightX();
  virtual long BottomRightY() { return y+(textheight/3); }
  virtual void rescale(double r);
};

#endif
