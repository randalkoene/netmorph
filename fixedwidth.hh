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
#ifndef __FIXEDWIDTH_HH
#define __FIXEDWIDTH_HH

template <unsigned r> struct Cb0 {
  static const unsigned Size = Cb0<r >> 1>::Size + 1;
};

template <> struct Cb0<0> {
  static const unsigned Size = 0;
};

struct Cb {
  static const unsigned Size = Cb0<(unsigned char) -1>::Size;
};

template <unsigned n> struct Int0 {
  typedef typename Int0<n + 1>::Type Type;
};

#define Int_Define(T, p) \
template <> struct Int0<sizeof(T) * 16 + p> { \
typedef T Type; \
};

// Lower number has priority for same sized types
Int_Define(signed char, 0)
Int_Define(short, 1)
Int_Define(int, 2)
Int_Define(long, 3)
  //Int_Define(long long, 4)
  
template <unsigned n> struct Int {
  typedef typename Int0<(n + Cb::Size - 1) / Cb::Size * 16>::Type Type;
};

typedef Int<8>::Type int8_t;
typedef Int<16>::Type int16_t;
typedef Int<32>::Type int32_t;
//typedef Int<64>::Type int64_t;

#endif
