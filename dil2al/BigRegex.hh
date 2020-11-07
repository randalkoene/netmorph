// This may look like C code, but it is really -*- C++ -*-
/* 
Copyright (C) 1988 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

Modified by Randal A. Koene, 20000228
    disabled #include <rx.h> due to unreliabilities
    POSIX.2 adaptation of Regex (!)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/


#ifndef _BigRegex_hh
#ifdef __GNUG__
#pragma interface
#endif
#define _BigRegex_hh 1
// block other possible redefinitions
#define _Regex_h 1

#undef OK

#include <sys/types.h>

/*
	Note: The rx.h supplied with Linux or with the libg++ library
	source code is fairly unreliable. To avoid these difficulties
	as much as possible, BigRegex adapts Regex to use POSIX.2
	regex functions.
*/

// if <regex.h> below does not point to a POSIX.2 regex implementation such
// as that of the GNU libc library, then obtain a POSIX.2 implementation and
// set _ALT_REGEX as well as _USE_ALT_REGEX
//#define _ALT_REGEX_H "regex.h"
//#define _USE_ALT_REGEX
#ifndef _CPP_REGEX
extern "C" {
#endif
#ifdef _USE_ALT_REGEX
	#include _ALT_REGEX_H
#else
	#ifdef SYSTEM_RX
		#include <regex.h>
	#else
		#include "rx.h"
	#endif
#endif
#ifndef _CPP_REGEX
}
#endif
// This slows things down a little, but guarantees that matches do not
// proceed beyond the length of a string
//#ifdef _BIGREGEX_SAFE_MATCHES
//	#include <string.h>
//#endif

//#include <iostream.h>
#include <iostream>

// Currently don't have error handler for this class yet
//#define _BIGREGEX_HAVE_ERROR_HANDLER

// The structure type regmatch_t is not very large, so it is quite possible
// to define a greater number of subexpression registers.
// Note that if the numberf of subexpressions exceeds the number of registers,
// searches are still correctly performed, but none of the registers are filled.
#define BIGREGEX_MAX_SUBEXPRESSIONS	16

class BigRegex
{
private:

                     BigRegex(const BigRegex&) { *re_assigned=false; }  // no X(X&)
  void               operator = (const BigRegex&) { *re_assigned=false; } // no assignment

protected:
	regex_t	*	re;
	regmatch_t *	rm;
	bool	*		re_assigned;

public:
                     BigRegex(const char* t, 
                           int fast = 0, 
                           int bufsize = 40, 
                           const char* transtable = 0);

                    ~BigRegex();

  int                match(const char* s, int len, int pos = 0) const;
  int                search(const char* s, int len, 
                            int& matchlen, int startpos = 0) const;
  int                match_info(int& start, int& length, int nth = 0) const;
  
  int                subpos(int nth = 0) const;
  int                sublen(int nth = 0) const;

  int                OK() const;  // representation invariant
};

// some built in regular expressions

extern const BigRegex BRXwhite;          // = "[ \n\t\r\v\f]+"
extern const BigRegex BRXint;            // = "-?[0-9]+"
extern const BigRegex BRXdouble;         // = "-?\\(\\([0-9]+\\.[0-9]*\\)\\|
                                     //    \\([0-9]+\\)\\|\\(\\.[0-9]+\\)\\)
                                     //    \\([eE][---+]?[0-9]+\\)?"
extern const BigRegex BRXalpha;          // = "[A-Za-z]+"
extern const BigRegex BRXlowercase;      // = "[a-z]+"
extern const BigRegex BRXuppercase;      // = "[A-Z]+"
extern const BigRegex BRXalphanum;       // = "[0-9A-Za-z]+"
extern const BigRegex BRXidentifier;     // = "[A-Za-z_][A-Za-z0-9_]*"

#endif
