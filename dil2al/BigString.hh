// This may look like C code, but it is really -*- C++ -*-
/* 
Copyright (C) 1988 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

Modified by Randal A. Koene, 20000224
    allocation types modified to enable big strings
    added ability to find index directly after a substring
    added constructors with istream, long and double
    added sub() function
Modified by Randal A. Koene, 20070622
    updated for new ISO C++ friend function scope rules

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


#ifndef _BigString_hh
#ifdef __GNUG__
#pragma interface
#endif
#define _BigString_hh 1
#define _String_hh 1

#include "BigRegex.hh"
//#include <iostream.h>
#include <iostream>
using namespace std;
#include <stdlib.h>

#undef OK

// Don't use named return values! It is deprecated!
#define _G_NO_NRV

//#define DEBUG_BIGSTRING

#ifdef DEBUG_BIGSTRING
extern bool globaldebugop;
#endif

#define BIGSTRING_MAXIMIZE_SPEED

#define BIGSTRING_ISTREAM_INITIALIZATION
#define BIGSTRING_LONGLONG_INITIALIZATION
#define BIGSTRING_BOOL_INITIALIZATION
#define BIGSTRING_DOUBLE_INITIALIZATION

#ifdef BIGSTRING_DOUBLE_INITIALIZATION
extern int untruncated_double_string_length; // provides the string length that would have been needed to represent a double if it is not truncated by snprintf() (30 or more indicates that the string was truncated)
#endif

#ifdef BIGSTRING_LONGLONG_INITIALIZATION
	#include <stdio.h>
#else
	#ifdef BIGSTRING_DOUBLE_INITIALIZATION
		#include <stdio.h>
	#endif
#endif

/*
#define _BIGSTRING_ALLOC		unsigned int	// was unsigned short
#define _BIGSTRING_LENGTH	unsigned int	// was unsigned short
#define _BIGSTRING_SIZE		int			// was int (-1=unknown)
#define _BIGSTRING_ITERATOR	int			// was int
#define _BIGSTRING_REVERSE	int			// was short
#define _BIGSTRING_ALLOCMANIP	unsigned long	// was unsigned int
*/

typedef unsigned int		_BIGSTRING_ALLOC;		// was unsigned short
typedef unsigned int		_BIGSTRING_LENGTH;		// was unsigned short
typedef int				_BIGSTRING_SIZE;		// was int (-1=unknown)
typedef int				_BIGSTRING_ITERATOR;	// was int
typedef int				_BIGSTRING_REVERSE;		// was short
typedef unsigned long		_BIGSTRING_ALLOCMANIP;	// was unsigned int

typedef char * BigString_char_ptr;

struct StrRep                     // internal String representations
{
  _BIGSTRING_LENGTH		len;    // string length 
  _BIGSTRING_ALLOC		sz;     // allocated space
  char				s[1];   // the string starts here 
                                 // (at least 1 char for trailing null)
                                 // allocated & expanded via non-public fcts
#ifdef DEBUG_BIGSTRING
	bool debugop;				// enables debugging of operations on this String
#endif
};

// primitive ops on StrReps -- nearly all String fns go through these.

StrRep*     Salloc(StrRep* old, const char* src, _BIGSTRING_SIZE srclen, _BIGSTRING_SIZE newlen);
StrRep*     Scopy(StrRep*, const StrRep*);
StrRep*     Scat(StrRep*, const char*, _BIGSTRING_SIZE, const char*, _BIGSTRING_SIZE);
StrRep*     Scat(StrRep*, const char*, _BIGSTRING_SIZE,const char*,_BIGSTRING_SIZE, const char*,_BIGSTRING_SIZE);
StrRep*     Sprepend(StrRep*, const char*, _BIGSTRING_SIZE);
StrRep*     Sreverse(const StrRep*, StrRep*);
StrRep*     Supcase(const StrRep*, StrRep*);
StrRep*     Sdowncase(const StrRep*, StrRep*);
StrRep*     Scapitalize(const StrRep*, StrRep*);
#ifdef BIGSTRING_ISTREAM_INITIALIZATION
StrRep*     Salloc_istream(StrRep* old, istream & src, _BIGSTRING_SIZE srclen, _BIGSTRING_SIZE newlen);
StrRep*     Scat_istream(StrRep*, const char*, _BIGSTRING_SIZE, istream&, _BIGSTRING_SIZE);
#endif

// These classes need to be defined in the order given

class String;
class SubString;

class SubString
{
  friend class      String;
protected:

  String&           S;        // The String I'm a substring of
  _BIGSTRING_LENGTH    pos;      // starting position in S's rep
  _BIGSTRING_LENGTH    len;      // length of substring

  void              assign(const StrRep*, const char*, _BIGSTRING_SIZE = -1);
                    SubString(String& x, _BIGSTRING_SIZE p, _BIGSTRING_SIZE l);
                    SubString(const SubString& x);

public:

// Note there are no public constructors. SubStrings are always
// created via String operations

                   ~SubString();

  SubString&        operator =  (const String&     y);
  SubString&        operator =  (const SubString&  y);
  SubString&        operator =  (const char* t);
  SubString&        operator =  (char        c);

// return 1 if target appears anywhere in SubString; else 0

  int               contains(char        c) const;
  int               contains(const String&     y) const;
  int               contains(const SubString&  y) const;
  int               contains(const char* t) const;
  int               contains(const BigRegex&       r) const;

// return 1 if target matches entire SubString

  int               matches(const BigRegex&  r) const;

// IO 

  friend ostream&   operator<<(ostream& s, const SubString& x);

// status

  _BIGSTRING_LENGTH      length() const;
  int               empty() const;
  const char*       chars() const;

  int               OK() const; 

};


class String
{
  friend class      SubString;

protected:
  StrRep*           rep;   // Strings are pointers to their representations

// some helper functions

  _BIGSTRING_SIZE   search(_BIGSTRING_SIZE, _BIGSTRING_SIZE, const char*, _BIGSTRING_SIZE = -1) const;
  _BIGSTRING_SIZE   search(_BIGSTRING_SIZE, _BIGSTRING_SIZE, char) const;
  _BIGSTRING_SIZE   match(_BIGSTRING_SIZE, _BIGSTRING_SIZE, int, const char*, _BIGSTRING_SIZE = -1) const;
  _BIGSTRING_SIZE   _gsub(const char*, _BIGSTRING_SIZE, const char* ,_BIGSTRING_SIZE);
  _BIGSTRING_SIZE   _gsub(const BigRegex&, const char*, _BIGSTRING_SIZE, char substrs = '\0');
  SubString         _substr(int, int);

public:
#ifdef DEBUG_BIGSTRING
	void debug_on() { rep->debugop=true; }
	void debug_off() { rep->debugop=false; }
	bool debug_state() { return rep->debugop; }
#endif

// constructors & assignment

                    String();
                    String(const String& x);
                    String(const SubString&  x);
                    String(const char* t);
                    String(const char* t, _BIGSTRING_SIZE len);
                    String(char c);
#ifdef BIGSTRING_ISTREAM_INITIALIZATION
				String(istream & ist);
#endif
#ifdef BIGSTRING_LONGLONG_INITIALIZATION
				String(long l);
#endif
#ifdef BIGSTRING_BOOL_INITIALIZATION
				String(bool b);
#endif
#ifdef BIGSTRING_DOUBLE_INITIALIZATION
				String(double d,const char * format);
#endif

                    ~String();

  String&           operator =  (const String&     y);
  String&           operator =  (const char* y);
  String&           operator =  (char        c);
  String&           operator =  (const SubString&  y);
#ifdef BIGSTRING_ISTREAM_INITIALIZATION
  String&           operator =  (istream&  y);
  String&	    readnoseek(istream& y); // useful for istreams such as cin
  String&           readtotag(istream& y, const char * t, char terminator = '\n'); // useful for data input from standard input
#endif

// concatenation

  String&           operator += (const String&     y); 
  String&           operator += (const SubString&  y);
  String&           operator += (const char* y);
  String&           operator += (char        c);
#ifdef BIGSTRING_ISTREAM_INITIALIZATION
  String&           operator += (istream&  y);
#endif

  void              prepend(const String&     y); 
  void              prepend(const SubString&  y);
  void              prepend(const char* t);
  void              prepend(char        c);


// procedural versions:
// concatenate first 2 args, store result in last arg

#ifdef BIGSTRING_ISTREM_INITIALIZE
  friend inline void	cat(const String& x, istream & y, String& r);
#endif
  friend inline void     cat(const String&, const String&, String&);
  friend inline void     cat(const String&, const SubString&, String&);
  friend inline void     cat(const String&, const char*, String&);
  friend inline void     cat(const String&, char, String&);

  friend inline void     cat(const SubString&, const String&, String&);
  friend inline void     cat(const SubString&, const SubString&, String&);
  friend inline void     cat(const SubString&, const char*, String&);
  friend inline void     cat(const SubString&, char, String&);

  friend inline void     cat(const char*, const String&, String&);
  friend inline void     cat(const char*, const SubString&, String&);
  friend inline void     cat(const char*, const char*, String&);
  friend inline void     cat(const char*, char, String&);

// double concatenation, by request. (yes, there are too many versions, 
// but if one is supported, then the others should be too...)
// Concatenate first 3 args, store in last arg

  friend inline void     cat(const String&,const String&, const String&,String&);
  friend inline void     cat(const String&,const String&,const SubString&,String&);
  friend inline void     cat(const String&,const String&, const char*, String&);
  friend inline void     cat(const String&,const String&, char, String&);
  friend inline void     cat(const String&,const SubString&,const String&,String&);
  inline friend void     cat(const String&,const SubString&,const SubString&,String&);
  friend inline void     cat(const String&,const SubString&, const char*, String&);
  friend inline void     cat(const String&,const SubString&, char, String&);
  friend inline void     cat(const String&,const char*, const String&,    String&);
  friend inline void     cat(const String&,const char*, const SubString&, String&);
  friend inline void     cat(const String&,const char*, const char*, String&);
  friend inline void     cat(const String&,const char*, char, String&);

  friend inline void     cat(const char*, const String&, const String&,String&);
  friend inline void     cat(const char*,const String&,const SubString&,String&);
  friend inline void     cat(const char*,const String&, const char*, String&);
  friend inline void     cat(const char*,const String&, char, String&);
  friend inline void     cat(const char*,const SubString&,const String&,String&);
  friend inline void     cat(const char*,const SubString&,const SubString&,String&);
  friend inline void     cat(const char*,const SubString&, const char*, String&);
  friend inline void     cat(const char*,const SubString&, char, String&);
  friend inline void     cat(const char*,const char*, const String&,    String&);
  friend inline void     cat(const char*,const char*, const SubString&, String&);
  friend inline void     cat(const char*,const char*, const char*, String&);
  friend inline void     cat(const char*,const char*, char, String&);


// searching & matching

// search for start of or position directly after match (default = start)

  enum BigString_Search_t { SEARCH_START, SEARCH_END };

// return position of target in string or -1 for failure

  _BIGSTRING_SIZE   index(char	   c, _BIGSTRING_SIZE startpos = 0) const;	 
  _BIGSTRING_SIZE   index(const String&     y, _BIGSTRING_SIZE startpos = 0) const;	  
  _BIGSTRING_SIZE   index(const SubString&  y, _BIGSTRING_SIZE startpos = 0) const;	  
  _BIGSTRING_SIZE   index(const char* t, _BIGSTRING_SIZE startpos = 0) const;  
  _BIGSTRING_SIZE   index(const BigRegex&      r, _BIGSTRING_SIZE startpos = 0) const;	   
  _BIGSTRING_SIZE   index(const BigRegex&      r, BigString_Search_t _st, _BIGSTRING_SIZE startpos = 0) const;	   

// return 1 if target appears anyhere in String; else 0

  int               contains(char        c) const;
  int               contains(const String&     y) const;
  int               contains(const SubString&  y) const;
  int               contains(const char* t) const;
  int               contains(const BigRegex&      r) const;

// return 1 if target appears anywhere after position pos 
// (or before, if pos is negative) in String; else 0

  int               contains(char        c, _BIGSTRING_SIZE pos) const;
  int               contains(const String&     y, _BIGSTRING_SIZE pos) const;
  int               contains(const SubString&  y, _BIGSTRING_SIZE pos) const;
  int               contains(const char* t, _BIGSTRING_SIZE pos) const;
  int               contains(const BigRegex&      r, _BIGSTRING_SIZE pos) const;

// return 1 if target appears at position pos in String; else 0

  int               matches(char        c, _BIGSTRING_SIZE pos = 0) const;
  int               matches(const String&     y, _BIGSTRING_SIZE pos = 0) const;
  int               matches(const SubString&  y, _BIGSTRING_SIZE pos = 0) const;
  int               matches(const char* t, _BIGSTRING_SIZE pos = 0) const;
  int               matches(const BigRegex&      r, _BIGSTRING_SIZE pos = 0) const;
  // BEWARE: An empty string matches anything, since the BigRegex match() function
  // returns the number of characters matched and that is compared with the length
  // of the string to produce the boolean result for String::matches()!

//  return number of occurences of target in String

  _BIGSTRING_SIZE   freq(char        c) const; 
  _BIGSTRING_SIZE   freq(const String&     y) const;
  _BIGSTRING_SIZE   freq(const SubString&  y) const;
  _BIGSTRING_SIZE   freq(const char* t) const;

// SubString extraction

// Note that you can't take a substring of a const String, since
// this leaves open the possiblility of indirectly modifying the
// String through the SubString

  SubString         at(_BIGSTRING_SIZE         pos, _BIGSTRING_SIZE len);
  SubString         operator () (_BIGSTRING_SIZE         pos, _BIGSTRING_SIZE len); // synonym for at

  SubString         at(const String&     x, _BIGSTRING_SIZE startpos = 0); 
  SubString         at(const SubString&  x, _BIGSTRING_SIZE startpos = 0); 
  SubString         at(const char* t, _BIGSTRING_SIZE startpos = 0);
  SubString         at(char        c, _BIGSTRING_SIZE startpos = 0);
  SubString         at(const BigRegex&      r, _BIGSTRING_SIZE startpos = 0); 

  SubString         before(_BIGSTRING_SIZE          pos);
  SubString         before(const String&      x, _BIGSTRING_SIZE startpos = 0);
  SubString         before(const SubString&   x, _BIGSTRING_SIZE startpos = 0);
  SubString         before(const char*  t, _BIGSTRING_SIZE startpos = 0);
  SubString         before(char         c, _BIGSTRING_SIZE startpos = 0);
  SubString         before(const BigRegex&       r, _BIGSTRING_SIZE startpos = 0);

  SubString         through(_BIGSTRING_SIZE          pos);
  SubString         through(const String&      x, _BIGSTRING_SIZE startpos = 0);
  SubString         through(const SubString&   x, _BIGSTRING_SIZE startpos = 0);
  SubString         through(const char*  t, _BIGSTRING_SIZE startpos = 0);
  SubString         through(char         c, _BIGSTRING_SIZE startpos = 0);
  SubString         through(const BigRegex&       r, _BIGSTRING_SIZE startpos = 0);

  SubString         from(_BIGSTRING_SIZE          pos);
  SubString         from(const String&      x, _BIGSTRING_SIZE startpos = 0);
  SubString         from(const SubString&   x, _BIGSTRING_SIZE startpos = 0);
  SubString         from(const char*  t, _BIGSTRING_SIZE startpos = 0);
  SubString         from(char         c, _BIGSTRING_SIZE startpos = 0);
  SubString         from(const BigRegex&       r, _BIGSTRING_SIZE startpos = 0);

  SubString         after(_BIGSTRING_SIZE         pos);
  SubString         after(const String&     x, _BIGSTRING_SIZE startpos = 0);
  SubString         after(const SubString&  x, _BIGSTRING_SIZE startpos = 0);
  SubString         after(const char* t, _BIGSTRING_SIZE startpos = 0);
  SubString         after(char        c, _BIGSTRING_SIZE startpos = 0);
  SubString         after(const BigRegex&      r, _BIGSTRING_SIZE startpos = 0);

// String copies of SubStrings indicated by search with BigRegex

  String            sub(const BigRegex & r,int nth);


// deletion

// delete len chars starting at pos
  void              del(_BIGSTRING_SIZE         pos, _BIGSTRING_SIZE len);

// delete the first occurrence of target after startpos

  void              del(const String&     y, _BIGSTRING_SIZE startpos = 0);
  void              del(const SubString&  y, _BIGSTRING_SIZE startpos = 0);
  void              del(const char* t, _BIGSTRING_SIZE startpos = 0);
  void              del(char        c, _BIGSTRING_SIZE startpos = 0);
  void              del(const BigRegex&      r, _BIGSTRING_SIZE startpos = 0);

// global substitution: substitute all occurrences of pat with repl

  _BIGSTRING_SIZE   gsub(const String&     pat, const String&	 repl);
  _BIGSTRING_SIZE   gsub(const SubString&  pat, const String&	 repl);
  _BIGSTRING_SIZE   gsub(const char* pat, const String&	repl);
  _BIGSTRING_SIZE   gsub(const char* pat, const char* repl);
  _BIGSTRING_SIZE   gsub(const BigRegex&	   pat, const String&	 repl, char substrs = '\0');

// friends & utilities

// split string into array res at separators; return number of elements

  friend _BIGSTRING_SIZE        split(const String& x, String res[], _BIGSTRING_SIZE maxn, 
                          const String& sep);
  friend _BIGSTRING_SIZE        split(const String& x, String res[], _BIGSTRING_SIZE maxn, 
                          const BigRegex&  sep);

  friend String     common_prefix(const String& x, const String& y, 
                                  _BIGSTRING_SIZE startpos = 0);
  friend String     common_suffix(const String& x, const String& y, 
                                  _BIGSTRING_SIZE startpos = -1);
  friend String     replicate(char        c, _BIGSTRING_SIZE n);
  friend String     replicate(const String&     y, _BIGSTRING_SIZE n);
  friend String     join(String src[], _BIGSTRING_SIZE n, const String& sep);

// simple builtin transformations

  friend inline String     reverse(const String& x);
  friend inline String     upcase(const String& x);
  friend inline String     downcase(const String& x);
  friend inline String     capitalize(const String& x);

// in-place versions of above

  void              reverse();
  void              upcase();
  void              downcase();
  void              capitalize();

// element extraction

  char&             operator [] (_BIGSTRING_ITERATOR i);
  const char&       operator [] (_BIGSTRING_ITERATOR i) const;
  char              elem(_BIGSTRING_ITERATOR i) const;
  char              firstchar() const;
  char              lastchar() const;

// conversion

                    operator const char*() const;
  const char*       chars() const;


// IO

  friend inline ostream&   operator<<(ostream& s, const String& x);
  friend ostream&   operator<<(ostream& s, const SubString& x);
  friend istream&   operator>>(istream& s, String& x);

  friend _BIGSTRING_SIZE   readline(istream& s, String& x, 
                             char terminator = '\n',
                             int discard_terminator = 1);

// status

  _BIGSTRING_LENGTH length() const;
  int               empty() const;

// preallocate some space for String
  void              alloc(_BIGSTRING_SIZE newsize);

// report current allocation (not length!)

  _BIGSTRING_ALLOC  allocation() const;


  void     error(const char* msg) const;

  int               OK() const;
};

// Randal A. Koene, update 20070622
// The new ISO C++ standard implemented in GCC 4.1 requires a
// definition/declaration of friend functions outside the class
// body in order to extend their scope beyond the class.

#ifdef BIGSTRING_ISTREM_INITIALIZE
inline void	cat(const String& x, istream & y, String& r);
#endif
inline void     cat(const String&, const String&, String&);
inline void     cat(const String&, const SubString&, String&);
inline void     cat(const String&, const char*, String&);
inline void     cat(const String&, char, String&);

inline void     cat(const SubString&, const String&, String&);
inline void     cat(const SubString&, const SubString&, String&);
inline void     cat(const SubString&, const char*, String&);
inline void     cat(const SubString&, char, String&);

inline void     cat(const char*, const String&, String&);
inline void     cat(const char*, const SubString&, String&);
inline void     cat(const char*, const char*, String&);
inline void     cat(const char*, char, String&);

// double concatenation, by request. (yes, there are too many versions, 
// but if one is supported, then the others should be too...)
// Concatenate first 3 args, store in last arg

inline void     cat(const String&,const String&, const String&,String&);
inline void     cat(const String&,const String&,const SubString&,String&);
inline void     cat(const String&,const String&, const char*, String&);
inline void     cat(const String&,const String&, char, String&);
inline void     cat(const String&,const SubString&,const String&,String&);
void     cat(const String&,const SubString&,const SubString&,String&);
inline void     cat(const String&,const SubString&, const char*, String&);
inline void     cat(const String&,const SubString&, char, String&);
inline void     cat(const String&,const char*, const String&,    String&);
inline void     cat(const String&,const char*, const SubString&, String&);
inline void     cat(const String&,const char*, const char*, String&);
inline void     cat(const String&,const char*, char, String&);

inline void     cat(const char*, const String&, const String&,String&);
inline void     cat(const char*,const String&,const SubString&,String&);
inline void     cat(const char*,const String&, const char*, String&);
inline void     cat(const char*,const String&, char, String&);
inline void     cat(const char*,const SubString&,const String&,String&);
inline void     cat(const char*,const SubString&,const SubString&,String&);
inline void     cat(const char*,const SubString&, const char*, String&);
inline void     cat(const char*,const SubString&, char, String&);
inline void     cat(const char*,const char*, const String&,    String&);
inline void     cat(const char*,const char*, const SubString&, String&);
inline void     cat(const char*,const char*, const char*, String&);
inline void     cat(const char*,const char*, char, String&);

_BIGSTRING_SIZE        split(const String& x, String res[], _BIGSTRING_SIZE maxn, 
			     const String& sep);
_BIGSTRING_SIZE        split(const String& x, String res[], _BIGSTRING_SIZE maxn, 
			     const BigRegex&  sep);

String     common_prefix(const String& x, const String& y, 
			 _BIGSTRING_SIZE startpos);
String     common_suffix(const String& x, const String& y, 
			 _BIGSTRING_SIZE startpos);
String     replicate(char        c, _BIGSTRING_SIZE n);
String replicate(const String&     y, _BIGSTRING_SIZE n);
String     join(String src[], _BIGSTRING_SIZE n, const String& sep);

// simple builtin transformations

inline String     reverse(const String& x);
inline String     upcase(const String& x);
inline String     downcase(const String& x);
inline String     capitalize(const String& x);
inline ostream&   operator<<(ostream& s, const String& x);
ostream&   operator<<(ostream& s, const SubString& x);
istream&   operator>>(istream& s, String& x);

_BIGSTRING_SIZE   readline(istream& s, String& x, 
			   char terminator,
			   int discard_terminator);

typedef String StrTmp; // for backward compatibility

// other externs

_BIGSTRING_SIZE        compare(const String&    x, const String&     y);
_BIGSTRING_SIZE        compare(const String&    x, const SubString&  y);
_BIGSTRING_SIZE        compare(const String&    x, const char* y);
_BIGSTRING_SIZE        compare(const SubString& x, const String&     y);
_BIGSTRING_SIZE        compare(const SubString& x, const SubString&  y);
_BIGSTRING_SIZE        compare(const SubString& x, const char* y);
_BIGSTRING_SIZE        fcompare(const String&   x, const String&     y); // ignore case

extern StrRep  _nilStrRep;
extern String _nilString;

// status reports, needed before defining other things

inline _BIGSTRING_LENGTH String::length() const {  return rep->len; }
inline int         String::empty() const { return rep->len == 0; }
inline const char* String::chars() const { return &(rep->s[0]); }
inline _BIGSTRING_ALLOC         String::allocation() const { return rep->sz; }

inline _BIGSTRING_LENGTH SubString::length() const { return len; }
inline int         SubString::empty() const { return len == 0; }
inline const char* SubString::chars() const { return &(S.rep->s[pos]); }


// constructors

//#define CONSTCHAR_INITIALIZE_TEST

inline String::String() 
  : rep(&_nilStrRep) {}
inline String::String(const String& x) 
  : rep(Scopy(0, x.rep)) {}
inline String::String(const char* t) 
  : rep(Salloc(0, t, -1, -1)) {
#ifdef CONSTCHAR_INITIALIZE_TEST
	cout << "dil2al: Debug CONSTCHAR_INITIALIZE_TEST: ";
	for (char * tp = (char *) t; *tp!='\0'; tp++) cout << (*tp) << '[' << (int) (*tp) << ']';
	cout << '\n';
	cout << '{' << (*this) << "}\n";
#endif
}
inline String::String(const char* t, _BIGSTRING_SIZE tlen)
  : rep(Salloc(0, t, tlen, tlen)) {}
inline String::String(const SubString& y)
  : rep(Salloc(0, y.chars(), y.length(), y.length())) {}
inline String::String(char c) 
  : rep(Salloc(0, &c, 1, 1)) {}
#ifdef BIGSTRING_ISTREAM_INITIALIZATION
inline String::String(istream & ist)
  : rep(Salloc_istream(0,ist,-1,-1)) {}
#endif
#ifdef BIGSTRING_LONGLONG_INITIALIZATION
inline String::String(long l) { char s[12]; snprintf(s,12,"%ld",l); rep = Salloc(0,s,-1,-1); }
#endif
#ifdef BIGSTRING_BOOL_INITIALIZATION
inline String::String(bool b) { if (b) rep = Salloc(0,"true",-1,-1); else rep = Salloc(0,"false",-1,-1); }
#endif
#ifdef BIGSTRING_DOUBLE_INITIALIZATION
inline String::String(double d,const char * format) { char s[30]; untruncated_double_string_length=snprintf(s,30,format,d); rep = Salloc(0,s,-1,-1); }
#endif

inline String::~String() { if (rep != &_nilStrRep) delete rep; }

inline SubString::SubString(const SubString& x)
  :S(x.S), pos(x.pos), len(x.len) {}
inline SubString::SubString(String& x, _BIGSTRING_SIZE first, _BIGSTRING_SIZE l)
  :S(x), pos(first), len(l) {}

inline SubString::~SubString() {}

// assignment

inline String& String::operator =  (const String& y)
{ 
  rep = Scopy(rep, y.rep);
  return *this;
}

inline String& String::operator=(const char* t)
{
  rep = Salloc(rep, t, -1, -1);
  return *this;
}

inline String& String::operator=(const SubString&  y)
{
  rep = Salloc(rep, y.chars(), y.length(), y.length());
  return *this;
}

inline String& String::operator=(char c)
{
  rep = Salloc(rep, &c, 1, 1);
  return *this;
}

#ifdef BIGSTRING_ISTREAM_INITIALIZATION
inline String & String::operator=(istream& y) {
	rep = Salloc_istream(rep, y, -1, -1);
	return *this;
}

inline String & String::readnoseek(istream& y) {
	char buf[4097]; int bufread;
	do {
		y.read(buf,4096);
		bufread = y.gcount();
		if (bufread>0) {
			buf[bufread] = '\0';
			(*this) += buf;
		}
	} while (bufread==4096);
	return *this;
}

inline String & String::readtotag(istream& y, const char * t, char terminator) {
  String buf;
  while (readline(y,buf,terminator,0)>0) {
    if (buf.contains(t)) break;
    else (*this) += buf;
  }
  return *this;
}
#endif

inline SubString& SubString::operator = (const char* ys)
{
  assign(0, ys);
  return *this;
}

inline SubString& SubString::operator = (char ch)
{
  assign(0, &ch, 1);
  return *this;
}

inline SubString& SubString::operator = (const String& y)
{
  assign(y.rep, y.chars(), y.length()); // INVOLVED in lfstr.at(...) = ...
  return *this;
}

inline SubString& SubString::operator = (const SubString& y)
{
  assign(y.S.rep, y.chars(), y.length());
  return *this;
}

// Zillions of cats...

#ifdef BIGSTRING_ISTREAM_INITIALIZE
inline void cat(const String& x, istream & y, String& r)
{
	r.rep = Scat_istream(r.rep, x.chars(), x.length(), y, -1);
}
#endif

inline void cat(const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y, -1);
}

inline void cat(const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), &y, 1);
}

inline void cat(const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y, -1);
}

inline void cat(const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y.chars(), y.length());
}

inline void cat(const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y.chars(), y.length());
}

inline void cat(const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y, -1);
}

inline void cat(const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, x, -1, &y, 1);
}

inline void cat(const String& a, const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y, -1);
}

inline void cat(const String& a, const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), &y, 1);
}

inline void cat(const String& a, const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y, -1);
}

inline void cat(const String& a, const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), &y, 1);
}

inline void cat(const String& a, const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y.chars(), y.length());
}

inline void cat(const String& a, const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y.chars(), y.length());
}

inline void cat(const String& a, const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y, -1);
}

inline void cat(const String& a, const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, &y, 1);
}


inline void cat(const char* a, const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y, -1);
}

inline void cat(const char* a, const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* a, const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y, -1);
}

inline void cat(const char* a, const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* a, const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y.chars(), y.length());
}

inline void cat(const char* a, const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y.chars(), y.length());
}

inline void cat(const char* a, const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y, -1);
}

inline void cat(const char* a, const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, &y, 1);
}

// operator versions

inline String& String::operator +=(const String& y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String::operator +=(const SubString& y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String::operator += (const char* y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String:: operator +=(char y)
{
  cat(*this, y, *this);
  return *this;
}

#ifdef BIGSTRING_ISTREAM_INITIALIZE
inline String& String::operator += (istream&  y) {
	cat(*this, y, *this);
	return *this;
}
#endif

// prepend

inline void String::prepend(const String& y)
{
  rep = Sprepend(rep, y.chars(), y.length());
}

inline void String::prepend(const char* y)
{
  rep = Sprepend(rep, y, -1); 
}

inline void String::prepend(char y)
{
  rep = Sprepend(rep, &y, 1); 
}

inline void String::prepend(const SubString& y)
{
  rep = Sprepend(rep, y.chars(), y.length());
}

// constructive concatenation

#if defined(_USE_NAMED_RETURN_VALUE_EXTENSION) && defined(__GNUG__) && !defined(_G_NO_NRV)

inline String operator + (const String& x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, const char* y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, char y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const char* y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, char y) return r;
{
  cat(x, y, r);
}

inline String operator + (const char* x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (char x, const String& y) return r;
{
  r = y;
  r.prepend(x);
}

inline String operator + (const char* x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String reverse(const String& x) return r;
{
  r.rep = Sreverse(x.rep, r.rep);
}

inline String upcase(const String& x) return r;
{
  r.rep = Supcase(x.rep, r.rep);
}

inline String downcase(const String& x) return r;
{
  r.rep = Sdowncase(x.rep, r.rep);
}

inline String capitalize(const String& x) return r;
{
  r.rep = Scapitalize(x.rep, r.rep);
}

#else /* NO_NRV */

inline String operator + (const String& x, const String& y)
{
  String r;  cat(x, y, r);  return r;
}

inline String operator + (const String& x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const String& x, const char* y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const String& x, char y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const String& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const char* y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, char y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const char* x, const String& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (char x, const String& y)
{
//cout << "<*>" ;
  String r;
  r = y;
//for (char *tp = (char *) r.chars(); *tp!='\0'; tp++) cout << (*tp) << '{' << (int) (*tp) << '}';
//cout << '\n';
  r.prepend(x);
//for (char *tp = (char *) r.chars(); *tp!='\0'; tp++) cout << (*tp) << '{' << (int) (*tp) << '}';
  return r;
}

inline String operator + (const char* x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String reverse(const String& x) 
{
  String r; r.rep = Sreverse(x.rep, r.rep); return r;
}

inline String upcase(const String& x) 
{
  String r; r.rep = Supcase(x.rep, r.rep); return r;
}

inline String downcase(const String& x) 
{
  String r; r.rep = Sdowncase(x.rep, r.rep); return r;
}

inline String capitalize(const String& x) 
{
  String r; r.rep = Scapitalize(x.rep, r.rep); return r;
}

#endif

// misc transformations


inline void String::reverse()
{
  rep = Sreverse(rep, rep);
}


inline void String::upcase()
{
  rep = Supcase(rep, rep);
}


inline void String::downcase()
{
  rep = Sdowncase(rep, rep);
}


inline void String::capitalize()
{
  rep = Scapitalize(rep, rep);
}

// element extraction

inline char&  String::operator [] (_BIGSTRING_ITERATOR i) 
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline const char&  String::operator [] (_BIGSTRING_ITERATOR i) const
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline char  String::elem (_BIGSTRING_ITERATOR i) const
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline char  String::firstchar() const
{ 
  return elem(0);
}

inline char  String::lastchar() const
{ 
  return elem(length() - 1);
}

// searching

inline _BIGSTRING_SIZE String::index(char c, _BIGSTRING_SIZE startpos) const
{
  return search(startpos, length(), c);
}

inline _BIGSTRING_SIZE String::index(const char* t, _BIGSTRING_SIZE startpos) const
{   
  return search(startpos, length(), t);
}

inline _BIGSTRING_SIZE String::index(const String& y, _BIGSTRING_SIZE startpos) const
{   
  return search(startpos, length(), y.chars(), y.length());
}

inline _BIGSTRING_SIZE String::index(const SubString& y, _BIGSTRING_SIZE startpos) const
{   
  return search(startpos, length(), y.chars(), y.length());
}

inline _BIGSTRING_SIZE String::index(const BigRegex& r, _BIGSTRING_SIZE startpos) const
{
  int unused;  return r.search(chars(), length(), unused, startpos);
}

inline _BIGSTRING_SIZE String::index(const BigRegex& r, BigString_Search_t _st, _BIGSTRING_SIZE startpos) const
{
  int matchlen;
  int matchstart = r.search(chars(), length(), matchlen, startpos);
  if ((matchstart<0) || (_st==SEARCH_START)) return matchstart;
  else return (matchstart+matchlen);
}

inline int String::contains(char c) const
{
  return search(0, length(), c) >= 0;
}

inline int String::contains(const char* t) const
{   
  return search(0, length(), t) >= 0;
}

inline int String::contains(const String& y) const
{   
  return search(0, length(), y.chars(), y.length()) >= 0;
}

inline int String::contains(const SubString& y) const
{   
  return search(0, length(), y.chars(), y.length()) >= 0;
}

inline int String::contains(char c, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 0, &c, 1) >= 0;
}

inline int String::contains(const char* t, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 0, t) >= 0;
}

inline int String::contains(const String& y, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 0, y.chars(), y.length()) >= 0;
}

inline int String::contains(const SubString& y, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 0, y.chars(), y.length()) >= 0;
}

inline int String::contains(const BigRegex& r) const
{
  int unused;  return r.search(chars(), length(), unused, 0) >= 0;
}

inline int String::contains(const BigRegex& r, _BIGSTRING_SIZE p) const
{
  return r.match(chars(), length(), p) >= 0;
}


inline int String::matches(const SubString& y, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 1, y.chars(), y.length()) >= 0;
}

inline int String::matches(const String& y, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 1, y.chars(), y.length()) >= 0;
}

inline int String::matches(const char* t, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 1, t) >= 0;
}

inline int String::matches(char c, _BIGSTRING_SIZE p) const
{
  return match(p, length(), 1, &c, 1) >= 0;
}

inline int String::matches(const BigRegex& r, _BIGSTRING_SIZE p) const
{
  _BIGSTRING_SIZE l = (p < 0)? -p : length() - p;
  return r.match(chars(), length(), p) == l;
}


inline int SubString::contains(const char* t) const
{   
  return S.search(pos, pos+len, t) >= 0;
}

inline int SubString::contains(const String& y) const
{   
  return S.search(pos, pos+len, y.chars(), y.length()) >= 0;
}

inline int SubString::contains(const SubString&  y) const
{   
  return S.search(pos, pos+len, y.chars(), y.length()) >= 0;
}

inline int SubString::contains(char c) const
{
  return S.search(pos, pos+len, c) >= 0;
}

inline int SubString::contains(const BigRegex& r) const
{
  int unused;  return r.search(chars(), len, unused, 0) >= 0;
}

inline int SubString::matches(const BigRegex& r) const
{
  // Note the comparison between int and converted unsigned int.
  // r.match can return -1 or -2 (internal error), which might in
  // rare cases be equal to a very large SubString len value!
  return r.match(chars(), len, 0) == (int) len;
}

inline _BIGSTRING_SIZE String::gsub(const String& pat, const String& r)
{
  return _gsub(pat.chars(), pat.length(), r.chars(), r.length());
}

inline _BIGSTRING_SIZE String::gsub(const SubString&  pat, const String& r)
{
  return _gsub(pat.chars(), pat.length(), r.chars(), r.length());
}

inline _BIGSTRING_SIZE String::gsub(const BigRegex& pat, const String& r, char substrs)
{
  return _gsub(pat, r.chars(), r.length(), substrs);
}

inline _BIGSTRING_SIZE String::gsub(const char* pat, const String& r)
{
  return _gsub(pat, -1, r.chars(), r.length());
}

inline _BIGSTRING_SIZE String::gsub(const char* pat, const char* r)
{
  return _gsub(pat, -1, r, -1);
}



inline  ostream& operator<<(ostream& s, const String& x)
{
   s << x.chars(); return s;
}

// a zillion comparison operators

inline int operator==(const String& x, const String& y) 
{
  return compare(x, y) == 0; 
}

inline int operator!=(const String& x, const String& y)
{
  return compare(x, y) != 0; 
}

inline int operator>(const String& x, const String& y)
{
  return compare(x, y) > 0; 
}

inline int operator>=(const String& x, const String& y)
{
  return compare(x, y) >= 0; 
}

inline int operator<(const String& x, const String& y)
{
  return compare(x, y) < 0; 
}

inline int operator<=(const String& x, const String& y)
{
  return compare(x, y) <= 0; 
}

inline int operator==(const String& x, const SubString&  y) 
{
  return compare(x, y) == 0; 
}

inline int operator!=(const String& x, const SubString&  y)
{
  return compare(x, y) != 0; 
}

inline int operator>(const String& x, const SubString&  y)      
{
  return compare(x, y) > 0; 
}

inline int operator>=(const String& x, const SubString&  y)
{
  return compare(x, y) >= 0; 
}

inline int operator<(const String& x, const SubString&  y) 
{
  return compare(x, y) < 0; 
}

inline int operator<=(const String& x, const SubString&  y)
{
  return compare(x, y) <= 0; 
}

inline int operator==(const String& x, const char* t) 
{
  return compare(x, t) == 0; 
}

inline int operator!=(const String& x, const char* t) 
{
  return compare(x, t) != 0; 
}

inline int operator>(const String& x, const char* t)  
{
  return compare(x, t) > 0; 
}

inline int operator>=(const String& x, const char* t) 
{
  return compare(x, t) >= 0; 
}

inline int operator<(const String& x, const char* t)  
{
  return compare(x, t) < 0; 
}

inline int operator<=(const String& x, const char* t) 
{
  return compare(x, t) <= 0; 
}

inline int operator==(const SubString& x, const String& y) 
{
  return compare(y, x) == 0; 
}

inline int operator!=(const SubString& x, const String& y)
{
  return compare(y, x) != 0;
}

inline int operator>(const SubString& x, const String& y)      
{
  return compare(y, x) < 0;
}

inline int operator>=(const SubString& x, const String& y)     
{
  return compare(y, x) <= 0;
}

inline int operator<(const SubString& x, const String& y)      
{
  return compare(y, x) > 0;
}

inline int operator<=(const SubString& x, const String& y)     
{
  return compare(y, x) >= 0;
}

inline int operator==(const SubString& x, const SubString&  y) 
{
  return compare(x, y) == 0; 
}

inline int operator!=(const SubString& x, const SubString&  y)
{
  return compare(x, y) != 0;
}

inline int operator>(const SubString& x, const SubString&  y)      
{
  return compare(x, y) > 0;
}

inline int operator>=(const SubString& x, const SubString&  y)
{
  return compare(x, y) >= 0;
}

inline int operator<(const SubString& x, const SubString&  y) 
{
  return compare(x, y) < 0;
}

inline int operator<=(const SubString& x, const SubString&  y)
{
  return compare(x, y) <= 0;
}

inline int operator==(const SubString& x, const char* t) 
{
  return compare(x, t) == 0; 
}

inline int operator!=(const SubString& x, const char* t) 
{
  return compare(x, t) != 0;
}

inline int operator>(const SubString& x, const char* t)  
{
  return compare(x, t) > 0; 
}

inline int operator>=(const SubString& x, const char* t) 
{
  return compare(x, t) >= 0; 
}

inline int operator<(const SubString& x, const char* t)  
{
  return compare(x, t) < 0; 
}

inline int operator<=(const SubString& x, const char* t) 
{
  return compare(x, t) <= 0; 
}


// a helper needed by at, before, etc.

inline SubString String::_substr(_BIGSTRING_SIZE first, _BIGSTRING_SIZE l)
{
  if (first < 0 || (unsigned)(first + l) > length() )
    return SubString(_nilString, 0, 0) ;
  else 
    return SubString(*this, first, l);
}

#endif
