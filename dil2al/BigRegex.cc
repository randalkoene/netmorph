/* 
Copyright (C) 1988 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

Modified by Randal A. Koene, 20000228
    disabled #include <rx.h> due to unreliabilities
    hacked problem of duplicate functions due to
    BigString/String and BigRegex/Regex versions
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

/* 
  BigRegex class implementation
 */

#ifdef __GNUG__
#pragma implementation
#endif

#include "BigRegex.hh"

// #define _BIGREGEX_PARANOID_POSIX

#define USE_FASTER
#ifdef USE_REGEX_GNU_ALLOC
	#define TALLOC(n, t) ((t *) malloc ((n) * sizeof (t)))
	#include <stdlib.h>
	#define TFREE(p) free(p)
#else
	#define TALLOC(n, t) new t[n]
	#define TFREE(p) delete[] p
#endif

// #define DEBUG
#ifdef DEBUG
	#include <iostream.h>
	#define DEBUG_OUT(dtxt) cerr << dtxt; cerr.flush();
#else
	#define DEBUG_OUT(dtxt)
#endif

BigRegex::~BigRegex() {
	if (re_assigned) if (*re_assigned) regfree(re);
	delete re;
	delete[] rm;
	delete re_assigned;
}

BigRegex::BigRegex(const char* t, int fast, int bufsize, const char* transtable) {
	re_assigned = new bool;
	if (t==0) *re_assigned = false;	// zero-length RE
	else {
		re = new regex_t;
//*** NOTE: bufsize doesn't really work yet, since it
//*** is not remembered anywhere and up to
//*** BIGREGEX_MAX_SUBEXPRESSIONS are filled in searches below
		 if (bufsize>BIGREGEX_MAX_SUBEXPRESSIONS) bufsize = BIGREGEX_MAX_SUBEXPRESSIONS;
		rm = new regmatch_t[bufsize];
		int res;
#ifdef __DIL2AL__
		// revert to Basic RE backslashing style as inherited from Regex
		// this reverses backslashed and non-backslashed meanings for: (){}|
		char * ttmp;
		int tlen = 0; while (t[tlen]) tlen++;
		ttmp = new char[(2*tlen)+1]; // allocate a character buffer twice as large in case everything has to be escaped with `\'
		int ttmplen = 0; int ii;
		for (int i=0; i<=tlen; i++) {
			ii = i+1;
			if ((t[i]=='\\') && ((t[ii]=='\\') || (t[ii]=='(') || (t[ii]==')') || (t[ii]=='{') || (t[ii]=='}') || (t[ii]=='|'))) {
				if (t[ii]=='\\') { // retain escaped backslash
					ttmp[ttmplen] = '\\';
					ttmplen++;
				}
				i++; // skip backslash
			} else if ((t[i]=='(') || (t[i]==')') || (t[i]=='{') || (t[i]=='}') || (t[i]=='|')) {
				ttmp[ttmplen] = '\\'; // add backslash
				ttmplen++;
			}
			ttmp[ttmplen] = t[i]; // copy character (possibly after skipping or adding a backslash)
			ttmplen++;
		}
		ttmplen--; // don't count '\0'
		res = regcomp(re,ttmp,REG_EXTENDED);
		delete[] ttmp;
#else
		res = regcomp(re,t,REG_EXTENDED);
#endif
		if (!res) *re_assigned = true;
		else {
			*re_assigned = false;
#ifdef _BIGREGEX_HAVE_ERROR_HANDLER
			char errbuf[256];
			if (regerror(res,re,errbuf,256)>0) (*lib_error_handler)("BigRegex",errbuf);
#endif			
		}
	}
}

int BigRegex::match_info(int& start, int& length, int nth) const {
	if ((unsigned)(nth) >= BIGREGEX_MAX_SUBEXPRESSIONS) return 0;
	start = rm[nth].rm_so;
	length = rm[nth].rm_eo - start;
	return start >= 0 && length >= 0;
}

int BigRegex::subpos(int nth) const {
	if ((unsigned)(nth) >= BIGREGEX_MAX_SUBEXPRESSIONS) return -1;
	return rm[nth].rm_so;
}

int BigRegex::sublen(int nth) const {
	if ((unsigned)(nth) >= BIGREGEX_MAX_SUBEXPRESSIONS) return -1;
	return rm[nth].rm_eo - rm[nth].rm_so;
}

#ifndef USE_FASTER

int BigRegex::search(const char* s, int len, int& matchlen, int startpos) const {
// Returns -1 if no match, otherwise index of match
// Negative startpos means scan from end (very slow with regular POSIX.2 regex)
	int res = -1;
	char * stmp;
	bool copied = true;
	int p = startpos;
	if (startpos<0) {
		p += len; // start searching at end
		len = p + 1;	// adjust string to avoid finding matches beyond p when searching backward
		if (p < 0) return -1;
	}
#ifndef _BIGREGEX_SAFE_MATCHES
	if (s[len]=='\0') { // assumes no '\0' prior to s[len]
		copied = false;
		stmp = s;
	} else
#endif
	stmp = new char[len+1];
	int i;
#ifdef _BIGREGEX_SAFE_MATCHES
	for (i=0; ((i<len) && (s[i]!='\0')); i++) stmp[i]=s[i];
#else
	for (i=0; i<len; i++) stmp[i]=s[i];
#endif
	stmp[i]='\0';
	while (p>=0) {
		res = regexec(re,&stmp[p],(size_t) BIGREGEX_MAX_SUBEXPRESSIONS,(regmatch_t *) rm,0);
		if ((startpos>=0) || (res==0)) break;
		p--;
	}
#ifndef _BIGREGEX_SAFE_MATCHES
	if (copied)
#endif
	delete[] stmp;
	if (res==0) {
		if (p>0) { // adjust rm to the offset p
			for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
				rm[i].rm_so += p;
				rm[i].rm_eo += p;
			}
		}
#if defined(_BIGREGEX_SAFE_MATCHES) && defined(_BIGREGEX_HAS_RM_SPEP)
		// change rm pointers back to original string
		for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
			rm[i].rm_sp = s + rm[i].rm_so;
			rm[i].rm_ep = s + rm[i].rm_eo;
		}
#endif
		matchlen = rm[0].rm_eo - rm[0].rm_so;
		return rm[0].rm_so; // start of match
	} else {
		matchlen = 0;
		if (res==REG_NOMATCH) return -1;
		else return -2;	// internal error
	}
}

#else

int BigRegex::search(const char* s, int len, int& matchlen, int startpos) const {
// Returns -1 if no match, otherwise index of match
// Negative startpos means scan from end (very slow with regular POSIX.2 regex)
	int res = -1;
	char * stmp;
#ifndef _BIGREGEX_SAFE_MATCHES
	bool copied = true;
#endif
	int p = startpos;
	if (startpos<0) {
DEBUG_OUT('~')
		len += p + 1;	// start searching at len-abs(startpos)
		if (len <= 0) return -1;
		p = 0;
	}
#ifndef _BIGREGEX_SAFE_MATCHES
	if (s[len]=='\0') { // assumes no '\0' prior to s[len]
		copied = false;
		stmp = s;
	} else
#endif
	stmp = new char[len+1];
	int i;
#ifdef _BIGREGEX_SAFE_MATCHES
	for (i=0; ((i<len) && (s[i]!='\0')); i++) stmp[i]=s[i];
#else
	for (i=0; i<len; i++) stmp[i]=s[i];
#endif
	stmp[i]='\0';
	if (startpos>=0) { // search forward
		res = regexec(re,&stmp[p],(size_t) BIGREGEX_MAX_SUBEXPRESSIONS,(regmatch_t *) rm,0);
	} else { // search backward
#ifdef _BIGREGEX_PARANOID_POSIX
		int tmpp = p, tmpres = 0;
		regmatch_t * tmprm = new regmatch_t[BIGREGEX_MAX_SUBEXPRESSIONS];
		res = REG_NOMATCH;
		while ((tmpp<len) && (tmpres==0)) {
			tmpres = regexec(re,&stmp[tmpp],(size_t) BIGREGEX_MAX_SUBEXPRESSIONS,(regmatch_t *) tmprm,0);
			if (tmpres!=REG_NOMATCH) { // furthest match found so far
				res = tmpres;
				p = tmpp;
				for (i=0; i<BIGREGEX_MAX_SUBEXPRESSIONS; i++) {
					rm[i].rm_so = tmprm[i].rm_so;
					rm[i].rm_eo = tmprm[i].rm_eo;
#if defined(_BIGREGEX_SAFE_MATCHES) && defined(_BIGREGEX_HAS_RM_SPEP)
					rm[i].rm_sp = tmprm[i].rm_sp;
					rm[i].rm_ep = tmprm[i].rm_ep;
#endif
				}
				tmpp=tmprm[0].rm_so+1;	// try next just after current match index
			}
		}
		delete[] tmprm;
#else
DEBUG_OUT('-')
		regex_t private_preg;
		struct re_registers regs;
		private_preg = *re;
		private_preg.not_bol = !!(0 & REG_NOTBOL);
		private_preg.not_eol = !!(0 & REG_NOTEOL);
		private_preg.regs_allocated = REGS_FIXED;
		regs.num_regs = BIGREGEX_MAX_SUBEXPRESSIONS;
		regs.start = TALLOC (BIGREGEX_MAX_SUBEXPRESSIONS, regoff_t);
		regs.end = TALLOC (BIGREGEX_MAX_SUBEXPRESSIONS, regoff_t);
		if (regs.start == NULL || regs.end == NULL) return (int) REG_NOMATCH;
		res = re_search (&private_preg, stmp, len, len-1, -(len-1), &regs);
	     if (res >= 0) {
DEBUG_OUT('!')
			for (unsigned r = 0; r < BIGREGEX_MAX_SUBEXPRESSIONS; r++) {
				rm[r].rm_so = regs.start[r];
				rm[r].rm_eo = regs.end[r];
			}
			res = 0;
		} else res = REG_NOMATCH;
		TFREE (regs.start);
		TFREE (regs.end);
#endif
	}
#ifndef _BIGREGEX_SAFE_MATCHES
	if (copied)
#endif
	delete[] stmp;
	if (res==0) {
		if (p>0) { // adjust rm to the offset p
			for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
				rm[i].rm_so += p;
				rm[i].rm_eo += p;
			}
		}
#if defined(_BIGREGEX_SAFE_MATCHES) && defined(_BIGREGEX_HAS_RM_SPEP)
		// change rm pointers back to original string
		for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
			rm[i].rm_sp = s + rm[i].rm_so;
			rm[i].rm_ep = s + rm[i].rm_eo;
		}
#endif
		matchlen = rm[0].rm_eo - rm[0].rm_so;
		return rm[0].rm_so; // start of match
	} else {
		matchlen = 0;
		if (res==REG_NOMATCH) return -1;
		else return -2;	// internal error
	}
}

#endif

int BigRegex::match(const char*s, int len, int p) const {
	// negative p is considered an offset from len backwards,
	// but note that the match is still scanned forward
	if (p < 0) {
		len += p + 1; // point len at len - abs(p) + 1
		if (len <= 0) return -1;
		p = 0;
	} else if (p > len) return -1;
	// match re to (char*) s, from position p to position len
	// store subexpression information
	// return the number of characters matched
	int res;
#ifndef _BIGREGEX_SAFE_MATCHES
	if (s[len]=='\0') res = regexec(re,&s[p],(size_t) BIGREGEX_MAX_SUBEXPRESSIONS,(regmatch_t *) rm,0);
	else {
#endif
		char * stmp; // copy up to len characters
		stmp = new char[len+1];
		int i;
#ifdef _BIGREGEX_SAFE_MATCHES
		for (i=0; ((i<len) && (s[i]!='\0')); i++) stmp[i]=s[i];
#else
		for (i=0; i<len; i++) stmp[i]=s[i];
#endif
		stmp[i]='\0';
		res = regexec(re,&stmp[p],(size_t) BIGREGEX_MAX_SUBEXPRESSIONS,(regmatch_t *) rm,0);
		delete[] stmp;
#ifndef _BIGREGEX_SAFE_MATCHES
	}
#endif
	if (res==0) {
		if (p>0) { // adjust rm to the offset p
			for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
				rm[i].rm_so += p;
				rm[i].rm_eo += p;
			}
		}
#if defined(_BIGREGEX_SAFE_MATCHES) && defined(_BIGREGEX_HAS_RM_SPEP)
		// change rm pointers back to original string
		for(i=0; ((i<BIGREGEX_MAX_SUBEXPRESSIONS) && (rm[i].rm_so>-1)); i++) {
			rm[i].rm_sp = s + rm[i].rm_so;
			rm[i].rm_ep = s + rm[i].rm_eo;
		}
#endif
		return rm[0].rm_eo - rm[0].rm_so; // number of characters matched
	} else {
		if (res==REG_NOMATCH) return 0;
		else return -2;	// internal error
	}
}

int BigRegex::OK() const {
#ifdef _BIGREGEX_HAVE_ERROR_HANDLER
	if (re_assigned) if (!(*re_assigned)) (*lib_error_handler)("BigRegex", "invariant failure");
#endif
	if (re_assigned) return *re_assigned;
	else return 0;
}

/*
 some built-in Regular expressions
*/

const BigRegex BRXwhite("[ \n\t\r\v\f]+", 1);
const BigRegex BRXint("-?[0-9]+", 1);
const BigRegex BRXdouble("-?\\(\\([0-9]+\\.[0-9]*\\)\\|\\([0-9]+\\)\\|\\(\\.[0-9]+\\)\\)\\([eE][---+]?[0-9]+\\)?", 1, 200);
const BigRegex BRXalpha("[A-Za-z]+", 1);
const BigRegex BRXlowercase("[a-z]+", 1);
const BigRegex BRXuppercase("[A-Z]+", 1);
const BigRegex BRXalphanum("[0-9A-Za-z]+", 1);
const BigRegex BRXidentifier("[A-Za-z_][A-Za-z0-9_]*", 1);
