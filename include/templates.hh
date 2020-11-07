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
// Useful C++ templates
//
// Randal A. Koene, 20000802
//
// Note: This version of templates.hh implements the version option
// from solution 2 of the file ~/src/include/cast-const.txt.
// Solution 1 is implemented in the file templates.classic.hh.
// The file templates.unfixed.hh is the implementation that does
// not compile properly with gcc 2.96 or greater.

#ifndef __TEMPLATES_HH
#define __TEMPLATES_HH

// Protected Linked List
//
// A linked list with protection against breakage and floating remnants,
// which is also type-safe by using define to create type-specific
// code instead of type-casting from generic code.
// The two classes that are created for each class that uses
// this Protected Linked List template are declared as friends
// of each other, to simplify intuitive processes for such tasks
// as linking and unlinking, without granting any private access to
// either the class using the Protected Linked List template or
// any classes derived from the two additional created classes,
// thereby effectively guarding the link structure against
// changes that could violate the protocol.
/* Schematic of the structure:

   +----+     +----+     +----+     +----+
   |obj1|     |obj2|     |obj3|     |obj4|
   |inh.|     |inh.|     |inh.|     |inh.|
   |PLL |     |PLL |     |PLL |     |PLL |
   +----+ next+----+     +----+     +----+
+->|PLL1|<--->|PLL2|<--->|PLL3|<--->|PLL4|<-+
|  +----+prev +----+     +----+     +----+  |
|     \          \         /          /     |   PLL and
|      \__________\_______/__________/      |   PLLRoot
|           root      |                     |   are
|                     v                     |   friends
|         head  +-----------+  tail         |
+---------------|  PLLRoot  |---------------+
                +-----------+
                     ### (local, global or new-allocated)
        +----------------------------+
        | class, function or program |
        +----------------------------+

*/

// a - Protected Linked List class and its Root class, add these
//     above the class header

#ifndef NULL
#define NULL 0
#endif

enum PLLHandle_Signal { PLL_NOSIGNAL, PLL_LIST_DESTRUCT, PLL_MULTI_LIST_DELEGATOR_DESTRUCTING };

template <class PLLType> class PLLRoot;
template <class PLLType>
// <A NAME="PLLHandle"></A>
class PLLHandle { 
  friend class PLLRoot<PLLType>;
  //template<class PLLType> friend PLLRoot; 
protected: 
  PLLRoot<PLLType> * pllroot; 
  PLLType * pllnext, * pllprev; 
  PLLHandle_Signal pllsignal; // most recent or ongoing operation
  PLLHandle(): pllroot(NULL), pllnext(NULL), pllprev(NULL), pllsignal(PLL_NOSIGNAL) {}
  virtual ~PLLHandle() { // virtual destructor calls actual complete most-derived object's destructor first
    //cerr << "PLLH=" << (int) this << " next=" << (int) pllnext << " s=" << (int) pllsignal << '\n'; cerr.flush();
    if (pllsignal==PLL_LIST_DESTRUCT) { // If pllroot calls deletion, delete from here onwards.
      if (pllnext) {
	pllnext->pllsignal=PLL_LIST_DESTRUCT;
	delete pllnext;
      }
    } else {
      // If un unexpected source calls deletion (e.g. if this element
      // was a local object that is destructed) then unlink and
      // delete.
      isolate();
    }
  }
  inline void isolate(); 
public: 
  PLLRoot<PLLType> * Root() const { return pllroot; } 
  PLLType * Next() const { return pllnext; } 
  PLLType * Prev() const { return pllprev; } 
  PLLType * head() const { return pllroot->head(); } 
  PLLType * tail() const { return pllroot->tail(); }
  PLLHandle_Signal signal() const { return pllsignal; }
  inline PLLType * el(int n) const; 
  inline PLLType * seek_head() const; 
  inline PLLType * seek_tail() const; 
  inline int length() const; 
  inline int fulllength() const; 
  inline bool unlink_before(PLLRoot<PLLType> & succ); 
  inline bool unlink_after(PLLRoot<PLLType> & pred); 
  inline void remove(); 
};

template <class PLLType>
class PLLRoot { 
	friend class PLLHandle<PLLType>; 
	//template<class PLLType> friend PLLHandle; 
	protected: 
		PLLType * pllhead; 
		PLLType * plltail; 
		inline void set_root(PLLType * h); 
	public: 
/* *** Deletion issue:
    Should perhaps enable deletion of a root without deletion of the list, if
    it is possible to have multiple roots pointing at the list.
*/
		~PLLRoot() {
		  //cerr << "PLLR=" << (int) this << " head=" << (int) pllhead << '\n'; cerr.flush();
		clear(); } // *** perhaps also make this virtual, but delete is not usually called on a PLLRoot<PLLType> object
		PLLRoot(): pllhead(NULL), plltail(NULL) {}
		PLLType * head() const { return pllhead; } 
		PLLType * tail() const { return plltail; } 
		PLLType * el(int n) const { if (pllhead) return pllhead->el(n); return NULL; } 
		PLLType * operator[](int n) const { return el(n); } 
		int length() const { if (pllhead) return pllhead->length(); return 0; } 
		inline void update(); 
		inline bool link_before(PLLRoot * succ); 
		inline bool link_before(PLLType * succ); 
		inline bool link_after(PLLRoot * pred); 
		inline bool link_after(PLLType * pred); 
		inline bool insert_after(PLLType * pred, PLLType * inshead); 
		inline bool insert_after(int n, PLLType * inshead); 
		inline bool insert_before(PLLType * succ, PLLType * inshead); 
		inline bool insert_before(int n, PLLType * inshead);
		// *** if the clear() function is for some reason too
		//     drastic, it can be rewritten as follows to require
		//     subsequent manual deletion:
		// PLLRoot<PLLType> * clear() { PLLRoot<PLLType> * clearroot = new PLLRoot<PLLType>; if (pllhead) pllhead->unlink_before(clearroot); return clearroot; }
		void clear() { if (pllhead) pllhead->pllsignal = PLL_LIST_DESTRUCT; delete pllhead; pllhead = NULL; plltail = NULL; }
};

// b - Protected Linked List inheritance, add this into the list
//     of inherited classes: public PLLHandle<PLLType>

// c - Protected Linked List inline functions
//     Note that the static_cast<>()s used below are safe, since
//     it is assured through the creation mechanisms of PLLHandle<PLLType>
//     and PLLRoot<PLLType> that all elements are of the derived class
//     PLLType with base PLLHandle<PLLType>. Also see corresponding
//     <A HREF="../../nnmodels/nnmodels.html">notes</A>.
//     Addendum: The use of static_cast<>() is now avoided as much as
//     possible, since it should not be necessary if the code is
//     properly programmed.

template <class PLLType>
inline bool PLLHandle<PLLType>::unlink_before(PLLRoot<PLLType> & succ) { 
// unlinks before this and gives this and after
// the new root succ
// succ is currently required to be an empty list
	if (succ.head()) return false; 
	if (pllprev) pllprev->pllnext = NULL; 
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	succ.pllhead = (PLLType *) this; 
	if (pllroot) { 
		succ.plltail = pllroot->tail(); 
		pllroot->plltail = pllprev; 
		if (!pllprev) pllroot->pllhead = NULL; 
	} else succ.plltail = seek_tail(); 
	pllprev = NULL; 
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	succ.set_root((PLLType *) this); 
	return true; 
} 

template <class PLLType>
inline bool PLLHandle<PLLType>::unlink_after(PLLRoot<PLLType> & pred) { 
// unlinks after this and gives this and before
// the new root pred
// pred is currently required to be an empty list
	if (pred.head()) return false; 
	if (pllnext) pllnext->pllprev = NULL; 
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	pred.plltail = (PLLType *) this; 
	if (pllroot) { 
		pred.pllhead = pllroot->head(); 
		pllroot->pllhead = pllnext; 
		if (!pllnext) pllroot->plltail = NULL; 
	} else pred.pllhead = seek_head(); 
	pllnext = NULL; 
	pred.set_root(pred.head()); 
	return true; 
} 

template <class PLLType>
inline void PLLHandle<PLLType>::isolate() { 
// unlink Protected Linked List element and update head and tail
	if (!pllroot) { 
		if (pllprev) pllprev->pllnext = pllnext; 
		if (pllnext) pllnext->pllprev = pllprev;
		pllprev = pllnext = NULL;
		return;
	} 
	PLLRoot<PLLType> * h = pllroot, r, t; 
	unlink_before(t); // t becomes root of this and after
	unlink_after(r); // r becomes root of this
	h->link_before(&t); // append t to h
	// avoid deletion as r is destructed
	r.pllhead = NULL; r.plltail = NULL; pllroot = NULL;
} 

template <class PLLType>
inline void PLLHandle<PLLType>::remove() { 
// unlink and delete Protected Linked List element
// and update head and tail
	isolate();
	pllsignal = PLL_LIST_DESTRUCT; // an orderly destruct
	delete this;
} 

template <class PLLType>
inline PLLType * PLLHandle<PLLType>::el(int n) const {
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	if (!n) return (PLLType *) this;
	if (n>0) { 
		if (pllnext) return pllnext->el(n-1); 
		else return NULL; 
	} else { 
		if (pllprev) return pllprev->el(n+1); 
		else return NULL; 
	} 
} 

template <class PLLType>
inline PLLType * PLLHandle<PLLType>::seek_head() const { 
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	if (!pllprev) return (PLLType *) this;
	return pllprev->seek_head(); 
} 

template <class PLLType>
inline PLLType * PLLHandle<PLLType>::seek_tail() const { 
 /* *** At this location the simple solution of PLLType * pllnext, etc. does not work.
        Using this method, it was necessary to resort to old-style C casting. */
	if (!pllnext) return (PLLType *) this; 
	return pllnext->seek_tail(); 
} 

template <class PLLType>
inline int PLLHandle<PLLType>::length() const { 
	if (Next()) return 1 + Next()->length(); 
	return 1; 
} 

template <class PLLType>
inline int PLLHandle<PLLType>::fulllength() const { 
	return head()->length(); 
} 

template <class PLLType>
inline void PLLRoot<PLLType>::set_root(PLLType * h) { 
	for (; (h); h = h->Next()) h->pllroot = this; 
} 

template <class PLLType>
inline void PLLRoot<PLLType>::update() { 
// clears head and tail if elements were linked to another root
// or fixes head and tail if either one pointed to an element with
// a different root
	if (pllhead) if (pllhead->Root()!=this) pllhead = NULL; 
	if (plltail) { 
		if (plltail->Root()!=this) { 
			if (!pllhead) plltail = NULL; 
			else plltail = pllhead->seek_tail(); 
		} else { 
			if (!pllhead) pllhead = plltail->seek_head(); 
		} 
	} 
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::link_before(PLLRoot<PLLType> * succ) {
// Links the tail of this to the head of succ, then makes this the
// root of elements previously belonging to succ
	if (!succ) return false; 
	if (!succ->head()) return false; 
	// *** could call update() here first to insure that pllhead and plltail
	//     are correct
	if (plltail) { 
		plltail->pllnext = succ->head(); 
		succ->head()->pllprev = plltail; 
	} else pllhead = succ->head(); 
	plltail = succ->tail(); 
	set_root(succ->head()); 
	succ->update(); 
	return true; 
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::link_before(PLLType * succ) { 
// when succ->Root()==NULL, this allows previously
// unprotected lists or newly created PLLType elements
// to be linked
	if (!succ) return false; 
	if (succ->Root()) return link_before(succ->Root()); 
	if (plltail) { 
		plltail->pllnext = succ;
		succ->pllprev = plltail; 
	} else pllhead = succ; 
	plltail = succ; // *** this may need to be succ->seek_tail(), unless lists without root cannot exist
	set_root(succ); 
	return true; 
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::link_after(PLLRoot<PLLType> * pred) { 
	if (!pred) return false; 
	if (!pred->tail()) return false; 
	// *** could call update() here first to insure that pllhead and plltail
	//     are correct
	if (pllhead) { 
		pllhead->pllprev = pred->tail(); 
		pred->tail()->pllnext = pllhead; 
	} else plltail = pred->tail(); 
	pllhead = pred->head(); 
	set_root(pred->head()); 
	pred->update(); 
	return true; 
}

template <class PLLType>
inline bool PLLRoot<PLLType>::link_after(PLLType * pred) { 
// when succ->Root()==NULL, this allows previously
// unprotected lists or newly created PLLType elements
// to be linked
	if (!pred) return false; 
	if (pred->Root()) return link_after(pred->Root()); 
	if (pllhead) { 
		pllhead->pllprev = pred; 
		pred->pllnext = pllhead; 
	} else plltail = pred; 
	pllhead = pred;  // *** this may need to be pred->seek_head(), unless lists without root cannot exist
	set_root(pred); 
	return true; 
}

template <class PLLType>
inline bool PLLRoot<PLLType>::insert_after(PLLType * pred, PLLType * inshead) {
// inserts inshead after element pred in this list
// if pred==NULL or pred->Root()!=this, inserts at tail of list
	if (!pred) return link_before(inshead); // link inshead to tail of this
	if (pred->Root()!=this) return link_before(inshead); // link inshead to tail of this
	PLLRoot<PLLType> h;
	pred->unlink_after(h); // get head of list up to pred into h
	bool res = h.link_before(inshead); // link inshead to tail of h
	link_after(&h); // link h to head of this
	return res;
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::insert_after(int n, PLLType * inshead) {
// inserts inshead after the nth element in this list
// if nth element does not exist, inserts at tail of list
	PLLType * pred = el(n);
	return insert_after(pred,inshead);
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::insert_before(PLLType * succ, PLLType * inshead) {
// inserts inshead before element succ in this list
// if succ==NULL or succ->Root()!=this, inserts at head of list
	if (!succ) return link_after(inshead); // link inshead to head of this
	if (succ->Root()!=this) return link_after(inshead); // link inshead to head of this
	PLLRoot<PLLType> t;
	succ->unlink_before(t); // get tail of list from succ into t
	bool res = t.link_after(inshead); // link inshead to head of t
	link_before(&t); // link t to tail of this
	return res;
} 

template <class PLLType>
inline bool PLLRoot<PLLType>::insert_before(int n, PLLType * inshead) {
// inserts inshead before the nth element in this list
// if nth element does not exist, inserts at head of list
	PLLType * succ = el(n);
	return insert_before(succ,inshead);
} 

// d - Protected Linked List macro that can be used to
//     access functions in linked elements at a relative
//     distance n from an element
#define Protected_Linked_List_Indexed_Access(PLLType,Here,There) \
	if (!n) return Here; \
	if (PLLType * e = el(n)) return e->There; else return NULL \

// e - Protected Linked List macro that loops forward through
//     linked elements contingent upon a test, providing the
//     pointer e to the current element
#define PLL_LOOP_FORWARD(PLLType,Init_El,Extra_Test) \
	for (PLLType * e = Init_El; ((e) && (Extra_Test)); e = e->Next())

// f - Protected Linked List macro that loops forward through
//     linked elements contingent upon a test, providing the
//     pointer to the current element in a specified variable
//     to enable nesting of PLL_LOOP_FORWARD_NESTED calls
#define PLL_LOOP_FORWARD_NESTED(PLLType,Init_El,Extra_Test,RefVar) \
	for (PLLType * RefVar = Init_El; ((RefVar) && (Extra_Test)); RefVar = RefVar->Next())

#endif
