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
// diagnostic.hh
// Randal A. Koene, 20050309
//
// Classes and definitions that aid in self-diagnostic tests for
// profiling and debugging.

#ifndef __DIAGNOSTIC_HH
#define __DIAGNOSTIC_HH

#ifdef EVALUATE_POSSIBLE_CONNECTION_PROFILING
extern long evaluate_possible_connection_calls;
#define DIAGNOSTIC_EVALUATE_POSSIBLE_CONNECTION_PROFILING evaluate_possible_connection_calls++
#else
#define DIAGNOSTIC_EVALUATE_POSSIBLE_CONNECTION_PROFILING
#endif

#ifndef MEMTESTS
#define DIAGNOSTIC_BEFORE_ALLOCATION(new_type)
#define DIAGNOSTIC_AFTER_ALLOCATION(new_type,n,sz)
#define DIAGNOSTIC_TEMP_VARIABLES(var_type,var_name,var_value)
#endif

#ifdef MEMTESTS

enum allocation_event {initial_program_size,new_fibre_segment,new_Fig_element,new_Spatial_Segment_Subset,new_synapse_structure};

// classess that collect statistics and store them in diagnostic files

class memory_allocation_diagnostic {
public:
  long num_allocations;
  long last_intended;
  long last_actual;
  long total_intended;
  long total_actual;
  memory_allocation_diagnostic(): num_allocations(0), last_intended(0), last_actual(0), total_intended(0), total_actual(0) {}
  void allocating(long intended, long actual);
};

class memory_statistics {
protected:
  void initialize();
public:
  long events;
  long t_ujiffies;
  long t_sjiffies;
  long actualused_before;
  long actualused_after;
  long actualused_intervening_num;
  long actualused_intervening;
  long allocated_program;
  long allocated_total;
  memory_allocation_diagnostic allocated_fibre_segment;
  memory_allocation_diagnostic allocated_Fig_element; 
  memory_allocation_diagnostic allocated_Spatial_Segment_Subset; 
  memory_allocation_diagnostic allocated_synapse_structure; 
  memory_statistics(): events(0), t_ujiffies(0), t_sjiffies(0), actualused_before(0), actualused_after(0), actualused_intervening_num(0), actualused_intervening(0), allocated_program(0), allocated_total(0) { initialize(); }
  ~memory_statistics();
  void before_allocation(allocation_event new_type);
  void after_allocation(allocation_event new_type, int n, long sz);
};

extern memory_statistics memstats;

// definitions that frame memory allocation and call functions that
// collect statistics

#define DIAGNOSTIC_BEFORE_ALLOCATION(new_type) memstats.before_allocation(new_type)
#define DIAGNOSTIC_AFTER_ALLOCATION(new_type,n,sz) memstats.after_allocation(new_type,n,sz)
#define DIAGNOSTIC_TEMP_VARIABLES(var_type,var_name,var_value) var_type var_name = var_value

#endif

struct debugging_data {
  unsigned long int tsemb_elongate_decreasing;
  unsigned long int tsemb_elongate_stopped;
  unsigned long int randomize_last_segment_elongation_negative;
  unsigned long int randomize_last_segment_elongation_zero_length;
  unsigned long int randomize_last_segment_elongation_fixedstepelongation_negative;
  unsigned long int randomize_last_segment_elongation_fixedstepelongation_zero;
  unsigned long int tension_direction_model_no_direction;
  debugging_data(): tsemb_elongate_decreasing(0), tsemb_elongate_stopped(0), randomize_last_segment_elongation_negative(0), randomize_last_segment_elongation_zero_length(0), randomize_last_segment_elongation_fixedstepelongation_negative(0), randomize_last_segment_elongation_fixedstepelongation_zero(0), tension_direction_model_no_direction(0) {}
};

extern debugging_data debugging;

void global_debug_report();

#endif
