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
// diagnostic.cc
// Randal A. Koene, 20050309

#include "diagnostic.hh"
#include <unistd.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
using namespace std;

#ifdef EVALUATE_POSSIBLE_CONNECTION_PROFILING
long evaluate_possible_connection_calls = 0;
#endif

#ifdef MEMTESTS

#include "file.hh"
#include "BigString.hh"
#include "BigRegex.hh"

// see man 5 proc for stat documentation
const BigRegex procrx("^\\([0-9]+\\) [^ ]+ [^ ] [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ \\([0-9]+\\) \\([0-9]+\\) [0-9]+ [0-9]+ [0-9]+ [0-9]+ 0 [0-9]+ \\([0-9]+\\) \\([0-9]+\\)");

memory_statistics memstats;

// statistics track actual memory usage in the system before and after
// a specific memory allocation, as well as the memory allocated to a
// specific type of object at a specific time

void memory_allocation_diagnostic::allocating(long intended, long actual) {
  if (intended>0) num_allocations++;
  if (intended<0) num_allocations--;
  last_intended = intended;
  last_actual = actual;
  total_intended += last_intended;
  total_actual += last_actual;
}

memory_statistics::~memory_statistics() {
  ofstream dfl("memory_statistics.ascii");
  if (!dfl) {
    warning("Warning: Unable to create memory_statistics.ascii in memory_statistics::~memory_statistics()\n");
  } else {
    dfl << "# Created by nibr diagnostic: memory statistics, Randal A. Koene\n";
    dfl << "# name: events\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(events) << '\n';
    dfl << "# name: jiffies\n# type: matrix\n# rows: 1\n# columns: 2\n";
    dfl << String(t_ujiffies) << ' ' << String(t_sjiffies) << '\n';
    dfl << "# name: actualused\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(actualused_after) << '\n';
    dfl << "# name: numintervening\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(actualused_intervening_num) << '\n';
    dfl << "# name: intervening\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(actualused_intervening) << '\n';
    dfl << "# name: total\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(allocated_total) << '\n';
    dfl << "# name: programsize\n# type: matrix\n# rows: 1\n# columns: 1\n";
    dfl << String(allocated_program) << '\n';
    dfl << "# name: fibre_segment\n# type: matrix\n# rows: 1\n# columns: 3\n";
    dfl << String(allocated_fibre_segment.num_allocations) << ' ' << String(allocated_fibre_segment.total_intended) << ' ' << String(allocated_fibre_segment.total_actual) << '\n';
    dfl << "# name: Fig_element\n# type: matrix\n# rows: 1\n# columns: 3\n";
    dfl << String(allocated_Fig_element.num_allocations) << ' ' << String(allocated_Fig_element.total_intended) << ' ' << String(allocated_Fig_element.total_actual) << '\n';
    dfl << "# name: Spatial_Segment_Subset\n# type: matrix\n# rows: 1\n# columns: 3\n";
    dfl << String(allocated_Spatial_Segment_Subset.num_allocations) << ' ' << String(allocated_Spatial_Segment_Subset.total_intended) << ' ' << String(allocated_Spatial_Segment_Subset.total_actual) << '\n';
    dfl << "# name: synapse_structure\n# type: matrix\n# rows: 1\n# columns: 3\n";
    dfl << String(allocated_synapse_structure.num_allocations) << ' ' << String(allocated_synapse_structure.total_intended) << ' ' << String(allocated_synapse_structure.total_actual) << '\n';
    dfl.close();
  }
  //  system("sed 's/^# rows: 1$/# rows: "+String(events)+"/g' memory_statistics.ascii_pre > memory_statistics.ascii");
  //  unlink("memory_statistics.ascii_pre");

  // This assumes that /usr/include/(asm|linux)/param.h defines HZ as 100
  // so that the time in ms is found by jiffies*1000/HZ=jiffies*10.
  String diagnosticscript("\
#! /usr/bin/octave\n# Randal A. Koene\n\
# nibr memory statistics diagnostic output\n\
\n\
LOADPATH=[LOADPATH,':',system('printf $HOME',1),'/octave/common-m//'];\n\
\n\
load -ascii memory_statistics.ascii\n\
\n\
");
  write_file_from_String("memory_statistics.m",diagnosticscript);
}

void memory_statistics::initialize() {
  before_allocation(initial_program_size);
  allocated_program = actualused_before;
  allocated_total = allocated_program;
  actualused_intervening = 0;
  actualused_intervening_num = 0;
  actualused_after = actualused_before;
}

void memory_statistics::before_allocation(allocation_event new_type) {
  int pid,ppid,pgrp,session,tty_nr,tpgid;
  unsigned long flags,minflt,cminflt,majflt,cmajflt,utime,stime;
  long cutime,cstime,priority,nice,placeholder,itrealvalue;
  unsigned long starttime,vsize;
  char comm[256], state;
  FILE * statstream = fopen("/proc/self/stat","r");
  if (fscanf(statstream,"%d %s %c %d %d %d %d %d %lu %lu %lu %lu %lu %lu %lu %ld %ld %ld %ld %ld %ld %lu %lu",&pid,comm,&state,&ppid,&pgrp,&session,&tty_nr,&tpgid,&flags,&minflt,&cminflt,&majflt,&cmajflt,&utime,&stime,&cutime,&cstime,&priority,&nice,&placeholder,&itrealvalue,&starttime,&vsize)==23) {
    t_ujiffies = utime;
    t_sjiffies = stime;
    actualused_before = vsize;
    if (actualused_before>actualused_after) {
      actualused_intervening_num++;
      actualused_intervening += (actualused_before - actualused_after);
    }
  } else {
    t_ujiffies = -1;
    t_sjiffies = -1;
    actualused_before = -1;
  }
  fclose(statstream);
}

void memory_statistics::after_allocation(allocation_event new_type, int n, long sz) {
  int pid,ppid,pgrp,session,tty_nr,tpgid;
  unsigned long flags,minflt,cminflt,majflt,cmajflt,utime,stime;
  long cutime,cstime,priority,nice,placeholder,itrealvalue;
  unsigned long starttime,vsize;
  char comm[256], state;
  FILE * statstream = fopen("/proc/self/stat","r");
  if (fscanf(statstream,"%d %s %c %d %d %d %d %d %lu %lu %lu %lu %lu %lu %lu %ld %ld %ld %ld %ld %ld %lu %lu",&pid,comm,&state,&ppid,&pgrp,&session,&tty_nr,&tpgid,&flags,&minflt,&cminflt,&majflt,&cmajflt,&utime,&stime,&cutime,&cstime,&priority,&nice,&placeholder,&itrealvalue,&starttime,&vsize)==23) {
    actualused_after = vsize;
    sz *= n;
    switch (new_type) {
    case initial_program_size: break;
    case new_Fig_element: allocated_Fig_element.allocating(sz,actualused_after-actualused_before); break;
    case new_Spatial_Segment_Subset: allocated_Spatial_Segment_Subset.allocating(sz,actualused_after-actualused_before); break;
    case new_synapse_structure: allocated_synapse_structure.allocating(sz,actualused_after-actualused_before); break;
    default: allocated_fibre_segment.allocating(sz,actualused_after-actualused_before);
    }
    allocated_total += sz;
    events++;
  } else {
    actualused_after = -1;
    allocated_total = -1;
  }
  fclose(statstream);
}

#endif

debugging_data debugging;

void global_debug_report() {
#ifdef DEBUGGING_ELONGATION
  cout << "ts_elongation_model_base::elongate decreasing *error*: " << debugging.tsemb_elongate_decreasing << '\n';
  cout << "ts_elongation_model_base::elongate stopped           : " << debugging.tsemb_elongate_stopped << '\n';
  cout << "randomize_last_segment_elongation negative *error*   : " << debugging.randomize_last_segment_elongation_negative << '\n';
  cout << "  fixedstepelongation negative *error*               : " << debugging.randomize_last_segment_elongation_fixedstepelongation_negative << '\n';
  cout << "randomize_last_segment_elongation zero length results: " << debugging.randomize_last_segment_elongation_zero_length << '\n';
  cout << "  fixedstepelongation zero                           : " << debugging.randomize_last_segment_elongation_fixedstepelongation_zero << '\n';
#endif
#ifdef DEBUGGING_DIRECTION
  cout << "tension_direction_model::direction no direction      : " << debugging.tension_direction_model_no_direction << '\n';
#endif
}
