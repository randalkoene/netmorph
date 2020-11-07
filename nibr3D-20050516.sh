#! /bin/sh
# © Copyright 2008 Randal A. Koene <randalk@netmorph.org>
# 
# With design assistance from J. van Pelt & A. van Ooyen, and support
# from the Netherlands Organization for Scientific Research (NWO)
# Program Computational Life Sciences grant CLS2003 (635.100.005) and
# from the EC Marie Curie Research and Training Network (RTN)
# NEURoVERS-it 019247.

# This file is part of NETMORPH.

# NETMORPH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# NETMORPH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.

# nibr3D-20050516.sh
# Randal A. Koene
#
# This script attempts to produce a 3D simulation with animation output that
# demonstrates network growth with a set of cubic polynomial growth functions.
./nibr dt=1000 shape=box electrodes=false fibreswithturns=true branchinsegment=true branchatinitlength=false growthfunction=polynomialO3 outattr_show_progress=true outattr_show_figure=true outattr_show_stats=true figattr_neurons=true figattr_connections=false figattr_presynaptic=true figattr_postsynaptic=true figattr_tsupd_visibly=true figattr_synapses=false figattr_partitions=true figattr_connection_eval=false sample_dt=43200 figureTeXwidth=3.5 figuresequence=true figattr_progress=none combinesequence=true autorotatesequence=true ROT_x=0 ROT_y=0 ROT_z=0 ROT_interval_x=0 ROT_interval_y=1.5 ROT_interval_z=0 figattr_use_color=true combinemagnification=0.5
