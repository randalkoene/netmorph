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

#
# van_Pelt.sh
#
# Randal A. Koene, 20041230

if [ "$1" = "" ]; then
    numneurons=100
else
    numneurons=$1
fi

./nibr "neurons=$numneurons" "sidelength=1168.0" "shape=rectangle" "connectivity=uniformrandom" "seconds=0.0" "figattr_presynaptic=false" "figattr_postsynaptic=false" "figattr_partitions=false" "figattr_connection_eval=false" "figattr_progress=none"

# net.uniform_random_connectivity(37.0,20,28,PRINCIPAL_NEURON,PRINCIPAL_NEURON); 
# net.uniform_random_connectivity(232.0,17,23,PRINCIPAL_NEURON,INTERNEURON,false);
# net.uniform_random_connectivity(82.0,60,88,INTERNEURON,PRINCIPAL_NEURON,false);
# net.uniform_random_connectivity(82.0,5,7,INTERNEURON,INTERNEURON,false);
