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
# hexagon.sh
#
# Randal A. Koene, 20041230

if [ "$1" = "" ]; then
    numneurons=21000
else
    numneurons=$1
fi

./nibr "neurons=$numneurons" "sidelength=292.0" "shape=hexagon" "electrodes=false" "seconds=0.0" "figattr_connections=false" "figattr_presynaptic=false" "figattr_postsynaptic=false" "figattr_partitions=false" "figattr_connection_eval=false" "figattr_progress=none"
