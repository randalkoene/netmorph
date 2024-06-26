#! /bin/bash
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
# rectangle_and_electrodes.sh
#
# Randal A. Koene, 20041230

curdir=`pwd`
if [ ! -f ../../nibr ]; then
	echo "This script should be run from its example directory"
	echo "examples/rectangle_and_electrodes/ and the nibr (netmorph) program needs"
	echo "to be two steps up from that."
	exit 1
fi

if [ "$1" = "" ]; then
    numneurons=21000
else
    numneurons=$1
fi

cd ../..

./nibr "neurons=$numneurons" "sidelength=1168.0" "shape=rectangle" "seconds=0.0" "figattr_connections=false" "figattr_presynaptic=false" "figattr_postsynaptic=false" "figattr_partitions=false" "figattr_connection_eval=false" "figattr_progress=none" "outattr_directory=examples/rectangle_and_electrodes/"

cd $curdir
