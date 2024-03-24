#! /bin/bash
# © Copyright 2024 Randal A. Koene <randalk@netmorph.org>
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
# fig8B-smaller.sh
#
# Randal A. Koene, 20240324

curdir=`pwd`
if [ ! -f ../../nibr ]; then
	echo "This script should be run from its example directory"
	echo "examples/fig8B/ and the nibr (netmorph) program needs"
	echo "to be two steps up from that."
	exit 1
fi

cd ../..

./nibr "include=examples/fig8B/fig8B-smaller" "outattr_directory=examples/fig8B/"

cd $curdir
