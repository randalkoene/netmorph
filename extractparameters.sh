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

# extractparameters.sh
# Randal A. Koene, 20070508
#
# This script extracts valid commands from the NETMORPH source code.

if [ "$1" = "-h" ]; then
    echo "extractparameters.sh [[-s] file-to-compare-with]"
    echo ""
    echo "  If file-to-compare-with is empty then all parameter/model commands"
    echo "  are printed. Otherwise, only those that are not found in the"
    echo "  designated file are printed."
    echo ""
    echo "  With -s, commands are shown side by side with the line within"
    echo "  the file, in which they appear."
    exit
fi

if [ "$1" = "" ]; then
    sed -n 's/^.*Specifies_Parameter(\([^)]*\)).*$/\1/p' *.cc *.hh | sort | uniq
    exit
fi

commandnames=`sed -n 's/^.*Specifies_Parameter(\([^)]*\)).*$/\1/p' *.cc *.hh | sort | uniq`

if [ "$1" = "-s" -a "$2" != "" ]; then
    for n in $commandnames; do
	acommand=`echo $n | sed 's/^[^"]*"\([^"]*\)".*/\1/g'`
	occurs=`grep "$acommand" $2`
	if [ $? -eq 0 ]; then
	    echo "$acommand [$n]:"
	    printf "$occurs\n"
	fi
    done
    exit
fi

if [ "$1" = "-s" -a "$2" == "" ]; then
    echo "Syntax error: Parameter -s requires a file name parameter."
    exit 1
fi

for n in $commandnames; do
    acommand=`echo $n | sed 's/^[^"]*"\([^"]*\)".*/\1/g'`
    occurs=`grep -c "$acommand" $1`
    if [ $occurs -eq 0 ]; then
	echo "$acommand [$n]"
    fi
done

