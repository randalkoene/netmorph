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

# package.sh
# Randal A. Koene, 20061030
#
# Create a versioned package of files needed to compile and use the program.

if [ "$1" = "-h" ]; then
    echo "Usage: ./package.sh [version-ID]"
    exit
else
    if [ "$1" != "" ]; then
	PACKAGEDATE="$1"
    else
	PACKAGEDATE=`date +"%Y%m%d.%H%M"`
  fi
fi

PACKAGENAME=netmorph
PACKAGE="$PACKAGENAME-$PACKAGEDATE"

CURDIR=`pwd`

echo "$PACKAGEDATE" > VERSION

mkdir $PACKAGE
mkdir $PACKAGE/mac
mkdir $PACKAGE/cygwin

cp README INSTALL VERSION LICENSE GPLHEADER Makefile Makefile.linux *.hh *.cc *.sh .nibr2Drc .nibrrc .nibrrc-form-defaults nibr.m evaluate-partitioning.present.m *.ascii $PACKAGE/
cp mac/* $PACKAGE/mac/
cp cygwin/* $PACKAGE/cygwin/
rm -f $PACKAGE/BigRegex.hh $PACKAGE/BigString.hh
cd $PACKAGE
mkdir common-m
cp ~/octave/common-m/plot_wait.m common-m/
mkdir include
cp ~/src/include/templates.hh include/
mkdir geometry
cp ~/src/geometry/Makefile ~/src/geometry/Readme ~/src/geometry/*.c ~/src/geometry/*.h ~/src/geometry/*.cc ~/src/geometry/*.hh geometry/
mkdir dil2al
cp ~/src/dil2al/Makefile ~/src/dil2al/BigRegex.cc ~/src/dil2al/BigRegex.hh ~/src/dil2al/BigString.cc ~/src/dil2al/BigString.hh ~/src/dil2al/regex-gnu.c ~/src/dil2al/regex-gnu.h ~/src/dil2al/rx.h dil2al/

cd $CURDIR
tar zcvf $PACKAGE.tar.gz $PACKAGE
rm -rf $PACKAGE
mv $PACKAGE.tar.gz ../
