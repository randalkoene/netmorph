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

# install.sh
# Randal A. Koene, 20061030, 20070208
#
# Install from a versioned package of files and compile the program.

which fig2dev
if [ $? -ne 0 ]; then
  echo "Installation cannot proceed:"
  echo "  To compile NETMORPH, the fig2dev program that is distributed with"
  echo "  the XFig package (see http://www.xfig.org) must be installed and"
  echo "  locatable by the command 'which fig2dev'."
  exit 1
fi

machinespec=`echo $MACHTYPE | grep -c "linux"`
if [ $machinespec -ne 1 ]; then
  echo "Your operating system does not appear to be Linux."
  echo "Please select one of the following:"
  echo "[1] Install for Linux."
  echo "[2] Install for Macintosh OSX."
  echo "[3] Install for Cygwin."
  echo "(Or cancel installation by entering a value other than 1-3.)"
  printf "Installation choice: "
  read machinespec
  if [ $machinespec -ne 1 -a $machinespec -ne 2 -a $machinespec -ne 3 ]; then
    echo "Installation cancelled."
    exit 1
  fi
fi

if [ $machinespec -eq 1 ]; then
  cp -f Makefile.linux Makefile
else
  if [ $machinespec -eq 2 ]; then
    cp -f Makefile.mac Makefile
  else
    if [ $machinespec -eq 3 ]; then
      cp -f Makefile.cygwin Makefile
      cd dil2al
      mv -f Makefile Makefile.linux
      sed 's/[-]Werror//g' Makefile.linux > Makefile
      cd ..
    else
      echo "Installation cancelled."
      exit
    fi
  fi
fi

PACKAGENAME=netmorph
PACKAGEDATE=`cat VERSION`
PACKAGE="$PACKAGENAME-$PACKAGEDATE"

mkdir -p ~/octave/common-m
mkdir -p ~/src/include
mkdir -p ~/src/geometry
mkdir -p ~/src/dil2al
mkdir -p ~/src/lib

mv -f common-m/* ~/octave/common-m/
mv -f include/* ~/src/include/
mv -f geometry/* ~/src/geometry/
mv -f dil2al/* ~/src/dil2al/

CURDIR=`pwd`

cd ~/src/dil2al

# in case of older versions
make clean
# new dependencies
ln -s regex-gnu.c regex.c
ln -s regex-gnu.h regex.h
make regex-gnu.o
make regex.o
make BigRegex.o
make BigString.o

cd ~/src/geometry

# in case of older versions
make clean
# new dependencies
make

cd ~/src/include

ln -s ../geometry/Point2D.hh
ln -s ../geometry/Point3D.hh
ln -s ../geometry/Vector2D.hh
ln -s ../geometry/Vector3D.hh

cd ~/src/lib

ln -s ../geometry/libgeometry.a

cd "$CURDIR"

ln -s ~/src/dil2al/BigRegex.hh
ln -s ~/src/dil2al/BigString.hh
ln -s ~/src/dil2al/regex-gnu.h
ln -s ~/src/dil2al/BigRegex.o
ln -s ~/src/dil2al/BigString.o
ln -s ~/src/dil2al/regex-gnu.o
ln -s evaluate-partitioning.present.m epp.m
make

if [ -f nibr ]; then
    echo ""
    echo "Compilation appears successful."
    echo "You may now run simulations with ./netmorph or ./netmorph2D."
    echo "For information about simulation parameters, see any of the following sources:"
    echo "  ./netmorph -h"
    echo "  http://rak.minduploading.org:8080/caspan/Members/randalk/commands-and-parameters/"
    echo "  http://rak.minduploading.org:8080/caspan/Members/randalk/tutorial/"
    echo ""
    echo "Happy modeling! Randal A. Koene (randalk@minduploading.org)"
else
    make > make.log 2>&1
    grep "replicate" make.log
    if [ $? -eq 0 ]; then
	echo ""
	echo "It appears that compilation may have run into a known problem on some platforms."
	echo "Attempting a workaround with the BIGSTRING_REPLICATE_PROBLEM option."
	echo ""

	NIBRCONDITIONAL=-DBIGSTRING_REPLICATE_PROBLEM; export NIBRCONDITIONAL

	make
	if [ -f nibr ]; then
	    echo ""
	    echo "Compilation appears successful."
	    echo "You may now run simulations with ./netmorph or ./netmorph2D."
	    echo "For information about simulation parameters, see any of the following sources:"
	    echo "  ./netmorph -h"
	    echo "  http://rak.minduploading.org:8080/caspan/Members/randalk/commands-and-parameters/"
	    echo "  http://rak.minduploading.org:8080/caspan/Members/randalk/tutorial/"
	    echo ""
	    echo "Happy modeling! Randal A. Koene (randalk@minduploading.org)"
	else
	    echo ""
	    echo "There appears to be an unresolved compilation error."
	    echo "Please send bug reports to randalk@minduploading.org (Randal A. Koene)."
	fi
    fi
fi
