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

echo "© Copyright 2008 Randal A. Koene <randalk@netmorph.org>"
echo ""
echo "With design assistance from J. van Pelt & A. van Ooyen, and support"
echo "from the Netherlands Organization for Scientific Research (NWO)"
echo "Program Computational Life Sciences grant CLS2003 (635.100.005) and"
echo "from the EC Marie Curie Research and Training Network (RTN)"
echo "NEURoVERS-it 019247."
echo ""
echo "This installation script is part of NETMORPH."
echo ""
echo "NETMORPH is free software: you can redistribute it and/or modify"
echo "it under the terms of the GNU General Public License as published by"
echo "the Free Software Foundation, either version 3 of the License, or"
echo "(at your option) any later version."
echo ""
echo "NETMORPH is distributed in the hope that it will be useful,"
echo "but WITHOUT ANY WARRANTY; without even the implied warranty of"
echo "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
echo "GNU General Public License for more details."
echo ""
echo "You should have received a copy of the GNU General Public License"
echo "along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>."
echo ""

includefig2dev=1

which sed
if [ $? -ne 0 ]; then
    echo "The sed utility is not available in the execution PATH."
    echo "The sed utility is distributed with all standard versions"
    echo "of the Linux/Unix operating systems."
    echo ""
    echo "Installation cancelled."
    exit
fi

which fig2dev
if [ $? -ne 0 ]; then
    echo "The fig2dev utility is not available in the execution PATH."
    echo "The fig2dev utility is distributed with the XFig package,"
    echo "available at http://www.xfig.org."
    echo ""
    echo "Would you like to compile NETMORPH without fig2dev capabilities?"
    printf "(Y/n): "
    read fig2devspec
    if [ "$fig2devspec" = "n" ]; then
	echo ""
	echo "Installation cancelled."
	exit
    fi
    includefig2dev=0
fi

usercfile=0
echo "Would you like to include the use of a Unix-style resource file,"
printf ".nibrrc, .nibr2Drc (y/N): "
read rcfilespec
if [ "$rcfilespec" = "y" ]; then
    echo "Use of a resource file will be compiled into NETMORPH."
    usercfile=1
fi

if [ "$1" = "-nodetect" ]; then
    machinespec=0
else 
    machinespec=`echo $MACHTYPE | grep -c "linux"`
fi

if [ $machinespec -ne 1 ]; then
  if [ "$1" = "-nodetect" ]; then
      echo "Installation invoked with the -nodetect option."
  else
      echo "Your operating system does not appear to be Linux."
  fi
  echo "Please select one of the following:"
  echo "[1] Install for Linux."
  echo "[2] Install for Linux64."
  echo "[3] Install for Macintosh OSX."
  echo "[4] Install for Cygwin."
  echo "(Or cancel installation by entering a value other than 1-4.)"
  printf "Installation choice: "
  read machinespec
  if [ $machinespec -ne 1 -a $machinespec -ne 2 -a $machinespec -ne 3 -a $machinespec -ne 4 ]; then
    echo "Installation cancelled."
    exit 1
  fi
fi

if [ $machinespec -eq 1 ]; then
# Linux
  cp -f Makefile.linux Makefile
else
  if [ $machinespec -eq 2 ]; then
  # Linux64
    sed 's/MACHSTR=x86/#MACHSTR=x86/g' Makefile.linux > Makefile 
  else
    if [ $machinespec -eq 3 ]; then
    # Mac
      sed 's/MACHSTR=x86/#MACHSTR=x86/g; s/^\(INCLUDES=.*\)$/\1 -I.\/mac/g' Makefile.linux > Makefile 
    else
      if [ $machinespec -eq 4 ]; then
      # Cygwin
	sed 's/^\(INCLUDES=.*\)$/\1 -I.\/cygwin/g' Makefile.linux > Makefile
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
fi

if [ $includefig2dev -eq 0 ]; then
    sed 's/^\(OUTPUTCHOICES=.*\)$/\1 -DCOMPILE_WITHOUT_FIG2DEV/g' Makefile > Makefile.nofig2dev
    mv -f Makefile.nofig2dev Makefile
fi

if [ $usercfile -ne 0 ]; then
    sed 's/^\(MODELCHOICES=.*\)$/\1 -DUSE_RCFILE/g' Makefile > Makefile.rcfile
    mv -f Makefile.rcfile Makefile
fi

PACKAGENAME=netmorph
PACKAGEDATE=`cat VERSION`
PACKAGE="$PACKAGENAME-$PACKAGEDATE"

CURDIR=`pwd`

echo "==> Preparing netmorph code..."

cd "$CURDIR"
ln -s evaluate-partitioning.present.m epp.m

make clean
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
