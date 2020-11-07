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

# Randal A. Koene, 20080321
# gplheader.sh
#
# Prepends each source file with the GPLHEADER.

# header and source files in the main source directory

for f in `ls *.cc *.hh`; do

    if [ "$f" != "BigString.hh" -a "$f" != "BigString.cc" -a "$f" != "BigRegex.hh" -a "$f" != "BigRegex.cc" -a "$f" != "mtprng.hh" -a "$f" != "mtprng.cc" ]; then
	hasgplheader=`grep -c 'This file is part of NETMORPH' $f`
	if [ $hasgplheader -eq 0 ]; then
	    cat GPLHEADER > $f.gplprepended
	    cat $f >> $f.gplprepended
	    mv -f $f.gplprepended $f
        fi
    fi

done

# shell scripts in the main source directory

sed '/^[/][*]/ d; s/^[*][/]//g; s/^.\(.*\)$/#\1/g' GPLHEADER > GPLHEADERforsh

for f in `ls *.sh`; do

    hasgplheader=`grep -c 'This file is part of NETMORPH' $f`
    if [ $hasgplheader -eq 0 ]; then
	sed -n '1 p' $f > $f.gplprepended
	cat GPLHEADERforsh >> $f.gplprepended
	sed '1 d' $f >> $f.gplprepended
	mv -f $f.gplprepended $f
	chmod 755 $f
    fi

done

# Makefiles

for f in `ls Makefile*`; do

    hasgplheader=`grep -c 'This file is part of NETMORPH' $f`
    if [ $hasgplheader -eq 0 ]; then
	cat GPLHEADERforsh > $f.gplprepended
	cat $f >> $f.gplprepended
	mv -f $f.gplprepended $f
    fi

done

rm -f GPLHEADERforsh

# templates.hh

hasgplheader=`grep -c 'This file is part of NETMORPH' ~/src/include/templates.hh`
if [ $hasgplheader -eq 0 ]; then
    cat GPLHEADER > ~/src/include/templates.hh.gplprepended
    cat ~/src/include/templates.hh >> ~/src/include/templates.hh.gplprepended
    mv -f ~/src/include/templates.hh.gplprepended ~/src/include/templates.hh
fi
