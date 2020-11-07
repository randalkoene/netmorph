# NETMORPH (netmorph.org)

A Framework for the Generation of Large Scale Neuronal Networks that is
based on the Simulated Development of Neuron Morphology and Synaptic
Connectivity by means of Neurite Growth and Synapse Formation Models

(N)NETMORPH

Randal A. Koene

========================================================================

## Distribution and Downloading

NETMORPH is free software under the terms of the GNU General Public License
and is made available via repository on Github at:

???

To obtain the software, you can clone the repository. Alternatively, you
can download a compressed tar archive of the latest version:

http://www.netmorph.org/Home/netmorph-pub1.0.tar.gz

For more information, to contact the author, or for additional network
development scripts and tools, please see the website at:

http://www.netmorph.org


## Installation


NETMORPH depends on the GNU system libraries and is written in standard
C++, although compilation during development is done only with the
g++/gcc GNU C++/C compilers.

1. Extract the package `netmorph-YYYYMMDD.tar.gz`.
   (YYYYMMDD specifies the creation date of the package.)

   A directory netmorph-YYYYMMDD will be created that contains
   the source files of NETMORPH.

2. Run the installation script:
   `./install.sh`

   Successful compilation results in the programs
   `netmorph` and `netmorph2D`.

At this point it is possible to use the netmorph programs to run
simulations of network development. In order to perform evaluation
of simulation results and to produce graphs of a resulting network
and its statistical data, the following additional tools need to
be installed:

`octave` (http://www.octave.org/)
`xfig` (http://www.xfig.org/)
`ImageMagick` (http://www.imagemagick.org/)

These are standard mathematics and graphics packages that are
available through the package managers of all major Linux
distributions.


  © Copyright 2008 Randal A. Koene <randalk@netmorph.org>
  
  With design assistance from J. van Pelt & A. van Ooyen, and support
  from the Netherlands Organization for Scientific Research (NWO)
  Program Computational Life Sciences grant CLS2003 (635.100.005) and
  from the EC Marie Curie Research and Training Network (RTN)
  NEURoVERS-it 019247.

  This file is part of NETMORPH.

  NETMORPH is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  NETMORPH is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with NETMORPH.  If not, see <http://www.gnu.org/licenses/>.