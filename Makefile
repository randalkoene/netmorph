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
# Makefile for NETMORPH
# (by Randal A. Koene, randalk@minduploading.org)

##########################################################
## Compiler
## --------
CC=gcc
CCPP=g++ 
##########################################################

INC_PATH= ./include
INCLUDES= -I$(INC_PATH)
LIB_PATH= ./lib
LIBRARIES= -L$(LIB_PATH)
OBJ_PATH= ./obj
BIN_PATH= ./bin

##########################################################
## Machine selection
## -----------------
## aurora (x86):
MACHSPEC=
#MACHSTR=x86
MACHOPT=
#MACHOPT=-march=pentium4 -mtune=pentium4 -mfpmath=sse
##
## kenji (SUN Ultra):
#MACHSPEC= -DSUN_RX
#MACHSTR=SUN
##
##########################################################
### Attempting a temporary fix for gcc 4.x compilation of dil2al
### (This should be replaced by actual source code modification.)
CPPTEMPWORKAROUND=-fno-access-control

##########################################################
## Regular Expression Library
## --------------------------
## (Suggestion: Using regex-gnu.h guarantees the greatest
## amount of compatibility and identical behaviour on all
## platforms, since it is included with the dil2al source
## code and has been adapted for reliable use with C++.)
##
## regex-gnu.h adapted for integration with C++:
GCC3FXTRA=$(if $(shell gcc -v 2>&1 | grep "^gcc version[^0-9]*[3-9][.]"),-Wno-unused-function)
ALT_REGEX=regex-gnu.o
ALT_REGEX_H=\"regex-gnu.h\"
MEMTESTS=
#MEMTESTS=-DMEMTESTS
#PROFILING=
PROFILING=-DEVALUATE_POSSIBLE_CONNECTION_PROFILING
#SYNAPSEINVENTORY=
SYNAPSEINVENTORY=-DSYNAPTOGENESIS_AND_LOSS_INVENTORY
#STATISTICSVERBOSITY=-DVERBOSE_EMPTY_DATA_SAMPLE_WARNINGS
#STATISTICS=-DSAMPLES_INCLUDE_NETWORK_GENERATED_STATISTICS $(STATISTICSVERBOSITY)
STATISTICS=-DSAMPLES_INCLUDE_NETWORK_STATISTICS_BASE $(STATISTICSVERBOSITY)
#STATISTICS=-DSAMPLES_INCLUDE_NETWORK_STATISTICS_BASE
MEMORYSAVING=
#MEMORYSAVING=-DREDUCED_MEMORY_PTSEM
#RNGCHOICE=-DFAST_RNG
RNGCHOICE=-DPRECISE_RNG
#MODELCHOICES=-DADM_PRIORITIZE_SPEED_OVER_SPACE -DTESTING_SIMPLE_ATTRACTION -DENABLE_APICAL_TUFTING_TRIGGER_DISTANCE $(RNGCHOICE) $(MEMORYSAVING) -DINCLUDE_SIMPLE_FIBER_DIAMETER
MODELCHOICES=-DADM_PRIORITIZE_SPEED_OVER_SPACE -DTESTING_SIMPLE_ATTRACTION $(RNGCHOICE) $(MEMORYSAVING) -DINCLUDE_SIMPLE_FIBER_DIAMETER -DUSE_RCFILE -DUSE_RCFILE
OUTPUTCHOICES=-DCHECK_XFIG_RANGES
#-DTESTING_SPIKING
#-DLEGACY_ELONGATION
SIMULATIONMETHODS=-DENABLE_FIXED_STEP_SIMULATION -DFIXED_STEP_SIMULATION_START_AT_DT -DTESTING_SPIKING
DEBUGGING=
#DEBUGGING=-DDEBUG_SPHERICAL_BOUNDARY
#DEBUGGING=-DTESTING_ELONGATION_TOTAL
#DEBUGGING=-DTEST_FOR_NAN
#DEBUGGING=-DINCLUDE_PDF_SAMPLING
#-DDEBUGGING_ELONGATION -DDEBUGGING_DIRECTION -DINCLUDE_PDF_SAMPLING
## Put conditional defines, such as -D__SPATIAL_SEGMENT_SUBSET_TEST into
## the variable NIBRCONDITIONAL and export before running make.
CFXTRA= -D__DIL2AL__ -DSTDC_HEADERS -pedantic -Wall -Wno-char-subscripts $(GCC3FXTRA) $(MEMTESTS) $(PROFILING) $(SYNAPSEINVENTORY) $(STATISTICS) $(MODELCHOICES) $(SIMULATIONMETHODS) $(DEBUGGING) $(OUTPUTCHOICES)
CPPXTRA= -D__DIL2AL__ -D_ALT_REGEX_H=$(ALT_REGEX_H) -D_USE_ALT_REGEX -D_CPP_REGEX $(NIBRCONDITIONAL) $(MACHOPT) -mieee-fp -ffast-math -pedantic -Wall $(GCC3FXTRA) $(MEMTESTS) $(PROFILING) $(SYNAPSEINVENTORY) $(STATISTICS) $(MODELCHOICES) $(SIMULATIONMETHODS) $(DEBUGGING) $(CPPTEMPWORKAROUND) $(OUTPUTCHOICES)
REGEXSTR=regex-gnu_for_C++
FIG2DEV=$(shell which fig2dev)
##
## rx.h:
#ALT_REGEX=rx.o
#ALT_REGEX_H=
#CFXTRA= -D__DIL2AL__
#CPPXTRA= -D__DIL2AL__ -D_BIGREGEX_HAS_RM_SPEP
#REGEXSTR=rx
##
## generic GNU regex (Linux/GNU specific):
#ALT_REGEX=regex.o
#ALT_REGEX_H=regex.h
#CFXTRA=
#CPPXTRA= -D_ALT_REGEX_H=$(ALT_REGEX_H) -D_USE_ALT_REGEX -D__DIL2AL__ -D_BIGREGEX_HAS_RM_SPEP
#REGEXSTR=regex
##
## system <regex.h> (SUN specific):
#ALT_REGEX=
#ALT_REGEX_H=
#CFXTRA= -DSYSTEM_RX
#CPPXTRA= -D__DIL2AL__ -DSYSTEM_RX
#REGEXSTR=system
##########################################################

##########################################################
## Safe Regular Expressions
## ------------------------
## assume '\0' within String length and rm structure has
## no sp/ep pointers:
#SAFEREGEX=
## make no assumptions:
SAFEREGEX= -D_BIGREGEX_SAFE_MATCHES
##########################################################

##########################################################
## Compiler Options
## ----------------
# Normal:
COMPOPT=-O3
OPTSTR=optimized
# Unspecified:
#COMPOPT=
#OPTSTR=unspecified
## debugging information:
# For gprof:
#GPROF=-pg
# For Valgrind:
#COMPOPT=-g -O0
#OPTSTR=Valgrind
# For gprof:
#COMPOPT=-pg -O3
#OPTSTR=Gprof
#COMPOPT= -g
#OPTSTR=debugging_info
##########################################################

##########################################################
## C++ Specific Compiler Options
## -----------------------------
## debugging information:
#CPPOPT=
## generate profile information for use with gprof:
#CPPOPT=
## optimized:
CPPOPT= -felide-constructors
##########################################################

CFLAGS= $(COMPOPT) $(MACHSPEC) $(CFXTRA)
CPPFLAGS= $(COMPOPT) $(CPPOPT) $(MACHSPEC) $(CPPXTRA) $(SAFEREGEX) $(INCLUDES)

LIB = libgeometry.a

NAME = nibr nibr2D wizard netmorph netmorph2D

all: SPECS nibr nibr2D wizard

SPECS:
	@echo '-------------------------------------------------------------------'
	@echo 'Compilation options: $(MACHSTR), $(OPTSTR), $(REGEXSTR)'
	@echo '-------------------------------------------------------------------'

$(OBJ_PATH)/diagnostic.o: diagnostic.cc diagnostic.hh file.hh BigString.hh BigRegex.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c diagnostic.cc -o $(OBJ_PATH)/diagnostic.o

$(OBJ_PATH)/state_storable.o: state_storable.cc state_storable.hh global.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c state_storable.cc -o $(OBJ_PATH)/state_storable.o

$(OBJ_PATH)/mtprng.o: mtprng.cc mtprng.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c mtprng.cc -o $(OBJ_PATH)/mtprng.o

$(OBJ_PATH)/global.o: global.cc global.hh file.hh Fig_Object.hh Command_Line_Parameters.hh mtprng.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c global.cc -o $(OBJ_PATH)/global.o

$(OBJ_PATH)/file.o: file.cc file.hh BigString.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c file.cc -o $(OBJ_PATH)/file.o

$(OBJ_PATH)/StringList.o: StringList.cc StringList.hh BigString.hh BigRegex.hh regex-gnu.h
	$(CCPP) $(CPPFLAGS) $(COLORS) -c StringList.cc -o $(OBJ_PATH)/StringList.o

$(OBJ_PATH)/event.o $(OBJ_PATH)/event2D.o: event.cc event.hh network.hh global.hh Command_Line_Parameters.hh state_storable.hh
	$(CCPP) $(CPPFLAGS) -DVECTOR3D -c event.cc -o $(OBJ_PATH)/event.o
	$(CCPP) $(CPPFLAGS) -DVECTOR2D -c event.cc -o $(OBJ_PATH)/event2D.o

$(OBJ_PATH)/spatial.o $(OBJ_PATH)/spatial2D.o: spatial.cc spatial.hh Spatial_Presentation.hh global.hh
	$(CCPP) $(CPPFLAGS) -DVECTOR3D -c spatial.cc -o $(OBJ_PATH)/spatial.o
	$(CCPP) $(CPPFLAGS) -DVECTOR2D -c spatial.cc -o $(OBJ_PATH)/spatial2D.o

$(OBJ_PATH)/network.o $(OBJ_PATH)/network2D.o: network.cc network.hh neuron.hh Network_Statistics.hh synapse.hh synapse_formation_model.hh global.hh file.hh Command_Line_Parameters.hh Network_Generated_Statistics.hh Sampled_Output.hh Network_Statistics.hh Spatial_Segment_Subset.hh Results.hh Connection_Statistics.hh dendritic_growth_model.hh neurite_diameter_model.hh spatial.hh event.hh state_storable.hh Txt_Object.hh VRML_Object.hh Catacomb_Object.hh slice.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c network.cc -o $(OBJ_PATH)/network.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c network.cc -o $(OBJ_PATH)/network2D.o

$(OBJ_PATH)/neuron.o $(OBJ_PATH)/neuron2D.o: neuron.cc neuron.hh network.hh prepost_structure.hh connection.hh fibre_elongation_model.hh axon_direction_model.hh branching_models.hh turning_models.hh Color_Table.hh spatial.hh event.hh global.hh Network_Generated_Statistics.hh Sampled_Output.hh Command_Line_Parameters.hh BigString.hh state_storable.hh Txt_Object.hh VRML_Object.hh Catacomb_Object.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c neuron.cc -o $(OBJ_PATH)/neuron.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c neuron.cc -o $(OBJ_PATH)/neuron2D.o

$(OBJ_PATH)/IFactivity.o $(OBJ_PATH)/IFactivity2D.o: IFactivity.cc IFactivity.hh neuron.hh connection.hh Sampled_Output.hh event.hh
	$(CCPP) $(CPPFLAGS) -DVECTOR3D -c IFactivity.cc -o $(OBJ_PATH)/IFactivity.o
	$(CCPP) $(CPPFLAGS) -DVECTOR2D -c IFactivity.cc -o $(OBJ_PATH)/IFactivity2D.o

$(OBJ_PATH)/synapse.o $(OBJ_PATH)/synapse2D.o: synapse.cc synapse.hh connection.hh synapse_structure.hh Fig_Object.hh Color_Table.hh global.hh state_storable.hh Txt_Object.hh VRML_Object.hh Catacomb_Object.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c synapse.cc -o $(OBJ_PATH)/synapse.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c synapse.cc -o $(OBJ_PATH)/synapse2D.o

$(OBJ_PATH)/connection.o $(OBJ_PATH)/connection2D.o: connection.cc connection.hh neuron.hh Fig_Object.hh Sampled_Output.hh state_storable.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c connection.cc -o $(OBJ_PATH)/connection.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c connection.cc -o $(OBJ_PATH)/connection2D.o

$(OBJ_PATH)/generic_synapse.o: generic_synapse.cc generic_synapse.hh nibr.hh BigString.hh BigRegex.hh regex-gnu.h
	$(CCPP) $(CPPFLAGS) $(COLORS) -c generic_synapse.cc -o $(OBJ_PATH)/generic_synapse.o

$(OBJ_PATH)/synapse_formation_model.o $(OBJ_PATH)/synapse_formation_model2D.o: synapse_formation_model.cc synapse_formation_model.hh fibre_structure.hh synapse_structure.hh synapse.hh Connection_Statistics.hh neuron.hh spatial.hh Command_Line_Parameters.hh diagnostic.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c synapse_formation_model.cc -o $(OBJ_PATH)/synapse_formation_model.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c synapse_formation_model.cc -o $(OBJ_PATH)/synapse_formation_model2D.o

$(OBJ_PATH)/branching_models.o $(OBJ_PATH)/branching_models2D.o: branching_models.cc branching_models.hh dendritic_growth_model.hh fibre_structure.hh neuron.hh spatial.hh Command_Line_Parameters.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c branching_models.cc -o $(OBJ_PATH)/branching_models.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c branching_models.cc -o $(OBJ_PATH)/branching_models2D.o

$(OBJ_PATH)/turning_models.o $(OBJ_PATH)/turning_models2D.o: turning_models.cc turning_models.hh dendritic_growth_model.hh fibre_structure.hh neuron.hh spatial.hh Command_Line_Parameters.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c turning_models.cc -o $(OBJ_PATH)/turning_models.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c turning_models.cc -o $(OBJ_PATH)/turning_models2D.o

$(OBJ_PATH)/axon_direction_model.o $(OBJ_PATH)/axon_direction_model2D.o: axon_direction_model.cc axon_direction_model.hh dendritic_growth_model.hh fibre_structure.hh neuron.hh spatial.hh Spatial_Presentation.hh Command_Line_Parameters.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c axon_direction_model.cc -o $(OBJ_PATH)/axon_direction_model.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c axon_direction_model.cc -o $(OBJ_PATH)/axon_direction_model2D.o

$(OBJ_PATH)/fibre_elongation_model.o $(OBJ_PATH)/fibre_elongation_model2D.o: fibre_elongation_model.cc fibre_elongation_model.hh dendritic_growth_model.hh axon_direction_model.hh fibre_structure.hh neuron.hh event.hh Command_Line_Parameters.hh global.hh Txt_Object.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c fibre_elongation_model.cc -o $(OBJ_PATH)/fibre_elongation_model.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c fibre_elongation_model.cc -o $(OBJ_PATH)/fibre_elongation_model2D.o

$(OBJ_PATH)/dendritic_growth_model.o $(OBJ_PATH)/dendritic_growth_model2D.o: dendritic_growth_model.cc dendritic_growth_model.hh axon_direction_model.hh fibre_structure.hh network.hh neuron.hh Connection_Statistics.hh Sampled_Output.hh Command_Line_Parameters.hh global.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c dendritic_growth_model.cc -o $(OBJ_PATH)/dendritic_growth_model.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c dendritic_growth_model.cc -o $(OBJ_PATH)/dendritic_growth_model2D.o

$(OBJ_PATH)/neurite_diameter_model.o $(OBJ_PATH)/neurite_diameter_model2D.o: neurite_diameter_model.cc neurite_diameter_model.hh fibre_structure.hh network.hh Command_Line_Parameters.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c neurite_diameter_model.cc -o $(OBJ_PATH)/neurite_diameter_model.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c neurite_diameter_model.cc -o $(OBJ_PATH)/neurite_diameter_model2D.o

$(OBJ_PATH)/environment_physics.o: environment_physics.cc environment_physics.hh neuron.hh network.hh event.hh Command_Line_Parameters.hh spatial.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c environment_physics.cc -o $(OBJ_PATH)/environment_physics.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c environment_physics.cc -o $(OBJ_PATH)/environment_physics2D.o

$(OBJ_PATH)/synapse_structure.o $(OBJ_PATH)/synapse_structure2D.o: synapse_structure.cc synapse_structure.hh fibre_structure.hh Fig_Object.hh Sampled_Output.hh spatial.hh global.hh state_storable.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c synapse_structure.cc -o $(OBJ_PATH)/synapse_structure.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c synapse_structure.cc -o $(OBJ_PATH)/synapse_structure2D.o

$(OBJ_PATH)/fibre_structure.o $(OBJ_PATH)/fibre_structure2D.o: fibre_structure.cc fibre_structure.hh dendritic_growth_model.hh axon_direction_model.hh branching_models.hh turning_models.hh fibre_elongation_model.hh environment_physics.hh neuron.hh event.hh Network_Generated_Statistics.hh Sampled_Output.hh Fig_Object.hh spatial.hh global.hh Color_Table.hh diagnostic.hh state_storable.hh Txt_Object.hh VRML_Object.hh Catacomb_Object.hh Command_Line_Parameters.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c fibre_structure.cc -o $(OBJ_PATH)/fibre_structure.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c fibre_structure.cc -o $(OBJ_PATH)/fibre_structure2D.o

$(OBJ_PATH)/prepost_structure.o $(OBJ_PATH)/prepost_structure2D.o: prepost_structure.cc prepost_structure.hh neuron.hh network.hh axon_direction_model.hh branching_models.hh turning_models.hh spatial.hh fibre_structure.hh Network_Generated_Statistics.hh global.hh state_storable.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c prepost_structure.cc -o $(OBJ_PATH)/prepost_structure.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c prepost_structure.cc -o $(OBJ_PATH)/prepost_structure2D.o

$(OBJ_PATH)/Connection_Statistics.o $(OBJ_PATH)/Connection_Statistics2D.o: Connection_Statistics.cc Connection_Statistics.hh neuron.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c Connection_Statistics.cc -o $(OBJ_PATH)/Connection_Statistics.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c Connection_Statistics.cc -o $(OBJ_PATH)/Connection_Statistics2D.o

$(OBJ_PATH)/Network_Statistics.o $(OBJ_PATH)/Network_Statistics2D.o: Network_Statistics.cc Network_Statistics.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c Network_Statistics.cc -o $(OBJ_PATH)/Network_Statistics.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c Network_Statistics.cc -o $(OBJ_PATH)/Network_Statistics2D.o

$(OBJ_PATH)/Spatial_Segment_Subset.o $(OBJ_PATH)/Spatial_Segment_Subset2D.o: Spatial_Segment_Subset.cc Spatial_Segment_Subset.hh fibre_structure.hh synapse_formation_model.hh Fig_Object.hh Color_Table.hh spatial.hh global.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c Spatial_Segment_Subset.cc -o $(OBJ_PATH)/Spatial_Segment_Subset.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c Spatial_Segment_Subset.cc -o $(OBJ_PATH)/Spatial_Segment_Subset2D.o

$(OBJ_PATH)/Command.o: Command.cc BigString.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Command.cc -o $(OBJ_PATH)/Command.o

$(OBJ_PATH)/Command_Line_Parameters.o: Command_Line_Parameters.cc Command_Line_Parameters.hh file.hh global.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Command_Line_Parameters.cc -o $(OBJ_PATH)/Command_Line_Parameters.o

$(OBJ_PATH)/Results.o: Results.cc Results.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Results.cc -o $(OBJ_PATH)/Results.o

$(OBJ_PATH)/Network_Generated_Statistics.o: Network_Generated_Statistics.cc Network_Generated_Statistics.hh Command_Line_Parameters.hh global.hh BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Network_Generated_Statistics.cc -o $(OBJ_PATH)/Network_Generated_Statistics.o

$(OBJ_PATH)/Spatial_Presentation.o $(OBJ_PATH)/Spatial_Presentation2D.o: Spatial_Presentation.cc Spatial_Presentation.hh global.hh Command_Line_Parameters.hh spatial.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -DFIG2DEV=\"$(FIG2DEV)\" -c Spatial_Presentation.cc -o $(OBJ_PATH)/Spatial_Presentation.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -DFIG2DEV=\"$(FIG2DEV)\" -c Spatial_Presentation.cc -o $(OBJ_PATH)/Spatial_Presentation2D.o

$(OBJ_PATH)/Sampled_Output.o $(OBJ_PATH)/Sampled_Output2D.o: Sampled_Output.cc Sampled_Output.hh Network_Generated_Statistics.hh Color_Table.hh Spatial_Presentation.hh Command_Line_Parameters.hh spatial.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -DFIG2DEV=\"$(FIG2DEV)\" -c Sampled_Output.cc -o $(OBJ_PATH)/Sampled_Output.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -DFIG2DEV=\"$(FIG2DEV)\" -c Sampled_Output.cc -o $(OBJ_PATH)/Sampled_Output2D.o

$(OBJ_PATH)/Color_Table.o: Color_Table.cc Color_Table.hh Fig_Object.hh Command_Line_Parameters.hh global.hh BigString.hh 
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Color_Table.cc -o $(OBJ_PATH)/Color_Table.o

$(OBJ_PATH)/Fig_Object.o: Fig_Object.cc Fig_Object.hh BigString.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Fig_Object.cc -o $(OBJ_PATH)/Fig_Object.o

$(OBJ_PATH)/Txt_Object.o: Txt_Object.cc Txt_Object.hh global.hh BigString.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c Txt_Object.cc -o $(OBJ_PATH)/Txt_Object.o

$(OBJ_PATH)/VRML_Object.o: VRML_Object.cc VRML_Object.hh global.hh spatial.hh BigString.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c VRML_Object.cc -o $(OBJ_PATH)/VRML_Object.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c VRML_Object.cc -o $(OBJ_PATH)/VRML_Object2D.o

$(OBJ_PATH)/Catacomb_Object.o: Catacomb_Object.cc Catacomb_Object.hh global.hh network.hh spatial.hh BigString.hh diagnostic.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c Catacomb_Object.cc -o $(OBJ_PATH)/Catacomb_Object.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c Catacomb_Object.cc -o $(OBJ_PATH)/Catacomb_Object2D.o

$(OBJ_PATH)/slice.o: slice.cc slice.hh spatial.hh Command_Line_Parameters.hh global.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -c slice.cc -o $(OBJ_PATH)/slice.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -c slice.cc -o $(OBJ_PATH)/slice2D.o

$(OBJ_PATH)/nibr.o $(OBJ_PATH)/nibr2D.o: nibr.cc nibr.hh neuron.hh network.hh prepost_structure.hh synapse.hh synapse_formation_model.hh axon_direction_model.hh environment_physics.hh Command_Line_Parameters.hh Sampled_Output.hh Txt_Object.hh Spatial_Presentation.hh file.hh spatial.hh event.hh slice.hh diagnostic.hh BigString.hh BigRegex.hh mtprng.hh regex-gnu.h
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR3D -DFIG2DEV=\"$(FIG2DEV)\" -c nibr.cc -o $(OBJ_PATH)/nibr.o
	$(CCPP) $(CPPFLAGS) $(COLORS) -DVECTOR2D -DFIG2DEV=\"$(FIG2DEV)\" -c nibr.cc -o $(OBJ_PATH)/nibr2D.o

$(OBJ_PATH)/wizard.o: wizard.cc wizard.hh StringList.hh BigString.hh BigRegex.hh regex-gnu.h
	$(CCPP) $(CPPFLAGS) $(COLORS) -c wizard.cc -o $(OBJ_PATH)/wizard.o

$(OBJ_PATH)/BigString.o: BigString.cc BigString.hh
	$(CCPP) $(CPPFLAGS) $(COLORS) -c BigString.cc -o $(OBJ_PATH)/BigString.o

$(OBJ_PATH)/BigRegex.o: BigRegex.cc BigRegex.hh rx.h regex-gnu.h
	$(CCPP) $(CPPFLAGS) $(COLORS) -c BigRegex.cc -o $(OBJ_PATH)/BigRegex.o

$(OBJ_PATH)/regex-gnu.o: regex-gnu.c regex-gnu.h
	$(CCPP) $(CFLAGS) -c regex-gnu.c -o $(OBJ_PATH)/regex-gnu.o

$(OBJ_PATH)/Point3D.o: Point3D.cc Point3D.hh Vector3D.hh
	$(CCPP) $(CPPFLAGS) -c Point3D.cc -o $(OBJ_PATH)/Point3D.o

$(OBJ_PATH)/Vector3D.o: Vector3D.cc Vector3D.hh Point3D.hh
	$(CCPP) $(CPPFLAGS) -c Vector3D.cc -o $(OBJ_PATH)/Vector3D.o

$(OBJ_PATH)/Point2D.o: Point2D.cc Point2D.hh Vector2D.hh
	$(CCPP) $(CPPFLAGS) -c Point2D.cc -o $(OBJ_PATH)/Point2D.o

$(OBJ_PATH)/Vector2D.o: Vector2D.cc Vector2D.hh Point2D.hh
	$(CCPP) $(CPPFLAGS) -c Vector2D.cc -o $(OBJ_PATH)/Vector2D.o

ptest: ptest.c $(LIB_PATH)/$(LIB)
	g++ -o ptest ptest.c $(LIB_PATH)/$(LIB)

vtest: vtest.c $(LIB_PATH)/$(LIB)
	g++ -o vtest vtest.c $(LIB_PATH)/$(LIB)

$(LIB_PATH)/$(LIB): $(OBJ_PATH)/Point3D.o $(OBJ_PATH)/Vector3D.o $(OBJ_PATH)/Point2D.o $(OBJ_PATH)/Vector2D.o
	mkdir -p ./lib
	ar rcsv $(LIB_PATH)/$(LIB) $(OBJ_PATH)/Point3D.o $(OBJ_PATH)/Vector3D.o $(OBJ_PATH)/Point2D.o $(OBJ_PATH)/Vector2D.o

nibr: $(OBJ_PATH)/nibr.o $(OBJ_PATH)/Txt_Object.o $(OBJ_PATH)/VRML_Object.o $(OBJ_PATH)/Catacomb_Object.o $(OBJ_PATH)/Fig_Object.o $(OBJ_PATH)/Color_Table.o $(OBJ_PATH)/Sampled_Output.o $(OBJ_PATH)/Spatial_Presentation.o $(OBJ_PATH)/Network_Generated_Statistics.o $(OBJ_PATH)/Results.o $(OBJ_PATH)/Command_Line_Parameters.o $(OBJ_PATH)/Command.o $(OBJ_PATH)/Spatial_Segment_Subset.o $(OBJ_PATH)/Network_Statistics.o $(OBJ_PATH)/Connection_Statistics.o $(OBJ_PATH)/prepost_structure.o $(OBJ_PATH)/fibre_structure.o $(OBJ_PATH)/synapse_structure.o $(OBJ_PATH)/generic_synapse.o $(OBJ_PATH)/dendritic_growth_model.o $(OBJ_PATH)/synapse_formation_model.o $(OBJ_PATH)/branching_models.o $(OBJ_PATH)/turning_models.o $(OBJ_PATH)/axon_direction_model.o $(OBJ_PATH)/fibre_elongation_model.o $(OBJ_PATH)/neurite_diameter_model.o $(OBJ_PATH)/environment_physics.o $(OBJ_PATH)/connection.o $(OBJ_PATH)/synapse.o $(OBJ_PATH)/IFactivity.o $(OBJ_PATH)/neuron.o $(OBJ_PATH)/network.o $(OBJ_PATH)/spatial.o $(OBJ_PATH)/event.o $(OBJ_PATH)/slice.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/file.o $(OBJ_PATH)/state_storable.o $(OBJ_PATH)/global.o $(OBJ_PATH)/diagnostic.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/mtprng.o $(OBJ_PATH)/regex-gnu.o $(LIB_PATH)/$(LIB)
	@echo "Design rule: Add only the physiology required by function."
	@echo "  ___connected neurons (connection)"
	@echo "  _____individual synapses (synapse)"
	@echo "  _______presynaptic and postsynaptic nodes and branches"
	@echo "  _________terminals, receptive zones/spines, soma, axon hillock"
	$(CCPP) $(CPPFLAGS) $(OBJ_PATH)/nibr.o $(OBJ_PATH)/Txt_Object.o $(OBJ_PATH)/VRML_Object.o $(OBJ_PATH)/Catacomb_Object.o $(OBJ_PATH)/Fig_Object.o $(OBJ_PATH)/Color_Table.o $(OBJ_PATH)/Sampled_Output.o $(OBJ_PATH)/Spatial_Presentation.o $(OBJ_PATH)/Network_Generated_Statistics.o $(OBJ_PATH)/Results.o $(OBJ_PATH)/Command_Line_Parameters.o $(OBJ_PATH)/Command.o $(OBJ_PATH)/Spatial_Segment_Subset.o $(OBJ_PATH)/Network_Statistics.o $(OBJ_PATH)/Connection_Statistics.o $(OBJ_PATH)/prepost_structure.o $(OBJ_PATH)/fibre_structure.o $(OBJ_PATH)/synapse_structure.o $(OBJ_PATH)/generic_synapse.o $(OBJ_PATH)/dendritic_growth_model.o $(OBJ_PATH)/synapse_formation_model.o $(OBJ_PATH)/branching_models.o $(OBJ_PATH)/turning_models.o $(OBJ_PATH)/axon_direction_model.o $(OBJ_PATH)/fibre_elongation_model.o $(OBJ_PATH)/neurite_diameter_model.o $(OBJ_PATH)/environment_physics.o $(OBJ_PATH)/connection.o $(OBJ_PATH)/synapse.o $(OBJ_PATH)/IFactivity.o $(OBJ_PATH)/neuron.o $(OBJ_PATH)/network.o $(OBJ_PATH)/spatial.o $(OBJ_PATH)/event.o $(OBJ_PATH)/slice.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/file.o $(OBJ_PATH)/state_storable.o $(OBJ_PATH)/global.o $(OBJ_PATH)/diagnostic.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/mtprng.o $(OBJ_PATH)/regex-gnu.o -lgeometry -o $(BIN_PATH)/netmorph $(LIBRARIES)
	cd $(BIN_PATH) && ln -f -s netmorph nibr

nibr2D: $(OBJ_PATH)/nibr2D.o $(OBJ_PATH)/Txt_Object.o $(OBJ_PATH)/VRML_Object2D.o $(OBJ_PATH)/Catacomb_Object2D.o $(OBJ_PATH)/Fig_Object.o $(OBJ_PATH)/Color_Table.o $(OBJ_PATH)/Sampled_Output2D.o $(OBJ_PATH)/Spatial_Presentation2D.o $(OBJ_PATH)/Network_Generated_Statistics.o $(OBJ_PATH)/Results.o $(OBJ_PATH)/Command_Line_Parameters.o $(OBJ_PATH)/Command.o $(OBJ_PATH)/Spatial_Segment_Subset2D.o $(OBJ_PATH)/Network_Statistics2D.o $(OBJ_PATH)/Connection_Statistics2D.o $(OBJ_PATH)/prepost_structure2D.o $(OBJ_PATH)/fibre_structure2D.o $(OBJ_PATH)/synapse_structure2D.o $(OBJ_PATH)/generic_synapse.o $(OBJ_PATH)/dendritic_growth_model2D.o $(OBJ_PATH)/synapse_formation_model2D.o $(OBJ_PATH)/branching_models2D.o $(OBJ_PATH)/turning_models2D.o $(OBJ_PATH)/axon_direction_model2D.o $(OBJ_PATH)/fibre_elongation_model2D.o $(OBJ_PATH)/neurite_diameter_model2D.o $(OBJ_PATH)/environment_physics2D.o $(OBJ_PATH)/connection2D.o $(OBJ_PATH)/synapse2D.o $(OBJ_PATH)/IFactivity2D.o $(OBJ_PATH)/neuron2D.o $(OBJ_PATH)/network2D.o $(OBJ_PATH)/spatial2D.o $(OBJ_PATH)/event2D.o $(OBJ_PATH)/slice2D.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/file.o $(OBJ_PATH)/state_storable.o $(OBJ_PATH)/global.o $(OBJ_PATH)/diagnostic.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/mtprng.o $(OBJ_PATH)/regex-gnu.o $(LIB_PATH)/$(LIB)
	@echo "Design rule: Add only the physiology required by function."
	@echo "  ___connected neurons (connection)"
	@echo "  _____individual synapses (synapse)"
	@echo "  _______presynaptic and postsynaptic nodes and branches"
	@echo "  _________terminals, receptive zones/spines, soma, axon hillock"
	$(CCPP) $(CPPFLAGS) $(OBJ_PATH)/nibr2D.o $(OBJ_PATH)/Txt_Object.o $(OBJ_PATH)/VRML_Object2D.o $(OBJ_PATH)/Catacomb_Object2D.o $(OBJ_PATH)/Fig_Object.o $(OBJ_PATH)/Color_Table.o $(OBJ_PATH)/Sampled_Output2D.o $(OBJ_PATH)/Spatial_Presentation2D.o $(OBJ_PATH)/Network_Generated_Statistics.o $(OBJ_PATH)/Results.o $(OBJ_PATH)/Command_Line_Parameters.o $(OBJ_PATH)/Command.o $(OBJ_PATH)/Spatial_Segment_Subset2D.o $(OBJ_PATH)/Network_Statistics2D.o $(OBJ_PATH)/Connection_Statistics2D.o $(OBJ_PATH)/prepost_structure2D.o $(OBJ_PATH)/fibre_structure2D.o $(OBJ_PATH)/synapse_structure2D.o $(OBJ_PATH)/generic_synapse.o $(OBJ_PATH)/dendritic_growth_model2D.o $(OBJ_PATH)/synapse_formation_model2D.o $(OBJ_PATH)/branching_models2D.o $(OBJ_PATH)/turning_models2D.o $(OBJ_PATH)/axon_direction_model2D.o $(OBJ_PATH)/fibre_elongation_model2D.o $(OBJ_PATH)/neurite_diameter_model2D.o $(OBJ_PATH)/environment_physics2D.o $(OBJ_PATH)/connection2D.o $(OBJ_PATH)/synapse2D.o $(OBJ_PATH)/IFactivity2D.o $(OBJ_PATH)/neuron2D.o $(OBJ_PATH)/network2D.o $(OBJ_PATH)/spatial2D.o $(OBJ_PATH)/event2D.o $(OBJ_PATH)/slice2D.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/file.o $(OBJ_PATH)/state_storable.o $(OBJ_PATH)/global.o $(OBJ_PATH)/diagnostic.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/mtprng.o $(OBJ_PATH)/regex-gnu.o -lgeometry -o $(BIN_PATH)/netmorph2D $(LIBRARIES)
	cd $(BIN_PATH) && ln -f -s netmorph2D nibr2D

wizard: $(OBJ_PATH)/wizard.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/regex-gnu.o
	$(CCPP) $(CPPFLAGS) $(OBJ_PATH)/wizard.o $(OBJ_PATH)/StringList.o $(OBJ_PATH)/BigString.o $(OBJ_PATH)/BigRegex.o $(OBJ_PATH)/regex-gnu.o -o $(BIN_PATH)/wizard $(LIBRARIES)

clean:
	rm -r -f $(BIN_PATH)/* $(OBJ_PATH)/*.o $(LIB_PATH)/$(LIB)
