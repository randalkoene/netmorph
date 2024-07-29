#!/bin/bash

EXTRA=
BLENDER=

if [ "$1" = "blend" ]; then
	BLENDER="$HOME/Downloads/blender-4.1.1-linux-x64/blender"
	if [ -f "$BLENDER" ]; then
		echo "Blender at $BLENDER"
	else
		BLENDER="$HOME/blender-4.1.1-linux-x64/blender"
		echo "Blender at $BLENDER"
	fi
	BLENDER="blender_exec_path=$BLENDER"
	EXTRA="include=../Source/examples/nesvbp/obj-and-blend"
fi

# Is this NetmorphCMake?
if [ -d ../../../Source ]; then
	cd ../../../Binaries
	./NetmorphCMake "include=../Source/examples/nesvbp/nesvbp-xor-res-sep-targets" "outattr_directory=../Source/examples/nesvbp/" $EXTRA $BLENDER
	cd ../Source/examples/nesvbp
else
	cd ../../bin
	./netmorph "include=../examples/nesvbp/nesvbp-xor-res-sep-targets" "outattr_directory=../examples/nesvbp/" $EXTRA $BLENDER
	cd ../examples/nesvbp
fi
