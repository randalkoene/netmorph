#!/bin/bash
# Is this NetmorphCMake?
if [ -d ../../../Source ]; then
	cd ../../../Binaries
	./NetmorphCMake "include=../Source/examples/nesvbp/nesvbp-autoassociative" "outattr_directory=../Source/examples/nesvbp/"
	cd ../Source/examples/nesvbp
else
	cd ../../bin
	./netmorph "include=../examples/nesvbp/nesvbp-autoassociative" "outattr_directory=../examples/nesvbp/"
	cd ../examples/nesvbp
fi
