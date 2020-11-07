# Makefile for geometry classes
# Randal A. Koene, 20050112

CC=gcc
CCPP=g++ 

#MACHOPT=
MACHOPT=-march=pentium4 -mtune=pentium4 -mfpmath=sse

GCC3FXTRA=$(if $(shell gcc -v 2>&1 | grep "^gcc version[^0-9]*[3-9][.]"),-Wno-unused-function)

CFXTRA= -DSTDC_HEADERS -pedantic -Wall -Werror -Wno-char-subscripts $(GCC3FXTRA)
CPPXTRA= -D_ALT_REGEX_H=$(ALT_REGEX_H) -D_USE_ALT_REGEX -D_CPP_REGEX $(MACHOPT) -mieee-fp -ffast-math -pedantic -Wall -Werror $(GCC3FXTRA)

COMPOPT= -O3
CPPOPT= -felide-constructors

CFLAGS= $(COMPOPT) $(CFXTRA)
CPPFLAGS= $(COMPOPT) $(CPPOPT) $(CPPXTRA)

#.c.o:
#	g++ -c -g $<

#HDR = common.h point.h vector.h
#SRC = point.c vector.c
#OBJ = point.o vector.o
LIB = libgeometry.a
TAR = geometry.tar

all: $(LIB)

Point3D.o: Point3D.cc Point3D.hh Vector3D.hh
	$(CCPP) $(CPPFLAGS) -c Point3D.cc -o Point3D.o

Vector3D.o: Vector3D.cc Vector3D.hh Point3D.hh
	$(CCPP) $(CPPFLAGS) -c Vector3D.cc -o Vector3D.o

Point2D.o: Point2D.cc Point2D.hh Vector2D.hh
	$(CCPP) $(CPPFLAGS) -c Point2D.cc -o Point2D.o

Vector2D.o: Vector2D.cc Vector2D.hh Point2D.hh
	$(CCPP) $(CPPFLAGS) -c Vector2D.cc -o Vector2D.o

ptest: ptest.c $(LIB)
	g++ -o ptest ptest.c $(LIB)

vtest: vtest.c $(LIB)
	g++ -o vtest vtest.c $(LIB)

$(LIB): Point3D.o Vector3D.o Point2D.o Vector2D.o
	ar rcsv $(LIB) Point3D.o Vector3D.o Point2D.o Vector2D.o

tar: $(TAR)
$(TAR): Readme Makefile *.hh *.cc
	tar cvf $(TAR) Readme Makefile *.hh *.cc

bare: clean
	rm -f $(LIB)
clean:
	rm -f *.o $(TAR)
