# USER CONFIGURATION
#
# CXX - c++ compiler. specify the full path if necessary
#
# GSLDIR - prefix directory for the gsl library, it is assumed
#          that the library is in GSLDIR/lib and the include
#          files in GSLDIR/include
#
# LHADIR - prefix directory for the Les Houches pdf library, it is
#          assumed that the library is in LHADIR/lib and the
#          include files are in LHADIR/include/LHAPDF

CXX = g++
GSLDIR = $(shell gsl-config --prefix)
LHADIR = $(shell lhapdf-config --prefix)

# DO NOT MODIFY

CXXFLAGS = -O -I$(GSLDIR)/include -I$(LHADIR)/include
LD = $(CXX)
LDFLAGS =
LIBS = -L$(GSLDIR)/lib -lgsl -lgslcblas -lm -L$(LHADIR)/lib -lLHAPDF
SRCS = Utilities.cpp FixedOrder.cpp \
  PartonicFlux.cpp Resummation.cpp \
  SubtrFlux.cpp lgamma.cpp top++.cpp psin.cpp
OBJS = $(patsubst %.cpp,%.o,$(SRCS))

all: top++

top++: $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f top++ *.o *.res *~

# DO NOT DELETE

Utilities.o: Utilities.h PartonicFlux.h SubtrFlux.h
Utilities.o: lgamma.h FixedOrder.h Resummation.h psin.h
FixedOrder.o: FixedOrder.h
PartonicFlux.o: PartonicFlux.h
Resummation.o: Resummation.h psin.h lgamma.h SubtrFlux.h PartonicFlux.h
SubtrFlux.o: SubtrFlux.h PartonicFlux.h lgamma.h
lgamma.o: lgamma.h
top++.o: PartonicFlux.h SubtrFlux.h lgamma.h FixedOrder.h Resummation.h psin.h
psin.o: psin.h
