# makefile for creating vector and matrix library
 CC 			= gcc
CXX 			= g++
FC			= gfortran

LIBS			= -lm -lgsl -lgslcblas

CXXFLAGS		= 

IFLAGS			= -I${GSL_HOME}/include/
LFLAGS			= -L${GSL_HOME}/lib/

#CPPFLAGS		= -W -Wall -O3 $(IFLAGS)
#CPPFLAGS		= -W -Wall -pg -O0 $(IFLAGS)
CPPFLAGS		= -W -Wall -g $(IFLAGS)

all:			runHLasso .depend

clean:		
			rm -f *.o
			rm runHLasso
			rm -f .depend

infiles			= runHLasso.o HLasso.o CLG.o logreg.o paraboliccylinder.o InputData.o StringSplitter.o StringConvertor.o rand.o

runHLasso: 	$(infiles)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o runHLasso $(infiles) $(LIBS)

dep .depend:
		$(CXX) $(CPPFLAGS) -MM *.cc >>.depend
