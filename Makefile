# DO NOT EDIT!!! Changes will be lost. Modify driver.pl instead
# This Makefile was written by driver.pl
# SPF v7 by G.A. and C.S.
# Platform: linux
# Model: PnictidesTwoOrbitals
#LDFLAGS = -L. /usr/lib64/liblapack.so.3.0 /usr/lib64/libblas.so.3.0 -lm -L../lib
LDFLAGS = -L.   -llapack -lblas -lm -L../lib
EXENAME = hf
CPPFLAGS = -DNDEBUG  -IPartialPsimag -IProgram
CXX = g++ -Werror -Wall -g3 -O3

$(EXENAME): clean main.o 
	$(CXX) -o $(EXENAME) main.o $(LDFLAGS) 


clean:
	rm -f core* $(EXENAME) *.o *.ii *.tt

######## End of Makefile ########

