CXX=/usr/bin/mpicxx
CFLAGS=-c -O2 -W -g
INCLUDE=-I ./include 
LIBPATH=-L ./lib 
OBJECTS=Data FftwToolkit BasicFunc Initialization\
		IntegrationToolkit MDESolver MemFree SCFTBaseAB
all:ab
./lib/lib%.so:./src/%.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(INCLUDE) -lfftw3_mpi -lfftw3 -I/usr/include -L/usr/lib
ab:$(addsuffix .so, $(addprefix ./lib/lib, $(OBJECTS))) ab.cpp
	$(CXX) -g $(INCLUDE) ab.cpp -o $@ $(LIBPATH) $(addprefix -l, $(OBJECTS)) -lfftw3_mpi -lfftw3 -I/usr/include -L/usr/lib
clean: 
	-rm ./lib/*.so
	-rm ab
