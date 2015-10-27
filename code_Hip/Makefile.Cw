#CXX        = /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc472p1_x86_64_slc6/slc6/x86_64-slc6-gcc47-opt/bin/g++
#CXX        = /usr/bin/g++
#Enable warnings and debugging info
CXX        = g++ -std=c++0x
#CXX        = g++
CXXFLAGS   = -Wall -g
INCS      = `root-config --cflags` -I.
LIBS      = `root-config --libs` 

CorrectionFactorsMain.exe: CorrectionFactorsMain.o EfficiencyCorrection.o
	$(CXX) $(CXXFLAGS) -o $@ $< EfficiencyCorrection.o $(INCS) $(LIBS)

CorrectionFactorsMain.o: CorrectionFactorsMain.cc EfficiencyCorrection.cc
	$(CXX) $(CXXFLAGS) -c $< $(INCS) $(LIBS)

EfficiencyCorrection.o : EfficiencyCorrection.cc EfficiencyCorrection.h 
	$(CXX) $(CXXFLAGS) -c $< $(INCS) $(LIBS)

Dict.o : Dict.cxx 
		$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS) 

Dict.cxx : EfficiencyCorrection.h 
		rootcint -f $@ -c $<+ 

clean : 
	rm -f *.o
