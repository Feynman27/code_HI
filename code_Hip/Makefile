#CXX        = /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc472p1_x86_64_slc6/slc6/x86_64-slc6-gcc47-opt/bin/g++
CXX        = g++ -std=c++0x
#Enable warnings and debugging info
CXXFLAGS   = -Wall -g
INCS      = `root-config --cflags` -I.
LIBS      = `root-config --libs` -lRooFit -lRooFitCore -lRooStats

WAnalysisHip : WAnalysisHip.o WAnalysisTool.o WBookKeeper.o 
	$(CXX) $(CXXFLAGS) -o $@ $< WAnalysisTool.o WBookKeeper.o $(INCS) $(LIBS)

WAnalysisHip.o : WAnalysisHip.cc WAnalysisTool.cc WBookKeeper.cc 
	$(CXX) $(CXXFLAGS) -c $< $(INCS) $(LIBS)

WBookKeeper.o : WBookKeeper.cc
	$(CXX) $(CXXFLAGS) -c $< $(INCS) $(LIBS)

WAnalysisTool.o : WAnalysisTool.cc WAnalysisTool.h
	$(CXX) $(CXXFLAGS) -c $< $(INCS) $(LIBS)

clean : 
	rm -f WAnalysisHip *.o
