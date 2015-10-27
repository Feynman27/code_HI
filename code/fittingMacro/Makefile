#CC        = /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64/bin/g++
CC        = /usr/bin/g++
CCFLAGS   = -Wall -g
INCS      = `root-config --cflags` -I.
LIBS      = `root-config --libs` -lRooFit -lRooFitCore -lRooStats

WAnalysis_DPD_2011Run : WAnalysis_DPD_2011Run.C flavourFitterHIDep_DPD_2011Run.o RooShiftedKeysPdf.o Dict.o
	$(CC) $(CCFLAGS) $< RooShiftedKeysPdf.o Dict.o -o $@ $(INCS) $(LIBS)

#hfDataDrivenFitter : hfDataDrivenFitter.C flavourFitterHIDep_DPD_2011Run.C RooShiftedKeysPdf.o Dict.o
#	$(CC) $(CCFLAGS) $< RooShiftedKeysPdf.o Dict.o -o $@ $(INCS) $(LIBS)

flavourFitterHIv2 : flavourFitterHIv2.C flavourFitterHIDep_DPD_2011Run.o RooShiftedKeysPdf.o Dict.o
	$(CC) $(CCFLAGS) $< RooShiftedKeysPdf.o Dict.o -o $@ $(INCS) $(LIBS)

flavourFitterHIDep_DPD_2011Run : flavourFitterHIDep_DPD_2011Run.C RooShiftedKeysPdf.o Dict.o
	$(CC) $(CCFLAGS) $< RooShiftedKeysPdf.o Dict.o -o $@ $(INCS) $(LIBS)


flavourFitterHIDep_DPD_2011Run.o : flavourFitterHIDep_DPD_2011Run.C
	$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS)

RooShiftedKeysPdf.o : RooShiftedKeysPdf.cxx RooShiftedKeysPdf.h
	$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS)

#hfDataDrivenFitter.o : hfDataDrivenFitter.C
#	$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS)

#AtlasUtils.o : ${HOME}/atlasstyle-00-03-04/AtlasUtils.C 
AtlasUtils.o : ${TestArea}/HeavyIonSingleMuons2011/code/AtlasUtils.C 
	$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS)

Dict.o : Dict.cxx
	$(CC) $(CCFLAGS) -c $< $(INCS) $(LIBS)

Dict.cxx : RooShiftedKeysPdf.h
	rootcint -f $@ -c $<+

clean : 
	rm -f flavourFitterHIDep_DPD_2011Run RooShiftedKeysPdf.o hfDataDrivenFitter Dict.o Dict.cxx Dict.h
	rm -f flavourFitterHIDep_DPD_2011Run RooShiftedKeysPdf.o Dict.o Dict.cxx Dict.h
# 	rm -f flavourFitterHIDep_DPD_2011Run RooShiftedKeysPdf.o Dict.o Dict.cxx Dict.h
