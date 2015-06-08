#LDFLAGS=`root-config --libs`
#CPPFLAGS= -Wall -Wno-long-long  -pthread -DCTHREAD_POSIX -D_THREAD_SAFE -D_REENTRANT -I$(ROOTSYS)/include 

CPPFLAGS=`root-config --cflags` -IRooUnfold-1.1.1/src/   -g 
#If running in CMSSW42...
#LDFLAGS = -L$(ROOTSYS)/lib -lNew -lRint -lTree -lTreePlayer -lCint -lThread -lGraf -lGraf3d -lHist -lHtml -lMatrix -lMinuit -lPostscript -lProof -lThread -lCore -lGX11 -lPhysics -lGpad -lGui -lTreeViewer -L/usr/X11R6/lib -lm -ldl -L/usr/lib -lpthread -rdynamic 

#if running in CMSSW53...
LDFLAGS =$(shell root-config --libs) RooUnfold-1.1.1/libRooUnfold.so


# FOR DATA
#MYCCOPT=-D DATA
#FOR OLDMC
#MYCCOPT2=-M OLDMC
# FOR mc
 MYCCOPT=


#wzDoGenAnalysisNewAutomatic1: wzDoGenAnalysisNewAutomatic1.C wzTools2.C WZ.C WZEvent.C
#	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzAnalysisData: wzAnalysisData.C wzToolsNew.C WZ2012Data.C WZEventMCOld.C HistogramFactory.C
	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzAnalysisDataWithMM: wzAnalysisDataWithMM.C wzToolsNew.C WZ2012Data.C UnfoldingHistogramFactory.C HistogramFactory.C WZEventMCOld.C
	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzMatrixMethod: wzMatrixMethod.C wzToolsNew.C WZ2012Data.C UnfoldingHistogramFactory.C WZEventMCOld.C
	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzAnalysisMCAll: wzAnalysisMCAll.C wzToolsNew.C WZ.C UnfoldingHistogramFactory.C HistogramFactory.C WZEventMCOld.C
	g++ -D OLDMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzBackground: wzBackground.C wzToolsNew.C WZ.C UnfoldingHistogramFactory.C HistogramFactory.C WZEventMCOld.C
	g++ -D OLDMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzAnalysisMC: wzAnalysisMC.C wzToolsNew.C WZGenEvent.C WZEvent.C UnfoldingHistogramFactory.C HistogramFactory.C
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzWZ: wzWZ.C wzToolsNew.C WZGenEvent.C WZEvent.C UnfoldingHistogramFactory.C
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzMCSignalAnalysis: wzMCSignalAnalysis.C wzToolsNew.C WZ.C WZGenEvent.C WZEvent.C UnfoldingAnalysis.C  UnfoldingAnalysis.h
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzMCUnfoldingAnalysis: wzMCUnfoldingAnalysis.C wzToolsNew.C WZGenEvent_v140710.C WZGenEvent.C WZEvent.C UnfoldingAnalysis.C WZAnalysis.C UnfoldingAnalysis.h UnfoldingHistogramFactory.C JetEnergyTool.C SystematicsManager.C MetSystematicsTool.C metsys.C
	g++ -D NEWMCPUFIX $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzExampleAnalysis: wzExampleAnalysis.C wzToolsNew.C WZGenEvent.C WZEvent.C WZAnalysis.C 
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

unfold: unfold.C 
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzDataUnfold: wzDataUnfold.C SystematicsManager.C 
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

BLUE_unfolding: BLUE_unfolding.C UnfoldingHistogramFactory.C 
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

Systematics: Systematics.C UnfoldingHistogramFactory.C 
	g++ -D NEWMC $(CPPFLAGS) $(LDFLAGS) -o $@ $^

mcfmPlots: mcfmPlots.C mcfmTree.C UnfoldingHistogramFactory.C 
	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^



#test: test.C wzTools2.C WZ2012Data.C
#	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^

#wzMatrixMethod: wzMatrixMethod.C wzTools2.C WZ2012Data.C
#	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^




