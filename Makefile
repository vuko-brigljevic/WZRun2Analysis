#LDFLAGS=`root-config --libs`
#CPPFLAGS= -Wall -Wno-long-long  -pthread -DCTHREAD_POSIX -D_THREAD_SAFE -D_REENTRANT -I$(ROOTSYS)/include 

CPPFLAGS=`root-config --cflags` -g
#-IRooUnfold-1.1.1/src/   
#If running in CMSSW42...
#LDFLAGS = -L$(ROOTSYS)/lib -lNew -lRint -lTree -lTreePlayer -lCint -lThread -lGraf -lGraf3d -lHist -lHtml -lMatrix -lMinuit -lPostscript -lProof -lThread -lCore -lGX11 -lPhysics -lGpad -lGui -lTreeViewer -L/usr/X11R6/lib -lm -ldl -L/usr/lib -lpthread -rdynamic 

#if running in CMSSW53...
LDFLAGS =$(shell root-config --libs) 
#RooUnfold-1.1.1/libRooUnfold.so


# FOR DATA
#MYCCOPT=-D DATA
#FOR OLDMC
#MYCCOPT2=-M OLDMC
# FOR mc
 MYCCOPT=


#wzDoGenAnalysisNewAutomatic1: wzDoGenAnalysisNewAutomatic1.C wzTools2.C WZ.C WZEvent.C
#	g++ $(CPPFLAGS) $(LDFLAGS) -o $@ $^

wzAnalysis: wzAnalysis.C WZEvent.C EventTree.C
	g++ -D DATA $(CPPFLAGS) $(LDFLAGS) -o $@ $^


#wzMCUnfoldingAnalysis: wzMCUnfoldingAnalysis.C wzToolsNew.C WZGenEvent_v140710.C WZGenEvent.C WZEvent.C UnfoldingAnalysis.C WZAnalysis.C UnfoldingAnalysis.h UnfoldingHistogramFactory.C JetEnergyTool.C SystematicsManager.C MetSystematicsTool.C metsys.C
#	g++ -D NEWMCPUFIX $(CPPFLAGS) $(LDFLAGS) -o $@ $^




