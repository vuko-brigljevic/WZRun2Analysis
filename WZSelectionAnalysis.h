#ifndef WZSelectionAnalysis_h
#define WZSelectionAnalysis_h

#include "GenericAnalysis.h"


class WZSelectionAnalysis : public GenericAnalysis
{
public:
  
  WZSelectionAnalysis(WZEvent *, TFile * fout);

  void Init();
  
  void EventAnalysis();
  
  void Finish();
  

protected:

  int nrAnalyzedEvents;

  float yieldsByChannel[4];
  std::ofstream eventLists[4];


};




#endif
