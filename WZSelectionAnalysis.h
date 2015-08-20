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

  float yieldsByChannel[4];


};




#endif
