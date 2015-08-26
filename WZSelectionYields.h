#ifndef WZSelectionYields_h
#define WZSelectionYields_h

#include "GenericAnalysis.h"


class WZSelectionYields : public GenericAnalysis
{

public:

  WZSelectionYields(WZEvent*, TFile* outputFile);

	unsigned int GetNAnalyzed() { return nAnalyzedEvents; }
	unsigned int GetNSelected() { return nSelectedEvents; }

  void Init();

  void EventAnalysis();

  void Finish();


protected:

  unsigned int nAnalyzedEvents;
	unsigned int nSelectedEvents;
  unsigned int yieldsByChannelPreselection[4];
  unsigned int yieldsByChannelZSelection[4];
  unsigned int yieldsByChannelWSelection[4];
  unsigned int yieldsByChannelFullSelection[4];

  std::ofstream eventLists[4];

};

#endif
