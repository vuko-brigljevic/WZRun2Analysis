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
  unsigned int yieldsByChannelPreselection[5];
  unsigned int yieldsByChannelZSelection[5];
  unsigned int yieldsByChannelWSelection[5];
  unsigned int yieldsByChannelFullSelection[5];

  std::ofstream eventLists[4];

  TH1D* hZmass[5];
  TH1D* hZpt[5];
  TH1D* hMET[5];
  TH1D* hMt[5];
  TH1D* hZl1pt[5];
  TH1D* hZl2pt[5];
  TH1D* hWlpt[5];
  TH1D* hNJets[5];
  TH1D* hNJetsNoMuIso[5];
  TH1D* hNJetsNoEleIso[5];
  TH1D* hNJetsNoIso[5];

};

#endif
