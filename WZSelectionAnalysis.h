#ifndef WZSelectionAnalysis_h
#define WZSelectionAnalysis_h

#include "GenericAnalysis.h"


class WZSelectionAnalysis : public GenericAnalysis
{

public:

  WZSelectionAnalysis(WZEvent *, TFile * fout);

  void Init();

  void EventAnalysis();

  //  void FindTrueWZLeptons();

  void Finish();


protected:

  int nrAnalyzedEvents;

  unsigned int yieldsByChannel[6];

  unsigned int signalYields[6];

  unsigned int genYieldsByChannel[6];

  // unsigned int yieldsByChannelPreselection[6];
  // unsigned int yieldsByChannelZSelection[6];
  // unsigned int yieldsByChannelWSelection[6];
  unsigned int yieldsByChannelFullSelection[6];

  std::ofstream eventLists[6];

  TH2D * hTrueDecayMode;

  TH2D * hRecoVsTrueChannel;

  GenParticle * trueWLepton;
  GenParticle * trueZLeptons[2];

};

#endif
