#ifndef WZEvent_h
#define WZEvent_h

#define WZBASECLASS EventTree_ggNtuplizer_V07_04_09_00

#include "EventTree_ggNtuplizer_V07_04_09_00.h"
#include "Leptons.h"

#include "TH1F.h"
#include "TRandom3.h"

#include <vector>
#include <iostream>


enum FinalState
{
  undefined,   // 0
  eee,         // 1
  eem,         // 2
  mme,         // 3
  mmm          // 4
};


enum SelectionLevel
{
  Undefined,            // 0
  FailsPreselection,    // 1
  Preselection,         // 2
  ZSelection,           // 3
  WSelection,           // 4
  FullSelection         // 5
};


class WZEvent : public WZBASECLASS
{

  friend class WZSelectionYields;
  friend class WZJetStudy;

public:

  WZEvent(TTree* tree);

  vector<bool> GetHLT25ns() { return fHLT25ns; }
  vector<bool> GetHLT50ns() { return fHLT50ns; }

  bool PassesPreselection();
  bool PassesZSelection();
  bool PassesWSelection();
  bool PassesFullSelection();

  FinalState GetFinalState() { return fFinalState; }
  SelectionLevel GetSelectionLevel() { return fSelectionLevel; }

  Lepton* GetWLepton() { return fLeptons.at(fWLeptonIndex); }
  pair<Lepton*, Lepton*> GetZLeptons();

  void ReadEvent();

  void DumpEvent(std::ostream& out, int verbosity=0);

  static TH1F* hScaleInEB;
  static TH1F* hScaleOutEB;
  static TH1F* hScaleEE;


protected:

  void Clear();

  // For various smearing functions
  TRandom3* fRandom;

  vector<bool> fHLT25ns;
// 0 - HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2         - bit  8
// 1 - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2               - bit 20
// 2 - HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2             - bit 11
// 3 - HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v2    - bit 41
// 4 - HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2   - bit 42
// 5 - HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v2             - bit  9
// 6 - HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v2                  - bit 43
// 7 - HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v2                   - bit 44
// 8 - HLT_TripleMu_12_10_5_v2                              - bit 28

  vector<bool> fHLT50ns;
// 0 - HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v2         - bit 16
// 1 - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v2               - bit 10
// 2 - HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v2             - bit 11
// 3 - HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v2    - bit 36
// 4 - HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v2   - bit 37
// 5 - HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v2             - bit 17
// 6 - HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v2                  - bit 38
// 7 - HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v2                   - bit 39
// 8 - HLT_TripleMu_12_10_5_v2                              - bit 35

  FinalState fFinalState;
  SelectionLevel fSelectionLevel;

  vector<Lepton*> fLeptons;
  vector<unsigned int> fTightLeptonsIndex;
  pair<unsigned int, unsigned int> fZLeptonsIndex;
  unsigned int fWLeptonIndex;

};

#endif
