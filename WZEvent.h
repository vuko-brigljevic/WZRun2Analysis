#ifndef WZEvent_h
#define WZEvent_h


#include "EventTree.h"
#include "Leptons.h"


#define WZBASECLASS EventTree


#include "TLorentzVector.h"
#include "TH2F.h"
#include "TRandom3.h"

#include <vector>

#include <iostream>


enum FinalState { undefined,
                  eee,
                  eem,
                  mme,
                  mmm };


enum PassedSelectionStep { selectionNotRun,
                           failsThreeLeptonFilter,
                           passesThreeLeptonFilter,
                           passesPreselection,
                           passesZSelection,
                           passesWSelection,
                           passesFullSelection};


class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree* tree);

  bool passesSelection();

  FinalState GetFinalState() { return final_state; }
  PassedSelectionStep GetSelectionLevel() { return selection_level; };
  bool passesFullSelection();

  void ReadEvent();

  void DumpEvent(std::ostream & out, int verbosity=0);

  static TH1F* hScaleInEB;
  static TH1F* hScaleOutEB;
  static TH1F* hScaleEE;

  unsigned int nZCand;
  std::vector<double> massZCand;

protected:

  void Cleanup();

  // For various smearing functions

  TRandom3 *random;

  FinalState final_state;
  PassedSelectionStep selection_level;

  std::vector<Lepton*> leptons;
  std::vector<Lepton*> leptonsNew;
  unsigned int nEleTightTwiki;
  unsigned int nMuTightTwiki;

  int zLeptonsIndex[2];
  int wLeptonIndex;

};


#endif // #ifdef WZ_cxx
