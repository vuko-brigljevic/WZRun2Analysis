#ifndef WZEvent_h
#define WZEvent_h


#include "EventTree.h"
#include "Leptons.h"


#define WZBASECLASS EventTree


#include "TLorentzVector.h"
#include "TH2F.h"
#include "TRandom3.h"

#include <vector>


enum FinalState { undefined,
                  eee,
                  eem,
                  mme,
                  mmm};


enum PassedSelectionStep { selectionNotRun,
                           passesThreeLeptonFilter,
                           passesPreselection,
                           
			   failsSelection,
			   passesZSelection,
			   passesWSelection,
			   passesFullSelection};


class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree *tree);

  bool passesSelection();
  FinalState GetFinalState() { return final_state; };
  PassedSelectionStep GetSelectionLevel() { return selection_level; };

  bool passesFullSelection();

  void ReadEvent();

  static TH1F* hScaleInEB;
  static TH1F* hScaleOutEB;
  static TH1F* hScaleEE;

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
};


#endif // #ifdef WZ_cxx
