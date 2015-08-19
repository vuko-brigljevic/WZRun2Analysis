#ifndef WZEvent_h
#define WZEvent_h


// #include "GGTree.h"
#include "EventTree.h"
#include "Leptons.h"



#define WZBASECLASS EventTree

//Lucija commented it
//#define WZBASECLASS WZGenEvent


#include "TLorentzVector.h"
#include "TH2F.h"
#include "TRandom3.h"

#include <vector>



// #include "WZGenNew.h"
//Lucija commented it
//#include "WZGenEvent.h"

enum FinalState { undefined,
		  eee,
		  eem,
		  mme,
		  mmm};

enum PassedSelectionStep { selectionNotRun,
			   failsSelection,
			   passesZSelection,
			   passesWSelection,
			   passesFullSelection};
		  



class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree *tree);

  bool passesSelection();

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
  PassedSelectionStep  selection_level;
  

  std::vector<Lepton* > leptons;


};


#endif // #ifdef WZ_cxx
