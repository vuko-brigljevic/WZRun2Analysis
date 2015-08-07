#ifndef WZEvent_h
#define WZEvent_h


// #include "GGTree.h"
#include "EventTree.h"

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
		  


class RecoLepton : public TLorentzVector {

public:

  RecoLepton(double pt, double eta, double phi,
	     float ch, float id ) {
    SetPtEtaPhiM(pt, eta, phi, 0.);
    charge = ch;
    pdgid  = id;
  };

  float PdgId() { return pdgid;};

protected:
  float pdgid;
  float charge;

};



class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree *tree);

  bool passesSelection();

  void ReadEvent();

  static TH1F* hScaleInEB;
  static TH1F* hScaleOutEB;
  static TH1F* hScaleEE;

protected:

  // For various smearing functions

  TRandom3 *random;

  FinalState final_state;
  PassedSelectionStep  selection_level;
  


};


#endif // #ifdef WZ_cxx
