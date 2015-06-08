#ifndef WZEvent_h
#define WZEvent_h


//Lucija added this to test working with data

#ifdef DATA
#define WZBASECLASS WZ2012Data
#include "WZ2012Data.h"
#endif
#ifdef OLDMC
#define WZBASECLASS WZ
#include "WZ.h"
#endif
#ifdef NEWMC
#define WZBASECLASS WZGenEvent
#include "WZGenEvent.h"
#endif
#ifdef NEWMCPUFIX
#define WZBASECLASS WZGenEvent_v140710
#include "WZGenEvent_v140710.h"
#endif

#define _ISWZMC_


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

  float GetScaleFactor();
  float LeadTriggerEff();
  float TrailTriggerEff();

  float PdgId() { return pdgid;};

protected:
  float pdgid;
  float charge;

  static TH2F * MuonSF;
  static TH2F * ElecSF;

  // Trigger efficiencies
  static TH2F * DoubleMuLeadEff;
  static TH2F * DoubleMuTrailEff;
  static TH2F * DoubleEleLeadEff;
  static TH2F * DoubleEleTrailEff;

};


class GenS1Lepton : public TLorentzVector {

public:

  GenS1Lepton(double pt, double eta, double phi,
	      float id=0, float oVpid=-9999., 
	      float imTau = -9999. ) {

    mass = 0;
    pdgid  = id;    
    if (abs(id)==11) {
      mass = 0.000511;
    } else if (abs(id) == 13) {
      mass = 0.1057;
    }

    SetPtEtaPhiM(pt, eta, phi, mass);

    motherBoson     = oVpid;
    isTauDescendent = imTau;

  };

  int Id() { return pdgid; };
  int MotherBoson() { return motherBoson; };
  int ComesFromTau() { return isTauDescendent; };


protected:
  int   pdgid;
  int   motherBoson;
  int   isTauDescendent;
  float mass;
};



class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree *tree);

  bool passesSelection();

  bool passesNewSelection();

  void ReadEvent();

  void DumpEvent(std::ostream & out, int verboseLevel = 0);

  void PrintSummary();

  bool PassesGenCuts();

  void SmearJets();

  float GetShiftedMET(int met_syst);

  void ApplyJESCorrection(double strength = 0.);
  void PrintNjetsLeadingJetPt();

  float LeptonPt(int i);
  float LeptonBDT(int i);
  float LeptonCharge(int i);
  float LeptonEta(int i);
  float LeptonPhi(int i);
  int   LeptonPass2012ICHEP(int i);
  float LeptonIso(int i);
  float LeptonIsomva(int i);

  float GetPileupWeight();
  float GetMCWeight();
  float GetTriggerEfficiency();

  float GetBrWeight();  
  FinalState GetFinalState(){return final_state;}
  PassedSelectionStep GetSelectionLevel() { return selection_level;}

  static TH1F* hScaleInEB;
  static TH1F* hScaleOutEB;
  static TH1F* hScaleEE;

#ifdef _ISWZMC_
  std::vector<GenS1Lepton> genLeptons;
  std::vector<int> genLeptonsIds;
  std::vector<TLorentzVector> genJets;
#endif
  std::vector<TLorentzVector> recoJets;

  std::vector<RecoLepton>  leptons;

  float NumZ() { return numZ;}
  float NumW() { return numW;}

  double SelectedZPt() {
    return selectedZPt;
  }

  TLorentzVector SelectedZP4() {
    return selectedZP4;
  }


  TLorentzVector * ZLepton (unsigned int index);
  TLorentzVector * WLepton();


protected:

  std::vector<float*> pt;
  std::vector<float*> bdt;
  std::vector<float*> pdgid;
  std::vector<float*> ch;
  std::vector<float*> eta;
  std::vector<float*> phi;
  std::vector<int*> pass2012ICHEP;
  std::vector<float*> iso;
  std::vector<float*> isomva;

  FinalState          final_state;
  PassedSelectionStep selection_level;

  double selectedZPt;
  TLorentzVector selectedZP4;

  int  wLeptonIndex;
  int  zLeptonsIndex[2];

  float numZ;
  float numW;
  float numMET;
  float num3e;
  float num2e1mu;
  float num1e2mu;
  float num3mu;
  float numMET3e;
  float numMET2e1mu;
  float numMET1e2mu;
  float numMET3mu;


  // For various smearing functions

  TRandom3 *random;


};


#endif // #ifdef WZ_cxx
