#include "WZEvent.h"

//#include "JetEnergyTool.h"
//#include "MetSystematicsTool.h"
//#include "SystematicsManager.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <math.h>



WZEvent::WZEvent(TTree * tree) :
  WZBASECLASS(tree)
  //
{

  // Set static pointer to Lepton class
  Lepton::SetWZEvent(this);

  
}


void WZEvent::Cleanup() {

  final_state     = undefined;
  selection_level = selectionNotRun;

  // Empty leptons
  std::vector<Lepton*>::iterator it;
  for ( it = leptons.begin(); it != leptons.end(); ) {
    delete * it;  
    it = leptons.erase(it);
  }

}

void WZEvent::ReadEvent()
{ // If you want to fill some own variables, do it here...

  Cleanup();


  //  selection_level = failsSelection;


  // Electrons
  for (int iele=0; iele<eleCharge->size(); iele++) {
    Electron* ele = new Electron(iele,
				 (*elePt)[iele],
				 (*eleEta)[iele],
				 (*elePhi)[iele],
				 (*eleCharge)[iele]);
    
    leptons.push_back( ele);
  }

  // Muons
  for (int imu=0; imu<muCharge->size(); imu++) {
    leptons.push_back( new Muon(imu,
     			    (*muPt)[imu],
     			    (*muEta)[imu],
     			    (*muPhi)[imu],
     			    (*muCharge)[imu]));
  }



}






bool WZEvent::passesSelection(){

  bool passed = false;

  // Do we have a Z decay ? 

  // Do we have a W candidate

  // MET cut

  if ( (nEle + nMu) == 3)
    passed = true;

  return passed;
  
}


bool WZEvent::passesFullSelection(){


  std::vector<int> tightLeptons;

  // Do we have exactly 3 tight leptons
  for (int ilep=0; ilep<leptons.size(); ilep++ ) {
    
    if (leptons[ilep]->IsTight()
	&& leptons[ilep]->Pt() > 10.) {
      tightLeptons.push_back(ilep);
    }

  }

  if (tightLeptons.size() != 3) return false;


  // Look for Z candidates

  int izlep1=-1;
  int izlep2=-1;
  int iwlep = -1;

  float dzmin = 1000.;

  //
  // Still 
  // 

  
  for (int id1=0; id1<tightLeptons.size(); id1++ ) {
    for (int id2=id1+1; id2<tightLeptons.size(); id2++ ) {
      // Same flavor, opposite charge
      int ilep1 = tightLeptons[id1];
      int ilep2 = tightLeptons[id2];

      if ( abs(leptons[ilep1]->PdgId()) != abs(leptons[ilep2]->PdgId()) ) continue;

      // 
      if ( leptons[ilep1]->Pt()<10. || leptons[ilep2]->Pt()<10.) continue;
      if( leptons[ilep1]->Pt()<20. && leptons[ilep2]->Pt()<20.) continue;
      

      float mcand = (*(leptons[ilep1]) + *(leptons[ilep2])).M();
      std::cout << "Z candidate mass = " << mcand << std::endl;
      if (fabs(mcand - 91.11) < dzmin ) {
	dzmin = fabs(mcand - 91.11);
	izlep1 = ilep1;
	izlep2 = ilep2;
      }
    }
  }


  // MISSING: Pt cuts on W & Z leptons
  // Opposite charge requirements

  // Third lepton 

  float ptlw = -10.;
  int nwlcand = 0;
  for (int id=0; id<tightLeptons.size(); id++ ) {
    int ilep = tightLeptons[id];
    if (ilep != izlep1 && ilep != izlep2) {
      float ptl = leptons[ilep]->Pt();
      if (ptl > 20.) {
	nwlcand++;
	if (ptl > ptlw) {
	  ptlw = ptl;
	  iwlep = ilep;
	}
      }
    }
  }


  if (nwlcand > 1) {
    std::cout << "SHOULD NOT BE: More than one W candidate lepton \n";
  }


  if (izlep1<0 || izlep2 < 0 || iwlep < 0) return false;

  // DR cut between Z and W leptons

  bool passesDRCut = false;

  if (  (*(leptons[izlep1])).DeltaR( (*(leptons[iwlep]))) < 0.1
	&&  (*(leptons[izlep2])).DeltaR( (*(leptons[iwlep]))) < 0.1 ) {

    passesDRCut = true;
  } 

  bool passesMET = false;


  // MET cut
  if (pfMET > 30) passesMET = true;

  // Which final state is it

  if (leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 11 ) 
    final_state = eee;
  if (leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 13 ) 
    final_state = eem;
  if (leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 11 ) 
    final_state = mme;
  if (leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 11 ) 
    final_state = mmm;


}


