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
{
  // Set static pointer to Lepton class
  Lepton::SetWZEvent(this);
}


void WZEvent::Cleanup()
{
  final_state     = undefined;
  selection_level = selectionNotRun;

  // Empty leptons
  std::vector<Lepton*>::iterator it;
  for ( it = leptons.begin(); it != leptons.end(); ) {
    delete * it;  
    it = leptons.erase(it);
  }

  for (vector<Lepton*>::const_iterator lIt = leptonsNew.begin(), end = leptonsNew.end();
       lIt != end; ++lIt) {
    delete *lIt;
    lIt = leptonsNew.erase(lIt);
  }
}

void WZEvent::ReadEvent()
{ // If you want to fill some own variables, do it here...

  Cleanup();

  selection_level = selectionNotRun;

  // Electrons
  for (unsigned int iele=0; iele < eleCharge->size(); iele++) {
    Electron* ele = new Electron(iele,
				 (*elePt)[iele],
				 (*eleEta)[iele],
				 (*elePhi)[iele],
				 (*eleCharge)[iele]);
    
    leptons.push_back( ele);
  }

  // Muons
  for (unsigned int imu=0; imu < muCharge->size(); imu++) {
    leptons.push_back( new Muon(imu,
     			    (*muPt)[imu],
     			    (*muEta)[imu],
     			    (*muPhi)[imu],
     			    (*muCharge)[imu]));
  }

  // Electrons - Sasa
  if (nEle != eleCharge->size()) {
    std::cout << "nEle != vector<ele>.size() !!!" << std::endl;
  } else {
    for (int inE = 0; inE < nEle; inE++) {
      Electron* ele = new Electron(inE, elePt->at(inE), eleEta->at(inE),
                                   elePhi->at(inE), eleCharge->at(inE));
      leptonsNew.push_back(ele);
    }
  }

  // Muons - Sasa
  if (nMu != muCharge->size()) {
    std::cout << "nMu != vector<ele>.size() !!!" << std::endl;
  } else {
    for (int inMu = 0; inMu < nMu; inMu++) {
      Muon* mu = new Muon(inMu, muPt->at(inMu), muEta->at(inMu), muPhi->at(inMu), muCharge->at(inMu));
      leptonsNew.push_back(mu);
    }
  }

}


bool WZEvent::passesSelection()
{
  bool passed = false;

  if (nEle + nMu < 3) {
    selection_level = failsThreeLeptonFilter;
    return passed;
  }

//  HLT triggers part missing


// Preselection

  std::vector<unsigned int> indexTightLeptons;
  nEleTightTwiki = 0;
  nMuTightTwiki = 0;

  std::vector<Lepton*>::iterator lIt;
  unsigned int i;
  for (lIt = leptonsNew.begin(), i = 0; lIt != leptonsNew.end(); ++lIt, i++) {
    if ((*lIt)->IsTightTwiki()) {
//      index = std::distance(leptonsNew.begin(), lIt);
      indexTightLeptons.push_back(i);
      if ((*lIt)->PdgId() == 11) {
        nEleTightTwiki++;
      }
      else if ((*lIt)->PdgId() == 13) {
        nMuTightTwiki++;
      }
    }
  }

  if (nMuTightTwiki + nEleTightTwiki != indexTightLeptons.size()) {
    std::cout << "ERROR: nTight != vector<index>.size() !!!" << std::endl;
    return passed;
  }

  if (nMuTightTwiki + nEleTightTwiki == indexTightLeptons.size() &&
      nMuTightTwiki + nEleTightTwiki == 3) {
    selection_level = passesPreselection;
    if (nMuTightTwiki == 0 && nEleTightTwiki == 3) {
      final_state = eee;
    }
    if (nMuTightTwiki == 1 && nEleTightTwiki == 2) {
      final_state = eem;
    }
    if (nMuTightTwiki == 2 && nEleTightTwiki == 1) {
      final_state = mme;
    }
    if (nMuTightTwiki == 3 && nEleTightTwiki == 0) {
      final_state = mmm;
    }
  }
  
  // Do we have a Z decay ? 

  // Do we have a W candidate

  // MET cut
  
}


bool WZEvent::passesFullSelection()
{

  bool passed = true;


  std::vector<int> tightLeptons;

  // Do we have exactly 3 tight leptons
  for (unsigned int ilep=0; ilep<leptons.size(); ilep++ ) {
    
    if (leptons[ilep]->IsTight()
	&& leptons[ilep]->Pt() > 10.) {
      tightLeptons.push_back(ilep);
    }

  }

  if (tightLeptons.size() != 3)
    return false;


  // Look for Z candidates

  int izlep1 = -1;
  int izlep2 = -1;
  int iwlep = -1;

  float dzmin = 1000.;

  //
  // Still 
  // 

  for (unsigned int id1=0; id1 < tightLeptons.size(); id1++ ) {
    for (unsigned int id2=id1+1; id2 < tightLeptons.size(); id2++ ) {
      // Same flavor, opposite charge
      int ilep1 = tightLeptons[id1];
      int ilep2 = tightLeptons[id2];

      if ( abs(leptons[ilep1]->PdgId()) != abs(leptons[ilep2]->PdgId()) )
        continue;

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


  if (izlep1 < 0 || izlep2 < 0 || iwlep < 0)
    return false;

  // DR cut between Z and W leptons

  bool passesDRCut = false;

  if ( (*(leptons[izlep1])).DeltaR( (*(leptons[iwlep]))) < 0.1
       &&  (*(leptons[izlep2])).DeltaR( (*(leptons[iwlep]))) < 0.1 ) {
    passesDRCut = true;
  } 

  bool passesMET = false;


  // MET cut
  if (pfMET > 30) passesMET = true;

  if (! (passesMET && passesDRCut) ) return false;


  // Which final state is it

  int zflavor = abs(leptons[izlep1]->PdgId());
  int wflavor = abs(leptons[iwlep]->PdgId());

  if (zflavor == 11) {
    if (wflavor == 11) {
      final_state = eee;
    }  else  if (wflavor == 13) {
      final_state = eem;
    }
  } else  if (zflavor == 13) {
    if (wflavor == 11) {
      final_state = mme;
    }  else  if (wflavor == 13) {
      final_state = mmm;
    }
  }

  return true;


}
