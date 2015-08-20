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

  for (vector<Lepton*>::iterator lIt = leptonsNew.begin(); lIt != leptonsNew.end(); ) {
    delete *lIt;
    lIt = leptonsNew.erase(lIt);
  }

  nEleTightTwiki = 0;
  nMuTightTwiki = 0;
  nZCand = 0;
  massZCand.clear();
}

void WZEvent::ReadEvent()
{ // If you want to fill some own variables, do it here...

  Cleanup();

  // Electrons
  for (unsigned int iele = 0; iele < eleCharge->size(); iele++) {
    Electron* ele = new Electron(iele,
				 (*elePt)[iele],
				 (*eleEta)[iele],
				 (*elePhi)[iele],
				 (*eleCharge)[iele]);
    
    leptons.push_back( ele);
  }

  // Muons
  for (unsigned int imu = 0; imu < muCharge->size(); imu++) {
    leptons.push_back( new Muon(imu,
     			    (*muPt)[imu],
     			    (*muEta)[imu],
     			    (*muPhi)[imu],
     			    (*muCharge)[imu]));
  }

  // Electrons - Sasa
  unsigned int nE = nEle;
  if (nE != eleCharge->size()) {
    std::cout << "nEle != vector<Ele>.size() !!!" << std::endl;
  } else {
    for (int inE = 0; inE < nEle; inE++) {
      Electron* ele = new Electron(inE, elePt->at(inE), eleEta->at(inE),
                                   elePhi->at(inE), eleCharge->at(inE));
      leptonsNew.push_back(ele);
    }
  }

  // Muons - Sasa
  unsigned int nM = nMu;
  if (nM != muCharge->size()) {
    std::cout << "nMu != vector<Mu>.size() !!!" << std::endl;
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

  if (!(nEle + nMu < 3)) {
    selection_level = passesThreeLeptonFilter;
  } else {
    selection_level = failsThreeLeptonFilter;
    return passed;
  }

//  HLT triggers part missing


// Preselection

  std::vector<unsigned int> indexTightLeptons;

  std::vector<Lepton*>::iterator lIt;
  unsigned int i;
  for (lIt = leptonsNew.begin(), i = 0; lIt != leptonsNew.end(); ++lIt, i++) {
    if ((*lIt)->IsTightTwiki()) {
//      index = std::distance(leptonsNew.begin(), lIt);
      indexTightLeptons.push_back(i);
      if ((*lIt)->GetPdgId() == 11) {
        nEleTightTwiki++;
      }
      else if ((*lIt)->GetPdgId() == 13) {
        nMuTightTwiki++;
      }
    }
  }

  if (nMuTightTwiki + nEleTightTwiki != indexTightLeptons.size()) {
    std::cout << "ERROR: nTight != vector<indexTight>.size() !!!" << std::endl;
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
  } else {
    return passed;
  }

// Z Selection
// massZCand window [60, 120] and leading lepton Pt > 20 GeV

  const double mZMin = 60.;
  const double mZMax = 120.;
  std::vector<pair<unsigned int, unsigned int> > indexZCand;
//  const double mZ = 91.11

  for (unsigned int ind1 = 0; ind1 < indexTightLeptons.size(); ind1++) {
    for (unsigned int ind2 = ind1+1; ind2 < indexTightLeptons.size(); ind2++) {
      unsigned int indL1 = indexTightLeptons.at(ind1);
      unsigned int indL2 = indexTightLeptons.at(ind2);

      if (abs(leptonsNew.at(indL1)->GetPdgId()) != abs(leptonsNew.at(indL2)->GetPdgId()) ||
          leptonsNew.at(indL1)->GetCharge() == leptonsNew.at(indL2)->GetCharge()) {
        continue;
      }

      const double mZCand = (*(leptonsNew.at(indL1)) + *(leptonsNew.at(indL2))).M();
      if (mZCand < mZMin || mZCand > mZMax ||
          !(TMath::Max(leptonsNew.at(indL1)->Pt(), leptonsNew.at(indL2)->Pt()) > 20.)) {
        continue;
      } else {
        nZCand++;
        massZCand.push_back(mZCand);
        indexZCand.push_back(make_pair(indL1, indL2));
        selection_level = passesZSelection;
      }

    }
  }

  if (nZCand != massZCand.size()) {
    std::cout << "ERROR: Number of Z candidates != vector<massZCand>.size() !!!" << std::endl;
    return passed;
  }

  if (nZCand > 2) {
    std::cout << "ERROR: More than 2 Z candidates in an event - go back to coding !!!" << std::endl;
    return passed;
  }

  if (nZCand == 0) {
    return passed;
  }

  unsigned int indexZl1, indexZl2;
  if (nZCand == 1) {
    indexZl1 = (indexZCand.at(0)).first;
    indexZl2 = (indexZCand.at(0)).second;
  }

  if (nZCand == 2) {
    if (massZCand.at(0) > massZCand.at(1)) {
      indexZl1 = (indexZCand.at(0)).first;
      indexZl2 = (indexZCand.at(0)).second;
    } else {
      indexZl1 = (indexZCand.at(1)).first;
      indexZl2 = (indexZCand.at(1)).second;
    }
  }

// W selection

  unsigned int indexWl;
  for (unsigned int iWl = 0; iWl < indexTightLeptons.size(); iWl++) {
    unsigned int tempIndexWl = indexTightLeptons.at(iWl);
    if (tempIndexWl != indexZl1 && tempIndexWl != indexZl2) {
      indexWl = tempIndexWl;
      if (leptonsNew.at(indexWl)->Pt() > 20. && pfMET > 30. &&
          leptonsNew.at(indexWl)->DeltaR(*(leptonsNew.at(indexZl1))) > 0.1 &&
          leptonsNew.at(indexWl)->DeltaR(*(leptonsNew.at(indexZl2))) > 0.1) {
        passed = true;
        selection_level = passesWSelection;
      } else {
      return passed;
      }
    }
  }

  return passed;
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

      if ( abs(leptons[ilep1]->GetPdgId()) != abs(leptons[ilep2]->GetPdgId()) )
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

  int zflavor = abs(leptons[izlep1]->GetPdgId());
  int wflavor = abs(leptons[iwlep]->GetPdgId());

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
