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
}

void WZEvent::ReadEvent()
{ // If you want to fill some own variables, do it here...

  Cleanup();

  selection_level = selectionNotRun;
  final_state = undefined;

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
    selection_level = selectionNotRun;
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









  // Do we have a W candidate

  // MET cut
  return passed;
}


bool WZEvent::passesFullSelection()
{
  std::vector<int> tightLeptons;

  // Do we have exactly 3 tight leptons
  for (unsigned int ilep=0; ilep < leptons.size(); ilep++ ) {
    if (leptons[ilep]->IsTight())
      tightLeptons.push_back(ilep);
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


  if (izlep1 < 0 || izlep2 < 0 || iwlep < 0)
    return false;

  // DR cut between Z and W leptons

  bool passesDRCut = false;

  if ( (*(leptons[izlep1])).DeltaR( (*(leptons[iwlep]))) < 0.1
       &&  (*(leptons[izlep2])).DeltaR( (*(leptons[iwlep]))) < 0.1 ) {
    passesDRCut = true;
  } 


  // MET cut


  // Which final state is it

  if ( leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 11 )
    final_state = eee;
  if ( leptons[izlep1]->PdgId() == 11 && leptons[iwlep]->PdgId() == 13 )
    final_state = eem;
  if ( leptons[izlep1]->PdgId() == 13 && leptons[iwlep]->PdgId() == 11 )
    final_state = mme;
  if ( leptons[izlep1]->PdgId() == 13 && leptons[iwlep]->PdgId() == 11 )
    final_state = mmm;

}
