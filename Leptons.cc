#include "Leptons.h"

#include "WZEvent.h"

#include "TMath.h"

#include <iostream>


// Initialize static members
WZEvent * Lepton::wztree = 0;


void Lepton::SetWZEvent(WZEvent * wzt) {
  wztree = wzt;
}

Lepton::Lepton(int ind, double pt, double eta, double phi,
	       float ch ) 
{
    SetPtEtaPhiM(pt, eta, phi, 0.);
    charge   = ch;
    index    = ind;
}




Electron::Electron(int ind, double pt, double eta, double phi,
		   float ch )
  : Lepton(ind,pt,eta,phi,ch)
{
  pdgid = 11;
}

bool Electron::IsTight() {

  if (wztree == 0 ) std::cout << "WZEvent pointer is ZERO!!!! \n";

  return true;


}


Muon::Muon(int ind, double pt, double eta, double phi,
	   float ch )
  : Lepton(ind,pt,eta,phi,ch)
{
  pdgid = 13;
}

bool Muon::IsTight() {

  if (wztree == 0 ){ 
    std::cout << "WZEvent pointer is ZERO!!!! \n";
    return false;
  }

  bool passesTight      = false;
  bool passesIsolation  = false;

  if ( (*(wztree->muChi2NDF))[index] < 10 
       && (*(wztree->muMuonHits))[index] > 0 
       && (*(wztree->muStations))[index] > 1 
       && (*(wztree->muD0))[index] < 0.2
       && (*(wztree->muDz))[index] < 0.5
       && (*(wztree->muPixelHits))[index] > 0
       && (*(wztree->muTrkLayers))[index] > 5 ) {

    passesTight = true;

  }

  // From UW Twiki
  // (mu.chargedHadronIso()
  // +max(mu.photonIso()+mu.neutralHadronIso()
  //       -0.5*mu.puChargedHadronIso,0.0))/mu.pt()

  if ( ( (*(wztree->muPFChIso))[index] + 
	 TMath::Max( (*(wztree->muPFPhoIso))[index] + (*(wztree->muPFNeuIso))[index] 
	      - 0.5 *  (*(wztree->muPFPUIso))[index] , 0.) ) 
       / Pt() < 0.12 ) {
    passesIsolation = true;
  }

  if (passesTight && passesIsolation) { 
    return true; 
  } else {
    return false;
  }

}
