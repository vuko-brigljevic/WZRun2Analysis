#include "Leptons.h"
#include "WZEvent.h"
#include "TMath.h"
#include <iostream>


// Initialize static members
WZEvent* Lepton::wztree = 0;


void Lepton::SetWZEvent(WZEvent* wzt) {
  wztree = wzt;
}

Lepton::Lepton(int ind, double pt, double eta, double phi, float ch) 
{
    SetPtEtaPhiM(pt, eta, phi, 0.);
    charge   = ch;
    index    = ind;
}


Electron::Electron(int ind, double pt, double eta, double phi, float ch)
  : Lepton(ind, pt, eta, phi, ch)
{
  pdgid = 11;
}

bool Electron::IsTight() {

  if (wztree == 0 ) {
    std::cout << "WZEvent pointer is ZERO!!!! \n";
    return false;
  }

  bool passesTight = false;
//  bool passesPtEtaCuts = false;
  
  if ( (wztree->eleIDbit)->at(index) == 2 ) {
    passesTight = true;
  }
  
//  if ( (wztree->elePt)->at(index) > 10 && (wztree->eleEta)->at(index) < 2.5 ) {
//    passesPtEtaCuts = true;
//  }

  if ( passesTight ) {
    return true;
  }
  else {
    return false;
  }
}


bool Electron::IsTightTwiki()
{
  if (wztree == 0 ) {
    std::cout << "WZEvent pointer is ZERO!!!! \n";
    return false;
  }

  bool passesTight = false;
  bool passesPtEtaCuts = false;
  
  if ( ((wztree->eleIDbit)->at(index))>>2&1 ) {
    passesTight = true;
  }
  
  if ( Pt() > 10 && Eta() < 2.5 ) {
    passesPtEtaCuts = true;
  }

  if ( passesTight && passesPtEtaCuts ) {
    return true;
  }
  else {
    return false;
  }
}


Muon::Muon(int ind, double pt, double eta, double phi, float ch)
  : Lepton(ind, pt, eta, phi, ch)
{
  pdgid = 13;
}

bool Muon::IsTight() {

  if (wztree == 0 ) { 
    std::cout << "WZEvent pointer is ZERO!!!! \n";
    return false;
  }

  bool passesTight      = false;
//  bool passesPtEtaCuts = false;
  bool passesIsolation  = false;

  if ( (wztree->muChi2NDF)->at(index) < 10 && (wztree->muMuonHits)->at(index) > 0 &&
       (wztree->muStations)->at(index) > 1 && (wztree->muD0)->at(index) < 0.2 &&
       (wztree->muDz)->at(index) < 0.5 && (wztree->muPixelHits)->at(index) > 0 &&
       (wztree->muTrkLayers)->at(index) > 5 ) {

    passesTight = true;
  }

//  if ( (wztree->muPt)->at(index) > 10 && (wztree->muEta)->at(index) < 2.4 ) {
//    passesPtEtaCuts = true;
//  }


  // From UW Twiki
  // (mu.chargedHadronIso()
  // +max(mu.photonIso()+mu.neutralHadronIso()
  //       -0.5*mu.puChargedHadronIso,0.0))/mu.pt()

  if ( ((wztree->muPFChIso)->at(index) +
        TMath::Max((wztree->muPFPhoIso)->at(index) + (wztree->muPFNeuIso)->at(index) -
                   0.5 * (wztree->muPFPUIso)->at(index) , 0.)) / Pt() < 0.12 ) {
    passesIsolation = true;
  }

  if (passesTight && passesIsolation) { 
    return true; 
  } else {
    return false;
  }
}


bool Muon::IsTightTwiki()
{
  if (wztree == 0 ) { 
    std::cout << "WZEvent pointer is ZERO!!!! \n";
    return false;
  }

  bool passesTight      = false;
  bool passesPtEtaCuts = false;
  bool passesIsolation  = false;

  if ( (wztree->muChi2NDF)->at(index) < 10 && (wztree->muMuonHits)->at(index) > 0 &&
       (wztree->muStations)->at(index) > 1 && (wztree->muD0)->at(index) < 0.2 &&
       (wztree->muDz)->at(index) < 0.5 && (wztree->muPixelHits)->at(index) > 0 &&
       (wztree->muTrkLayers)->at(index) > 5 ) {
    passesTight = true;
  }

  if ( Pt() > 10 && Eta() < 2.4 ) {
    passesPtEtaCuts = true;
  }

  if ( (((wztree->muPFChIso)->at(index) +
        TMath::Max((wztree->muPFPhoIso)->at(index) + (wztree->muPFNeuIso)->at(index) -
                   0.5 * (wztree->muPFPUIso)->at(index) , 0.0)) / Pt()) < 0.12 ) {
    passesIsolation = true;
  }

  if (passesTight && passesPtEtaCuts && passesIsolation) { 
    return true; 
  } else {
    return false;
  }
}

