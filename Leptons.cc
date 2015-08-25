#include "Leptons.h"
#include "WZEvent.h"
#include "Constants.h"

#include <iostream>


using namespace std;


// Initialize static members
WZEvent* Lepton::fWZTree = 0;


void Lepton::SetWZEvent(WZEvent* wzt)
{
  fWZTree = wzt;
}


Lepton::Lepton(unsigned int index, double pt, double eta, double phi, double charge) 
{
  SetPtEtaPhiM(pt, eta, phi, 0.);
  fCharge   = charge;
  fIndex    = index;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 11;
}


pair<bool, bool> Electron::IsLooseTight()
{
  pair<bool, bool> id(false, false);

  if (fWZTree == 0) {
    cout << "WZEvent pointer is ZERO!!!! \n";
    return id;
  }

  if(!(fWZTree->nEle)) {
  return id;
  }

  if (fWZTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT&1) {
    id.first = true;
  }

  if (fWZTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT&1) {
    id.second = true;
  }

  return id;
}


bool Electron::PassesPtMinCut()
{
  bool passPt = false;

  if (Pt() > ELE_PTMIN) {
  passPt = true;
  }

  return passPt;
}


bool Electron::PassesEtaMaxCut()
{
  bool passEta = false;

  if (abs(Eta()) < ELE_ETAMAX) {
  passEta = true;
  }

  return passEta;
}


Muon::Muon(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 13;
}


pair<bool, bool> Muon::IsLooseTight()
{
  pair<bool, bool> id(false, false);

  if (fWZTree == 0) {
    cout << "WZEvent pointer is ZERO!!!! \n";
    return id;
  }

  if (!(fWZTree->nMu)) {
    return id;
  }

// From UW Twiki
// (mu.chargedHadronIso() + max(mu.photonIso()+mu.neutralHadronIso()-0.5*mu.puChargedHadronIso,0.0)) / mu.pt()
  const double relIso = (fWZTree->muPFChIso->at(fIndex) + max(fWZTree->muPFPhoIso->at(fIndex) +
                         fWZTree->muPFNeuIso->at(fIndex) - 0.5*fWZTree->muPFPUIso->at(fIndex), 0.))
                        / Pt();

  if ((fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 ||
      fWZTree->muType->at(fIndex)>>TRACKERMUON_BIT&1) &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1 &&
      relIso < MULOOSE_RELISOMIN) {
    id.first = true;
  }

  if (fWZTree->muChi2NDF->at(fIndex) < 10. && fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      fWZTree->muD0->at(fIndex) < 0.2 && fWZTree->muDz->at(fIndex) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 && fWZTree->muTrkLayers->at(fIndex) > 5 &&
      relIso < MUTIGHT_RELISOMIN &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1 ) {
    id.second = true;
  }

  return id;
}


bool Muon::PassesPtMinCut()
{
  bool passPt = false;

  if (Pt() > MU_PTMIN) {
  passPt = true;
  }

  return passPt;
}


bool Muon::PassesEtaMaxCut()
{
  bool passEta = false;

  if (abs(Eta()) < MU_ETAMAX) {
  passEta = true;
  }

  return passEta;
}

