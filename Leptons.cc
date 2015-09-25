#include "WZEvent.h"
#include "Constants.h"


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


Lepton::Lepton(unsigned int index, double pt, double eta, double phi, double charge, double relIso)
  : Lepton(index, pt, eta, phi, charge) 
{
  fRelIso = relIso;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 11;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge, double relIso)
  : Lepton(index, pt, eta, phi, charge, relIso)
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

  if (fWZTree->nEle == 0)  return id;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE &&
        fRelIso < RELISO_BARREL_LOOSE &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM &&
        fRelIso < RELISO_BARREL_MEDIUM &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.second = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE &&
        fRelIso < RELISO_ENDCAP_LOOSE &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM &&
        fRelIso < RELISO_ENDCAP_MEDIUM &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.second = true;
  }

  pair<bool, bool> idVID(false, false);
  if (fWZTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT&1)   idVID.first = true;
  if (fWZTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT&1)  idVID.second = true;

  if (id.first != idVID.first)    cout << "Error: VID different from Cut Based for LOOSE !!!" << endl;
  if (id.second != idVID.second)  cout << "Error: VID different from Cut Based for MEDIUM !!!" << endl;

  return id;
}


bool Electron::PassesPtMinCut()
{
  bool passPt = false;

  if (Pt() > ELE_PTMIN)  passPt = true;

  return passPt;
}


bool Electron::PassesEtaMaxCut()
{
  bool passEta = false;

  if (abs(Eta()) < ELE_ETAMAX)  passEta = true;

  return passEta;
}


Muon::Muon(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 13;
}


Muon::Muon(unsigned int index, double pt, double eta, double phi, double charge, double relIso)
  : Lepton(index, pt, eta, phi, charge, relIso)
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

  if (!(fWZTree->nMu))  return id;

  if ((fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 ||
      fWZTree->muType->at(fIndex)>>TRACKERMUON_BIT&1) &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1 &&
      fRelIso < MULOOSE_RELISOMIN)
    id.first = true;

  if (fWZTree->muChi2NDF->at(fIndex) < 10. && fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 && abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 && fWZTree->muTrkLayers->at(fIndex) > 5 &&
      fRelIso < MUTIGHT_RELISOMIN &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1)
    id.second = true;

  return id;
}


bool Muon::PassesPtMinCut()
{
  bool passPt = false;

  if (Pt() > MU_PTMIN)  passPt = true;

  return passPt;
}


bool Muon::PassesEtaMaxCut()
{
  bool passEta = false;

  if (abs(Eta()) < MU_ETAMAX)  passEta = true;

  return passEta;
}
