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


Lepton::Lepton(unsigned int index, double pt, double eta, double phi, double charge,
               double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns)
  : Lepton(index, pt, eta, phi, charge) 
{
  fRelIsoDeltaB = relIsoDeltaB;
  fRelIsoEffArea25ns = relIsoEffA25ns;
  fRelIsoEffArea50ns = relIsoEffA50ns;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 11;
}


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge,
               double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns)
  : Lepton(index, pt, eta, phi, charge, relIsoDeltaB, relIsoEffA25ns, relIsoEffA50ns)
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

  if (fWZTree->eleIDbit->at(fIndex)>>ELELOOSE_BIT&1)   id.first = true;
  if (fWZTree->eleIDbit->at(fIndex)>>ELEMEDIUM_BIT&1)  id.second = true;

  return id;
}


pair<bool, bool> Electron::IsLooseTightCutBased25ns()
{
  pair<bool, bool> id(false, false);

  if (fWZTree == 0) {
    cout << "WZEvent pointer is ZERO!!!! \n";
    return id;
  }

  if (fWZTree->nEle == 0)  return id;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE_25ns &&
        fRelIsoEffArea25ns < RELISOWITHEA_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM_25ns &&
        fRelIsoEffArea25ns < RELISOWITHEA_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.second = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE_25ns &&
        fRelIsoEffArea25ns < RELISOWITHEA_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM_25ns &&
        fRelIsoEffArea25ns < RELISOWITHEA_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.second = true;
  }

  return id;
}


pair<bool, bool> Electron::IsLooseTightCutBased50ns()
{
  pair<bool, bool> id(false, false);

  if (fWZTree == 0) {
    cout << "WZEvent pointer is ZERO!!!! \n";
    return id;
  }

  if (!(fWZTree->nEle))  return id;

  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));
  if (!(absEleSCEta > ETASCBARREL)) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE_50ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE_50ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE_50ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE_50ns &&
        fRelIsoEffArea50ns < RELISOWITHEA_BARREL_LOOSE_50ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE_50ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE_50ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE_50ns &&
        !(fWZTree->eleMissHits->at(fIndex) > EXPMISSINNERHITS_BARREL) &&
        fWZTree->eleConvVeto->at(fIndex) == PASSCONVERSIONVETO)
      id.first = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE_50ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE_50ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE_50ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE_50ns &&
        fRelIsoEffArea50ns < RELISOWITHEA_ENDCAP_LOOSE_50ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE_50ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE_50ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE_50ns &&
        !(fWZTree->eleMissHits->at(fIndex) > EXPMISSINNERHITS_ENDCAP) &&
        fWZTree->eleConvVeto->at(fIndex) == PASSCONVERSIONVETO)
      id.first = true;
  }

  if (!(absEleSCEta > ETASCBARREL)) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM_50ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM_50ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM_50ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM_50ns &&
        fRelIsoEffArea50ns < RELISOWITHEA_BARREL_MEDIUM_50ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM_50ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM_50ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM_50ns &&
        !(fWZTree->eleMissHits->at(fIndex) > EXPMISSINNERHITS_BARREL) &&
        fWZTree->eleConvVeto->at(fIndex) == PASSCONVERSIONVETO)
      id.second = true;
  } else if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM_50ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM_50ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM_50ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM_50ns &&
        fRelIsoEffArea50ns < RELISOWITHEA_ENDCAP_MEDIUM_50ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_MEDIUM_50ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_MEDIUM_50ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_MEDIUM_50ns &&
        !(fWZTree->eleMissHits->at(fIndex) > EXPMISSINNERHITS_ENDCAP) &&
        fWZTree->eleConvVeto->at(fIndex) == PASSCONVERSIONVETO)
      id.second = true;
  }

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


Muon::Muon(unsigned int index, double pt, double eta, double phi, double charge,
               double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns)
  : Lepton(index, pt, eta, phi, charge, relIsoDeltaB, relIsoEffA25ns, relIsoEffA50ns)
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
      fRelIsoDeltaB < MULOOSE_RELISOMIN)
    id.first = true;

  if (fWZTree->muChi2NDF->at(fIndex) < 10. && fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 && abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 && fWZTree->muTrkLayers->at(fIndex) > 5 &&
      fRelIsoDeltaB < MUTIGHT_RELISOMIN &&
      fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1)
    id.second = true;

  return id;
}


pair<bool, bool> Muon::IsLooseTightCutBased25ns()
{
  return IsLooseTight();
}


pair<bool, bool> Muon::IsLooseTightCutBased50ns()
{
  return IsLooseTight();
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
