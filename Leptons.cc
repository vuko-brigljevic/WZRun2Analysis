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


Electron::Electron(unsigned int index, double pt, double eta, double phi, double charge)
  : Lepton(index, pt, eta, phi, charge)
{
  fPdgId = 11;
}


double Electron::EffA25ns()
{
  double effA = 0;
  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));

// Check Slide 4 in https://indico.cern.ch/event/370507/contribution/1/attachments/1140657/1633761/Rami_eleCB_ID_25ns.pdf
// and Slide 12 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// and Line 407 for abs(eleSCEta) in https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
  if (absEleSCEta >= 0.0 && absEleSCEta < 1.0)    effA = 0.1752;
  if (absEleSCEta >= 1.0 && absEleSCEta < 1.479)  effA = 0.1862;
  if (absEleSCEta >= 1.479 && absEleSCEta < 2.0)  effA = 0.1411;
  if (absEleSCEta >= 2.0 && absEleSCEta < 2.2)    effA = 0.1534;
  if (absEleSCEta >= 2.2 && absEleSCEta < 2.3)    effA = 0.1903;
  if (absEleSCEta >= 2.3 && absEleSCEta < 2.4)    effA = 0.2243;
  if (absEleSCEta >= 2.4 && absEleSCEta < 2.5)    effA = 0.2687;

  return effA;
}


double Electron::EffA50ns()
{
  double effA = 0;
  const double absEleSCEta = abs(fWZTree->eleSCEta->at(fIndex));

// check Slide 4 in https://indico.cern.ch/event/369239/contribution/6/attachments/1134836/1623383/Rami_eleCB_ID_50ns.pdf
// and Slide 8 in https://indico.cern.ch/event/369235/contribution/4/attachments/734635/1007867/Rami_EffAreas.pdf
// and Line 407 for abs(EleSCEta) in https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
  if      (absEleSCEta >= 0.0 && absEleSCEta < 0.8)  effA = 0.0973;
  else if (absEleSCEta >= 0.8 && absEleSCEta < 1.3)  effA = 0.0954;
  else if (absEleSCEta >= 1.3 && absEleSCEta < 2.0)  effA = 0.0632;
  else if (absEleSCEta >= 2.0 && absEleSCEta < 2.2)  effA = 0.0727;
  else if (absEleSCEta >= 2.2 && absEleSCEta < 2.5)  effA = 0.1337;

  return effA;
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
  const double effA = EffA25ns();

// Check Slide 2 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                        max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
  const double relIsoWithEA = (fWZTree->elePFChIso->at(fIndex) +
    max(0.0, fWZTree->elePFPhoIso->at(fIndex) - fWZTree->rho * effA +
             fWZTree->elePFNeuIso->at(fIndex))) / Pt();

  fEoverPinv = fWZTree->eleEoverPInv->at(fIndex);

  /*
  if( el->ecalEnergy() == 0 ){
        printf("Electron energy is zero!\n");
        ooEmooP_.push_back( 1e30 );
      }else if( !std::isfinite(el->ecalEnergy())){
        printf("Electron energy is not finite!\n");
        ooEmooP_.push_back( 1e30 );
      }else{
        ooEmooP_.push_back( fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
      }
*/


  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE_25ns &&
        relIsoWithEA < RELISOWITHEA_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_LOOSE_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_LOOSE_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
  }
  if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_LOOSE_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_LOOSE_25ns &&
        relIsoWithEA < RELISOWITHEA_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_ENDCAP_LOOSE_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_ENDCAP_LOOSE_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_ENDCAP &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.first = true;
  }

  if (absEleSCEta <= ETASCBARREL) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_MEDIUM_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_MEDIUM_25ns &&
        relIsoWithEA < RELISOWITHEA_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleEoverPInv->at(fIndex)) < OOEMOOP_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleD0->at(fIndex)) < D0_BARREL_MEDIUM_25ns &&
        abs(fWZTree->eleDz->at(fIndex)) < DZ_BARREL_MEDIUM_25ns &&
        fWZTree->eleMissHits->at(fIndex) <= EXPMISSINNERHITS_BARREL &&
        fWZTree->eleConvVeto->at(fIndex) == true)
      id.second = true;
  }
  if (absEleSCEta > ETASCBARREL && absEleSCEta < ETASCENDCAP) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_ENDCAP_MEDIUM_25ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_ENDCAP_MEDIUM_25ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_ENDCAP_MEDIUM_25ns &&
        relIsoWithEA < RELISOWITHEA_ENDCAP_MEDIUM_25ns &&
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
  const double effA = EffA50ns();
// Check Slide 8 in https://indico.cern.ch/event/369235/contribution/4/attachments/734635/1007867/Rami_EffAreas.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                       max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
  const double relIsoWithEA = (fWZTree->elePFChIso->at(fIndex) +
    max(0.0, fWZTree->elePFPhoIso->at(fIndex) - fWZTree->rho * effA +
             fWZTree->elePFNeuIso->at(fIndex))) / Pt();

  if (!(absEleSCEta > ETASCBARREL)) {
    if (fWZTree->eleSigmaIEtaIEtaFull5x5->at(fIndex) < FULL5x5_SIGMAIETAIETA_BARREL_LOOSE_50ns &&
        abs(fWZTree->eledEtaAtVtx->at(fIndex)) < DETAIN_BARREL_LOOSE_50ns &&
        abs(fWZTree->eledPhiAtVtx->at(fIndex)) < DPHIIN_BARREL_LOOSE_50ns &&
        fWZTree->eleHoverE->at(fIndex) < HOVERE_BARREL_LOOSE_50ns &&
        relIsoWithEA < RELISOWITHEA_BARREL_LOOSE_50ns &&
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
        relIsoWithEA < RELISOWITHEA_ENDCAP_LOOSE_50ns &&
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
        relIsoWithEA < RELISOWITHEA_BARREL_MEDIUM_50ns &&
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
        relIsoWithEA < RELISOWITHEA_ENDCAP_MEDIUM_50ns &&
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


pair<bool, bool> Muon::IsLooseTight()
{
  pair<bool, bool> id(false, false);

  if (fWZTree == 0) {
    cout << "WZEvent pointer is ZERO!!!! \n";
    return id;
  }

  if (!(fWZTree->nMu))  return id;

// From UW Twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZ13TeV
// (mu.chargedHadronIso() +
//  max(mu.photonIso() + mu.neutralHadronIso() - 0.5 * mu.puChargedHadronIso,0.0))
// / mu.pt()
  const double relIso = (fWZTree->muPFChIso->at(fIndex) + max(fWZTree->muPFPhoIso->at(fIndex) +
                         fWZTree->muPFNeuIso->at(fIndex) - 0.5*fWZTree->muPFPUIso->at(fIndex), 0.))
                        / Pt();

  if ((fWZTree->muType->at(fIndex)>>GLOBALMUON_BIT&1 ||
      fWZTree->muType->at(fIndex)>>TRACKERMUON_BIT&1) &&
      fWZTree->muType->at(fIndex)>>PFMUON_BIT&1 &&
      relIso < MULOOSE_RELISOMIN)
    id.first = true;

  if (fWZTree->muChi2NDF->at(fIndex) < 10. && fWZTree->muMuonHits->at(fIndex) > 0 &&
      fWZTree->muStations->at(fIndex) > 1 &&
      abs(fWZTree->muD0->at(fIndex)) < 0.2 && abs(fWZTree->muDz->at(fIndex)) < 0.5 &&
      fWZTree->muPixelHits->at(fIndex) > 0 && fWZTree->muTrkLayers->at(fIndex) > 5 &&
      relIso < MUTIGHT_RELISOMIN &&
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
