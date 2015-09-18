#include "WZEvent.h"
#include "Constants.h"

#include <iomanip>


using namespace std;


double EffArea25ns(double absEleSCEta)
{
  double effA25ns = 0;
// Slide 4 in https://indico.cern.ch/event/370507/contribution/1/attachments/1140657/1633761/Rami_eleCB_ID_25ns.pdf
// and Slide 12 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// and Line 407 for abs(eleSCEta) in https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
  if (absEleSCEta >= 0.0 && absEleSCEta < 1.0)         effA25ns = 0.1752;
  else if (absEleSCEta >= 1.0 && absEleSCEta < 1.479)  effA25ns = 0.1862;
  else if (absEleSCEta >= 1.479 && absEleSCEta < 2.0)  effA25ns = 0.1411;
  else if (absEleSCEta >= 2.0 && absEleSCEta < 2.2)    effA25ns = 0.1534;
  else if (absEleSCEta >= 2.2 && absEleSCEta < 2.3)    effA25ns = 0.1903;
  else if (absEleSCEta >= 2.3 && absEleSCEta < 2.4)    effA25ns = 0.2243;
  else if (absEleSCEta >= 2.4 && absEleSCEta < 2.5)    effA25ns = 0.2687;
  else  effA25ns = 0;

  return effA25ns;
}


double EffArea50ns(double absEleSCEta)
{
  double effA50ns = 0;
// Slide 4 in https://indico.cern.ch/event/369239/contribution/6/attachments/1134836/1623383/Rami_eleCB_ID_50ns.pdf
// and Slide 8 in https://indico.cern.ch/event/369235/contribution/4/attachments/734635/1007867/Rami_EffAreas.pdf
// and Line 407 for abs(EleSCEta) in https://github.com/ikrav/EgammaWork/blob/ntupler_and_VID_demos_747/ElectronNtupler/plugins/SimpleElectronNtupler.cc
  if      (absEleSCEta >= 0.0 && absEleSCEta < 0.8)  effA50ns = 0.0973;
  else if (absEleSCEta >= 0.8 && absEleSCEta < 1.3)  effA50ns = 0.0954;
  else if (absEleSCEta >= 1.3 && absEleSCEta < 2.0)  effA50ns = 0.0632;
  else if (absEleSCEta >= 2.0 && absEleSCEta < 2.2)  effA50ns = 0.0727;
  else if (absEleSCEta >= 2.2 && absEleSCEta < 2.5)  effA50ns = 0.1337;
  else  effA50ns = 0;

  return effA50ns;
}


WZEvent::WZEvent(TTree* tree) :
  WZBASECLASS(tree)
{
// Set static pointer to Lepton class
  Lepton::SetWZEvent(this);
  Particle::SetWZEvent(this);
}


void WZEvent::Clear()
{
  fHLT25ns.clear();
  fHLT50ns.clear();

  fFinalState = undefined;
  fSelectionLevel = Undefined;

// Empty leptons
  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ) {
    delete *lIt;  
    lIt = fLeptons.erase(lIt);
  }

  fTightLeptonsIndex.clear();
  fZLeptonsIndex = make_pair(999, 999);
  fWLeptonIndex = 999;

  fGenWLepton         = NULL;
  fGenZLeptons        = make_pair( (GenParticle*) NULL, (GenParticle*) NULL);
  fGenWDecayFlavor    = 0;
  fGenZDecayFlavor    = 0;

  // Empty Gen Particles

  for (vector<GenParticle*>::iterator gIt = fGenParticles.begin(); gIt != fGenParticles.end(); ) {
    delete *gIt;  
    gIt = fGenParticles.erase(gIt);
  }

}


void WZEvent::ReadEvent()
{
  Clear();

  if (HLT50ns) {
    const vector<unsigned int> hlt50nsBits { 16, 10, 11, 36, 37, 38, 39, 35 };
    for (vector<unsigned int>::const_iterator bIt = hlt50nsBits.begin();
         bIt != hlt50nsBits.end(); ++bIt) {
      (HLT50ns>>(*bIt)&1)  ?  fHLT50ns.push_back(true)  :  fHLT50ns.push_back(false);
    }
  }

// in ggNtuplizer versions V07-04-09+
  if (HLTEleMuX) {
    const vector<unsigned int> hlt25nsBits { 8, 20, 21, 41, 42, 9, 43, 44, 28 };
        for (vector<unsigned int>::const_iterator bIt = hlt25nsBits.begin();
         bIt != hlt25nsBits.end(); ++bIt) {
      (HLTEleMuX>>(*bIt)&1)  ?  fHLT25ns.push_back(true)  :  fHLT25ns.push_back(false);
    }
  }

// Electrons
  for (unsigned int indexEle = 0; indexEle < eleCharge->size(); indexEle++) {
// From UW Twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZ13TeV
// (mu.chargedHadronIso() +
//  max(mu.photonIso() + mu.neutralHadronIso() - 0.5 * mu.puChargedHadronIso,0.0))
// / mu.pt()
    const double relIsoDeltaB = (elePFChIso->at(indexEle) + max(elePFPhoIso->at(indexEle) +
                                 elePFNeuIso->at(indexEle) - 0.5 * elePFPUIso->at(indexEle), 0.0))
                                / elePt->at(indexEle);

    const double absEleSCEta = abs(eleSCEta->at(indexEle));
// Slide 2 in https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                        max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
    const double relIsoEffA25ns = (elePFChIso->at(indexEle) + max(elePFPhoIso->at(indexEle) +
                                   elePFNeuIso->at(indexEle) - rho * EffArea25ns(absEleSCEta), 0.0))
                                  / elePt->at(indexEle);    
// Slide 8 in https://indico.cern.ch/event/369235/contribution/4/attachments/734635/1007867/Rami_EffAreas.pdf
// relIsoWithEA = 1/pt * (pfIso.sumChargedHadronPt +
//                       max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffArea)
    const double relIsoEffA50ns = (elePFChIso->at(indexEle) + max(elePFPhoIso->at(indexEle) +
                                   elePFNeuIso->at(indexEle) - rho * EffArea50ns(absEleSCEta), 0.0))
                                  / elePt->at(indexEle);    

    Electron* ele = new Electron(indexEle, elePt->at(indexEle), eleEta->at(indexEle),
                                 elePhi->at(indexEle), eleCharge->at(indexEle),
                                 relIsoDeltaB, relIsoEffA25ns, relIsoEffA50ns);
    fLeptons.push_back(ele);
  }

// Muons
  for (unsigned int indexMu = 0; indexMu < muCharge->size(); indexMu++) {
// From UW Twiki https://twiki.cern.ch/twiki/bin/viewauth/CMS/WZ13TeV
// (mu.chargedHadronIso() +
//  max(mu.photonIso() + mu.neutralHadronIso() - 0.5 * mu.puChargedHadronIso,0.0))
// / mu.pt()
    const double relIsoDeltaB = (muPFChIso->at(indexMu) + max(muPFPhoIso->at(indexMu) +
                                 muPFNeuIso->at(indexMu) - 0.5 * muPFPUIso->at(indexMu), 0.0))
                                / muPt->at(indexMu);

    Muon* mu = new Muon(indexMu, muPt->at(indexMu), muEta->at(indexMu),
                        muPhi->at(indexMu), muCharge->at(indexMu), relIsoDeltaB, 0.0, 0.0);
    fLeptons.push_back(mu);
  }

  if (!isData) {
    ReadGenEvent();
  }


}




void WZEvent::ReadGenEvent() {

  // Gen Particles
  for (unsigned int indexGen = 0; indexGen < mcPID->size(); indexGen++) {

    GenParticle * gp = new GenParticle(indexGen,
				       mcPt->at(indexGen),
				       mcEta->at(indexGen),
				       mcPhi->at(indexGen));
    fGenParticles.push_back(gp);

  }

  GetGenWZFinalState();

}

void WZEvent::GetGenWZFinalState() {



  // Find out W and Z decay channels

  int trueWDecayMode = -1;
  int trueZDecayMode = -1;

  for (int igen=0; igen<fGenParticles.size(); igen++) {
    GenParticle * gp = fGenParticles.at(igen);
    int pdgId = gp->PdgId();
    if (abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15) {
      if (mcMomPID->at(igen) == 23) {
	trueZDecayMode = abs(pdgId);
      }
      if (abs(mcMomPID->at(igen)) == 24) {
	trueWDecayMode = abs(pdgId);
      }
    } else if (abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) {
      // Look also at neutrinos for W
      if (abs(mcMomPID->at(igen)) == 24) {
	trueWDecayMode = abs(pdgId)-1;
      }
    }
  }

  fGenWDecayFlavor = trueWDecayMode;
  fGenZDecayFlavor = trueZDecayMode;

  // Now look for the 3 leptons if it's eee,eemu,mumue,mumumu

  if (abs(trueZDecayMode) != 11 && abs(trueZDecayMode) != 13) return;
  if (abs(trueWDecayMode) != 11 && abs(trueWDecayMode) != 13) return;


  std::vector<GenParticle * > stablePromptLeptons;

  std::vector<GenParticle *> trueWLeptons;
  std::vector<GenParticle *> trueZLeptons;

  for (int igen=0; igen<fGenParticles.size(); igen++) {
    GenParticle * gp = fGenParticles.at(igen);
    int pdgId = gp->PdgId();
    if (abs(pdgId) == 11 || abs(pdgId) == 13) {
      int index = gp->GetIndex();
      if (mcStatus->at(index) == 1) {
	stablePromptLeptons.push_back(gp);
      }

      int momId     =  mcMomPID->at(index);
      int granMomId =  mcGMomPID->at(index);

      if (abs(momId) == 23 || abs(granMomId) == 23) {
	trueZLeptons.push_back(gp);
      }
      if (abs(momId) == 24 || abs(granMomId) == 24) {
	trueWLeptons.push_back(gp);	
      }
    }
  }

  std::cout << "Stable Leptons " << stablePromptLeptons.size() << std::endl;

  if (stablePromptLeptons.size() != 3) return;

  if (trueWLeptons.size() != 1 
      || trueZLeptons.size() != 2) {

    std::cout << "WEIRDODS NUMBERS OF W & Z LEPTONS: W " 
	      << trueWLeptons.size()
	      << "\t Z : " << trueZLeptons.size() << std::endl;
  } else {
    fGenWLepton  = trueWLeptons[0];
    fGenZLeptons = make_pair(trueZLeptons[0],trueZLeptons[1]);

  }


}

bool WZEvent::IsInGenXSPhaseSpace() {

  if (  abs(GetGenWDecayFlavor()) != 11 && abs(GetGenWDecayFlavor()) != 13) return false;

  if (  abs(GetGenZDecayFlavor()) != 11 && abs(GetGenZDecayFlavor()) != 13) return false;

  if (  !fGenZLeptons.first ||  !fGenZLeptons.second ) return false;


  std::cout << "Is in XS phase space: z lepton pointers : " 
	    <<  fGenZLeptons.first << "\t" <<  fGenZLeptons.second << std::endl;

  TLorentzVector zp4 = *fGenZLeptons.first + *fGenZLeptons.second;

  double zmass = zp4.M();

  if (zmass>MZ_MIN && zmass<MZ_MAX) {
    return true;
  } else {
    return false;
  }

}


bool WZEvent::IsInGenFiducialPhaseSpace() {

  // MZ cut
  bool passed = IsInGenXSPhaseSpace();
  
  if (!passed) return false;

  // All 3 leptons: |eta|<2.5
  if ( abs(fGenZLeptons.first->Eta())>2.5 ||  abs(fGenZLeptons.second->Eta())>2.5 || 
       abs(fGenWLepton->Eta())>2.5 ) 
    passed = false;

  // All 3 leptons: Pt>10
  if ( abs(fGenZLeptons.first->Pt())<10 ||  abs(fGenZLeptons.second->Pt())<10. || 
       abs(fGenWLepton->Pt())<10)
    passed = false;

  // At least one out of 3: Pt>20
  if ( abs(fGenZLeptons.first->Pt())<20 &&  abs(fGenZLeptons.second->Pt())<20. &&
       abs(fGenWLepton->Pt())<20)
    passed = false;

  return passed;

}


bool WZEvent::PassesPreselection()
{
  if      (!(fSelectionLevel < Preselection))     return true;
  else if (fSelectionLevel == FailsPreselection)  return false;

  bool passed = false;
  fSelectionLevel = FailsPreselection;

  if((nEle + nMu) < 3)  return passed;

  unsigned int nEleTight = 0;
  unsigned int nMuTight = 0;

  for (vector<Lepton*>::iterator lIt = fLeptons.begin(); lIt != fLeptons.end(); ++lIt) {
    if ((*lIt)->PassesPtMinCut() && (*lIt)->PassesEtaMaxCut() &&
        (*lIt)->IsLooseTight().second) {
      unsigned int index = distance(fLeptons.begin(), lIt);
      fTightLeptonsIndex.push_back(index);
      if ((*lIt)->GetPdgId() == 11)       nEleTight++;
      else if ((*lIt)->GetPdgId() == 13)  nMuTight++;
    }
  }

  if (nMuTight + nEleTight != fTightLeptonsIndex.size()) {
    cout << "ERROR: Number of Tight leptons DIFFERS from number of Tight lepton indices !!!" << endl;
    cout << nMuTight << " + " << nEleTight << " != " << fTightLeptonsIndex.size() << endl;
    return passed;
  }

// ### put the condition on number of tight leptons here!!! ###
  if (fTightLeptonsIndex.size() == N_TIGHTLEPTONS)  passed = true;

  if (passed) {
    fSelectionLevel = Preselection;
    if      (nEleTight == 3 && nMuTight == 0)  fFinalState = eee;
    else if (nEleTight == 2 && nMuTight == 1)  fFinalState = eem;
    else if (nEleTight == 1 && nMuTight == 2)  fFinalState = mme;
    else if (nEleTight == 0 && nMuTight == 3)  fFinalState = mmm;
  }

  return passed;
}


bool WZEvent::PassesZSelection()
{
  if (!(fSelectionLevel < ZSelection))  return true;

  bool passed = false;
  if (!PassesPreselection())  return passed;

  vector<pair<unsigned int, unsigned int> > iZl;
  vector<double> dMassZ;

  for (vector<unsigned int>::const_iterator iIt1 = fTightLeptonsIndex.begin();
       iIt1 != fTightLeptonsIndex.end(); ++iIt1) {
    for (vector<unsigned int>::const_iterator iIt2 = fTightLeptonsIndex.begin()+1;
         iIt2 != fTightLeptonsIndex.end(); ++iIt2) {
      unsigned int il1 = *iIt1;
      unsigned int il2 = *iIt2;

      if (abs(fLeptons.at(il1)->GetPdgId()) != abs(fLeptons.at(il2)->GetPdgId()) ||
          fLeptons.at(il1)->GetCharge() == fLeptons.at(il2)->GetCharge())
        continue;

      const double mZCand = (*(fLeptons.at(il1)) + *(fLeptons.at(il2))).M();
      if (mZCand < MZ_MIN || mZCand > MZ_MAX ||
          !(max(fLeptons.at(il1)->Pt(), fLeptons.at(il2)->Pt()) > ZLEADINGLEPTON_PTMIN))
        continue;
      else {
        passed = true;
        iZl.push_back(make_pair(il1, il2));
        const double dMZ = abs(mZCand - PDG_ZMASS);
        dMassZ.push_back(dMZ);
      }
    }
  }

  if (iZl.size() == 0)  return false;
  else {
    vector<double>::iterator bestMin = min_element(dMassZ.begin(), dMassZ.end());
    unsigned int iMin = distance(dMassZ.begin(), bestMin);
    fZLeptonsIndex = iZl.at(iMin);
  }

  if (fZLeptonsIndex.first != fZLeptonsIndex.second) {
    if (fLeptons.at(fZLeptonsIndex.first)->Pt() < fLeptons.at(fZLeptonsIndex.second)->Pt())
      swap(fZLeptonsIndex.first, fZLeptonsIndex.second);
  }

  if (passed)  fSelectionLevel = ZSelection;

  return passed;
}


bool WZEvent::PassesWSelection()
{
  if (!(fSelectionLevel < WSelection))  return true;

  bool passed = false;
  if (!PassesPreselection() || !PassesZSelection())  return passed;

  vector<unsigned int> iWl;
  vector<double> ptWl;

  for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
       iIt != fTightLeptonsIndex.end(); ++iIt) {
    if (*iIt != fZLeptonsIndex.first && *iIt != fZLeptonsIndex.second) {
      const double wlPt = fLeptons.at(*iIt)->Pt();
      const double deltaR1 = fLeptons.at(*iIt)->DeltaR(*(fLeptons.at(fZLeptonsIndex.first)));
      const double deltaR2 = fLeptons.at(*iIt)->DeltaR(*(fLeptons.at(fZLeptonsIndex.second)));

      if (wlPt > WLEPTON_PTMIN && deltaR1 > WZ_DELTARMIN && deltaR2 > WZ_DELTARMIN) {
/*        
        if (fLeptons.at(*iIt)->GetPdgId() == 11) {
          const unsigned int index = fLeptons.at(*iIt)->GetIndex();
          if (!(eleIDbit->at(index)>>ELETIGHT_BIT&1))  continue;
        }
*/
        passed = true;
        iWl.push_back(*iIt);
        ptWl.push_back(wlPt);
      } else continue;
    }
  }

  if (iWl.size() > 1)  cout << "WARNING: More than 1 W lepton - highest Pt criterium !!!" << endl;

  if (iWl.size() == 0)  return false;
  else {
    vector<double>::iterator bestMax = max_element(ptWl.begin(), ptWl.end());
    unsigned int iMax = distance(ptWl.begin(), bestMax);
    fWLeptonIndex = iWl.at(iMax);
  }

  if (passed)  fSelectionLevel = WSelection;

  return passed;
}


bool WZEvent::PassesFullSelection()
{
  if (!(fSelectionLevel < FullSelection))  return true;

  bool passed = false;
  if (!PassesPreselection() || !PassesZSelection() || !PassesWSelection())  return passed;

  if (pfMET > METMIN &&
      fZLeptonsIndex.first != fZLeptonsIndex.second &&
      fZLeptonsIndex.first != fWLeptonIndex)
    passed = true;

  if (passed) {
    fSelectionLevel = FullSelection;
    unsigned int zFlavor = abs(fLeptons.at(fZLeptonsIndex.first)->GetPdgId());
    unsigned int wFlavor = abs(fLeptons.at(fWLeptonIndex)->GetPdgId());

    if (zFlavor == 11) {
      if      (wFlavor == 11)  fFinalState = eee;
      else if (wFlavor == 13)  fFinalState = eem;
    } else if (zFlavor == 13) {
      if      (wFlavor == 11)  fFinalState = mme;
      else if (wFlavor == 13)  fFinalState = mmm;
    }
  }

  return passed;
}


pair<Lepton*, Lepton*> WZEvent::GetZLeptons()
{
  Lepton* zL1 = fLeptons.at(fZLeptonsIndex.first);
  Lepton* zL2 = fLeptons.at(fZLeptonsIndex.second);
  pair<Lepton*, Lepton*> zL = make_pair(zL1, zL2);
  return zL;
}


bool WZEvent::PassesFinalSelection()
{
  if (!(fSelectionLevel < FinalSelection))  return true;

  bool passed = false;
  if (!PassesPreselection() || !PassesZSelection() || !PassesWSelection() || !PassesFullSelection())
    return passed;

  const double mass3L = (*(GetZLeptons().first) + *(GetZLeptons().second) + *(GetWLepton())).M();
  (mass3L > MASS3LMIN) ? passed = true : passed = false;

  if (passed)  fSelectionLevel = FinalSelection;

  return passed;
}


void WZEvent::DumpEvent(ostream& out, int verbosity)
{
  out << run << ":" << lumis << ":" << event;

// Preselection (verbosity = 5) format for event listing:
// l1_pt:l1_eta:l1_phi:l1_iso:l2_pt:l2_eta:l2_phi:l2_iso:l3_pt:l3_eta:l3_phi:l3_iso:
// deltaR_l1_l2:deltaR_l1_l3:deltaR_l2_l3:0:MET:MET_phi:3l_mass

  if (verbosity == 5) {
    const unsigned int nTight = fTightLeptonsIndex.size();
    if (nTight == N_TIGHTLEPTONS) {
      vector<Lepton*> tightLeptons;
      for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
           iIt != fTightLeptonsIndex.end(); ++iIt)
        if (*iIt < fLeptons.size())  tightLeptons.push_back(fLeptons.at(*iIt));

      sort(tightLeptons.begin(), tightLeptons.end(), HigherPt);

      TLorentzVector total3L;
      for (vector<Lepton*>::const_iterator lIt = tightLeptons.begin();
           lIt != tightLeptons.end(); ++lIt) {
        out << fixed << setprecision(4)
            << ":" << (*lIt)->Pt() << ":" << (*lIt)->Eta() << ":" << (*lIt)->Phi()
            << ":" << (*lIt)->GetRelIsoDeltaB();
        total3L += *(*lIt);
      }

      for (vector<Lepton*>::const_iterator lIt1 = tightLeptons.begin();
           lIt1 != tightLeptons.end(); ++lIt1) {
        for (vector<Lepton*>::const_iterator lIt2 = tightLeptons.begin() + 1;
             lIt2 != tightLeptons.end(); ++lIt2)
          out << fixed << setprecision(4) << ":" << (*lIt1)->DeltaR(*(*lIt2));
      }

      const double mass3L = total3L.M();
      out << fixed << setprecision(4)
          << ":" << 0 << ":" << pfMET << ":" << pfMETPhi << ":" << mass3L;

      for (vector<Lepton*>::iterator lIt = tightLeptons.begin(); lIt != tightLeptons.end(); ) {
        delete *lIt;  
        lIt = tightLeptons.erase(lIt);
      }

    }
  }

// Event listing format for Z Selection (verbosity = 7) and W/Full/Final Selection (verbosity >= 10):
//  Zl1_pt:Zl1_eta:Zl1_phi:Zl1_iso:Zl2_pt:Zl2_eta:Zl2_phi:Zl2_iso:Wl_pt:Wl_eta:Wl_phi:Wl_iso:
//  deltaR_Zl1_Zl2:deltaR_Zl1_Wl:deltaR_Zl2_Wl:Z_mass:MET:MET_phi:3l_mass

  if (verbosity == 7) {
    if (fZLeptonsIndex.first < fLeptons.size() && fZLeptonsIndex.second < fLeptons.size()) {
      out << fixed << setprecision(4)
          << ":" << GetZLeptons().first->Pt()
          << ":" << GetZLeptons().first->Eta()
          << ":" << GetZLeptons().first->Phi()
          << ":" << GetZLeptons().first->GetRelIsoDeltaB()
          << ":" << GetZLeptons().second->Pt()
          << ":" << GetZLeptons().second->Eta()
          << ":" << GetZLeptons().second->Phi()
          << ":" << GetZLeptons().second->GetRelIsoDeltaB();

      for (vector<unsigned int>::const_iterator iIt = fTightLeptonsIndex.begin();
           iIt != fTightLeptonsIndex.end(); ++iIt) {
        if (*iIt != fZLeptonsIndex.first && *iIt != fZLeptonsIndex.second && *iIt < fLeptons.size()) {
          out << ":" << fLeptons.at(*iIt)->Pt()
              << ":" << fLeptons.at(*iIt)->Eta()
              << ":" << fLeptons.at(*iIt)->Phi()
              << ":" << fLeptons.at(*iIt)->GetRelIsoDeltaB()
              << ":" << GetZLeptons().first->DeltaR(*(GetZLeptons().second))
              << ":" << GetZLeptons().first->DeltaR(*(fLeptons.at(*iIt)))
              << ":" << GetZLeptons().second->DeltaR(*(fLeptons.at(*iIt)))
              << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second)).M()
              << ":" << pfMET << ":" << pfMETPhi
              << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second) + *(fLeptons.at(*iIt))).M();
        }
      }
    }
  }

  if (!(verbosity < 10)) {
    vector<unsigned int> indexWZLeptons =
      { fZLeptonsIndex.first, fZLeptonsIndex.second, fWLeptonIndex };
    for (vector<unsigned int>::iterator iIt = indexWZLeptons.begin();
         iIt != indexWZLeptons.end(); ++iIt) {
      if (*iIt < fLeptons.size()) {
        out << fixed << setprecision(4)
            << ":" << fLeptons.at(*iIt)->Pt()
            << ":" << fLeptons.at(*iIt)->Eta()
            << ":" << fLeptons.at(*iIt)->Phi()
            << ":" << fLeptons.at(*iIt)->GetRelIsoDeltaB();
      }
    }
    out << fixed << setprecision(4)
        << ":" << GetZLeptons().first->DeltaR(*(GetZLeptons().second))
        << ":" << GetZLeptons().first->DeltaR(*(GetWLepton()))
        << ":" << GetZLeptons().second->DeltaR(*(GetWLepton()))
        << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second)).M()
        << ":" << pfMET << ":" << pfMETPhi
        << ":" << (*(GetZLeptons().first) + *(GetZLeptons().second) + *(GetWLepton())).M();
  }

  out << endl;
}




void WZEvent::DumpGenEvent(ostream& out) 
{

  out << "MC Tree : process ID = " << processID << std::endl;
  out << "===============\n";
  for (int igen = 0; igen<nMC; igen++) {
    //    char  pyName[20];
    //    TPythia::Pyname((mcPID)->at(igen),pyName);

    unsigned short statusFlag = (mcStatusFlag)->at(igen);
    unsigned short isPrompt = (statusFlag>>1 & 1);

    out << "Gen Particle: " << (mcPID)->at(igen)
	<< "\t status   : " << (mcStatus)->at(igen)
	<< "\t mother ID: " << (mcMomPID)->at(igen)
	<< "\t GrandMa  ID: " << (mcGMomPID)->at(igen)
	<< "\t from hard process: " << (statusFlag & 1)
	<< "\t Is Prompt: " << isPrompt
	 << "\t Pt = " << (mcPt)->at(igen)
	<< "\t Eta = " << (mcEta)->at(igen)
	<< endl;
  }



}
