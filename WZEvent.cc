#include "WZEvent.h"
#include "Constants.h"


using namespace std;


WZEvent::WZEvent(TTree* tree) :
  WZBASECLASS(tree)
{
  // Set static pointer to Lepton class
  Lepton::SetWZEvent(this);
}


void WZEvent::Clear()
{
  fHLT25ns.clear();
  fHLT50ns.clear();

  fFinalState = undefined;
  fSelectionLevel = Undefined;

  // Empty leptons
  vector<Lepton*>::iterator lIt;
  for (lIt = fLeptons.begin(); lIt != fLeptons.end(); ) {
    delete *lIt;  
    lIt = fLeptons.erase(lIt);
  }

  fTightLeptonsIndex.clear();
  fZLeptonsIndex = make_pair(999, 999);
  fWLeptonIndex = 999;
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
    Electron* ele = new Electron(indexEle, elePt->at(indexEle), eleEta->at(indexEle),
                                 elePhi->at(indexEle), eleCharge->at(indexEle));
    fLeptons.push_back(ele);
  }

  // Muons
  for (unsigned int indexMu = 0; indexMu < muCharge->size(); indexMu++) {
    Muon* mu = new Muon(indexMu, muPt->at(indexMu), muEta->at(indexMu),
                        muPhi->at(indexMu), muCharge->at(indexMu));
    fLeptons.push_back(mu);
  }
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


void WZEvent::DumpEvent(ostream& out, int verbosity)
{
  out << run << ":" << lumis << ":" << event;

  vector<unsigned int> indexWZLeptons =
  { fZLeptonsIndex.first, fZLeptonsIndex.second, fWLeptonIndex };

  if (verbosity > 1) {
    for (vector<unsigned int>::iterator iIt = indexWZLeptons.begin();
         iIt != indexWZLeptons.end(); ++iIt) {
      if (*iIt < fLeptons.size()) {
        out << " Lepton : " << distance(indexWZLeptons.begin(), iIt)
            << " Pt = " << fLeptons.at(*iIt)->Pt()
            << " Eta = " << fLeptons.at(*iIt)->Eta();
      }
    }
  }

  out << endl;
}
