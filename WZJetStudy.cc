#include "WZJetStudy.h"
#include "Constants.h"
#include "TLatex.h"


using namespace std;


WZJetStudy::WZJetStudy(WZEvent* e, TFile* outputFile)
 : GenericAnalysis(e, outputFile)
{

}


void WZJetStudy::Init()
{
  nAnalyzedEvents = 0;
  nSelectedEvents = 0;
  nGoodJets = 0;

  hNJets = bookTH1D("hNJets", "All Jets", 31, -0.5, 30.5);
  hNGoodJets = bookTH1D("hNGoodJets", "Good Jets", 11, -0.5, 10.5);

  hNJetsWL = bookTH1D("hNJetsWL", "Jets with Tight Leptons", 11, -0.5, 10.5);
  hNJetsWMu = bookTH1D("hNJetsWMu", "Jets with Tight Mu", 11, -0.5, 10.5);
  hNJetsWEle = bookTH1D("hNJetsWElE", "Jets with Tight Ele", 11, -0.5, 10.5);

  hNJetsVsNL = bookTH2D("hNJetsVsNL", "nJets vs. nLeptons", 11, -0.5, 10.5, 6, -0.5, 5.5);
  hNJetsVsNMu = bookTH2D("hNJetsVsNMu", "nJets vs. nMu", 11, -0.5, 10.5, 5, -0.5, 4.5);
  hNJetsVsNEle = bookTH2D("hNJetsVsNEle", "nJets vs. nEle", 11, -0.5, 10.5, 5, -0.5, 4.5);

  hDeltaRL = bookTH1D("hDeltaRL", "#deltaR_{min} (Jet, Lepton)", 310, -0.1, 6.1);
  hDeltaRMu = bookTH1D("hDeltaRMu", "#deltaR_{min} (Jet, Mu)", 310, -0.1, 6.1);
  hDeltaREle = bookTH1D("hDeltaREle", "#deltaR_{min} (Jet, Ele)", 310, -0.1, 6.1);
}


void WZJetStudy::EventAnalysis()
{
  nAnalyzedEvents++;
  fWZEvent->PassesFullSelection();
  hNJets->Fill(fWZEvent->nJet);

  unsigned int nEle = 0;
  unsigned int nMu = 0;

  for (vector<unsigned int>::const_iterator iIt = fWZEvent->fTightLeptonsIndex.begin();
       iIt != fWZEvent->fTightLeptonsIndex.end(); ++iIt) {
    if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 11)       nEle++;
    else if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 13)  nMu++;
  }

  unsigned int nJ = 0;
  unsigned int nJL = 0;
  unsigned int nJMu = 0;
  unsigned int nJEle = 0;
  vector<TLorentzVector> jets;
  for (unsigned int i = 0; i < fWZEvent->jetPt->size(); i++) {

//  in ggNtuplizer V07-04-05 use explicit cuts
/*    if (!(fWZEvent->jetNHF->at(i) < 0.99) || !(fWZEvent->jetNEF->at(i) < 0.99) ||
        !(fWZEvent->jetCEF->at(i) < 0.99) || !(fWZEvent->jetNConstituents->at(i) > 1) ||
        fWZEvent->jetCHF->at(i) == 0 || fWZEvent->jetNCH->at(i) == 0)
      continue;
    */

//  while in version V07-04-09 use bool for 'PURE09' and 'LOOSE' (defined with cuts above)
    for (vector<bool>::const_iterator bIt = fWZEvent->jetPFLooseId->begin();
         bIt != fWZEvent->jetPFLooseId->end(); ++bIt)
      if (!(*bIt))  continue;

    const double ptJet = fWZEvent->jetPt->at(i);
    const double etaJet = fWZEvent->jetEta->at(i);
    if (!(ptJet > JET_PTMIN) || !(abs(etaJet) < JET_ETAMAX))  continue;
    nJ++;
    if (fWZEvent->fTightLeptonsIndex.size() == 0)  continue;
    const double phiJet = fWZEvent->jetPhi->at(i);
//    const double eJet = fWZEvent->jetEn->at(i); // present only in V07-04-09+, otherwise use M=0
    TLorentzVector lJet;
//    lJet.SetPtEtaPhiE(ptJet, etaJet, phiJet, eJet);
    lJet.SetPtEtaPhiM(ptJet, etaJet, phiJet, 0);
    jets.push_back(lJet);

    bool hasMu = false;
    bool hasEle = false;
    for (vector<unsigned int>::const_iterator iIt = fWZEvent->fTightLeptonsIndex.begin();
         iIt != fWZEvent->fTightLeptonsIndex.end(); ++iIt) {
      if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 11)       hasEle = true;
      else if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 13)  hasMu = true;
    }
    if (hasMu || hasEle)  nJL++;
    if (hasMu) nJMu++;
    if (hasEle) nJEle++;
  }
  nGoodJets += jets.size();
  hNGoodJets->Fill(nJ);
  if (jets.size() && fWZEvent->fTightLeptonsIndex.size()) {
    nSelectedEvents++;
    hNJetsWL->Fill(nJL);
    hNJetsWMu->Fill(nJMu);
    hNJetsWEle->Fill(nJEle);
    hNJetsVsNL->Fill(jets.size(), fWZEvent->fTightLeptonsIndex.size());
    hNJetsVsNMu->Fill(jets.size(), nMu);
    hNJetsVsNEle->Fill(jets.size(), nEle);
  }

  vector<double> dRJetL;
  vector<double> dRJetMu;
  vector<double> dRJetEle;
  for (vector<TLorentzVector>::const_iterator lIt = jets.begin(); lIt != jets.end(); ++lIt) {
    for (vector<unsigned int>::const_iterator iIt = fWZEvent->fTightLeptonsIndex.begin();
         iIt != fWZEvent->fTightLeptonsIndex.end(); ++iIt) {
      const double dRjl = fWZEvent->fLeptons.at(*iIt)->DeltaR(*lIt);
      dRJetL.push_back(dRjl);
      if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 11)       dRJetEle.push_back(dRjl);
      else if (fWZEvent->fLeptons.at(*iIt)->GetPdgId() == 13)  dRJetMu.push_back(dRjl);
    }

    if (dRJetL.size() > 0)    hDeltaRL->Fill(*(min_element(dRJetL.begin(), dRJetL.end())));
    if (dRJetMu.size() > 0)   hDeltaRMu->Fill(*(min_element(dRJetMu.begin(), dRJetMu.end())));
    if (dRJetEle.size() > 0)  hDeltaREle->Fill(*(min_element(dRJetEle.begin(), dRJetEle.end())));

    dRJetL.clear();
    dRJetMu.clear();
    dRJetEle.clear();
  }

  jets.clear();
}


void WZJetStudy::Finish()
{
  cout << endl << "DONE." << endl;
}
