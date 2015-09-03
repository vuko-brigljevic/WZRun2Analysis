#include "WZSelectionYields.h"
#include "Constants.h"

#include <cmath>
#include <sstream>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


WZSelectionYields::WZSelectionYields(WZEvent* e, TFile* outputFile)
 : GenericAnalysis(e, outputFile)
{

}


void WZSelectionYields::Init()
{
  nAnalyzedEvents = 0;
  nSelectedEvents = 0;

  for (unsigned int i = 0; i <= 4; i++) {
    hZmass[i] = bookTH1D(("hZmass_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z mass", 64, 58, 122);
    hZpt[i] = bookTH1D(("hZpt_" + boost::lexical_cast<string>(i)).c_str(),
                       "Z pt", 100, 0, 200);
    hMET[i] = bookTH1D(("hMET_" + boost::lexical_cast<string>(i)).c_str(),
                       "Missing Et", 100, 0, 200);
    hMt[i] = bookTH1D(("hMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "Mt", 100, 0, 200);
    hZl1pt[i] = bookTH1D(("hZl1pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton1 pt", 100, 0, 200);
    hZl2pt[i] = bookTH1D(("hZl2pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton2 pt", 100, 0, 200);
    hWlpt[i] = bookTH1D(("hWlpt_" + boost::lexical_cast<string>(i)).c_str(),
                        "W Lepton pt", 100, 0, 200);
    hNJets[i] = bookTH1D(("hNJets_" + boost::lexical_cast<string>(i)).c_str(),
                         "No. Jets", 11, -0.5, 10.5);
    hNJetsNoMuIso[i] = bookTH1D(("hNJetsNoMuIso_" + boost::lexical_cast<string>(i)).c_str(),
                                "No. Jets", 11, -0.5, 10.5);
    hNJetsNoEleIso[i] = bookTH1D(("hNJetsNoEleIso_" + boost::lexical_cast<string>(i)).c_str(),
                                 "No. Jets", 11, -0.5, 10.5);
    hNJetsNoIso[i] = bookTH1D(("hNJetsNoIso_" + boost::lexical_cast<string>(i)).c_str(),
                              "No. Jets", 11, -0.5, 10.5);
  }

  for (int i = 0; i <= 4; i++) {
    yieldsByChannelPreselection[i] = 0;
    yieldsByChannelZSelection[i] = 0;
    yieldsByChannelWSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
  }

  // Setup selected event lists 
  for (int i = 1; i <= 4; i++) {
    ostringstream outputFileName;
    outputFileName << "output/yields/test_" << i << ".txt";
    cout << "File name : " << outputFileName.str() << endl;
    eventLists[i-1].open(outputFileName.str().c_str());
  }
}


void WZSelectionYields::EventAnalysis()
{
  nAnalyzedEvents++;

  if (fWZEvent->PassesFullSelection()) {
    yieldsByChannelFullSelection[fWZEvent->GetFinalState()-1]++;
    yieldsByChannelFullSelection[4]++;
  }
  if (fWZEvent->PassesWSelection()) {
    yieldsByChannelWSelection[fWZEvent->GetFinalState()-1]++;
    yieldsByChannelWSelection[4]++;
  }
  if (fWZEvent->PassesZSelection()){
    yieldsByChannelZSelection[fWZEvent->GetFinalState()-1]++;
    yieldsByChannelZSelection[4]++;
  }
  if (fWZEvent->PassesPreselection()) {
    yieldsByChannelPreselection[fWZEvent->GetFinalState()-1]++;
    yieldsByChannelPreselection[4]++;
  }

  if (!(fWZEvent->PassesFullSelection()))  return;

  nSelectedEvents++;
  const double massZ = (*(fWZEvent->GetZLeptons().first) + *(fWZEvent->GetZLeptons().second)).M();
  const double ptZ = (*(fWZEvent->GetZLeptons().first) + *(fWZEvent->GetZLeptons().second)).Pt();
  const double ptZl1 = fWZEvent->GetZLeptons().first->Pt();
  const double ptZl2 = fWZEvent->GetZLeptons().second->Pt();
  const double ptWl = fWZEvent->GetWLepton()->Pt();
  const double met = fWZEvent->pfMET;
  const double phiMET= fWZEvent->pfMETPhi;
  const double pxMET = met * cos(phiMET);
  const double pyMET = met * sin(phiMET);
  TLorentzVector lMET(pxMET, pyMET, 0., met);
  const double mt = sqrt(2 * met * ptWl * (1 - cos(fWZEvent->GetWLepton()->DeltaPhi(lMET))));

// Counting accompanying jets
  unsigned int nSelectedJets = 0;
  unsigned int nSelectedJetsNoMuIso = 0;
  unsigned int nSelectedJetsNoEleIso = 0;
  unsigned int nSelectedJetsNoIso = 0;
  for (unsigned int i = 0; i < fWZEvent->jetPt->size(); i++) {
    /*
    if (!(fWZEvent->jetNHF->at(i) < 0.99) || !(fWZEvent->jetNEF->at(i) < 0.99) ||
        !(fWZEvent->jetCEF->at(i) < 0.99) || !(fWZEvent->jetNConstituents->at(i) > 1) ||
        fWZEvent->jetCHF->at(i) == 0 || fWZEvent->jetNCH->at(i) == 0)
      continue;
//  above jet variables were removed from ggNtuple with V07-04-09-00 */

    const double ptJet = fWZEvent->jetPt->at(i);
    const double etaJet = fWZEvent->jetEta->at(i);
    if (!(ptJet > JET_PTMIN) || !(abs(etaJet) < JET_ETAMAX))  continue;

    nSelectedJetsNoIso++;
    const double phiJet = fWZEvent->jetPhi->at(i);
    const double eJet = fWZEvent->jetEn->at(i);  // present only in V07-04-09+, use M=0 instead
    TLorentzVector lJet;
    lJet.SetPtEtaPhiE(ptJet, etaJet, phiJet, eJet);
    const double deltaRJetWl = fWZEvent->GetWLepton()->DeltaR(lJet);
    const double deltaRJetZl1 = fWZEvent->GetZLeptons().first->DeltaR(lJet);
    const double deltaRJetZl2 = fWZEvent->GetZLeptons().second->DeltaR(lJet);

    if (deltaRJetWl > LEPTONJET_DELTARMIN &&
        deltaRJetZl1 > LEPTONJET_DELTARMIN && deltaRJetZl2 > LEPTONJET_DELTARMIN)
      nSelectedJets++;

    if (fWZEvent->GetFinalState() == mmm) {
      nSelectedJetsNoMuIso++;
      if (deltaRJetWl > ELEJET_DELTARMIN &&
          deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoEleIso++;
    }
    else if (fWZEvent->GetFinalState() == mme) {
      if (deltaRJetWl > ELEJET_DELTARMIN)  nSelectedJetsNoMuIso++;
      else if (deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoEleIso++;
    }
    else if (fWZEvent->GetFinalState() == eem) {
      if (deltaRJetWl > ELEJET_DELTARMIN)  nSelectedJetsNoEleIso++;
      else if (deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoMuIso++;
    }
    else if (fWZEvent->GetFinalState() == eee) {
      nSelectedJetsNoEleIso++;
      if (deltaRJetWl > ELEJET_DELTARMIN &&
          deltaRJetZl1 > ELEJET_DELTARMIN && deltaRJetZl2 > ELEJET_DELTARMIN)
        nSelectedJetsNoMuIso++;
    }
    else  continue;
  }

  hZmass[4]->Fill(massZ);
  hZpt[4]->Fill(ptZ);
  hMET[4]->Fill(met);
  hMt[4]->Fill(mt);
  hZl1pt[4]->Fill(ptZl1);
  hZl2pt[4]->Fill(ptZl2);
  hWlpt[4]->Fill(ptWl);
  hNJets[4]->Fill(nSelectedJets);
  hNJetsNoMuIso[4]->Fill(nSelectedJetsNoMuIso);
  hNJetsNoEleIso[4]->Fill(nSelectedJetsNoEleIso);
  hNJetsNoIso[4]->Fill(nSelectedJetsNoIso);

  hZmass[fWZEvent->GetFinalState()-1]->Fill(massZ);
  hZpt[fWZEvent->GetFinalState()-1]->Fill(ptZ);
  hMET[fWZEvent->GetFinalState()-1]->Fill(met);
  hMt[fWZEvent->GetFinalState()-1]->Fill(mt);
  hZl1pt[fWZEvent->GetFinalState()-1]->Fill(ptZl1);
  hZl2pt[fWZEvent->GetFinalState()-1]->Fill(ptZl2);
  hWlpt[fWZEvent->GetFinalState()-1]->Fill(ptWl);
  hNJets[fWZEvent->GetFinalState()-1]->Fill(nSelectedJets);
  hNJetsNoMuIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoMuIso);
  hNJetsNoEleIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoEleIso);
  hNJetsNoIso[fWZEvent->GetFinalState()-1]->Fill(nSelectedJetsNoIso);

  fWZEvent->DumpEvent(eventLists[fWZEvent->GetFinalState()-1], 1);
}


void WZSelectionYields::Finish()
{
  cout << "CHANNEL \tPreselection \tZ Selection \tW Selection \tFull Selection" << "\n";
  for (int i = 1; i <= 5; i++) {
    cout << i << "\t\t" << yieldsByChannelPreselection[i-1]
              << "\t\t" << yieldsByChannelZSelection[i-1]
              << "\t\t" << yieldsByChannelWSelection[i-1]
              << "\t\t" << yieldsByChannelFullSelection[i-1] << "\n";
  }
  cout << endl;
}
