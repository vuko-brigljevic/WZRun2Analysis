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
    hMETWMt[i] = bookTH1D(("hMETWMt_" + boost::lexical_cast<string>(i)).c_str(),
                      "MET+W Mt", 100, 0, 200);
    hMETMt[i] = bookTH1D(("hMETMt_" + boost::lexical_cast<string>(i)).c_str(),
                       "MET Mt", 100, 0, 200);
    hWMt[i] = bookTH1D(("hWMt_" + boost::lexical_cast<string>(i)).c_str(),
                        "W Mt", 100, 0, 200);
    hZl1pt[i] = bookTH1D(("hZl1pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton1 pt", 100, 0, 200);
    hZl2pt[i] = bookTH1D(("hZl2pt_" + boost::lexical_cast<string>(i)).c_str(),
                         "Z Lepton2 pt", 100, 0, 200);
    hWlpt[i] = bookTH1D(("hWlpt_" + boost::lexical_cast<string>(i)).c_str(),
                        "W Lepton pt", 100, 0, 200);
    hNJets[i] = bookTH1D(("hNJets_" + boost::lexical_cast<string>(i)).c_str(),
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
    outputFileName << "output/yields/selectedEvents_" << i << "_testNew2.txt";
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

  const double pxMET = fWZEvent->pfMET * cos(fWZEvent->pfMETPhi);
  const double pyMET = fWZEvent->pfMET * sin(fWZEvent->pfMETPhi);
  TLorentzVector lMET(pxMET, pyMET, 0., fWZEvent->pfMET);
  TLorentzVector lMET1;
  lMET1.SetPtEtaPhiM(fWZEvent->pfMET, 0, fWZEvent->pfMETPhi, 0);
  const double mtMETW = (lMET + *(fWZEvent->GetWLepton())).Mt();
  const double mtMET = lMET.Mt();
  const double mtW = fWZEvent->GetWLepton()->Mt();
  
  const double mt = sqrt(2 * fWZEvent->pfMET * fWZEvent->GetWLepton()->Pt() *
                             (1 - cos(fWZEvent->GetWLepton()->DeltaPhi(lMET))));

  unsigned int nSelectedJets = 0;
  for (unsigned int i = 0; i < fWZEvent->jetPt->size(); i++) {
    if ((fWZEvent->jetPt->at(i) > JET_PTMIN) &&
        (abs(fWZEvent->jetEta->at(i)) < JET_ETAMAX) &&
        (fWZEvent->jetNConstituents->at(i) > 10))
      nSelectedJets++;
  }

  hZmass[4]->Fill(massZ);
  hZpt[4]->Fill(ptZ);
  hMET[4]->Fill(fWZEvent->pfMET);
  hMt[4]->Fill(mt);
  hMETWMt[4]->Fill(mtMETW);
  hMETMt[4]->Fill(mtMET);
  hWMt[4]->Fill(mtW);
  hZl1pt[4]->Fill(fWZEvent->GetZLeptons().first->Pt());
  hZl2pt[4]->Fill(fWZEvent->GetZLeptons().second->Pt());
  hWlpt[4]->Fill(fWZEvent->GetWLepton()->Pt());
  hNJets[4]->Fill(nSelectedJets);

  hZmass[fWZEvent->GetFinalState()-1]->Fill(massZ);
  hZpt[fWZEvent->GetFinalState()-1]->Fill(ptZ);
  hMET[fWZEvent->GetFinalState()-1]->Fill(fWZEvent->pfMET);
  hMt[fWZEvent->GetFinalState()-1]->Fill(mt);
  hMETWMt[fWZEvent->GetFinalState()-1]->Fill(mtMETW);
  hMETMt[fWZEvent->GetFinalState()-1]->Fill(mtMET);
  hWMt[fWZEvent->GetFinalState()-1]->Fill(mtW);
  hZl1pt[fWZEvent->GetFinalState()-1]->Fill((fWZEvent->GetZLeptons().first)->Pt());
  hZl2pt[fWZEvent->GetFinalState()-1]->Fill((fWZEvent->GetZLeptons().second)->Pt());
  hWlpt[fWZEvent->GetFinalState()-1]->Fill((fWZEvent->GetWLepton())->Pt());
  hNJets[fWZEvent->GetFinalState()-1]->Fill(nSelectedJets);

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

