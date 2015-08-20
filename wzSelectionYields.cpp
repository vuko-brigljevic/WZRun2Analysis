#define DEBUG  false

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "WZEvent.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iterator>
#include <set>


// function declaration goes here:


int main(int argc, char **argv)
{
  using namespace std;

//  bool debug=DEBUG;

  char* outputName(0);
  char* inputFileName(0);
//  char * binningFileName(0);

  char* fileList(0);
  bool useInputList = false;

  bool gotInput  = false;
  bool gotOutput = false;
//  bool gotHistoBinning = false;
//  bool gotSystematicsConfig = false;
//  char * systConfigFile(0);
  char c;

  while ((c = getopt (argc, argv, "i:o:l:")) != -1)
    switch (c)
      {
      case 'o':
	gotOutput = true;
	outputName = new char[strlen(optarg)+1];
	strcpy(outputName,optarg);
	break;
      case 'i':
	gotInput = true;
	inputFileName = new char[strlen(optarg)+1];
	strcpy(inputFileName,optarg);
	break;
      case 'l':
	useInputList = true;
	fileList = new char[strlen(optarg)+1];
	strcpy(fileList,optarg);
	break;
      default:
	std::cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
	abort ();
      }

  // OUTPUT ROOT FILE

  TFile * fout;
  if (gotOutput) {
    fout = new TFile(outputName, "RECREATE");    
  } else {
    fout = new TFile("wzJets-test.root", "RECREATE");
  }

  TH1F* hNele = new TH1F ("hNele", "Number of electrons", 15, -0.5, 14.5);
  TH1F* hNmu = new TH1F ("hNmu", "Number of muons", 15, -0.5, 14.5);
  TH1F* hElePt = new TH1F ("hElePt", "Electron Pt", 100, 0., 100.);
  TH1F* hMuPt = new TH1F ("hMuPt", "Muon Pt", 100, 0., 200.);
  TH1F* hMet = new TH1F ("hMet", "Missing Et", 100, 0., 200.);
//  TH1D * hZmassMu1         = new TH1D ("hZmassMu1", "hZmassMu1", 100, 60, 120);
//  TH1D * hZmassEl1         = new TH1D ("hZmassEl1", "hZmassEl1", 100, 60, 120);

  // INPUT TREES

  std::vector<TString> inputName;
  TChain wz("ggNtuplizer/EventTree");

  if (useInputList) {
    ifstream list(fileList);
    TString name;
    while (list>>name) {
      inputName.push_back(name);
    }
  } else  if (gotInput) {
    inputName.push_back(inputFileName);
  } else {
    std::cout << "Got no input ROOT file: quit \n";
    return 1;
  }

  for (unsigned int input = 0; input < inputName.size(); input++) {
    std::cout << "Adding: " << input << std::endl;
    wz.Add(inputName[input]);
    std::cout << "Added \n";
  }

  TTree* wz_tTree=(TTree*)&wz;
  WZEvent* cWZ= new WZEvent(wz_tTree);
  Int_t events = wz_tTree->GetEntries();
  
  std::cout << "Total number of events: " << events << std::endl;

  unsigned int nNoSelection = 0;
  unsigned int n3Lepton = 0;
  unsigned int nPreselected = 0;

  unsigned int nUndefined = 0;
  unsigned int nEEE = 0;
  unsigned int nEEM = 0;
  unsigned int nEMM = 0;
  unsigned int nMMM = 0;

  unsigned int testPassed = 0;
  unsigned int testPassed1 = 0;
  unsigned int testPassed2 = 0;

  unsigned int nSelected =0;

  // Event loop

  int nEvents = 0;

  for  (Int_t k = 0; k < events && k < 200000; k++) {
    std::cerr <<  "  " << int(100 * 100 * (k+1.) / events + 0.5) / 100. << " %              \r" << flush;
    nEvents++;

    if ( !(k % 50000) ) {
      std::cout << "Processed " << k << " events \n";
    }

    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

    hMet->Fill(cWZ->pfMET);
    hNele->Fill(cWZ->nEle);
    hNmu->Fill(cWZ->nMu);

    for (unsigned int iEle = 0; iEle < cWZ->elePt->size(); iEle++) {
      float ePt = cWZ->elePt->at(iEle);
      hElePt->Fill(ePt);
    }

    for (unsigned int iMu = 0; iMu < cWZ->muPt->size(); iMu++) {
      float mPt = cWZ->muPt->at(iMu);
      hMuPt->Fill(mPt);
    }

    if  (cWZ->nEle + cWZ->nMu == 3) {
      testPassed++;
    }

    if  (cWZ->nEle + cWZ->nMu < 3) {
      testPassed1++;
    }

    if  (cWZ->nEle + cWZ->nMu > 3) {
      testPassed2++;
    }

    bool selected = cWZ->passesSelection();

    if (selected) {
      nSelected++;
    }

    if (cWZ->GetSelectionLevel() == selectionNotRun) {
      nNoSelection++;
    }

    if (cWZ->GetSelectionLevel() == passesThreeLeptonFilter) {
      n3Lepton++;
    }

    if (cWZ->GetSelectionLevel() == passesPreselection) {
      nPreselected++;
    }

    if (cWZ->GetFinalState() == undefined) {
      nUndefined++;
    }
    if (cWZ->GetFinalState() == eee) {
      nEEE++;
    }
    if (cWZ->GetFinalState() == eem) {
      nEEM++;
    }
    if (cWZ->GetFinalState() == mme) {
      nEMM++;
    }
    if (cWZ->GetFinalState() == mmm) {
      nMMM++;
    }
  }

  cout << "Total number of processed events : " << nEvents << "\n"
       << "Events with 3 leptons: " << testPassed << "\n"
       << "Events with more than 3 leptons: " << testPassed2 << "\n"
       << "Events with less than 3 leptons: " << testPassed1 << endl << endl;

  cout << "Passed FULL selection: " << nSelected << "\n"
       << "No selection run: " << nNoSelection << "\n"
       << "3LeptonFilter passed: " << n3Lepton << "\n"
       << "Preselected: " << nPreselected << "\n"
       << "Total : " << nNoSelection + n3Lepton + nPreselected << endl << endl;

  cout << "undefined: " << nUndefined << "\n"
       << "eee: " << nEEE << "\n"
       << "eem: " << nEEM << "\n"
       << "mme: " << nEMM << "\n"
       << "mmm: " << nMMM << endl;

  fout->cd();

  hNele->Write();
  hNmu->Write();
  hElePt->Write();
  hMuPt->Write();
  hMet->Write();

  fout->Close();
}
