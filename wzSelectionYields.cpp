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

  TH1F* hNele = new TH1F ("hNele", "Number electrons", 15, -0.5, 14.5);
  TH1F* hNmu = new TH1F ("hNmu", "Number muons", 15, -0.5, 14.5);
  TH1F* hElePt = new TH1F ("hElePt", "Electron Pt", 100, 0., 100.);
  TH1F* hMuPt = new TH1F ("hMuPt", "Muon Pt", 100, 0., 200.);
  TH1F* hMet = new TH1F ("hMet", "Missing Et", 100, 0., 200.);

  TH1F* hNZCand = new TH1F ("hNZCand", "Number Z candidates", 6, -1.5, 4.5);
  TH1F* hMassZCand = new TH1F ("hMassZCand", "Z Candidate Mass", 80, 50., 130.);
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
  unsigned int nNo3Lepton = 0;
  unsigned int n3Lepton = 0;
  unsigned int nPreselected = 0;
  unsigned int nZSelected = 0;
  unsigned int nWSelected = 0;

  unsigned int nZSelectedCheck = 0;
  unsigned int nZ1 = 0;
  unsigned int nZ2 = 0;
  unsigned int nZ3 = 0;

  unsigned int nUndefined0 = 0;

  unsigned int nUndefined1 = 0;
  unsigned int nEEE1 = 0;
  unsigned int nEEM1 = 0;
  unsigned int nEMM1 = 0;
  unsigned int nMMM1 = 0;

  unsigned int nUndefined2 = 0;
  unsigned int nEEE2 = 0;
  unsigned int nEEM2 = 0;
  unsigned int nEMM2 = 0;
  unsigned int nMMM2 = 0;

  unsigned int nUndefined3 = 0;
  unsigned int nEEE3 = 0;
  unsigned int nEEM3 = 0;
  unsigned int nEMM3 = 0;
  unsigned int nMMM3 = 0;

  unsigned int testPassed = 0;
  unsigned int testPassed1 = 0;
  unsigned int testPassed2 = 0;

  unsigned int nSelected =0;

  // Event loop

  int nEvents = 0;

  for  (Int_t k = 0; k < events /* && k < 200000 */; k++) {
    std::cerr <<  "  " << int(100 * 100 * (k+1.) / events + 0.5) / 100. << " %              \r" << flush;
    nEvents++;

    if ( !(k % 100000) ) {
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

    if (cWZ->GetSelectionLevel() == failsThreeLeptonFilter) {
      nNo3Lepton++;
    }

    if (cWZ->GetSelectionLevel() == passesThreeLeptonFilter) {
      n3Lepton++;
    }

    if (cWZ->GetFinalState() == undefined) {
      nUndefined0++;
    }

    if (cWZ->GetSelectionLevel() == passesPreselection) {
      nPreselected++;
      if (cWZ->GetFinalState() == undefined) {
        nUndefined1++;
      }
      if (cWZ->GetFinalState() == eee) {
        nEEE1++;
      }
      if (cWZ->GetFinalState() == eem) {
        nEEM1++;
      }
      if (cWZ->GetFinalState() == mme) {
      nEMM1++;
      }
      if (cWZ->GetFinalState() == mmm) {
        nMMM1++;
      }
    }

    hNZCand->Fill(cWZ->nZCand);
    if (cWZ->nZCand) {
      nZSelectedCheck++;
      for (unsigned int i = cWZ->nZCand-1; i < cWZ->massZCand.size(); i++) {
        hMassZCand->Fill(cWZ->massZCand.at(i));
      }
      if (cWZ->nZCand == 1) {
        nZ1++;
      } else if (cWZ->nZCand == 2) {
        nZ2++;
      } else {
        nZ3++;
      }
    }

    if (cWZ->GetSelectionLevel() == passesZSelection) {
      nZSelected++;
      if (cWZ->GetFinalState() == undefined) {
        nUndefined2++;
      }
      if (cWZ->GetFinalState() == eee) {
        nEEE2++;
      }
      if (cWZ->GetFinalState() == eem) {
        nEEM2++;
      }
      if (cWZ->GetFinalState() == mme) {
      nEMM2++;
      }
      if (cWZ->GetFinalState() == mmm) {
        nMMM2++;
      }
    }

    if (cWZ->GetSelectionLevel() == passesWSelection) {
      nZSelected++;
      if (cWZ->GetFinalState() == undefined) {
        nUndefined3++;
      }
      if (cWZ->GetFinalState() == eee) {
        nEEE3++;
      }
      if (cWZ->GetFinalState() == eem) {
        nEEM3++;
      }
      if (cWZ->GetFinalState() == mme) {
      nEMM3++;
      }
      if (cWZ->GetFinalState() == mmm) {
        nMMM3++;
      }
    }

  }



  cout << "Total number of processed events : " << nEvents << "\n"
       << "Events with 3 leptons: " << testPassed << "\n"
       << "Events with more than 3 leptons: " << testPassed2 << "\n"
       << "Events with less than 3 leptons: " << testPassed1 << endl << endl;

  cout << "No selection run: " << nNoSelection << "\n"
       << "Fails 3+ lepton filter: " << nNo3Lepton << "\n"
       << "3+ lepton filter passed ONLY: " << n3Lepton << "\n"
       << "Preselected: " << nPreselected + nZSelected + nWSelected << "\n"
       << "Z Selected: " << nZSelected + nWSelected << "\n"
       << "Number of events with Z candidate: " << nZSelectedCheck << "\n"
       << "out of which: " << nZ1 << " with 1, " << nZ2 << " with 2 and " << nZ3
       << " with 3 Z candidates." << "\n"
       << "W Selected: " << nWSelected << "\n"
       << "Passed FULL selection: " << nSelected << "\n"
       << "Total : " << nNoSelection + nNo3Lepton + n3Lepton +
                        nPreselected + nZSelected + nWSelected << endl << endl;

  cout << "Total undefined: " << nUndefined0 << endl << endl;

  cout << "Preselection : " << "\n"
       << "undefined: " << nUndefined1 + nUndefined2 + nUndefined3 << "\n"
       << "eee: " << nEEE1 << " + " << nEEE2 << " + " << nEEE3 << " = "
       << nEEE1 + nEEE2 + nEEE3 << "\n"
       << "eem: " << nEEM1 << " + " << nEEM2 << " + " << nEEM3 << " = "
       << nEEM1 + nEEM2 + nEEM3 << "\n"
       << "mme: " << nEMM1 << " + " << nEMM2 << " + " << nEMM3 << " = "
       << nEMM1 + nEMM2 + nEMM3 << "\n"
       << "mmm: " << nMMM1 << " + " << nMMM2 << " + " << nMMM3 << " = "
       << nMMM1 + nMMM2 + nMMM3 << endl << endl;

  cout << "Z Selection : " << "\n"
       << "undefined: " << nUndefined2 + nUndefined3 << "\n"
       << "eee: " << nEEE2 << " + " << nEEE3 << " = " << nEEE2 + nEEE3 << "\n"
       << "eem: " << nEEM2 << " + " << nEEM3 << " = " << nEEM2 + nEEM3 << "\n"
       << "mme: " << nEMM2 << " + " << nEMM3 << " = " << nEMM2 + nEMM3 << "\n"
       << "mmm: " << nMMM2 << " + " << nMMM3 << " = " << nMMM2 + nMMM3 << endl << endl;

  cout << "W Selection : " << "\n"
       << "undefined: " << nUndefined3 << "\n"
       << "eee: " << nEEE3 << "\n"
       << "eem: " << nEEM3 << "\n"
       << "mme: " << nEMM3 << "\n"
       << "mmm: " << nMMM3 << endl << endl;

  fout->cd();

  hNele->Write();
  hNmu->Write();
  hElePt->Write();
  hMuPt->Write();
  hMet->Write();
  hNZCand->Write();
  hMassZCand->Write();

  fout->Close();
}
