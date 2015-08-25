#define DEBUG  false

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"

//#include "TH2D.h"

#include "WZEvent.h"

#include <iostream>
#include <fstream>
#include <string>


using namespace std;


int main(int argc, char **argv)
{

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
	cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
	abort ();
      }

  // OUTPUT ROOT FILE

  TFile * fout;
  if (gotOutput) {
    fout = new TFile(outputName, "RECREATE");    
  } else {
    fout = new TFile("wzJets-test.root", "RECREATE");
  }

  TH1D* hMet = new TH1D ("hMet", "Missing Et", 100, 0., 200.);

//  TH1D* hNZCand = new TH1D ("hNZCand", "Number Z candidates", 6, -1.5, 4.5);
  TH1D* hMassZCand = new TH1D ("hMassZCand", "Z Candidate Mass", 80, 50., 130.);
//  TH1D * hZmassMu1         = new TH1D ("hZmassMu1", "hZmassMu1", 100, 60, 120);
//  TH1D * hZmassEl1         = new TH1D ("hZmassEl1", "hZmassEl1", 100, 60, 120);

  // INPUT TREES

  vector<TString> inputName;
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
    cout << "Got no input ROOT file: quit \n";
    return 1;
  }

  for (unsigned int input = 0; input < inputName.size(); input++) {
    cout << "Adding: " << input << endl;
    wz.Add(inputName[input]);
    cout << "Added \n";
  }

  TTree* wz_tTree=(TTree*)&wz;
  WZEvent* cWZ= new WZEvent(wz_tTree);
  Int_t events = wz_tTree->GetEntries();
  
  cout << "Total number of events: " << events << endl;

  unsigned int nFailsPreselection = 0;
  unsigned int nUndefined = 0;
  unsigned int nPreselection = 0;
  unsigned int nZSelection = 0;
  unsigned int nWSelection = 0;
  unsigned int nFullSelection = 0;

  // Event loop

  int nEvents = 0;

  for  (Int_t k = 0; k < events /* && k < 200000 */; k++) {
    cerr <<  "  " << int(100 * 100 * (k+1.) / events + 0.5) / 100. << " %              \r" << flush;
    nEvents++;

    if ( !(k % 100000) ) {
      cout << "Processed " << k << " events \n";
    }

    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

/*    hNele->Fill(cWZ->nEle);
    hNmu->Fill(cWZ->nMu);

    for (unsigned int iEle = 0; iEle < cWZ->elePt->size(); iEle++) {
      float ePt = cWZ->elePt->at(iEle);
      hElePt->Fill(ePt);
    }

    for (unsigned int iMu = 0; iMu < cWZ->muPt->size(); iMu++) {
      float mPt = cWZ->muPt->at(iMu);
      hMuPt->Fill(mPt);
    }
*/

    if (cWZ->PassesFullSelection())  nFullSelection++;
    if (cWZ->GetSelectionLevel() == Undefined)          nUndefined++;
    if (cWZ->GetSelectionLevel() == FailsPreselection)  nFailsPreselection++;
    if (cWZ->GetSelectionLevel() == Preselection)       nPreselection++;
    if (cWZ->GetSelectionLevel() == ZSelection)         nZSelection++;
    if (cWZ->GetSelectionLevel() == WSelection)         nWSelection++;

/*
    const double massZCand = (*(cWZ->fLeptons.at(fZLeptonsIndex.first)) +
                              *(cWZ->fLeptons.at(fZLeptonsIndex.second)).M();
    hMassZCand->Fill(massZcand);
  }
*/
    hMet->Fill(cWZ->pfMET);

//    hNZCand->Fill(cWZ->nZCand);
//    if (cWZ->nZCand) {
//      nZSelectedCheck++;
//      for (unsigned int i = cWZ->nZCand-1; i < cWZ->massZCand.size(); i++) {
//        hMassZCand->Fill(cWZ->massZCand.at(i));
//      }
//      if (cWZ->nZCand == 1) {
//        nZ1++;
//      } else if (cWZ->nZCand == 2) {
//        nZ2++;
//      } else {
//        nZ3++;
//      }
//    }
  }

  cout << " Selection Level" << endl
       << " Undefined : " << nUndefined << endl
       << " FailsPreselection : " << nFailsPreselection << endl
       << " Preselection : " << nPreselection << endl
       << " Z Selection : " << nZSelection << endl
       << " W Selection : " << nWSelection << endl
       << " Complete Selection : " << nFullSelection << endl;

  fout->cd();

/*  hNele->Write();
  hNmu->Write();
  hElePt->Write();
  hMuPt->Write(); */
  hMet->Write();
//  hNZCand->Write();
  hMassZCand->Write();

  fout->Close();
}
