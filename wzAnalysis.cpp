#define DEBUG  false

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

// #include "RooUnfoldResponse.h"

// Replace this with the new tree
#include "WZEvent.h"

#include "WZSelectionAnalysis.h"

//#include "WZ.h"
//#include "WZ2012Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>


// TH1D * createJetPtHisto(std::string key, std::string title) {

//   TH1D * h = new TH1D(key.c_str(),title.c_str(), 10,0., 500.);
//   return h;

// }


//function declaration goes here:
//
//
//

int main(int argc, char **argv)
{
  using namespace std;


//  bool debug=DEBUG;

  char * outputName(0);
  char * inputFileName(0);
//  char * binningFileName(0);

  char * fileList;
  bool useInputList=false;

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


  TH1F * hMet              = new TH1F ("hMet","Missing ET",100,0.,200.);
  TH1F * hNele             = new TH1F ("hNele","Nr of electrons",15,-0.5,14.5);
  TH1F * hNmu              = new TH1F ("hNmu","Nr of muons",15,-0.5,14.5);

//  TH1D * hZmassMu1         = new TH1D ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
//  TH1D * hZmassEl1         = new TH1D ("hZmassEl1", "hZmassEl1", 100, 60, 120);  

//  TH2D * hRecoVsGenChannel = new TH2D ("hRecoVsGenChannel", "Reco vs Gen Channel", 5, -0.5, 4.5, 5, -0.5, 4.5);


  // INPUT TREES

  std::vector<TString>inputName;
  TChain wz("ggNtuplizer/EventTree");

  if (useInputList) {
    ifstream list(fileList);
    TString name;
    while (list>>name) {
      //      cout << "adding name: " << name << endl;
      inputName.push_back(name);
    }
  } else   if (gotInput) {
    inputName.push_back(inputFileName);
  } else {
    std::cout << "Got no input ROOT file: quit \n";
    return 1;
  }
  for (unsigned int input=0; input < inputName.size(); input++){
    std::cout << "Adding: " << input << std::endl;
    wz.Add(inputName[input]);
    std::cout << "Added \n";
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZEvent *cWZ= new WZEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;

  //  WZAnalysis   genAnalysis(cWZ);
  //  genAnalysis.Init();

  WZSelectionAnalysis * myAnalysis = new WZSelectionAnalysis(cWZ,fout);
  myAnalysis->Init();

  int nselected = 0;


  //
  // Event loop
  // 
  int nevents=0;

  for  (Int_t k = 0; k<events ;k++) {
  //  for  (Int_t k = 0; k<events; k++) {

    nevents++;

    if ( !(k%100000)  ) std::cout << "Processed " << k << " events \n";
    
    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

    hMet->Fill(cWZ->pfMET);
    hNele->Fill(cWZ->nEle);
    hNmu->Fill(cWZ->nMu);

    myAnalysis->EventAnalysis();


  }

  std::cout << "Events : " << nevents
	    << "\t Total Selected = " << nselected << std::endl;


  myAnalysis->Finish();

  fout->cd();
  hMet->Write();

  hMet->Write();
  hNele->Write();
  hNmu->Write();

  fout->Close();
}
