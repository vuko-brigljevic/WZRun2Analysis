#define DEBUG  false

#include "WZSelectionAnalysis.h"

#include <sstream>
#include <string>
#include <set>


using namespace std;


int main(int argc, char **argv)
{
//  bool debug=DEBUG;

  char * outputName(0);
  char * inputFileName(0);
//  char * binningFileName(0);

  char * fileList;
  bool useInputList=false;

  bool gotInput  = false;
  bool gotOutput = false;
  int  maxEvents = -1;
//  bool gotHistoBinning = false;
//  bool gotSystematicsConfig = false;
//  char * systConfigFile(0);
  char c;

  while ((c = getopt (argc, argv, "i:o:l:M:")) != -1)
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
      case 'M':
	maxEvents = atoi(optarg);
	break;
      default:
	cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
	abort ();
      }

  // OUTPUT ROOT FILE

  std::cout << "Max events to process : " << maxEvents << std::endl;

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

  vector<TString>inputName;
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
    cout << "Got no input ROOT file: quit \n";
    return 1;
  }
  for (unsigned int input=0; input < inputName.size(); input++){
    cout << "Adding: " << input << endl;
    wz.Add(inputName[input]);
    cout << "Added \n";
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZEvent *cWZ= new WZEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  cout<<"number of events: "<<events << endl;

  //  WZAnalysis   genAnalysis(cWZ);
  //  genAnalysis.Init();

  WZSelectionAnalysis* myAnalysis = new WZSelectionAnalysis(cWZ,fout);
  myAnalysis->Init();

  int nselected = 0;
  
  //
  // Event loop
  // 
  int nevents=0;

  for  (Int_t k = 0; k<events && (maxEvents<0 || k<maxEvents) ;k++) {
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

  cout << "Events : " << nevents
	    << "\t Total Selected = " << nselected << endl;


  myAnalysis->Finish();
  delete myAnalysis;

  fout->cd();
  hMet->Write();

  hMet->Write();
  hNele->Write();
  hNmu->Write();

  fout->Close();
}
