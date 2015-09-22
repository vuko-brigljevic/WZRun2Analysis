#include "WZSelectionAnalysis.h"

#include <sstream>

//#include "TPythia6.h"


using namespace std;


WZSelectionAnalysis::WZSelectionAnalysis(WZEvent * e, TFile * fout)
 : GenericAnalysis(e, fout)
{

}


void WZSelectionAnalysis::Init()
{
  nrAnalyzedEvents = 0;

  for (int i=0; i<4; i++) {
    yieldsByChannel[i] = 0;
    // yieldsByChannelPreselection[i] = 0;
    // yieldsByChannelZSelection[i] = 0;
    // yieldsByChannelWSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
  }

  // Setup selected event lists 
  for (int i=0; i<4; i++) {
    ostringstream fileName;
    fileName << "passedEvents_" << i+1;
    cout << "file name : " << fileName.str() << endl;
    eventLists[i].open(fileName.str().c_str());
  }
  
  // Book some histos

  hTrueDecayMode = bookTH2D("hTrueDecayMode","W Decay vs Z decay",20,-1.5,18.5,20,-1.5,18.5);
  hRecoVsTrueChannel = bookTH2D("hRecoVsTrueChannel","WZ Reco vs True channel"
				,41,-.5,40.5,10,-0.5,9.5);


}


void WZSelectionAnalysis::EventAnalysis()
{

  std::cout << "EVENT ANALYSIS \n ===================== \n";

  //  FindTrueWZLeptons();

  std::cout << "True W Decay: " << fWZEvent->GetGenWDecayFlavor() << std::endl;  
  std::cout << "True Z Decay: " << fWZEvent->GetGenZDecayFlavor() << std::endl;  

  int genWDecayMode = fWZEvent->GetGenWDecayFlavor();
  int genZDecayMode = fWZEvent->GetGenZDecayFlavor();
  
  int wzGenDecayMode = 0;
  std::cout << "DETERMINE WZDECYAMODE \n";

  if (abs(genZDecayMode) == 11) { 
    wzGenDecayMode += 10;
  }   if (abs(genZDecayMode) == 13) { 
    wzGenDecayMode += 20;
  }   if (abs(genZDecayMode) == 15) { 
    wzGenDecayMode += 30;
  }
  if (abs(genWDecayMode) == 11) { 
    wzGenDecayMode += 1;
  }   if (abs(genWDecayMode) == 13) { 
    wzGenDecayMode += 2;
  }   if (abs(genWDecayMode) == 15) { 
    wzGenDecayMode += 3;
  }

  FinalState genFinalState=undefined;

  if (abs(genZDecayMode) == 11) {   
    if (abs(genWDecayMode) == 11) { 
      genFinalState = eee; 
    } else if (abs(genWDecayMode) == 13) { 
      genFinalState = eem; 
    }
  }
  if (abs(genZDecayMode) == 13) {   
    if (abs(genWDecayMode) == 11) { 
      genFinalState = mme;
    } else if (abs(genWDecayMode) == 13) { 
      genFinalState = mmm;
    }
  }


  hTrueDecayMode->Fill(fWZEvent->GetGenZDecayFlavor(),
		       fWZEvent->GetGenWDecayFlavor());

  // MC Printout

  fWZEvent->DumpGenEvent(cout);

  nrAnalyzedEvents++;


  if (fWZEvent->PassesFullSelection()) {
    yieldsByChannelFullSelection[fWZEvent->GetFinalState()-1]++;

  }

  if (fWZEvent->IsInGenXSPhaseSpace() ) {
      genYieldsByChannel[genFinalState]++;
  }


  // PLOT RECO VS TRUE DECAY MODE


  if (! fWZEvent->PassesFullSelection() ) return;

  yieldsByChannel[fWZEvent->GetFinalState()-1]++;  

  if (fWZEvent->GetFinalState() == genFinalState) {
    signalYields[genFinalState]++;
  }


  fWZEvent->DumpEvent(eventLists[fWZEvent->GetFinalState()-1],1);

  hRecoVsTrueChannel->Fill(wzGenDecayMode, fWZEvent->GetFinalState());



}


void WZSelectionAnalysis::Finish()
{
  cout << "Total number of analyzed events : " <<   nrAnalyzedEvents << endl << endl;

  for (int i=0; i<4; i++)
    cout << "Yield for Complete Selection " << i << "\t:\t" << yieldsByChannelFullSelection[i] << "\n";
  cout << endl;

  for (int i=0; i<4; i++) {
    cout << "Yield for final state " << i << "\t:\t" << yieldsByChannel[i] << endl;
  }

  std::cout << "A*eff \n\n";
  for (int i=1; i<5; i++) {
    double Aeff = (double) signalYields[i] / (double) genYieldsByChannel[i];
    double dAeff = sqrt( Aeff*(1-Aeff)/(double) genYieldsByChannel[i]);
    std::cout << "Channel "  << i 
	      << "\t selected : " <<  signalYields[i]
	      << "\t generated : " <<  genYieldsByChannel[i]
	      << "\t A*eff = " << Aeff 
	      << " +/- " << dAeff << std::endl;
  }
}

