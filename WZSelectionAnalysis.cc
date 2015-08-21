#include "WZSelectionAnalysis.h"

#include <iostream>
#include <fstream>
#include <sstream>

//#include "TPythia6.h"


  
WZSelectionAnalysis::WZSelectionAnalysis(WZEvent * e, TFile * fout) 

  :   GenericAnalysis(e, fout)

{

}



void WZSelectionAnalysis::Init() {


  nrAnalyzedEvents = 0;

  for (int i=0; i<4; i++) {
    yieldsByChannel[i] = 0;
  }

  // Setup selected event lists 
  for (int i=0; i<4; i++) {
    std::ostringstream fileName;
    fileName << "passedEvents_" << i+1;
    std::cout << "file name : "  << fileName.str() << std::endl;
    eventLists[i].open(fileName.str().c_str());
  }


}


void WZSelectionAnalysis::EventAnalysis() {

  nrAnalyzedEvents++;

  if (! wzevt->passesFullSelection() ) return;

  yieldsByChannel[wzevt->GetFinalState()-1]++;  

  wzevt->DumpEvent(eventLists[wzevt->GetFinalState()-1],2);


  // MC Analysis

  std::cout << "MC Tree \n";
  std::cout << "===============\n";
  for (int igen = 0; igen<wzevt->nMC; igen++) {
    //    char  pyName[20];
    //    TPythia::Pyname((wzevt->mcPID)->at(igen),pyName);

    std::cout << "Gen Particle: " << (wzevt->mcPID)->at(igen)
	      << "\t status   : " << (wzevt->mcStatus)->at(igen)
      //	      << "  " << pyName
	      << "\t mother ID: " << (wzevt->mcMomPID)->at(igen)
	      << "\t GrandMa  ID: " << (wzevt->mcGMomPID)->at(igen)
	      << std::endl;


  }


}


void WZSelectionAnalysis::Finish() {

  std::cout << "Total number of analyzed events : " <<   nrAnalyzedEvents << std::endl << std::endl;


  for (int i=0; i<4; i++) {
    std::cout << "Yield for final state " << i << "\t:\t" << yieldsByChannel[i] << std::endl;
  }

}
  

