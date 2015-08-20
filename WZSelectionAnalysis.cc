#include "WZSelectionAnalysis.h"

#include <iostream>
#include <fstream>
#include <sstream>


  
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

  wzevt->DumpEvent(eventLists[wzevt->GetFinalState()-1]);

}


void WZSelectionAnalysis::Finish() {

  std::cout << "Total number of analyzed events : " <<   nrAnalyzedEvents << std::endl << std::endl;


  for (int i=0; i<4; i++) {
    std::cout << "Yield for final state " << i << "\t:\t" << yieldsByChannel[i] << std::endl;
  }

}
  

