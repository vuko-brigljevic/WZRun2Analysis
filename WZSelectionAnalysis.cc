#include "WZSelectionAnalysis.h"

#include <iostream>

  
WZSelectionAnalysis::WZSelectionAnalysis(WZEvent * e, TFile * fout) 

  :   GenericAnalysis(e, fout)

{



}



void WZSelectionAnalysis::Init() {


  for (int i=0; i<4; i++) {
    yieldsByChannel[4] = 0;
  }

}


void WZSelectionAnalysis::EventAnalysis() {


  if (! wzevt->passesFullSelection() ) return;

  yieldsByChannel[wzevt->GetFinalState()-1]++;  


}


void WZSelectionAnalysis::Finish() {


  for (int i=0; i<4; i++) {
    std::cout << "Yield for final state " << i << "\t:\t" << yieldsByChannel[i];

      }

}
  

