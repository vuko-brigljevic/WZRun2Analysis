#include "WZEvent.h"

//#include "JetEnergyTool.h"
//#include "MetSystematicsTool.h"
//#include "SystematicsManager.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <math.h>



WZEvent::WZEvent(TTree * tree) :
  WZBASECLASS(tree)
  //
{

}

void WZEvent::ReadEvent()
{ // If you want to fill some own variables, do it here...

  final_state     = undefined;
  selection_level = selectionNotRun;
  //  selection_level = failsSelection;


}



bool WZEvent::passesSelection(){

  bool passed = false;

  // Do we have a Z decay ? 

  // Do we have a W candidate

  // MET cut

  if ( (nEle + nMu) == 3)
    passed = true;

  
  return passed;
  
}


