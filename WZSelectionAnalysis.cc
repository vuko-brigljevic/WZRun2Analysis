#include "WZSelectionAnalysis.h"

#include <iostream>
#include <fstream>
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
    yieldsByChannelPreselection[i] = 0;
    yieldsByChannelZSelection[i] = 0;
    yieldsByChannelWSelection[i] = 0;
    yieldsByChannelFullSelection[i] = 0;
  }

  // Setup selected event lists 
  for (int i=0; i<4; i++) {
    ostringstream fileName;
    fileName << "selectedEvents_" << i+1;
    cout << "file name : " << fileName.str() << endl;
    eventLists[i].open(fileName.str().c_str());
  }
}


void WZSelectionAnalysis::EventAnalysis()
{
  nrAnalyzedEvents++;

  if (fWZEvent->PassesPreselection()) {
    yieldsByChannelPreselection[fWZEvent->GetFinalState()-1]++;
  }

  if (fWZEvent->PassesZSelection()) {
    yieldsByChannelZSelection[fWZEvent->GetFinalState()-1]++;
  }
 
    if (fWZEvent->PassesWSelection()) {
    yieldsByChannelWSelection[fWZEvent->GetFinalState()-1]++;
  }

    if (fWZEvent->PassesFullSelection()) {
    yieldsByChannelFullSelection[fWZEvent->GetFinalState()-1]++;
  }

  if (! fWZEvent->PassesFullSelection() ) return;

  yieldsByChannel[fWZEvent->GetFinalState()-1]++;  

  fWZEvent->DumpEvent(eventLists[fWZEvent->GetFinalState()-1],1);


  // MC Analysis
/*
  cout << "MC Tree \n";
  out << "===============\n";
  for (int igen = 0; igen<wzevt->nMC; igen++) {
    //    char  pyName[20];
    //    TPythia::Pyname((wzevt->mcPID)->at(igen),pyName);

    cout << "Gen Particle: " << (wzevt->mcPID)->at(igen)
         << "\t status   : " << (wzevt->mcStatus)->at(igen)
      //	      << "  " << pyName
	      << "\t mother ID: " << (wzevt->mcMomPID)->at(igen)
	      << "\t GrandMa  ID: " << (wzevt->mcGMomPID)->at(igen)
	      << endl;
  }
*/
}


void WZSelectionAnalysis::Finish()
{
  cout << "Total number of analyzed events : " <<   nrAnalyzedEvents << endl << endl;

  for (int i=0; i<4; i++)
    cout << "Yield for Preselection " << i << "\t:\t" << yieldsByChannelPreselection[i] << "\n";
  cout << endl;
  for (int i=0; i<4; i++)
    cout << "Yield for Z Selection " << i << "\t:\t" << yieldsByChannelZSelection[i] << "\n";
  cout << endl;
  for (int i=0; i<4; i++)
    cout << "Yield for W Selection " << i << "\t:\t" << yieldsByChannelWSelection[i] << "\n";
  cout << endl;
  for (int i=0; i<4; i++)
    cout << "Yield for Complete Selection " << i << "\t:\t" << yieldsByChannelFullSelection[i] << "\n";
  cout << endl;

  for (int i=0; i<4; i++)
    cout << "Yield for final state " << i << "\t:\t" << yieldsByChannel[i] << endl;
}

