#include "WZSelectionYields.h"

#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;


WZSelectionYields::WZSelectionYields(WZEvent* e, TFile* outputFile)
 : GenericAnalysis(e, outputFile)
{

}


void WZSelectionYields::Init()
{
  nAnalyzedEvents = 0;
	nSelectedEvents = 0;

  for (int i = 1; i <= 4; i++) {
    yieldsByChannelPreselection[i-1] = 0;
    yieldsByChannelZSelection[i-1] = 0;
    yieldsByChannelWSelection[i-1] = 0;
    yieldsByChannelFullSelection[i-1] = 0;
  }

  // Setup selected event lists 
  for (int i = 1; i <= 4; i++) {
    ostringstream outputFileName;
    outputFileName << "output/yields/selectedEvents_" << i << "_test.txt";
    cout << "File name : " << outputFileName.str() << endl;
    eventLists[i-1].open(outputFileName.str().c_str());
  }
}


void WZSelectionYields::EventAnalysis()
{
  nAnalyzedEvents++;

  if (fWZEvent->PassesFullSelection())  yieldsByChannelFullSelection[fWZEvent->GetFinalState()-1]++;
  if (fWZEvent->PassesWSelection())     yieldsByChannelWSelection[fWZEvent->GetFinalState()-1]++;
  if (fWZEvent->PassesZSelection())     yieldsByChannelZSelection[fWZEvent->GetFinalState()-1]++;
  if (fWZEvent->PassesPreselection())   yieldsByChannelPreselection[fWZEvent->GetFinalState()-1]++;

  if (!(fWZEvent->PassesFullSelection))  return;
  nSelected++;

  fWZEvent->DumpEvent(eventLists[fWZEvent->GetFinalState()-1], 1);
}


void WZSelectionYields::Finish()
{
  cout << "CHANNEL \t|\t Preselection \t|\t Z Selection \t|\t W Selection \t|\t Full Selection" << "\n";
  for (int i = 1; i <= 4; i++) {
    cout << i << "\t|\t" << yieldsByChannelPreselection[i-1]
              << "\t|\t" << yieldsByChannelZSelection[i-1]
              << "\t|\t" << yieldsByChannelWSelection[i-1]
              << "\t|\t" << yieldsByChannelFullSelection[i-1] << "\n";
  }

  cout << endl;
}
