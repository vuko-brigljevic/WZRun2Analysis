#include "GenericAnalysis.h"

#include <ios>
#include <iostream>
#include <fstream>


GenericAnalysis::GenericAnalysis(WZEvent * e, TFile * fout) {

  wzevt = e;
  outputRootFile = fout;

}

TH1F * GenericAnalysis::bookTH1F(TString key, 
				 TString title, 
				 int nbins, float min, float max) {

  TH1F * h = new TH1F(key, title, nbins, min, max);
  listOfHistos.push_back(h);
  return h;


}


TH2F * GenericAnalysis::bookTH2F(TString key, TString title, 
		int nbinsx , float xmin, float xmax,
		int nbinsy , float ymin, float ymax) {


  TH2F * h = new TH2F(key, title, 
		      nbinsx , xmin, xmax,
		      nbinsy , ymin, ymax);
  
  listOfHistos.push_back(h);
  return h;


}

void GenericAnalysis::Init() {


}

void GenericAnalysis::EventAnalysis() {


}

GenericAnalysis::~GenericAnalysis() {

  if (outputRootFile) {
    outputRootFile->cd();
    for (int i=0; i<listOfHistos.size(); i++ ) {
      listOfHistos[i]->Write();
    }
  }
}
