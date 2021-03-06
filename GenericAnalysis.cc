#include "GenericAnalysis.h"

#include <ios>


GenericAnalysis::GenericAnalysis(WZEvent* wzevt, TFile* fout)
{
  fWZEvent = wzevt;
  outputRootFile = fout;
}


TH1D* GenericAnalysis::bookTH1D(TString key, TString title, unsigned int nbins, double min, double max)
{
  TH1D* h = new TH1D(key, title, nbins, min, max);
  h->SetFillColor(kCyan);
  fListOfHistos.push_back(h);
  return h;
}


TH2D* GenericAnalysis::bookTH2D(TString key, TString title,
                                 unsigned int nbinsx, double xmin, double xmax,
                                 unsigned int nbinsy, double ymin, double ymax)
{
  TH2D* h = new TH2D(key, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax);
  h->SetOption("COLZ");
  fListOfHistos.push_back(h);
  return h;
}


void GenericAnalysis::Init()
{

}


void GenericAnalysis::EventAnalysis()
{

}


GenericAnalysis::~GenericAnalysis()
{
  std::cout << "Closing file: writing to ROOT.... \n";

  if (outputRootFile) {
    outputRootFile->cd();

    std::cout << "# objects to write : " << fListOfHistos.size() << std::endl;
    for (unsigned int i = 0; i < fListOfHistos.size(); i++) {
      if (typeid(*(fListOfHistos.at(i))) == typeid(TH2D))  fListOfHistos.at(i)->Draw();
      else if (typeid(*(fListOfHistos.at(i))) == typeid(TH1D))  fListOfHistos.at(i)->Draw();

      std::cout << "Writing object : " << fListOfHistos.at(i)->GetName() << std::endl;
      fListOfHistos.at(i)->Write();
    }
  }
}
