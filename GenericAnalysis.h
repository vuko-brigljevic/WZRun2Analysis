#ifndef GenericAnalysis_h
#define GenericAnalysis_h

#include "WZEvent.h"
#include "TH1D.h"
#include "TH2D.h"

#include <fstream>


class GenericAnalysis
{

public:

  GenericAnalysis(WZEvent* e, TFile* fout = 0);
  ~GenericAnalysis();

  virtual void Init();
  virtual void EventAnalysis();

  TH1D* bookTH1D(TString key, TString title, unsigned int nbins, double min, double max);
  TH2D* bookTH2D(TString key, TString title,
                 unsigned int nbinsx, double xmin, double xmax,
                 unsigned int nbinsy, double ymin, double ymax);


protected:

  WZEvent* fWZEvent;
  
  std::vector<TObject*> fListOfHistos;
  TFile* outputRootFile;

};

#endif