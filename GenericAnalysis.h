#ifndef GenericAnalysis_h
#define GenericAnalysis_h

#include "WZEvent.h"
#include <fstream>

#include <vector>


class GenericAnalysis {

public:

  GenericAnalysis(WZEvent * e, TFile * fout=0);

  ~GenericAnalysis();

  virtual void Init();

  virtual void EventAnalysis();

  //  virtual void Finish(TFile * fout =0);

  TH1F * bookTH1F(TString key, TString title, int nbins, float min, float max);
  TH2F * bookTH2F(TString key, TString title, 
		  int nbinsx , float xmin, float xmax,
		  int nbinsy , float ymin, float ymax);



protected:

  WZEvent * wzevt;
  
  std::vector<TObject * > listOfHistos;
  TFile * outputRootFile;


};

#endif
