#ifndef Leptons_h
#define Leptons_h



#include "TLorentzVector.h"

class WZEvent;


class Lepton : public TLorentzVector {

public:
  Lepton(int ind, double pt, double eta, double phi,
	 float ch );
  float PdgId() { return pdgid;};
  virtual bool IsTight()=0;
  int Index() { return index;}
  static void SetWZEvent(WZEvent * wzt);


protected:
  float pdgid;
  float charge;
  int index;

  static WZEvent * wztree;

};


class Electron : public Lepton {

public:
  Electron(int index, double pt, double eta, double phi,
	   float ch );
  bool IsTight();

protected:

  //  float a[1000];
  

};

class Muon : public Lepton {

public:
  Muon(int index, double pt, double eta, double phi,
	   float ch );
  bool IsTight();

protected:

  //  float amu[1000];

};



#endif // #ifdef 
