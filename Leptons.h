#ifndef Leptons_h
#define Leptons_h

#include "TLorentzVector.h"


class WZEvent;


class Lepton : public TLorentzVector
{

public:

  Lepton(unsigned int ind, double pt, double eta, double phi, double ch);
  int GetPdgId() { return fPdgId; }
  double GetCharge() { return fCharge; }
  int GetIndex() { return fIndex; };

  virtual std::pair<bool, bool> IsLooseTight() = 0;
  virtual std::pair<bool, bool> IsLooseTightCutBased25ns() = 0;
  virtual std::pair<bool, bool> IsLooseTightCutBased50ns() = 0;
  virtual bool PassesPtMinCut() = 0;
  virtual bool PassesEtaMaxCut() = 0;

  static void SetWZEvent(WZEvent* wzt);


protected:

  int fPdgId;
  double fCharge;
  int fIndex;

  static WZEvent* fWZTree;

};


class Electron : public Lepton
{

public:

  Electron(unsigned int index, double pt, double eta, double phi, double charge);
  std::pair<bool, bool> IsLooseTight();
  std::pair<bool, bool> IsLooseTightCutBased25ns();
  std::pair<bool, bool> IsLooseTightCutBased50ns();
  bool PassesPtMinCut();
  bool PassesEtaMaxCut();


protected:

  double EffA25ns();
  double EffA50ns();

};


class Muon : public Lepton
{

public:

  Muon(unsigned int index, double pt, double eta, double phi, double charge);
  std::pair<bool, bool> IsLooseTight();
  std::pair<bool, bool> IsLooseTightCutBased25ns();
  std::pair<bool, bool> IsLooseTightCutBased50ns();
  bool PassesPtMinCut();
  bool PassesEtaMaxCut();

};

#endif
