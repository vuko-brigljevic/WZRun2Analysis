#ifndef Leptons_h
#define Leptons_h

#include "TLorentzVector.h"


class WZEvent;


class Lepton : public TLorentzVector
{

public:

  Lepton(unsigned int ind, double pt, double eta, double phi, double ch);
  Lepton(unsigned int ind, double pt, double eta, double phi, double ch,
         double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns);
  int GetPdgId() { return fPdgId; }
  double GetCharge() { return fCharge; }
  int GetIndex() { return fIndex; }
  double GetRelIsoDeltaB() { return fRelIsoDeltaB; }
  double GetRelIsoEffArea25ns() { return fRelIsoEffArea25ns; }
  double GetRelIsoEffArea50ns() { return fRelIsoEffArea50ns; }

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

  double fRelIsoDeltaB;
  double fRelIsoEffArea25ns;
  double fRelIsoEffArea50ns;

  static WZEvent* fWZTree;

};


class Electron : public Lepton
{

public:

  Electron(unsigned int index, double pt, double eta, double phi, double charge);
  Electron(unsigned int ind, double pt, double eta, double phi, double ch,
           double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns);
  std::pair<bool, bool> IsLooseTight();
  std::pair<bool, bool> IsLooseTightCutBased25ns();
  std::pair<bool, bool> IsLooseTightCutBased50ns();
  bool PassesPtMinCut();
  bool PassesEtaMaxCut();

};


class Muon : public Lepton
{

public:

  Muon(unsigned int index, double pt, double eta, double phi, double charge);
  Muon(unsigned int ind, double pt, double eta, double phi, double ch,
       double relIsoDeltaB, double relIsoEffA25ns, double relIsoEffA50ns);
  std::pair<bool, bool> IsLooseTight();
  std::pair<bool, bool> IsLooseTightCutBased25ns();
  std::pair<bool, bool> IsLooseTightCutBased50ns();
  bool PassesPtMinCut();
  bool PassesEtaMaxCut();

};


inline
bool HigherPt(const Lepton* l1, const Lepton* l2)
{
  return l1->Pt() > l2->Pt();
}


#endif
