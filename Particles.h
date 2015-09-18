#ifndef Particles_h
#define Particles_h

#include "TLorentzVector.h"






class WZEvent;


class Particle : public TLorentzVector
{

public:

  Particle(unsigned int ind, double pt, double eta, double phi);

  int GetIndex() { return fIndex; }

  static void SetWZEvent(WZEvent* wzt);


  


protected:

  int fIndex;

  static WZEvent* fWZTree;

};

class GenParticle : public Particle
{

public:

  GenParticle (unsigned int index, double pt, double eta, double phi);

  int PdgId() { return fPdgId;}

protected:

  int fPdgId;

};


class GenDressedLepton : public TLorentzVector
{

public:

  GenDressedLepton(GenParticle * originalLepton, TLorentzVector p4);
  GenDressedLepton(GenParticle * originalLepton, WZEvent * evt);

  GenParticle * GetOriginalLepton();

protected:

  GenParticle * fUndressedLepton;

};



#endif
