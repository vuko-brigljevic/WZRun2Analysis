#include "WZEvent.h"

#include "JetEnergyTool.h"
#include "MetSystematicsTool.h"
#include "SystematicsManager.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include <iostream>
#include <math.h>


//
//  Exteran function declarations
//

bool Z_muons(WZBASECLASS *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, TLorentzVector* analysisLepton, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, TLorentzVector * analysisLepton);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

float trigger3sameLeptons(float* eL, float* eT);

float trigger2sameLeptons(float* eL, float* eT);

float triggerDifferentLeptons(float* eL, float* eT);

double ReturnBranchingWeight(int type);

TLorentzVector GetMET(float metModule, float metPhi);

TH1F* GetHistogramFromGraph(TString hname, TString gname);
// Initialize static data members
TH2F * RecoLepton::MuonSF = 0;
TH2F * RecoLepton::ElecSF = 0;
TH2F * RecoLepton::DoubleMuLeadEff    = 0;
TH2F * RecoLepton::DoubleMuTrailEff   = 0;
TH2F * RecoLepton::DoubleEleLeadEff   = 0;
TH2F * RecoLepton::DoubleEleTrailEff  = 0;
TH1F* WZEvent::hScaleInEB  =0;
TH1F* WZEvent::hScaleOutEB =0;
TH1F* WZEvent::hScaleEE    =0;


float RecoLepton::LeadTriggerEff() {

  TH2F * hlead;

  float eff = 1.;

  if (fabs(pdgid)==11) {
    if (DoubleEleLeadEff == 0) {
        DoubleEleLeadEff=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElLead", "DoubleElLead");
    }
    hlead = DoubleEleLeadEff;
  } else if (fabs(pdgid)==13) {
    if (DoubleMuLeadEff == 0) {
        DoubleMuLeadEff =LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuLead", "DoubleMuLead");
    }
    hlead = DoubleMuLeadEff;
  } else {
    std::cout << "RecoLepton::LeadTriggerEff> UNKNOWN LEPTON ID : " 
	      << pdgid << std::endl;
    return -1;
  }
  if (hlead) {
    //    int bin = hlead->FindBin(Pt(),Eta());
    //    eff = hlead->GetBinContent(bin);
    float factor= GetFactor(hlead, Pt(), Eta());
    //     std::cout << "eff = " << eff << "\t" << factor << std::endl;
    eff = factor;
  } else {
    std::cout << "Missing eff histo for " << pdgid << std::endl;
  }

  return eff;

}

float RecoLepton::TrailTriggerEff() {

  TH2F * htrail;

  if (fabs(pdgid)==11) {
    if (DoubleEleTrailEff == 0) {
        DoubleEleTrailEff=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElTrail", "DoubleElTrail");
    }
    htrail = DoubleEleTrailEff;
  } else if (fabs(pdgid)==13) {
    if (DoubleMuTrailEff == 0) {
      DoubleMuTrailEff =LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuTrail", "DoubleMuTrail");
    }
    htrail = DoubleMuTrailEff;
  }
  //  int bin = htrail->FindBin(Pt(),Eta());
  //  float eff = htrail->GetBinContent(bin);

  float eff= GetFactor(htrail, Pt(), Eta());

  return eff;
}

float WZEvent::GetBrWeight(){
  int genType=-999;
  
  double Wel(false), Wmu(false), Wtau(false), Zel(false), Zmu(false), Ztau(false);
  int numW(0), numZ(0), indexW(-999), indexZ1(-999), indexZ2(-999);

  for (int igl=0; igl<genLeptons.size() ; igl++) {
    int bosonId = genLeptons[igl].MotherBoson();
  //W
    if (abs(bosonId) == 24) {
      numW++;
      indexW=igl;
    }
    //Z
    if (abs(bosonId) == 23) {
      numZ++;
      if (numZ==2) indexZ2=igl;
      if (numZ==1) indexZ1=igl;
    }
  }
  if (numW > 1) return ReturnBranchingWeight(genType);
  if (numZ > 2) return ReturnBranchingWeight(genType);
  /*
  int Zid1=abs(cWZ->genLeptons[indexZ1].Id());
  int Zid2=abs(cWZ->genLeptons[indexZ2].Id());
  int wid=abs(cWZ->genLeptons[indexW].Id());
  */

  //hadronic decays
  if (numW==0) Wtau=true;
  if (numZ<2) Ztau=true;

  //W lepton
  if (numW>0){
    if (((genLeptons[indexW].ComesFromTau()))) Wtau=true;
    else {
      if ((abs(genLeptons[indexW].Id()))==11) Wel=true;
      if ((abs(genLeptons[indexW].Id()))==13) Wmu=true;
    }
  }
  //Z lepton
  if (numZ>1){
    if ((genLeptons[indexZ1].ComesFromTau()) && (genLeptons[indexZ2].ComesFromTau())) Ztau=true;
    else {
      if (((abs(genLeptons[indexZ1].Id()))==11) && ((abs(genLeptons[indexZ2].Id()))==11)) Zel=true;
      if (((abs(genLeptons[indexZ1].Id()))==13) && ((abs(genLeptons[indexZ2].Id()))==13)) Zmu=true;
    }
  }
  if (Zel && Wel) genType=0;
  if (Zel && Wmu) genType=1;
  if (Zmu && Wel) genType=2;
  if (Zmu && Wmu) genType=3;
  if (Ztau && Wel) genType=4;
  if (Ztau && Wmu) genType=5;
  if (Ztau && Wtau) genType=6;
  if (Zel && Wtau) genType=7;
  if (Zmu && Wtau) genType=8;
  

  float weight;
  weight=ReturnBranchingWeight(genType);
  return weight;
}



float RecoLepton::GetScaleFactor() {

  if (MuonSF) { 
    //std::cout << "eff map defined \n";
  }  else {
    //std::cout << "defining eff map \n";
    MuonSF=LoadHistogram("auxiliaryFiles/MuSF_2012.root", "h2inverted", "MuonSF");
  }
  if (ElecSF) { 
    // >>>>>>> b9fc1514f52dbc4f21337b6468cc3b7ba153e7c7
    //    std::cout << "eff map defined \n";
  }  else {
    //std::cout << "defining eff map \n";
    ElecSF=LoadHistogram("auxiliaryFiles/EleSF_2012.root", "h2inverted", "ElecSF");
  }

  float leptonEta=Eta();
  float leptonPt= Pt();
  float factor(0);
  if (fabs(pdgid)==11)
    factor= GetFactor(ElecSF, leptonPt, leptonEta);
  else if (fabs(pdgid)==13)
    factor= GetFactor(MuonSF, leptonPt, leptonEta);
  // >>>>>>> b9fc1514f52dbc4f21337b6468cc3b7ba153e7c7
  return factor;

}

WZEvent::WZEvent(TTree * tree) :
  WZBASECLASS(tree),
  //
  numZ(0),
  numW(0),
  numMET(0),
  num3e(0),
  num2e1mu(0),
  num1e2mu(0),
  num3mu(0),
  numMET3e(0),
  numMET2e1mu(0),
  numMET1e2mu(0),
  numMET3mu(0)
{
  pt.push_back(&pt1);
  pt.push_back(&pt2);
  pt.push_back(&pt3);
  pt.push_back(&pt4);
  //
  bdt.push_back(&bdt1);
  bdt.push_back(&bdt2);
  bdt.push_back(&bdt3);
  bdt.push_back(&bdt4);
  // 
  ch.push_back(&ch1);
  ch.push_back(&ch2);
  ch.push_back(&ch3);
  ch.push_back(&ch4);
  //
  eta.push_back(&eta1);
  eta.push_back(&eta2);
  eta.push_back(&eta3);
  eta.push_back(&eta4);
  // 
  phi.push_back(&phi1);
  phi.push_back(&phi2);
  phi.push_back(&phi3);
  phi.push_back(&phi4);
  // 
  pdgid.push_back(&pdgid1);
  pdgid.push_back(&pdgid2);
  pdgid.push_back(&pdgid3);
  pdgid.push_back(&pdgid4);
  // 
  pass2012ICHEP.push_back(&pass2012ICHEP1);
  pass2012ICHEP.push_back(&pass2012ICHEP2);
  pass2012ICHEP.push_back(&pass2012ICHEP3);
  pass2012ICHEP.push_back(&pass2012ICHEP4);
  //
  iso.push_back(&iso1);
  iso.push_back(&iso2);
  iso.push_back(&iso3);
  iso.push_back(&iso4);
  // 
  isomva.push_back(&isomva1);
  isomva.push_back(&isomva2);
  isomva.push_back(&isomva3);
  isomva.push_back(&isomva4);
}

void WZEvent::ReadEvent()
{

  final_state     = undefined;
  selection_level = selectionNotRun;
  //  selection_level = failsSelection;

  wLeptonIndex     = -1;
  zLeptonsIndex[0] = -1;
  zLeptonsIndex[1] = -1;


  // Red Reco leptons

  leptons.clear();

  float lepton_pt[4]  = {pt1,pt2,pt3,pt4};
  float lepton_eta[4] = {eta1,eta2,eta3,eta4};
  float lepton_phi[4] = {phi1,phi2,phi3,phi4};
  float lepton_ch[4]  = {ch1,ch2,ch3,ch4};
  float lepton_id[4]  = {pdgid1,pdgid2,pdgid3,pdgid4};

  for (int il=0; il<4; il++) {
    if (lepton_pt[il]>1.) {
      leptons.push_back(RecoLepton(lepton_pt[il],
				   lepton_eta[il],
				   lepton_phi[il],
				   lepton_ch[il],
				   lepton_id[il]));
    }

  }

 
  // Read generated leptons

  genLeptons.clear();
  TLorentzVector gl1,gl2,gl3; // ,gl;
  //  gl1.SetPtEtaPhiM(leptonGenpt1,leptonGeneta1,leptonGenphi1,0.);
  //  gl2.SetPtEtaPhiM(leptonGenpt2,leptonGeneta2,leptonGenphi2,0.);
  //  gl3.SetPtEtaPhiM(leptonGenpt3,leptonGeneta3,leptonGenphi3,0.);

  float genLeptons_phi[3] = {genVV_S1lepton1_phi,genVV_S1lepton2_phi,genVV_S1lepton3_phi};
  float genLeptons_eta[3] = {genVV_S1lepton1_eta,genVV_S1lepton2_eta,genVV_S1lepton3_eta};
  float genLeptons_pt[3]  = {genVV_S1lepton1_pt,genVV_S1lepton2_pt,genVV_S1lepton3_pt};
  float genLeptons_pid[3]  = {genVV_S1lepton1_pid,genVV_S1lepton2_pid,genVV_S1lepton3_pid};
  float genLeptons_oVpid[3]  = {genVV_S1lepton1_oVpid,genVV_S1lepton2_oVpid,genVV_S1lepton3_oVpid};
  float genLeptons_imTau[3]  = {genVV_S1lepton1_imTau,genVV_S1lepton2_imTau,genVV_S1lepton3_imTau};

  //  gl1.SetPtEtaPhiM(genVV_S1lepton1_pt,genVV_S1lepton1_eta,genVV_S1lepton1_phi,0);
  //  gl2.SetPtEtaPhiM(genVV_S1lepton2_pt,genVV_S1lepton2_eta,genVV_S1lepton2_phi,0);
  //  gl3.SetPtEtaPhiM(genVV_S1lepton3_pt,genVV_S1lepton3_eta,genVV_S1lepton3_phi,0);
  //  genLeptons.push_back(gl1);
  //  genLeptons.push_back(gl2);
  //  genLeptons.push_back(gl3);

  for (int i=0; i<3; i++) {
    if (genLeptons_pt[i] > -1.) {
      GenS1Lepton gl(genLeptons_pt[i],genLeptons_eta[i],genLeptons_phi[i],
		     genLeptons_pid[i],genLeptons_oVpid[i],genLeptons_imTau[i]);

      //      gl.SetPtEtaPhiM(genLeptons_pt[i],genLeptons_phi[i],genLeptons_eta[i],0);
      genLeptons.push_back(gl);
    } 
  }


  for (int i=0; i<genLeptons.size(); i++) {
    float ptl = genLeptons[i].Pt();
    //    if (ptl>900.) {
      //      std::cout << "Gen lepton " << i << " \t Pt = " << ptl
      //	<< "\t channel: " << WZchan << std::endl;
    //}
  }

  float genJets_phi[5] = {genVV_jet1_phi, genVV_jet2_phi, genVV_jet3_phi, genVV_jet4_phi, genVV_jet5_phi};
  float genJets_eta[5] = {genVV_jet1_eta, genVV_jet2_eta, genVV_jet3_eta, genVV_jet4_eta, genVV_jet5_eta};
  float genJets_pt[5]  = {genVV_jet1_pt, genVV_jet2_pt, genVV_jet3_pt, genVV_jet4_pt, genVV_jet5_pt};


  // Read generated jets
  TLorentzVector gj;  
  genJets.clear();
  for (int i=0; i<3; i++) {
    if (genJets_pt[i]>0.) {
      gj.SetPtEtaPhiM(genJets_pt[i],genJets_eta[i],genJets_phi[i],0);
      genJets.push_back(gj);
    }
  }

  
  // Read recoJets  
  recoJets.clear();

  float recoJets_pt[6]  = {jetpt1,jetpt2,jetpt3,jetpt4,jetpt5,jetpt6};
  float recoJets_eta[6] = {jeteta1,jeteta2,jeteta3,jeteta4,jeteta5,jeteta6};
  float recoJets_phi[6] = {jetphi1,jetphi2,jetphi3,jetphi4,jetphi5,jetphi6};
  float recoJets_m[6]   = {jetmass1,jetmass2,jetmass3,jetmass4,jetmass5,jetmass6};

  TLorentzVector rj;  

  for (int i=0; i<6; i++) {
    if (recoJets_pt[i]>0.) {
      rj.SetPtEtaPhiM(recoJets_pt[i],recoJets_eta[i],recoJets_phi[i],recoJets_m[i]);
      recoJets.push_back(rj);
    }
    // else {
    //      std::cout << "RECO JET WITH NEGATIVE PT: " << i << "\t"
    //		<< recoJets_pt[i] << " PHI = " 
    //		<< recoJets_phi[i] << std::endl;
    //    }
  }
  //  for (int i=0; i<genJets.size(); i++) {
  //    std::cout << "ZERO GEN JET: " << i << " " ;
  //    genJets[i].Print();
  //  }


  SystematicsManager * sysManager = SystematicsManager::GetInstance();

  float jes_strength = sysManager->GetValue("JES");

  SmearJets();


  if (jes_strength!=0.) 
    {
      ApplyJESCorrection(jes_strength);
    }

  // MET Systematics

  int met_syst = sysManager->GetValue("MET");

  if (met_syst != 0) {

    float met_new = GetShiftedMET(met_syst);

//    std::cout << "Old MET = " << pfmetTypeI
//	      << "\t new = " << met_new << std::endl;


    pfmetTypeI = met_new;

  }


}


void WZEvent::DumpEvent(std::ostream & out, int verbosity) {

  out << run << "\t" 
      << event ;

  float totalWeight = GetPileupWeight()*GetBrWeight()*GetMCWeight();

  if (verbosity>0) {
    out << "  " << totalWeight
	<< "  " << GetPileupWeight()
	<< "  " << GetBrWeight()
	<< "  " << GetMCWeight()
	<< "  " << GetTriggerEfficiency();

  }


  int leptonIndices[3] = { zLeptonsIndex[0],
			   zLeptonsIndex[1],
			   wLeptonIndex };
  //
  if (verbosity>1) {
    for (int i=0; i<3; i++) {
      if (leptonIndices[i]>=0 
	  && leptonIndices[i]<leptons.size()) {
	int index = leptonIndices[i];
	out << " Lepton : " << i
	    << " Pt = " << leptons[index].Pt()
	    << " Eta = " << leptons[index].Eta()
	    << " SF = " << leptons[index].GetScaleFactor()
	    << " Effl = " << leptons[index].LeadTriggerEff()
	    << " effT = " << leptons[index].TrailTriggerEff();
      }
    }
  }

  out << std::endl;

}

void WZEvent::PrintSummary()
{

  std::cout<<" SUMMARY " << std::endl << std::endl;
  std::cout<<" ==================== \n \n";

}

bool WZEvent::passesSelection(){

  // Check if we've run already on this event
  // If so, no need to rerun the full selection  
  if (selection_level != selectionNotRun) {
    if (selection_level == passesFullSelection) {
      return true;
    } else 
      return false;
  }

  selection_level = failsSelection;

  bool passed = false;


  selectedZPt=-88888.;
  selectedZP4.SetPtEtaPhiM(-9999.,-9999.,0.,-9999.);
  //rejecting run 201191

  if (run==201191) return false;
  if (!trigger)  return false;

  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);

  float pileUpWeight=this->GetPileupWeight();
    
  //rejecting run 201191
  if (this->run==201191) return passed;
  
  //  if (!(this->trigger)) return passed;
      float pts[leptonNumber]={pt1, pt2, pt3, pt4};
  float charges[leptonNumber]={ch1, ch2, ch3, ch4};
  float phis[leptonNumber]={phi1, phi2, phi3, phi4};
  float etas[leptonNumber]={eta1, eta2, eta3, eta4};
  TLorentzVector analysisLepton[leptonNumber];
  TLorentzVector analysisLeptonOld[leptonNumber];
  TLorentzVector EventMET;
  float pfmet= this->pfmetTypeI;
  float pfmetphi=this->pfmetTypeIphi;
  EventMET= GetMET(pfmet, pfmetphi); 
  /*
  for (int i1=0; i1<leptonNumber; i1++){
     if ((fabs(*pdgid[i1])==11)&& (*pt[i1]>0))
       analysisLepton[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], electronMass);
     if ((fabs(*pdgid[i1])==13) && (*pt[i1]>0))
       analysisLepton[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], muonMass);
  }
  */
  //some things for electron scale systematics
  
  
  if (hScaleInEB==0)
    hScaleInEB  = GetHistogramFromGraph("hScaleInEB",  "gScaleInEB");
  if (hScaleOutEB==0)
    hScaleOutEB = GetHistogramFromGraph("hScaleOutEB", "gScaleOutEB");
  if (hScaleEE==0)
    hScaleEE    = GetHistogramFromGraph("hScaleEE",    "gScaleEE");
  
  SystematicsManager * sysManager = SystematicsManager::GetInstance();

  int scaleSyst_mu= sysManager->GetValue("mu_scale_syst");
  int scaleSyst_el= sysManager->GetValue("ele_scale_syst");

  /*
  bool muScaleSyst(false);
  bool elScaleSyst(false);
  
  double elScale(1.0);
  */

  double muScale(0.002);  
  for (int i1=0; i1<leptonNumber; i1++){
    if ((fabs(*pdgid[i1])==13)&& (*pt[i1]>0)){
      analysisLepton[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], muonMass);
      analysisLeptonOld[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], muonMass);
      if (abs(scaleSyst_mu)){
	double spt=(*pt[i1])+ (*pt[i1])*muScale*scaleSyst_mu;
	double factScale=(*pt[i1])/spt;
	analysisLepton[i1]*=factScale;
	EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
      }
    }
  

    if ((fabs(*pdgid[i1])==11) && (*pt[i1])>0){
      analysisLepton[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], electronMass);
      analysisLeptonOld[i1].SetPtEtaPhiM(*pt[i1], *eta[i1], *phi[i1], electronMass);
      if (abs(scaleSyst_el)){
	double scale;
	const Float_t InEBMax  = hScaleInEB ->GetXaxis()->GetBinCenter(hScaleInEB ->GetNbinsX());
	const Float_t OutEBMax = hScaleOutEB->GetXaxis()->GetBinCenter(hScaleOutEB->GetNbinsX());
	const Float_t EEMax    = hScaleEE   ->GetXaxis()->GetBinCenter(hScaleEE   ->GetNbinsX());
	const Float_t scaleInEB  = hScaleInEB ->GetBinContent(hScaleInEB ->FindBin(std::min((*pt[i1]), InEBMax)));
	const Float_t scaleOutEB = hScaleOutEB->GetBinContent(hScaleOutEB->FindBin(std::min((*pt[i1]), OutEBMax)));
	const Float_t scaleEE    = hScaleEE   ->GetBinContent(hScaleEE   ->FindBin(std::min((*pt[i1]), EEMax)));


	const Float_t aeta = fabs(*eta[i1]);

	if (aeta < 0.8)
	  {
	    scale = scaleInEB;
	  }
	else if (aeta >= 0.8 && aeta < 1.479)
	  {
	    scale = scaleOutEB;
	  }
	else
	  {
	    scale = scaleEE;
	  }
	double spt=*pt[i1]+ (*pt[i1])*scale*scaleSyst_el;
	double factScale=(*pt[i1])/spt;
	analysisLepton[i1]*=factScale;
      }
      EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
    }
    
  }
  
  //find Z boson, save the index of it
  std::vector<int> good_muons;
  std::vector<int> good_electrons;
  TLorentzVector v_nizEl[9];
  TLorentzVector v_nizMu[9];
  TLorentzVector v_3Lepton(0.,0.,0.,0.);
  
  int WZcandidates[3]; 
    
  //here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
  int lepNum(0);
  for (int i=0; i<leptonNumber; i++){
    if ((*pt[i]>10) && (*pt[i]!=-9999))
      lepNum++;
    //      if ((*bdt[i]>-10) && (*pass2012ICHEP[i]) && (*pt[i]>10)){
    if (( *pdgid[i]==-11 || *pdgid[i]==11)   && (*pt[i]>10)){
      good_electrons.push_back(i);
      v_nizEl[i].SetPtEtaPhiM(*pt[i],*eta[i], *phi[i], electronMass);
      v_3Lepton=v_3Lepton+v_nizEl[i];
    }
    //      if ((*bdt[i]<-10) && (*pass2012ICHEP[i]) && (*pt[i]>10)){
    if (( *pdgid[i]==-13 || *pdgid[i]==13)   && (*pt[i]>10)){
      good_muons.push_back(i);
      v_nizMu[i].SetPtEtaPhiM(*pt[i],*eta[i], *phi[i], muonMass);
      v_3Lepton=v_3Lepton+v_nizMu[i];
    }
  }
  
  if (lepNum!=3) return false;
  
  if (v_3Lepton.M()<100) return false;
  
  bool foundZel(false), foundZmu(false);
  double massMu(-999), massEl(0), Zpt(0);
  
  foundZmu=  Z_muons(this, &good_muons, WZcandidates, v_nizMu, analysisLepton, charges, massMu, Zpt);
  
  foundZel= Z_muons(this, &good_electrons, WZcandidates, v_nizEl, analysisLepton, charges, massEl, Zpt);
  
  // reject all events without Z boson--and candidates with double Z boson
  
  ///////////////////////////rejecting Zel&Zmu
  if ((foundZel) && (foundZmu)){
    return false;
  }
  
  /////////////////////////////rejecting two independent Zels
  if (Z_independent(charges, &good_muons, WZcandidates, v_nizMu)) 
    return false;
  
  
  ////////////////////////////rejecting two independent Zmus
  if (Z_independent(charges, &good_electrons, WZcandidates, v_nizEl)) 
    return false;
  
  ///////////////////////////rejecting when found no Zel& no Zmu
  if ((!foundZel) && (!foundZmu))
    return false;
  
  selection_level = passesZSelection;
  
  ///////////////////////FILLING HISTOGRAMS/////////////////////////////////////
  //later...
  
  //jet number...
  
  
  //     std::cout << "Run: " << run << "\t Event: " << event 
  // 	      << " Found Z: " << foundZel << foundZmu 
  // 	      << "\t Good lepton numbers: " 
  // 	      << good_electrons.size() 
  // 	      << "\t " << good_muons.size() << std::endl;
  
  
  //////////////////////FIND W BOSON///////////////////////////////////////////
  int ZmuWelCounter(0),ZmuWmuCounter(0), ZelWelCounter(0), ZelWmuCounter(0); 
  int ZmuWel_index(0), ZmuWmu_index(0), ZelWel_index(0), ZelWmu_index(0);
  bool testEl1(false),testMu1(false), testEl2(false), testMu2(false);
  bool ZmuWel(false), ZmuWmu(false), ZelWel(false),ZelWmu(false);
  bool ev3mu(false), ev1e2mu(false), ev2e1mu(false), ev3e(false);
  
  //*****Z->mumu*************
  if (foundZmu){
    //      std::cout << "Found Zmumu : " << ZmuWel << ZmuWmu << std::endl;
    for (int iel1=0; iel1< (good_electrons.size()); iel1++){
      int elIndex1= good_electrons[iel1];
      if (*pt[elIndex1]>20){
	ZmuWelCounter++;
	ZmuWel_index=elIndex1; 
	testEl1=true;  
      }
    }
    if (ZmuWelCounter==1) ZmuWel=true;
    int ZmuWmuCounter(0);
    
    for (int imu1=0; imu1< (good_muons.size());imu1++){
      int muIndex1=good_muons[imu1];
      if ((muIndex1==WZcandidates[0]) || (muIndex1==WZcandidates[1])) continue;
      
      if (*pt[muIndex1]>20){
	ZmuWmuCounter++;
	ZmuWmu_index=muIndex1;
	testMu1=true;
      }
    }
    if (ZmuWmuCounter==1) ZmuWmu=true;
    
    if ((ZmuWmu) && (ZmuWel))   return false;
    
    if ((!ZmuWmu) && (!ZmuWel)) return false;
    
    if ((ZmuWmu) && (!ZmuWel))
      {
	WZcandidates[2]=ZmuWmu_index;
	ev3mu=true;
      }
    
    if ((ZmuWel) && (!ZmuWmu))
      {
	WZcandidates[2]=ZmuWel_index;
	ev1e2mu=true;
      }
  }
  
  //**********Z->ee***************
  if (foundZel) {
    int ZelWelCouter(0);
    for (int iel2=0; iel2< (good_electrons.size()); iel2++){
      int elIndex2=good_electrons[iel2];
      if ((elIndex2==WZcandidates[0]) || (elIndex2==WZcandidates[1])) continue;
      if (*pt[elIndex2]>20)
	{
	  ZelWelCounter++;
	  ZelWel_index=elIndex2;
	  testEl2=true;
	}
    }
    if (ZelWelCounter==1)
      { 
	ZelWel=true;
      }
    int ZelWmuCounter(0);
    for (int imu2=0; imu2<(good_muons.size()); imu2++){
      int muIndex2=good_muons[imu2];
      if (*pt[muIndex2]>20)
	{
	  ZelWmuCounter++;
	  ZelWmu_index=muIndex2;
	  testMu2=true;
	}
    }
    if (ZelWmuCounter==1) ZelWmu=true;
    
    if ((ZelWel) && (ZelWmu)){
      return false;
    }
    if ((!ZelWel) && (!ZelWmu)) return false;
    
    if ((ZelWmu) && (!ZelWel))
      {
	WZcandidates[2]=ZelWmu_index;
	ev2e1mu=true;
      }
    if ((ZelWel) && (!ZelWmu))
      {
	WZcandidates[2]=ZelWel_index;
	ev3e=true;
      }
  }
  
  //    std::cout << ev3e << ev3mu << ev1e2mu << ev2e1mu << std::endl;    
  
  if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) return false;
  
  //deltaR condition
  if (!passDeltaRWleptZlept(WZcandidates, analysisLepton)) return false;
  numW+=pileUpWeight;  
  
  selection_level = passesWSelection;
  
  if (ev3e){
    final_state = eee;
  }
  if (ev2e1mu){
    final_state = eem;
  }
  if (ev1e2mu){
    final_state = mme;
  }
  if (ev3mu){
    final_state = mmm;
  }
  //////////////////////////////////////MET CUT//////////////////////////
  
  //  if ((this->pfmet)<30)  return false;
  //  if ((this->pfmetTypeI)<30)  return false;
  if (EventMET.Et()<30) return false;
  //    if ((this->pfmetTypeI)<30) continue;   ///CHANGE THIS
  
  selection_level = passesFullSelection;
  
  //    std::cout << "PASSED FINAL SELECTION: " << state << std::endl;





  TLorentzVector zl1,zl2,zp4;;
  

  zl1.SetPtEtaPhiM(*pt[WZcandidates[0]],*eta[WZcandidates[0]],
		   *phi[WZcandidates[0]], muonMass);

  zl2.SetPtEtaPhiM(*pt[WZcandidates[1]],*eta[WZcandidates[1]],
		   *phi[WZcandidates[1]], muonMass);

  zLeptonsIndex[0] =   WZcandidates[0]; 
  zLeptonsIndex[1] =   WZcandidates[1]; 
  wLeptonIndex     =   WZcandidates[2]; 

  zp4 = zl1+zl2;

  selectedZPt = zp4.Pt();

  selectedZP4 = zp4;
  
  return true;
  
}


TLorentzVector * WZEvent::ZLepton(unsigned int i) {

  if (passesSelection()
      && i<=1) {
      return &leptons[zLeptonsIndex[i]];
  } else {
      return 0;
  }
}


TLorentzVector * WZEvent::WLepton() {

  if (passesSelection()) {
    return &leptons[wLeptonIndex];
  } else 
    return 0;
}
  



bool WZEvent::PassesGenCuts(){

  return (MZ>71.1876 && MZ<111.1876);

}


float WZEvent::GetPileupWeight() {
  
  SystematicsManager * sysManager = SystematicsManager::GetInstance();    
  int pu_sys  = sysManager->GetValue("pu_syst"); 
  
#ifdef NEWMCPUFIX
  if (pu_sys==0){
    return puW_new;}
  else if (pu_sys==-1){
    return puW_down;}
  else if (pu_sys==1){
    return puW_up;}
 #else
  return puW;
#endif
}


float  WZEvent::GetMCWeight() {

  // The following makes sense only for accepted reco events
  //   MC weight is zero otherwise

  if (! passesSelection() ) {
    return 0.;
  }
  SystematicsManager * sysManager = SystematicsManager::GetInstance();

  int systEl= sysManager->GetValue("ele_SF");
  int systMu= sysManager->GetValue("mu_SF");
  
  float trigEff = GetTriggerEfficiency();

  // Get lepton scale factors
  double leptonSF = 1.;

  int leptonIndices[3] = { zLeptonsIndex[0],
			   zLeptonsIndex[1],
			   wLeptonIndex };
  //
  //  std::cout<<"SYST: "<<systMu<<std::endl;
  
  for (int i=0; i<3; i++) {
    if (leptonIndices[i]>=0 
	&& leptonIndices[i]<leptons.size()) {
      double lsf = leptons[leptonIndices[i]].GetScaleFactor();
      if (fabs(leptons[i].PdgId())==11){
	lsf += systEl*0.01;
      }
      if (fabs(leptons[i].PdgId())==13){
	lsf += systMu*0.01;
      }

      leptonSF *= lsf;
    } else {
      std::cout << "Lepton Index out of bounds!!! \n";
    }
  }
  //  std::cout<<"Mu number: "<<muNumb<<std::endl;
  //  float weight = (leptonSF+syst*0.01)*trigEff;
  float weight = (leptonSF)*trigEff;

  return weight;

}

float  WZEvent::GetTriggerEfficiency(){

  double triggerEff=1.;

  float zpid = fabs(leptons[zLeptonsIndex[0]].PdgId());
  float wpid = fabs(leptons[wLeptonIndex].PdgId());

  int leptonIndices[3] = { zLeptonsIndex[0],
			   zLeptonsIndex[1],
			   wLeptonIndex };

  float effLead[3],effTrail[3];

  for (int i=0; i<3; i++) { 
    effLead[i]  = leptons[leptonIndices[i]].LeadTriggerEff();
    effTrail[i] = leptons[leptonIndices[i]].TrailTriggerEff();
  }

  if (zpid==wpid 
      || zpid!=wpid   /// JONATAN'S WAY: CORRECT AFTER ALL :)
      ) {
    triggerEff = trigger3sameLeptons(effLead, effTrail);
  } else { // The real thing
    // Here is another way to write it down for the x-channels, 
    // actually completely equivalent

    std::cout << "Trigger eff: alternate formulation - We should not be running this... \n";

    double eff_LL      = trigger2sameLeptons(effLead, effTrail);
    double eff_XTrig1  = triggerDifferentLeptons(effLead, effTrail);
    double eff_XTrig2  = triggerDifferentLeptons(effTrail, effLead);

    double PLLandXT1 = eff_LL*effTrail[2];
    double PLLandXT2 = eff_LL*effLead[2];
    double PXT1_and_XT2 = triggerDifferentLeptons(effLead, effLead);

    double PLLandXT = eff_LL*effLead[2];

    triggerEff = eff_LL + eff_XTrig1 + eff_XTrig2
      - PLLandXT1 - PLLandXT2 - PXT1_and_XT2 + PLLandXT;

  }

  // Sanity check
  if (triggerEff <0. || triggerEff > 1.) {
    std::cout << "NONSENSE trigger eff value: " << triggerEff
	      << "\t W: " << wpid
	      << "\t Z: " << wpid << std::endl;
  }

  return triggerEff;

}


void WZEvent::SmearJets()
{
  // Apply additional jet smearing to account for data/MC
  // differences in resolution, following 2nd way described at:
  // 
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures

  if (!random) random = new TRandom3();


  float sigmaJetResolution[5] = { 0.0342, 0.0363, 0.0504, 0.0534, 0.0898 };

  //  from WW up
  float sigmaJetResolutionUp[5] = { 0.0516, 0.0521, 0.0663, 0.0712, 0.1220 };
  //  from WW down
  float sigmaJetResolutionDown[5] = { 0.0000, 0.0047, 0.0286, 0.0293, 0.0477 };


  SystematicsManager * sysManager = SystematicsManager::GetInstance();

  float jer_syst = sysManager->GetValue("JER");

  float * resolutionCorrections;

  if (jer_syst == 0.) {
    resolutionCorrections = sigmaJetResolution;
  } else if (jer_syst == 1.) {
    resolutionCorrections = sigmaJetResolutionUp;
  } else if (jer_syst == -1.) {
    resolutionCorrections = sigmaJetResolutionDown;
  }

  float etaBoundaries[6] = { 0.0, 0.5, 1.1, 1.7, 2.3, 5.0};

  for (int ijet=0; ijet<recoJets.size(); ijet++) {
    float pt    = recoJets[ijet].Pt();
    float eta   = recoJets[ijet].Eta();
    float phi   = recoJets[ijet].Phi();
    float mass  = recoJets[ijet].M();

    float sigma = -1;

    float abseta = fabs(eta);

    for (int ie=0; ie<5; ie++) {
      if (abseta  >= etaBoundaries[ie]
	  && abseta  < etaBoundaries[ie+1]) {
	sigma = resolutionCorrections[ie];
	break;
      }
    }


    float scale = random->Gaus(1, sigma);
    double newpt   = scale * pt;
    double newmass = mass * pt;

    recoJets[ijet].SetPtEtaPhiM(newpt,eta,phi,newmass);
  }

}


void WZEvent::ApplyJESCorrection(double strength)
{



  // Apply jet energy scale correction for systematic studies
  // see:
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures


  if (!random) random = new TRandom3();

  float etaBoundaries[6] = { 0.0, 0.5, 1.1, 1.7, 2.3, 5.0};

  JetEnergyTool * jesTool = JetEnergyTool::GetInstance();

  for (int ijet=0; ijet<recoJets.size(); ijet++) {
    float pt    = recoJets[ijet].Pt();
    float eta   = recoJets[ijet].Eta();
    float phi   = recoJets[ijet].Phi();
    float mass  = recoJets[ijet].M();
    
    float scale = jesTool->GetJetEnergyScale(pt,eta)*strength;
    //    float xxx = jesTool->GetJetEnergyScale(pt,eta)*strength;

    //    std::cout << "   s=" << scale << " new=" << xxx;

    double newpt   = (1+scale)* pt;
    double newmass = (1+scale)*mass;

    //    std::cout << "\t Jet JES: eta " << eta << " pt : "  << pt << "\t scale = " 
    //	      << scale << " newPt =" << newpt << std::endl;

    recoJets[ijet].SetPtEtaPhiM(newpt,eta,phi,newmass);
    
  }
}

float WZEvent::GetShiftedMET(int met_syst){

  float newmet = -999.;

  MetSystematicsTool * metTool = MetSystematicsTool::GetInstance();

  newmet = metTool->GetMETValue(run,event,met_syst);

  return newmet;

}

void WZEvent::PrintNjetsLeadingJetPt()
{
  int nRecoJets = 0;
  float leadingRecoJetPt = -9999.;
  int leadingRecojet = -1;
    for (int i=0; i<recoJets.size(); i++) {

      if (recoJets[i].Pt() > 30 && fabs(recoJets[i].Eta()) < 2.5) {

	bool closeToLepton = false;
	float drMin = 3.;
	for (int il=0; il<leptons.size(); il++) {
	  if (recoJets[i].DeltaR(leptons[il])<0.5) {
	    closeToLepton = true;
	  }
	}
	if (closeToLepton) continue;
	
	nRecoJets++;
	if (recoJets[i].Pt() > leadingRecoJetPt) {

	  leadingRecoJetPt = recoJets[i].Pt();
	  leadingRecojet = i;
	}
      }
      
    }
    std::cout<<"Njets: "<<nRecoJets<<std::endl;
    std::cout<<"LeadingJet Pt: "<<leadingRecoJetPt<<std::endl;
}
