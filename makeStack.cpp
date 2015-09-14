#include "XSections.h"
#include "MyStyle.h"

#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TMath.h"

#include <iomanip>
#include <cmath>
#include <string>
#include <boost/lexical_cast.hpp>


using namespace std;


int
MakeStack(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50, TFile* f_TT,
          TString histKey, unsigned int channel, TString xAxisTitle = "", double binWidth = 1.)
{
  TH1D* h_WZ = (TH1D*) (f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*) (f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*) (f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*) (f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*) (f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_TT = (TH1D*) (f_TT->Get(histKey))->Clone(histKey + "_TT");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_TT->Scale(scale_TT);

  THStack* stack = new THStack(histKey, histKey);
  stack->Add(h_WW);
  stack->Add(h_ZZ2L2Q);  
  stack->Add(h_TT);    
  stack->Add(h_ZZ4L);
  stack->Add(h_DYM50);
  stack->Add(h_WZ);   

  h_WZ->SetFillColor(kOrange-2);
  h_ZZ4L->SetFillColor(kRed+1);
  h_ZZ2L2Q->SetFillColor(kRed+2);
  h_WW->SetFillColor(kYellow-2);
  h_DYM50->SetFillColor(kGray+1);
  h_TT->SetFillColor(kAzure);

  stack->SetMinimum(0.0);
  stack->SetMaximum(stack->GetMaximum() * 1.33);

 	stack->Draw();
  stack->GetXaxis()->SetTitle(xAxisTitle);
  if (binWidth > 0.)  stack->GetYaxis()->SetTitle(Form("Events / %1.1f GeV", binWidth));
  else  stack->GetYaxis()->SetTitle("Events");

  stack->GetXaxis()->SetLabelFont(132);
  stack->GetYaxis()->SetLabelFont(132);
  stack->GetXaxis()->SetLabelOffset(0.007);
  stack->GetYaxis()->SetLabelOffset(0.007);
  stack->GetXaxis()->SetLabelSize(0.03);
  stack->GetYaxis()->SetLabelSize(0.03);
  stack->GetXaxis()->SetTitleFont(132);
  stack->GetYaxis()->SetTitleFont(132);
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.04);
  stack->GetXaxis()->SetTitleOffset(1.);
  stack->GetYaxis()->SetTitleOffset(1.5);
  stack->GetXaxis()->SetNdivisions(510);
  stack->GetYaxis()->SetNdivisions(510);

  TLegend* leg = new TLegend(0.7, 0.75, 0.9, 0.95);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_DYM50, "Drell Yan (Z + jets)", "f");
  leg->AddEntry(h_ZZ4L, "ZZ #rightarrow 4l", "f");
  leg->AddEntry(h_TT, "TTJets", "f");
  leg->AddEntry(h_ZZ2L2Q, "ZZ #rightarrow 2l2q", "f");
  leg->AddEntry(h_WW, "WW", "f");

	leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
	latexLabel.SetTextAlign(12);
  latexLabel.SetTextSize(0.04);
  latexLabel.DrawLatex(0.19, 0.9, Form("#font[132]{#intL dt = %1.1f fb^{-1}}", 1.));

  latexLabel.SetTextSize(0.05);
  if (channel == 0) latexLabel.DrawLatex(0.46, 0.95, "#font[132]{Channel eee}");
  if (channel == 1) latexLabel.DrawLatex(0.46, 0.95, "#font[132]{Channel ee#mu}");
  if (channel == 2) latexLabel.DrawLatex(0.46, 0.95, "#font[132]{Channel #mu#mue}");
  if (channel == 3) latexLabel.DrawLatex(0.46, 0.95, "#font[132]{Channel #mu#mu#mu}");
  if (channel == 4) latexLabel.DrawLatex(0.46, 0.95, "All Channels");

  return 1;
}


int
MakeStackDeltaR(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50, TFile* f_TT,
            		TString histKey, TString xAxisTitle, double max)
{
  TH1D* h_WZ = (TH1D*) (f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*) (f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*) (f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*) (f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*) (f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_TT = (TH1D*) (f_TT->Get(histKey))->Clone(histKey + "_TT");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_TT->Scale(scale_TT);

  THStack* stack = new THStack(histKey, histKey);
  stack->Add(h_WW);
  stack->Add(h_ZZ2L2Q);  
  stack->Add(h_TT);    
  stack->Add(h_ZZ4L);
  stack->Add(h_DYM50);
  stack->Add(h_WZ);   

  h_WZ->SetFillColor(kOrange-2);
  h_ZZ4L->SetFillColor(kRed+1);
  h_ZZ2L2Q->SetFillColor(kRed+2);
  h_WW->SetFillColor(kYellow-2);
  h_DYM50->SetFillColor(kGray+1);
  h_TT->SetFillColor(kAzure);

  stack->SetMinimum(0.0);
  stack->SetMaximum(max);

 	stack->Draw();
  stack->GetXaxis()->SetTitle(xAxisTitle);
  stack->GetYaxis()->SetTitle("Jets / 0.02");

  stack->GetXaxis()->SetLabelFont(132);
  stack->GetYaxis()->SetLabelFont(132);
  stack->GetXaxis()->SetLabelOffset(0.007);
  stack->GetYaxis()->SetLabelOffset(0.007);
  stack->GetXaxis()->SetLabelSize(0.03);
  stack->GetYaxis()->SetLabelSize(0.03);
  stack->GetXaxis()->SetTitleFont(132);
  stack->GetYaxis()->SetTitleFont(132);
  stack->GetXaxis()->SetTitleSize(0.045);
  stack->GetYaxis()->SetTitleSize(0.04);
  stack->GetXaxis()->SetTitleOffset(1.);
  stack->GetYaxis()->SetTitleOffset(1.5);
  stack->GetXaxis()->SetNdivisions(510);
  stack->GetYaxis()->SetNdivisions(510);

  TLegend* leg = new TLegend(0.7, 0.75, 0.9, 0.95);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(h_WZ, "WZ", "f");
  leg->AddEntry(h_DYM50, "Drell Yan (Z + jets)", "f");
  leg->AddEntry(h_ZZ4L, "ZZ #rightarrow 4l", "f");
  leg->AddEntry(h_TT, "TTJets", "f");
  leg->AddEntry(h_ZZ2L2Q, "ZZ #rightarrow 2l2q", "f");
  leg->AddEntry(h_WW, "WW", "f");

	leg->Draw();

  TLatex latexLabel;
  latexLabel.SetNDC();
	latexLabel.SetTextAlign(12);
  latexLabel.SetTextSize(0.04);
  latexLabel.DrawLatex(0.19, 0.9, Form("#font[132]{#intL dt = %1.1f fb^{-1}}", 1.));

  return 1;
}


pair<double, double>
SignalAboveCutBin(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50, TFile* f_TT,
            		TString histKey, TFile* f_Signal, double scaleSignal, unsigned int cutBin)
{
  TH1D* h_WZ = (TH1D*) (f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*) (f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*) (f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*) (f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*) (f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_TT = (TH1D*) (f_TT->Get(histKey))->Clone(histKey + "_TT");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_TT->Scale(scale_TT);

  TH1D* h_All = (TH1D*) h_WZ->Clone("h_All");
  h_All->Add(h_WW);
  h_All->Add(h_ZZ2L2Q);  
  h_All->Add(h_TT);    
  h_All->Add(h_ZZ4L);
  h_All->Add(h_DYM50);

  TH1D* h_Signal = (TH1D*) (f_Signal->Get(histKey))->Clone("h_Signal");
  h_Signal->Scale(scaleSignal);

  double signal = h_Signal->Integral(cutBin, 300);
  double total = h_All->Integral(cutBin, 300);
  double signalFraction = 0;
  double significance = 0;
  if (total) {
    signalFraction = signal / total;
    significance = signal / sqrt(total);
  } else {
    signalFraction = 0;
    significance = 0;
  }

  pair<double, double> result = make_pair(signalFraction, significance);
  return result;
}


pair<double, double>
TotalSignalAboveCutBin(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50, TFile* f_TT,
            		       TString histKey, unsigned int cutBin)
{
  TH1D* h_WZ = (TH1D*) (f_WZ->Get(histKey))->Clone(histKey + "_WZ");
  TH1D* h_ZZ4L = (TH1D*) (f_ZZ4L->Get(histKey))->Clone(histKey + "_ZZ4L");
  TH1D* h_ZZ2L2Q = (TH1D*) (f_ZZ2L2Q->Get(histKey))->Clone(histKey + "_ZZ2L2Q");
  TH1D* h_WW = (TH1D*) (f_WW->Get(histKey))->Clone(histKey + "_WW");
  TH1D* h_DYM50 = (TH1D*) (f_DYM50->Get(histKey))->Clone(histKey + "_DYM50");
  TH1D* h_TT = (TH1D*) (f_TT->Get(histKey))->Clone(histKey + "_TT");

  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_ZZ2L2Q = EXPCMS_LUMINOSITY /  (NG_ZZTo2L2Q / XS_ZZTo2L2Q);
  const double scale_WW = EXPCMS_LUMINOSITY /  (NG_WWTo2L2Nu / XS_WWTo2L2Nu);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);

  h_WZ->Scale(scale_WZ);
  h_ZZ4L->Scale(scale_ZZ4L);
  h_ZZ2L2Q->Scale(scale_ZZ2L2Q);
  h_WW->Scale(scale_WW);
  h_DYM50->Scale(scale_DYM50);
  h_TT->Scale(scale_TT);

  TH1D* h_All = (TH1D*) h_WZ->Clone("h_All");
  h_All->Add(h_WW);
  h_All->Add(h_ZZ2L2Q);  
  h_All->Add(h_TT);    
  h_All->Add(h_ZZ4L);
  h_All->Add(h_DYM50);

  double signal = h_All->Integral(cutBin, 300);
  double total = h_All->Integral(0, 300);
  double signalFraction = 0;
  double significance = 0;
  if (total) {
    signalFraction = signal / total;
    significance = signal / sqrt(total);
  } else {
    signalFraction = 0;
    significance = 0;
  }

  pair<double, double> result = make_pair(signalFraction, significance);
  return result;
}


int
MakeSignificanceGraphs(TFile* f_WZ, TFile* f_ZZ4L, TFile* f_ZZ2L2Q, TFile* f_WW, TFile* f_DYM50, TFile* f_TT,
                TString histKey, unsigned int channel, TCanvas* canvas, TString xAxisTitle,
                double start, double end, double binWidth)
{
  const double scale_WZ = EXPCMS_LUMINOSITY / (NG_WZTo3LNu / XS_WZTo3LNu);
  const double scale_ZZ4L = EXPCMS_LUMINOSITY /  (NG_ZZTo4L / XS_ZZTo4L);
  const double scale_DYM50 = EXPCMS_LUMINOSITY /  (NG_DYJetsToLL_M50 / XS_DYJetsToLL_M50);
  const double scale_TT = EXPCMS_LUMINOSITY /  (NG_TTJets / XS_TTJets);

  vector<double> xValues;
  vector<double> totalSignalFractions, totalSignificances;
  vector<double> signalFractions_WZ, signalFractions_DYM50, signalFractions_ZZ4L, signalFractions_TT;
  vector<double> significances_WZ, significances_DYM50, significances_ZZ4L, significances_TT;
  const unsigned int endBin = (unsigned int)((end - start) / binWidth);
  for (unsigned int i = 0; i < endBin+1; i++) {
    const double xValue = start + i * binWidth;
    xValues.push_back(xValue);
    const double signalFraction_WZ =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_WZ, scale_WZ, i).first;
    signalFractions_WZ.push_back(signalFraction_WZ);
    const double significance_WZ =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_WZ, scale_WZ, i).second;
    significances_WZ.push_back(significance_WZ);
    const double signalFraction_DYM50 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_DYM50, scale_DYM50, i).first;
    signalFractions_DYM50.push_back(signalFraction_DYM50);
    const double significance_DYM50 =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_DYM50, scale_DYM50, i).second;
    significances_DYM50.push_back(significance_DYM50);
    const double signalFraction_ZZ4L =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_ZZ4L, scale_ZZ4L, i).first;
    signalFractions_ZZ4L.push_back(signalFraction_ZZ4L);
    const double significance_ZZ4L =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_ZZ4L, scale_ZZ4L, i).second;
    significances_ZZ4L.push_back(significance_ZZ4L);
    const double signalFraction_TT =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_TT, scale_TT, i).first;
    signalFractions_TT.push_back(signalFraction_TT);
    const double significance_TT =
      SignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
            		        histKey, f_TT, scale_TT, i).second;
    significances_TT.push_back(significance_TT);
    const double signalFraction_Total =
      TotalSignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, histKey, i).first;
    totalSignalFractions.push_back(signalFraction_Total);
    const double significance_Total =
      TotalSignalAboveCutBin(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT, histKey, i).second;
    totalSignificances.push_back(significance_Total);
  }

  canvas->Divide(1, 2);
  
  canvas->cd(1);
  TGraph* g_sF_WZ = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_WZ.at(0));
  TGraph* g_sF_DYM50 = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_DYM50.at(0));
  TGraph* g_sF_ZZ4L = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_ZZ4L.at(0));
  TGraph* g_sF_TT = new TGraph(xValues.size(), &xValues.at(0), &signalFractions_TT.at(0));
  TGraph* g_sF_Total = new TGraph(xValues.size(), &xValues.at(0), &totalSignalFractions.at(0));

  g_sF_WZ->SetMarkerColor(kOrange-2);
  g_sF_WZ->SetMarkerSize(0.6);
  g_sF_WZ->SetMaximum(1.44);
  g_sF_WZ->SetMinimum(0.);
  g_sF_WZ->Draw("AP");
  g_sF_WZ->GetYaxis()->SetTitle("S/(S+B)");
  g_sF_WZ->GetXaxis()->SetTitle(xAxisTitle);

  g_sF_WZ->GetXaxis()->SetLabelFont(132);
  g_sF_WZ->GetYaxis()->SetLabelFont(132);
  g_sF_WZ->GetXaxis()->SetLabelOffset(0.007);
  g_sF_WZ->GetYaxis()->SetLabelOffset(0.007);
  g_sF_WZ->GetXaxis()->SetLabelSize(0.03);
  g_sF_WZ->GetYaxis()->SetLabelSize(0.03);
  g_sF_WZ->GetXaxis()->SetTitleFont(132);
  g_sF_WZ->GetYaxis()->SetTitleFont(132);
  g_sF_WZ->GetXaxis()->SetTitleSize(0.045);
  g_sF_WZ->GetYaxis()->SetTitleSize(0.04);
  g_sF_WZ->GetXaxis()->SetTitleOffset(1.);
  g_sF_WZ->GetYaxis()->SetTitleOffset(1.5);
  g_sF_WZ->GetXaxis()->SetNdivisions(510);
  g_sF_WZ->GetYaxis()->SetNdivisions(510);

  g_sF_DYM50->SetMarkerColor(kGray+1);
  g_sF_DYM50->SetMarkerSize(0.6);
  g_sF_DYM50->SetMaximum(1.44);
  g_sF_DYM50->SetMinimum(0.);
  g_sF_DYM50->Draw("P");

  g_sF_ZZ4L->SetMarkerColor(kRed+1);
  g_sF_ZZ4L->SetMarkerSize(0.6);
  g_sF_ZZ4L->SetMaximum(1.44);
  g_sF_ZZ4L->SetMinimum(0.);
  g_sF_ZZ4L->Draw("P");

  g_sF_TT->SetMarkerColor(kAzure);
  g_sF_TT->SetMarkerSize(0.6);
  g_sF_TT->SetMaximum(1.44);
  g_sF_TT->SetMinimum(0.);
  g_sF_TT->Draw("P");

  g_sF_Total->SetMarkerColor(kBlack);
  g_sF_Total->SetMarkerSize(0.4);
  g_sF_Total->SetMaximum(1.44);
  g_sF_Total->SetMinimum(0.);
  g_sF_Total->Draw("*");

  TLegend* leg = new TLegend(0.18, 0.75, 0.38, 0.95);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetShadowColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(g_sF_WZ, "WZ", "p");
  leg->AddEntry(g_sF_DYM50, "Drell Yan (Z + jets)", "p");
  leg->AddEntry(g_sF_ZZ4L, "ZZ #rightarrow 4l", "p");
  leg->AddEntry(g_sF_TT, "TTJets", "p");
  leg->AddEntry(g_sF_Total, "Total", "p");

  leg->Draw();

  TLatex latexLabel_sF;
  latexLabel_sF.SetNDC();
  latexLabel_sF.SetTextAlign(12);
  latexLabel_sF.SetTextSize(0.05);
  latexLabel_sF.DrawLatex(0.8, 0.9, "Signal Fraction");

  if (channel == 0) latexLabel_sF.DrawLatex(0.5, 0.95, "#font[132]{Channel eee}");
  if (channel == 1) latexLabel_sF.DrawLatex(0.5, 0.95, "#font[132]{Channel ee#mu}");
  if (channel == 2) latexLabel_sF.DrawLatex(0.5, 0.95, "#font[132]{Channel #mu#mue}");
  if (channel == 3) latexLabel_sF.DrawLatex(0.5, 0.95, "#font[132]{Channel #mu#mu#mu}");
  if (channel == 4) latexLabel_sF.DrawLatex(0.5, 0.95, "All Channels");

  canvas->cd(2);
  TGraph* g_significance_WZ = new TGraph(xValues.size(), &xValues.at(0), &significances_WZ.at(0));
  TGraph* g_significance_DYM50 = new TGraph(xValues.size(), &xValues.at(0), &significances_DYM50.at(0));
  TGraph* g_significance_ZZ4L = new TGraph(xValues.size(), &xValues.at(0), &significances_ZZ4L.at(0));
  TGraph* g_significance_TT = new TGraph(xValues.size(), &xValues.at(0), &significances_TT.at(0));
  TGraph* g_significance_Total = new TGraph(xValues.size(), &xValues.at(0), &totalSignificances.at(0));

  vector<double> max_significances =
    { TMath::MaxElement(g_significance_WZ->GetN(), g_significance_WZ->GetY()),
      TMath::MaxElement(g_significance_DYM50->GetN(), g_significance_DYM50->GetY()),
      TMath::MaxElement(g_significance_ZZ4L->GetN(), g_significance_ZZ4L->GetY()),
      TMath::MaxElement(g_significance_TT->GetN(), g_significance_TT->GetY()),
      TMath::MaxElement(g_significance_Total->GetN(), g_significance_Total->GetY()) };
  double max_significance = *(max_element(max_significances.begin(), max_significances.end()));

  g_significance_WZ->SetMarkerColor(kOrange-2);
  g_significance_WZ->SetMarkerSize(0.6);
  g_significance_WZ->SetMaximum(max_significance * 1.44);
  g_significance_WZ->SetMinimum(0.);
  g_significance_WZ->Draw("AP");
  g_significance_WZ->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  g_significance_WZ->GetXaxis()->SetTitle(xAxisTitle);

  g_significance_WZ->GetXaxis()->SetLabelFont(132);
  g_significance_WZ->GetYaxis()->SetLabelFont(132);
  g_significance_WZ->GetXaxis()->SetLabelOffset(0.007);
  g_significance_WZ->GetYaxis()->SetLabelOffset(0.007);
  g_significance_WZ->GetXaxis()->SetLabelSize(0.03);
  g_significance_WZ->GetYaxis()->SetLabelSize(0.03);
  g_significance_WZ->GetXaxis()->SetTitleFont(132);
  g_significance_WZ->GetYaxis()->SetTitleFont(132);
  g_significance_WZ->GetXaxis()->SetTitleSize(0.045);
  g_significance_WZ->GetYaxis()->SetTitleSize(0.04);
  g_significance_WZ->GetXaxis()->SetTitleOffset(1.);
  g_significance_WZ->GetYaxis()->SetTitleOffset(1.5);
  g_significance_WZ->GetXaxis()->SetNdivisions(510);
  g_significance_WZ->GetYaxis()->SetNdivisions(510);

  g_significance_DYM50->SetMarkerColor(kGray+1);
  g_significance_DYM50->SetMarkerSize(0.6);
  g_significance_DYM50->SetMaximum(max_significance * 1.4);
  g_significance_DYM50->SetMinimum(0.);
  g_significance_DYM50->Draw("P");

  g_significance_ZZ4L->SetMarkerColor(kRed+1);
  g_significance_ZZ4L->SetMarkerSize(0.6);
  g_significance_ZZ4L->SetMaximum(max_significance * 1.44);
  g_significance_ZZ4L->SetMinimum(0.);
  g_significance_ZZ4L->Draw("P");

  g_significance_TT->SetMarkerColor(kAzure);
  g_significance_TT->SetMarkerSize(0.6);
  g_significance_TT->SetMaximum(max_significance * 1.44);
  g_significance_TT->SetMinimum(0.);
  g_significance_TT->Draw("P");

  g_significance_Total->SetMarkerColor(kBlack);
  g_significance_Total->SetMarkerSize(0.4);
  g_significance_Total->SetMaximum(max_significance * 1.44);
  g_significance_Total->SetMinimum(0.);
  g_significance_Total->Draw("*");

  leg->Draw();

  TLatex latexLabel_significance;
  latexLabel_significance.SetNDC();
  latexLabel_significance.SetTextAlign(12);
  latexLabel_significance.SetTextSize(0.05);
  latexLabel_significance.DrawLatex(0.82, 0.9, "Significance");

  if (channel == 0) latexLabel_significance.DrawLatex(0.5, 0.95, "#font[132]{Channel eee}");
  if (channel == 1) latexLabel_significance.DrawLatex(0.5, 0.95, "#font[132]{Channel ee#mu}");
  if (channel == 2) latexLabel_significance.DrawLatex(0.5, 0.95, "#font[132]{Channel #mu#mue}");
  if (channel == 3) latexLabel_significance.DrawLatex(0.5, 0.95, "#font[132]{Channel #mu#mu#mu}");
  if (channel == 4) latexLabel_significance.DrawLatex(0.5, 0.95, "All Channels");

  return 1;
}


int
main()
{
  const MyStyle rootStyle(600);

  TFile* f_WZ   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/WZ_eleTight.root");
  TFile* f_ZZ4L   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/ZZ4L_eleTight.root");
  TFile* f_ZZ2L2Q   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/ZZ2L2Q_eleTight.root");
  TFile* f_WW   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/WW_eleTight.root");
  TFile* f_DYM50   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/DYM50_eleTight.root");
  TFile* f_TT   = TFile::Open("/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/TT_eleTight.root");

  unsigned int n = 5;

  vector<string> histoName = { "Zmass", "Zpt", "MET", "Mt", "3LMass", "Zl1pt", "Zl2pt", "Wlpt" };
  vector<string> jetsName = { "NJets", "NJetsNoMuIso", "NJetsNoEleIso", "NJetsNoIso" };
  const string path = "/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/plots/stack/";
  vector<string> xAxisHisto =
    { "M_{Z} [GeV]", "Z p_{t} [GeV]","missing E_{t} [GeV]", "M_{t} [GeV]", "M_{3l} [GeV]",
      "Zl_{lead} p_{t} [GeV]", "Zl_{trail} p_{t} [GeV]", "Wl p_{t} [GeV]" };
  const string xAxisJets = "Number of Jets";
  vector<double> binWidthHisto = { 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 };

  for (unsigned int histo = 0; histo < histoName.size(); histo++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((histoName.at(histo) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              histoName.at(histo).c_str());
      canvas[i]->cd();
      MakeStack(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
                ("h" + histoName.at(histo) + "_" + boost::lexical_cast<string>(i)).c_str(),
                i, xAxisHisto.at(histo).c_str(), binWidthHisto.at(histo));
      canvas[i]->SaveAs((path + histoName.at(histo) + "_" +
                        boost::lexical_cast<string>(i) + "_eleTight.pdf").c_str());
      delete canvas[i];
    }
  }

  for (unsigned int jets = 0; jets < jetsName.size(); jets++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((jetsName.at(jets) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              jetsName.at(jets).c_str());
      canvas[i]->cd();
      MakeStack(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
                ("h" + jetsName.at(jets) + "_" + boost::lexical_cast<string>(i)).c_str(),
                i, xAxisJets.c_str(), 0.);
      canvas[i]->SaveAs((path + jetsName.at(jets) + "_" +
                        boost::lexical_cast<string>(i) + "_eleTight.pdf").c_str());
      delete canvas[i];
    }
  }

  vector<string> deltaRJets = { "DeltaRL", "DeltaRMu", "DeltaREle" };
  vector<string> xAxisDeltaR =
    { "#deltaR_{min}(jet, lepton)", "#deltaR_{min}(jet, muon)", "#deltaR_{min}(jet, electron)" };
  for (unsigned int i = 0; i < deltaRJets.size(); i++) {
    TCanvas* canvas = new TCanvas(deltaRJets.at(i).c_str(), deltaRJets.at(i).c_str());
    canvas->cd();
    double max = 10000.;
    i ? max = 6000. : max = 12000.;
    MakeStackDeltaR(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
                    ("h" + deltaRJets.at(i)).c_str(), xAxisDeltaR.at(i).c_str(), max);
    canvas->SaveAs((path + deltaRJets.at(i) + "_eleTight.pdf").c_str());
    delete canvas;
  }

  const string pathGraph = "/users/msasa/work/cms/wz/ggAna/code/WZRun2Analysis/output/selection/WSelection/plots/cuts/";
  vector<string> graphName = { "3LMass", "MET" };
  vector<string> xAxisGraph = { "M_{3l} [GeV]","missing E_{t} [GeV]" };
  vector<double> startGraph = { 50., 0. };
  vector<double> endGraph = { 350., 200. };
  vector<double> binWidthGraph = { 2., 2. };
  for (unsigned int graph = 0; graph < graphName.size(); graph++) {
    TCanvas* canvas[n];
    for (unsigned int i = 0; i < n; i++) {
      canvas[i] = new TCanvas((graphName.at(graph) + "_" + boost::lexical_cast<string>(i)).c_str(),
                              graphName.at(graph).c_str());
      canvas[i]->cd();
      MakeSignificanceGraphs(f_WZ, f_ZZ4L, f_ZZ2L2Q, f_WW, f_DYM50, f_TT,
                             ("h" + graphName.at(graph) + "_" + boost::lexical_cast<string>(i)).c_str(),
                             i, canvas[i], xAxisGraph.at(graph).c_str(),
                             startGraph.at(graph), endGraph.at(graph), binWidthGraph.at(graph));
      canvas[i]->SaveAs((pathGraph + graphName.at(graph) + "_" +
                        boost::lexical_cast<string>(i) + "_eleTight.pdf").c_str());
      delete canvas[i];
    }
  }

  return 1;
}

