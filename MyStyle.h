#include <TStyle.h>
#include <TROOT.h>


class MyStyle {

  static const int fN = 2;

public:
  MyStyle(const int width, const double aspect = 4/3.)
  {
    fStyle[0] = new TStyle("myRectangularStyle", "foo");
    fStyle[1] = new TStyle("myQuadraticStyle", "bar");

    for (int i = 0; i < fN; ++i) {
      fStyle[i]->SetCanvasDefH(i ? width : int(width / aspect));
      fStyle[i]->SetCanvasDefW(width);
      fStyle[i]->SetFrameBorderMode(0);
      fStyle[i]->SetCanvasBorderMode(0);
      fStyle[i]->SetPadBorderMode(0);
      fStyle[i]->SetPadColor(0);
      fStyle[i]->SetCanvasColor(0);
      fStyle[i]->SetTitleFillColor(0);
      fStyle[i]->SetTitleBorderSize(1);
      fStyle[i]->SetStatColor(0);
      fStyle[i]->SetStatBorderSize(1);
      fStyle[i]->SetOptTitle(0);
      fStyle[i]->SetOptStat(1);
      fStyle[i]->SetOptFit(0);
      fStyle[i]->SetPalette(1,0);
      fStyle[i]->SetTitleBorderSize(0); // border size of Title PavelLabel
      fStyle[i]->SetTitleX(0.1f);
      fStyle[i]->SetTitleW(0.8f);
      fStyle[i]->SetStatY(0.975f);
      fStyle[i]->SetStatX(0.95f);
      fStyle[i]->SetStatW(0.2f);
      fStyle[i]->SetStatH(0.15f);
      fStyle[i]->SetTitleXOffset(1);
      fStyle[i]->SetTitleYOffset(i ? 1.4f : 1.f);
      fStyle[i]->SetMarkerStyle(20);
      fStyle[i]->SetMarkerSize(0.5);
      fStyle[i]->SetFuncWidth(1.);
      fStyle[i]->SetPadTopMargin(0.025f);
      fStyle[i]->SetPadBottomMargin(0.175f);
      fStyle[i]->SetPadLeftMargin(i ? 0.18f : 0.15f);
      fStyle[i]->SetPadRightMargin(0.05f);
      fStyle[i]->SetTitleSize(0.06f, "X");
      fStyle[i]->SetTitleSize(0.06f, "Y");
      fStyle[i]->SetErrorX(0.);
      fStyle[i]->SetNumberContours(99);
    }

    SetRectangular();
  }

  ~MyStyle()
  {
    for (int i = 0; i < fN; ++i)
      delete fStyle[i];
  }

  void
  SetRectangular()
    const
  {
    gROOT->SetStyle("myRectangularStyle");
    gROOT->ForceStyle();
  }

  void
  SetQuadratic()
    const
  {
    gROOT->SetStyle("myQuadraticStyle");
    gROOT->ForceStyle();
  }

private:
  MyStyle(const MyStyle&);
  MyStyle& operator=(const MyStyle&);

  TStyle* fStyle[fN];
};
