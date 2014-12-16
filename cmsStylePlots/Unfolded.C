#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMath.h"
#include "TColor.h"


const Font_t _cmsTextFont   =  61;
const Font_t _extraTextFont =  52;
const Font_t _lumiTextFont  =  42;
const Size_t _msize         = 1.1;


void Unfolded(TString variable = "p_{T,max}^{#font[12]{l}}",
	      TString filename = "XSLeadingPt_AN")
{  
  gInterpreter->ExecuteMacro("GoodStyle.C");
 
  TFile* file = new TFile(filename + ".root", "read");
 
  TH1F* xsValue          = (TH1F*)xsValue->Clone();
  TH1F* xsValue_Powheg   = (TH1F*)xsValue_Powheg->Clone();
  TH1F* xsValue_Madgraph = (TH1F*)xsValue_Madgraph->Clone();
  TH1F* xsValue_MCnlo    = (TH1F*)xsValue_MCnlo->Clone();
  TH1F* systHisto        = (TH1F*)systHisto->Clone();


  // Set the data errors
  //----------------------------------------------------------------------------
  for (UInt_t i=1; i<xsValue->GetNbinsX(); i++) {
  
    float err_total  = systHisto->GetBinError(i); // included already both, statistical and systematics
    
    xsValue->SetBinError(i, err_total);
  }


  // Data cosmetics
  //----------------------------------------------------------------------------
  xsValue->SetLineWidth(1);
  xsValue->SetMarkerSize(_msize);
  xsValue->SetMarkerStyle(kFullCircle);


  // Powheg cosmetics
  //----------------------------------------------------------------------------
  xsValue_Powheg->SetFillColor(kGreen-9);
  xsValue_Powheg->SetFillStyle(1001);
  xsValue_Powheg->SetLineColor(kGreen+3);
  xsValue_Powheg->SetLineWidth(1);
  xsValue_Powheg->SetMarkerColor(kGreen+3);
  xsValue_Powheg->SetMarkerSize(_msize);
  xsValue_Powheg->SetMarkerStyle(22);


  // Madgraph cosmetics
  //----------------------------------------------------------------------------

  // Only includes statistical error. 

  xsValue_Madgraph->SetFillColor(kMagenta-10);
  xsValue_Madgraph->SetFillStyle(1001);
  xsValue_Madgraph->SetLineColor(kPink+10);
  xsValue_Madgraph->SetLineWidth(1);
  xsValue_Madgraph->SetMarkerColor(kPink+10);
  xsValue_Madgraph->SetMarkerSize(_msize);
  xsValue_Madgraph->SetMarkerStyle(21);


  // MCNLO cosmetics
  //----------------------------------------------------------------------------
  xsValue_MCnlo->SetFillColor(kOrange);
  xsValue_MCnlo->SetFillStyle(1001);
  xsValue_MCnlo->SetLineColor(kRed);
  xsValue_MCnlo->SetLineWidth(1);
  xsValue_MCnlo->SetMarkerColor(kRed);
  xsValue_MCnlo->SetMarkerSize(_msize);
  xsValue_MCnlo->SetMarkerStyle(24);


  // Set the canvas and pads
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas("wwxs", "wwxs", 600, 850);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.55, 1, 1.000);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.39, 1, 0.552);
  TPad* pad3 = new TPad("pad3", "pad3", 0, 0.23, 1, 0.392);
  TPad* pad4 = new TPad("pad4", "pad4", 0, 0.00, 1, 0.232);
 
  pad1->SetTopMargin(0.09);
  pad2->SetTopMargin(0);
  pad3->SetTopMargin(0);
  pad4->SetTopMargin(0);

  pad1->SetBottomMargin(0);
  pad2->SetBottomMargin(0);
  pad3->SetBottomMargin(0);
  pad4->SetBottomMargin(0.35);

  pad1->SetLeftMargin(0.16);
  pad2->SetLeftMargin(0.16);
  pad3->SetLeftMargin(0.16);
  pad4->SetLeftMargin(0.16);


  // pad1
  //----------------------------------------------------------------------------
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();


  // Draw
  //----------------------------------------------------------------------------
  AxisFonts(xsValue->GetXaxis(), variable + " (GeV)");
  AxisFonts(xsValue->GetYaxis(), "#frac{1}{#sigma} d#sigma(WW#rightarrow#mu#nue#nu + <1 jet)/d" + variable);
 
  TH1F* hpowError = (TH1F*)xsValue_Powheg->Clone();
  TH1F* hmadError = (TH1F*)xsValue_Madgraph->Clone();
  TH1F* hmcError  = (TH1F*)xsValue_MCnlo->Clone();

  xsValue         ->Draw("pz");
  hmadError       ->Draw("e2,same"); 
  hmcError        ->Draw("e2,same");
  hpowError       ->Draw("e2,same");
  xsValue_Madgraph->Draw("pz,same");
  xsValue_MCnlo   ->Draw("pz,same");
  xsValue_Powheg  ->Draw("pz,same");
  xsValue         ->Draw("pz,same");
 
  if ( filename == "XSLeadingPt_AN")  xsValue->SetMinimum(1.1e-4);
  if ( filename == "XSDileptonPt_AN")  xsValue->SetMinimum(5.0e-5);
  if ( filename == "XSMll_AN")  xsValue->SetMinimum(9.0e-4);
  if ( filename == "XSDeltaPhi_AN")  { 
    xsValue->SetMinimum(9.0e-2);
    xsValue->SetMaximum(1.1);
  }


 
  // Legend
  //----------------------------------------------------------------------------

  if ( filename != "XSDeltaPhi_AN") { 
    DrawLegend(0.718, 0.80, xsValue,   " Data",     "lp");
    DrawLegend(0.718, 0.74, hmadError, " Madgraph", "flp");  
    DrawLegend(0.718, 0.68, hmcError,  " MC@NLO",   "flp");
    DrawLegend(0.718, 0.62, hpowError, " Powheg",   "flp");
  } else {
    DrawLegend(0.218, 0.80, xsValue,   " Data",     "lp");
    DrawLegend(0.218, 0.74, hmadError, " Madgraph", "flp");  
    DrawLegend(0.218, 0.68, hmcError,  " MC@NLO",   "flp");
    DrawLegend(0.218, 0.62, hpowError, " Powheg",   "flp");
  }




  // Draw text 
  //----------------------------------------------------------------------------
  DrawLatex(_cmsTextFont,   0.173, 0.935, 0.065, 11, "CMS");
  DrawLatex(_extraTextFont, 0.268, 0.935, 0.035, 11, "Preliminary");
  DrawLatex(_lumiTextFont,  0.940, 0.935, 0.050, 31, "19.4 fb^{-1} (8 TeV)");


  // Prepare the ratios
  //----------------------------------------------------------------------------
  TH1F* ratio_pow    = xsValue_Powheg->Clone("ratio");
  TH1F* ratio_mad    = xsValue_Madgraph->Clone("ratio");
  TH1F* ratio_mcnlo  = xsValue_MCnlo->Clone("ratio");
  TH1F* hratio_pow   = xsValue_Powheg->Clone("ratio");
  TH1F* hratio_mad   = xsValue_Madgraph->Clone("ratio");
  TH1F* hratio_mcnlo = xsValue_MCnlo->Clone("ratio");
  TH1F* ratioErr     = xsValue->Clone("ratio");


  ratioErr->SetFillColor  (kGray+2);
  ratioErr->SetFillStyle  (   3004);
  ratioErr->SetLineColor  (kGray+2);
  ratioErr->SetMarkerColor(kGray+2);
  ratioErr->SetMarkerSize (      0);


  // Set the bin content
  //----------------------------------------------------------------------------
  for (UInt_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {
   
    Double_t powValue = xsValue_Powheg->GetBinContent(ibin);
    Double_t powError = xsValue_Powheg->GetBinError  (ibin);
   
    Double_t madValue = xsValue_Madgraph->GetBinContent(ibin);
    Double_t madError = xsValue_Madgraph->GetBinError  (ibin);
   
    Double_t mcnloValue = xsValue_MCnlo->GetBinContent(ibin);
    Double_t mcnloError = xsValue_MCnlo->GetBinError  (ibin);
   
    Double_t dataValue = xsValue->GetBinContent(ibin);
    Double_t statError = xsValue->GetBinError  (ibin);
    Double_t systError = systHisto->GetBinError(ibin);
   
    Double_t dataError = systError;
   
    Double_t ratioValue_pow = (powValue > 0) ? powValue / dataValue : 0.0;
    Double_t ratioError_pow = (powValue > 0) ? powError / dataValue : 0.0;
   
    Double_t ratioValue_mad = (madValue > 0) ? madValue / dataValue : 0.0;
    Double_t ratioError_mad = (madValue > 0) ? madError / dataValue : 0.0;
   
    Double_t ratioValue_mcnlo = (mcnloValue > 0) ? mcnloValue / dataValue : 0.0;
    Double_t ratioError_mcnlo = (mcnloValue > 0) ? mcnloError / dataValue : 0.0;
   
    Double_t uncertaintyError = (dataValue > 0) ? dataError / dataValue : 0.0;
   
    ratio_pow ->SetBinContent(ibin, ratioValue_pow);
    hratio_pow->SetBinContent(ibin, ratioValue_pow);
    hratio_pow->SetBinError  (ibin, ratioError_pow);
   
    ratio_mad ->SetBinContent(ibin, ratioValue_mad);
    hratio_mad->SetBinContent(ibin, ratioValue_mad);
    hratio_mad->SetBinError  (ibin, ratioError_mad);
   
    ratio_mcnlo ->SetBinContent(ibin, ratioValue_mcnlo);
    hratio_mcnlo->SetBinContent(ibin, ratioValue_mcnlo);
    hratio_mcnlo->SetBinError  (ibin, ratioError_mcnlo);
   
    ratioErr->SetBinContent(ibin, 1.0);
    ratioErr->SetBinError  (ibin, uncertaintyError);
  }


  AxisFontsRatio(ratioErr->GetYaxis(), "y", "Theory / Data");
  AxisFontsRatio(ratioErr->GetXaxis(), "x", variable + " (GeV)");


  // Draw pad2
  //----------------------------------------------------------------------------
  canvas->cd();
  pad2->Draw();
  pad2->cd();
  
  ratioErr  ->Draw("e2");
  hratio_mad->Draw("e2,same");
  ratio_mad ->Draw("pz,same");

  ratioErr->GetYaxis()->SetRangeUser(0.35, 1.65);  

  pad2->Modified();
  
  DrawLatex(43, 0.2, 0.12, 15.0, 11, "Madgraph+Pythia normalized to #sigma_{NNLO}");
  

  // Draw pad3
  //----------------------------------------------------------------------------
  canvas->cd();
  pad3->Draw();
  pad3->cd();
  
  ratioErr    ->Draw("e2");
  hratio_mcnlo->Draw("e2,same");
  ratio_mcnlo ->Draw("pz,same");

  pad3->Modified();
  
  DrawLatex(43, 0.2, 0.12, 15.0, 11, "MC@NLO+Herwig normalized to #sigma_{NNLO}");
  
    
  // Draw pad4
  //----------------------------------------------------------------------------
  canvas->cd();
  pad4->Draw();
  pad4->cd();
  
  ratioErr  ->Draw("e2");
  hratio_pow->Draw("e2,same");
  ratio_pow ->Draw("pz,same");

  pad4->Modified();
  
  DrawLatex(43, 0.2, 0.44, 15.0, 11, "Powheg+Pythia normalized to #sigma_{NNLO}");


  // Save
  //----------------------------------------------------------------------------
  pad1->cd(); pad1->GetFrame()->DrawClone();
  pad2->cd(); pad2->GetFrame()->DrawClone();
  pad3->cd(); pad3->GetFrame()->DrawClone();
  pad4->cd(); pad4->GetFrame()->DrawClone();

  canvas->cd();

  canvas->SaveAs("Unfolded.pdf");
  canvas->SaveAs("Unfolded.png");
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis* axis, TString title)
{
  axis->SetLabelFont  (   43);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize  (   22);
  axis->SetNdivisions (  110);
  axis->SetTitle      (title);
  axis->SetTitleFont  (   43);
  axis->SetTitleOffset(  2.5);
  axis->SetTitleSize  (   22);
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFontsRatio(TAxis*  axis,
		    TString which,
		    TString title)
{
  if (which.Contains("x"))
    {
      axis->SetLabelFont  (   43);
      axis->SetLabelOffset(0.022);
      axis->SetLabelSize  (   22);
      axis->SetNdivisions (  510);
      axis->SetTitle      (title);
      axis->SetTitleFont  (   43);
      axis->SetTitleOffset( 5.85);
      axis->SetTitleSize  (   22);
    }
  else
    {
      axis->CenterTitle   ( true);
      axis->SetLabelFont  (   43);
      axis->SetLabelOffset(0.012);
      axis->SetLabelSize  (   15);
      axis->SetNdivisions (  505);
      axis->SetTitle      (title);
      axis->SetTitleFont  (   43);
      axis->SetTitleOffset(  3.4);
      axis->SetTitleSize  (   15);
    }
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
		Float_t y1,
		TH1F*   hist,
		TString label,
		TString option)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + 0.22,
				y1 + 0.05);

  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (0.045);

  legend->AddEntry (hist, label.Data(), option.Data());
  legend->Draw();
}


//------------------------------------------------------------------------------
// DrawLatex
//------------------------------------------------------------------------------
void DrawLatex(Font_t      tfont,
	       Double_t    x,
	       Double_t    y,
	       Double_t    tsize,
	       Short_t     align,
	       const char* text,
	       Bool_t      setndc = true)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC      (setndc);
  tl->SetTextAlign( align);
  tl->SetTextFont ( tfont);
  tl->SetTextSize ( tsize);

  tl->Draw("same");
}
