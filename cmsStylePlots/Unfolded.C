// .x finalPlot.C+(0,5,"E_{T}^{miss}","GeV","histo.root","histo_zhpresel_met",1,125,24.4, doDataMCRatio, signalZoom)

#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TColor.h"


// const Float_t _tsize = 0.033;
const Float_t _tsize = 0.049;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;




void 
Unfolded  (
     bool drawRatio = 1,
     int differential = 0,
     int nsel = 0,
     int ReBin = 1,
     TString XTitle = "p_{T,max}^{l} (GeV)",
     TString units = "", 
     TString plotName = "XSLeadingPt_AN.root", 
     TString outputName = "WW_LeadingPt_final",
     bool isLogY = false,
     double lumi = 19.5
 )
{  
 
 gInterpreter->ExecuteMacro("GoodStyle.C");
 
 std::cout << "reading " << plotName << std::endl ;
 TFile* file = new TFile(plotName, "read");
 
 //---- prepare the object that is making the plots
 //---- ---- ---- ---- ---- ---- ---- ---- ---- ----
 
 std::cout << "setting up the plot object " << std::endl ;

 
 TH1F* xsValue = (TH1F*) xsValue->Clone();
 TH1F* xsValue_Powheg =  (TH1F*) xsValue_Powheg->Clone();
 TH1F* xsValue_Madgraph = (TH1F*) xsValue_Madgraph->Clone();
 TH1F* xsValue_MCnlo = (TH1F*) xsValue_MCnlo->Clone();
 TH1F* systHisto = (TH1F*) systHisto->Clone();
 
 
 TCanvas* canvas ;
 TPad *pad1, *pad2, *pad3, *pad4;
 
 if (drawRatio) {
  canvas = new TCanvas("wwxs", "wwxs", 600, 850);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  canvas->SetHighLightColor(2);
  canvas->Range(0,0,1,1);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->SetLeftMargin(0.16);
  canvas->SetRightMargin(0.02);
  canvas->SetTopMargin(0.05);
  canvas->SetBottomMargin(0.13);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  
  
  
  pad1 = new TPad("pad1", "pad1",0.01,0.55,0.99,0.99);
  pad1->Draw();
  pad1->cd();
  pad1->Range(-147.2973,-5.811723,960.8108,2.535539);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetBorderSize(2);
  pad1->SetLogy();
  pad1->SetTickx(1);
  pad1->SetTicky(1);
  pad1->SetLeftMargin(0.16);
  pad1->SetBottomMargin(0);
  pad1->SetFrameFillStyle(0);
  pad1->SetFrameBorderMode(0);
  pad1->SetFrameFillStyle(0);
  pad1->SetFrameBorderMode(0);
   
 }
 else { 
  canvas = new TCanvas("wwxs", "wwxs", 550, 550);
 }
 
 if (drawRatio) pad1->cd();
 
 std::cout << " now plot " << std::endl;
 
 //Plot Data
 xsValue->SetLineWidth(1);
 xsValue->SetMarkerSize(1.0);
 
 int NBins = xsValue->GetNbinsX();
 
 for(int i=1; i <NBins; i++) {
  
  float err_stat = xsValue->GetBinError(i);
  float err_syst = systHisto->GetBinError(i);
  float err_total = sqrt(err_stat*err_stat + err_syst*err_syst);
  
  xsValue->SetBinError(i, err_total);
 }
 
 
 //-- Plot Powheg
 
 TH1F *hpowError  = (TH1F*) xsValue_Powheg->Clone();
 
 xsValue_Powheg->SetMarkerColor(kGreen+3);
 xsValue_Powheg->SetLineWidth(1);
 xsValue_Powheg->SetLineColor(kGreen+3);
 xsValue_Powheg->SetMarkerStyle(22);
 xsValue_Powheg->SetMarkerSize(1.2);
 
 
 hpowError->SetLineWidth(0);
 hpowError->SetMarkerSize (      0);  
 hpowError->SetFillColor  (kGreen-9);
 


//-- Plot Madgraph
 
 TH1F *hmadError  = (TH1F*) xsValue_Madgraph->Clone();
 
 xsValue_Madgraph->SetMarkerColor(kPink+10);
 xsValue_Madgraph->SetLineWidth(1);
 xsValue_Madgraph->SetLineStyle(1);
 xsValue_Madgraph->SetMarkerStyle(21);
 xsValue_Madgraph->SetMarkerSize(1.0);
 
 hmadError->SetLineWidth(0);
 hmadError->SetMarkerSize (      0); 
 hmadError->SetFillColor (kMagenta-10);

 
 //-- Plot MCNLO
 
 TH1F *hmcError  = (TH1F*) xsValue_MCnlo->Clone();
 
 xsValue_MCnlo->SetMarkerColor(kRed);
 xsValue_MCnlo->SetLineColor(kRed);
 xsValue_MCnlo->SetLineWidth(1);
 xsValue_MCnlo->SetLineStyle(1);
 xsValue_MCnlo->SetMarkerStyle(24);
 xsValue_MCnlo->SetMarkerSize(1.2);
 
 hmcError->SetLineWidth(0);
 hmcError->SetMarkerSize (      0); 
 hmcError->SetFillColor  (kOrange);
 
 
 
 
 //-- Plot Data
 
 xsValue->SetMarkerStyle(kFullCircle);
 
 if (differential == 0) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} d#sigma(WW#rightarrow#mu#nue#nu + < 1 jet)/dp_{T,max}^{l}");
 if (differential == 1) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} d#sigma(WW#rightarrow#mu#nue#nu + < 1 jet)/dp_{T}(ll)}");
 if (differential == 2) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} d#sigma(WW#rightarrow#mu#nue#nu + < 1 jet)/dm_{#font[12]{ll}}");
 if (differential == 3) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} d#sigma(WW#rightarrow#mu#nue#nu + < 1 jet)/d#Delta#phi_{ll}");

 
//  if (differential == 0) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dp_{T,max}^{l}}");
//  if (differential == 1) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dp_{T}(ll)}");
//  if (differential == 2) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dm_{#font[12]{ll}}}");
//  if (differential == 3) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{d#Delta#phi_{ll}}");

 AxisFonts (xsValue->GetXaxis(), XTitle);
 
 
 
 xsValue->Draw("p");
 hmadError->Draw("e2,same"); 
 xsValue_Madgraph->Draw("pe1,same");
 hmcError->Draw("e2,same");
 xsValue_MCnlo->Draw("pe1,same");
 hpowError->Draw("e2,same");
 xsValue_Powheg->Draw("pe1,same");
 //systHisto->Draw("e2, same");
 xsValue->Draw("pe1,same");
 
 // Legend
 //----------------------------------------------------------------------------
 
 DrawLegend (0.65, 0.85, xsValue, "Data", "P");
 DrawLegend (0.65, 0.80, hpowError,   "", "F");
 DrawLegend (0.65, 0.80, xsValue_Madgraph,   "Madgraph", "PL");  
 DrawLegend (0.65, 0.75, hmadError,   "", "F");
 DrawLegend (0.65, 0.75, xsValue_MCnlo,   "MCNLO", "LP");
 DrawLegend (0.65, 0.70, hmcError,   "", "F");
 DrawLegend (0.65, 0.70, xsValue_Powheg,   "Powheg", "PL");
 
 canvas->GetFrame()->DrawClone();
 
 
 
 // Draw text 
 //----------------------------------------------------------------------------
 TLatex *   tex = new TLatex(0.17,0.96,"CMS            
                             19.4 fb^{-1} (8TeV)");
 tex->SetNDC();
 tex->SetTextAlign(12);
 tex->SetTextFont(42);
 tex->SetTextSize(0.06);
 tex->SetLineWidth(2);
 tex->Draw();
 
//  TLatex * CMSLabel = new TLatex (0.18, 0.96, "#bf{CMS}");
//  CMSLabel->SetNDC ();
//  CMSLabel->SetTextAlign (10);
//  CMSLabel->SetTextFont (42);
//  CMSLabel->SetTextSize (_tsize);
//  CMSLabel->Draw ("same") ;
//  
//  
//  TLatex * _lumiLabel = new TLatex (0.95, 0.96, "19.4fb#lower[0.3]{^{-1}} (8 TeV)");
//  _lumiLabel->SetNDC ();
//  _lumiLabel->SetTextAlign (30);
//  _lumiLabel->SetTextFont (42);
//  _lumiLabel->SetTextSize (_tsize);
//  _lumiLabel->Draw ("same") ;
 
 
 // Draw also ratio
 //----------------------------------------------------------------------------
 if (drawRatio) {
  
  //---- prepare the distributions

  TH1F* ratio_pow       = xsValue_Powheg->Clone("ratio");
  TH1F* ratio_mad       = xsValue_Madgraph->Clone("ratio");
  TH1F* ratio_mcnlo     = xsValue_MCnlo->Clone("ratio");
  TH1F* hratio_pow      = xsValue_Powheg->Clone("ratio");
  TH1F* hratio_mad      = xsValue_Madgraph->Clone("ratio");
  TH1F* hratio_mcnlo    = xsValue_MCnlo->Clone("ratio");
  TH1F* ratioErr        = xsValue->Clone("ratio");
  
  
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
   
   Double_t ratioValue_pow           = (powValue > 0) ? powValue/dataValue : 0.0;
   Double_t ratioError_pow           = (powValue > 0) ? powError / dataValue : 0.0;
   
   Double_t ratioValue_mad           = (madValue > 0) ? madValue/dataValue : 0.0;
   Double_t ratioError_mad           = (madValue > 0) ? madError/dataValue : 0.0;
   
   Double_t ratioValue_mcnlo         = (mcnloValue > 0) ? mcnloValue/dataValue : 0.0;
   Double_t ratioError_mcnlo         = (mcnloValue > 0) ? mcnloError/dataValue : 0.0;
   
   Double_t uncertaintyError         = (dataValue > 0) ? dataError/dataValue : 0.0;
   
   
   //dataError/dataValue 
   ratio_pow->SetBinContent(ibin, ratioValue_pow);
   hratio_pow->SetBinContent(ibin, ratioValue_pow);
   hratio_pow->SetBinError  (ibin, ratioError_pow);
   
   ratio_mad->SetBinContent(ibin, ratioValue_mad);
   hratio_mad->SetBinContent(ibin, ratioValue_mad);
   hratio_mad->SetBinError  (ibin, ratioError_mad);
   
   ratio_mcnlo->SetBinContent(ibin, ratioValue_mcnlo);
   hratio_mcnlo->SetBinContent(ibin, ratioValue_mcnlo);
   hratio_mcnlo->SetBinError  (ibin, ratioError_mcnlo);
   
   ratioErr->SetBinContent(ibin, 1.0);
   ratioErr->SetBinError  (ibin, uncertaintyError);
  }
  
  ratioErr->SetTitle("");

  ratioErr->SetFillColor  (kGray+2);
  ratioErr->SetFillStyle  (   3004);
//   ratioErr->SetFillStyle  (   3345);
  ratioErr->SetLineColor  (kGray+2);
  ratioErr->SetMarkerColor(kGray+2);
  ratioErr->SetMarkerSize (      0);
  
  
  
  //---- now draw
  
  canvas->cd();
  pad2 = new TPad("pad2", "pad2",0.01,0.39,0.99,0.55);
  pad2->Draw();
  pad2->cd();
  pad2->Range(-147.2973,0.2,960.8108,1.8);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetBorderSize(2);
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  pad2->SetLeftMargin(0.16);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->SetFrameFillStyle(0);
  pad2->SetFrameBorderMode(0);
  pad2->SetFrameFillStyle(0);
  pad2->SetFrameBorderMode(0);
  
  
  ratioErr  ->Draw("e2");
  ratio_mad      ->SetLineColor(kPink+10);
  ratio_mad      ->SetMarkerSize(1.0);
  ratio_mad      ->SetLineWidth(1);
  ratio_mad      ->SetMarkerStyle(21);
  hratio_mad     ->SetLineWidth(0);
  hratio_mad     ->SetMarkerSize (      0);  
  hratio_mad     ->SetFillColor  (kMagenta-10);
  hratio_mad     ->SetFillStyle  (1001);
  hratio_mad     ->Draw("e2,same");
  ratio_mad      ->Draw("e1p,same");
  AxisFontsRatio (ratioErr->GetYaxis(), "y", "Theory/data");
  AxisFontsRatio (ratioErr->GetXaxis(), "x", XTitle);
  ratioErr->Draw("sameaxis");
  ratioErr->GetYaxis()->SetRangeUser(0.4, 1.6);  
  pad2->Modified();
  
  TLatex* tex_mad = new TLatex(0.2,0.1,"Madgraph+Pythia, normalized to #sigma_{NNLO}");
  tex_mad->SetNDC();
  tex_mad->SetTextAlign(12);
  tex_mad->SetTextFont(42);
  tex_mad->SetTextSize(0.1);
  tex_mad->SetLineWidth(2);
  tex_mad->Draw();
  
  
  
  canvas->cd();
  pad3 = new TPad("pad3", "pad3",0.01,0.23,0.99,0.39);
  pad3->Draw();
  pad3->cd();
  pad3->Range(-147.2973,0.2,960.8108,1.8);
  pad3->SetFillColor(0);
  pad3->SetBorderMode(0);
  pad3->SetBorderSize(2);
  pad3->SetTickx(1);
  pad3->SetTicky(1);
  pad3->SetLeftMargin(0.16);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0);
  pad3->SetFrameFillStyle(0);
  pad3->SetFrameBorderMode(0);
  pad3->SetFrameFillStyle(0);
  pad3->SetFrameBorderMode(0);
  
  ratioErr  ->Draw("e2");
  ratio_mcnlo     ->SetLineColor(kRed);
  ratio_mcnlo     ->SetMarkerSize(1.2);
  ratio_mcnlo      ->SetLineWidth(1);
  ratio_mcnlo     ->SetMarkerStyle(24);
  hratio_mcnlo    ->SetLineWidth(0);
  hratio_mcnlo    ->SetMarkerSize (      0);  
  hratio_mcnlo    ->SetFillColor  (kOrange);
  hratio_mcnlo     ->SetFillStyle  (1001);
  hratio_mcnlo    ->Draw("e2,same");
  ratio_mcnlo     ->Draw("ep,same");
//   AxisFontsRatio (ratioErr->GetYaxis(), "y", "Theory/data");
//   AxisFontsRatio (ratioErr->GetXaxis(), "x", XTitle);
  ratioErr->Draw("sameaxis");
  ratioErr->GetYaxis()->SetRangeUser(0.4, 1.6);  
  pad3->Modified();
  
  TLatex* tex_mcnlo = new TLatex(0.2,0.1,"MC@NLO+Herwig, normalized to #sigma_{NNLO}");
  tex_mcnlo->SetNDC();
  tex_mcnlo->SetTextAlign(12);
  tex_mcnlo->SetTextFont(42);
  tex_mcnlo->SetTextSize(0.1);
  tex_mcnlo->SetLineWidth(2);
  tex_mcnlo->Draw();
  
  
  
  
  canvas->cd();
  pad4 = new TPad("pad4", "pad4",0.01,0.01,0.99,0.23);
  pad4->Draw();
  pad4->cd();
  pad4->Range(-147.2973,-0.4857143,960.8108,1.8);
  pad4->SetFillColor(0);
  pad4->SetBorderMode(0);
  pad4->SetBorderSize(2);
  pad4->SetTickx(1);
  pad4->SetTicky(1);
  pad4->SetLeftMargin(0.16);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.3);
  pad4->SetFrameFillStyle(0);
  pad4->SetFrameBorderMode(0);
  pad4->SetFrameFillStyle(0);
  pad4->SetFrameBorderMode(0);
  
  ratioErr  ->Draw("e2");
  ratio_pow      ->SetLineColor(kGreen+3);
  ratio_pow      ->SetMarkerSize(1.2);
  ratio_pow      ->SetLineWidth(1);
  ratio_pow      ->SetMarkerStyle(22);
  hratio_pow     ->SetLineWidth(0);
  hratio_pow     ->SetMarkerSize (      0);  
  hratio_pow     ->SetFillColor  (kGreen-9);
  hratio_pow     ->SetFillStyle  (1001);
  hratio_pow     ->Draw("e2,same");
  ratio_pow      ->Draw("e1p,same");
//   AxisFontsRatio (ratioErr->GetYaxis(), "y", "Theory/data", 1);
//   AxisFontsRatio (ratioErr->GetXaxis(), "x", XTitle);
  ratioErr->Draw("sameaxis");
  ratioErr->GetYaxis()->SetRangeUser(0.4, 1.6);  
  pad4->Modified();
  
  TLatex* tex_pow = new TLatex(0.2,0.39,"Powheg+Pythia, normalized to #sigma_{NNLO}");
  tex_pow->SetNDC();
  tex_pow->SetTextAlign(12);
  tex_pow->SetTextFont(42);
  tex_pow->SetTextSize(0.1);
  tex_pow->SetLineWidth(2);
  tex_pow->Draw();
  
  canvas->cd();
  
  
 }
 
 
}




//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts (TAxis* axis,
                TString title)
{
 axis->SetNdivisions(110);
 axis->SetLabelFont(42);
 axis->SetLabelOffset(0.007);
 axis->SetLabelSize(0.12/2);
 axis->SetTitleSize(0.13/2);
 axis->SetTitleFont(42);
 axis->SetTitleOffset ( 1.0);
 axis->SetTitle (title);
 
//  axis->SetLabelFont ( 42);
//  axis->SetLabelOffset (0.015);
//  axis->SetLabelSize (0.040);
//  axis->SetNdivisions ( 505);
//  axis->SetTitleFont ( 42);
//  axis->SetTitleOffset ( 1.9);
//  axis->SetTitleSize (0.040);
//  axis->SetTitle (title);
}

//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFontsRatio (TAxis* axis, TString which = "x", 
                     TString title, int small = 0)
{
 
 if (which == "x" ) {
  axis->SetNdivisions(110);
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(0.12);
  axis->SetTitleSize(0.13);
  axis->SetTitleFont(42);
  axis->SetTitleOffset ( 1.0);
  axis->SetTitle (title);
 }
 
 if (which == "y" ) {
  axis->SetTitle("Theory/Data");
  axis->CenterTitle(true);
  axis->SetNdivisions(505);
  axis->SetLabelFont(42);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize(0.12);
  axis->SetTitleSize(0.12);
  axis->SetTitleOffset(0.56);
  axis->SetTitleFont(42);
  axis->SetTitle (title);
 }
 
}



//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend (Float_t x1,
                 Float_t y1,
                 TH1F* hist,
                 TString label,
                 TString option)
{
 TLegend* legend = new TLegend (x1,y1,x1 + _xoffset,y1 + _yoffset);
 legend->SetBorderSize (0) ;
 legend->SetFillColor (0) ;
 legend->SetTextAlign (12) ;
 legend->SetTextFont (40) ;
 legend->SetTextSize (_tsize) ;
 legend->AddEntry (hist, label.Data (), option.Data ());
 legend->Draw (); 
 
}

