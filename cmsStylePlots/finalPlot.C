#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "StandardPlot.C"
#include "TSystem.h"
#include "CMS_lumi.C"


const Float_t _tsize = 0.033;
const Float_t _xoffset = 0.20;
const Float_t _yoffset = 0.05;


void finalPlot(bool drawRatio = 1, int differential = 0, int nsel = 0, int ReBin = 1, TString XTitle = "p_{T,max}^{l} (GeV)", TString units = "", TString plotName = "XSLeadingPt_AN.root", TString outputName = "WW_LeadingPt_final", bool isLogY = false,  double lumi = 19.5) {

  gInterpreter->ExecuteMacro("GoodStyle.C");

  TFile* file = new TFile(plotName.Data(), "read");
 
  TH1F* xsValue = (TH1F*) xsValue->Clone();
  TH1F* xsValue_Powheg =  (TH1F*) xsValue_Powheg->Clone();
  TH1F* xsValue_Madgraph = (TH1F*) xsValue_Madgraph->Clone();
  TH1F* xsValue_MCnlo = (TH1F*) xsValue_MCnlo->Clone();
  TH1F* systHisto = (TH1F*) systHisto->Clone();

 
  TCanvas* canvas ;
  TPad *pad1, *pad2;
 
   if (drawRatio) {
    canvas = new TCanvas("wwxs", "wwxs", 550, 1.2*600);
    
    pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetTopMargin   (0.05);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    
    pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.31); 
    pad2->SetTopMargin   (0.08);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
  
  }       else { 
    canvas = new TCanvas("wwxs", "wwxs", 550, 550);
    }

   if (drawRatio) pad1->cd();


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

   xsValue_Powheg->SetMarkerColor(kAzure-3);
   xsValue_Powheg->SetLineWidth(1);
   xsValue_Powheg->SetLineColor(kBlue+2);
   xsValue_Powheg->SetMarkerStyle(22);
   xsValue_Powheg->SetMarkerSize(1.2);
  

   hpowError->SetLineWidth(0);
   hpowError->SetMarkerSize (      0);  
   hpowError->SetFillColor  (kAzure-9);


  //-- Plot Madgraph

  TH1F *hmadError  = (TH1F*) xsValue_Madgraph->Clone();

  xsValue_Madgraph->SetMarkerColor(kPink-9);
  xsValue_Madgraph->SetLineWidth(1);
  xsValue_Madgraph->SetLineStyle(1);
  xsValue_Madgraph->SetMarkerStyle(21);
  xsValue_Madgraph->SetMarkerSize(1.0);

  hmadError->SetLineWidth(0);
  hmadError->SetMarkerSize (      0); 
  hmadError->SetFillColor  (kPink+1);


  //-- Plot MCNLO

  TH1F *hmcError  = (TH1F*) xsValue_MCnlo->Clone();

  xsValue_MCnlo->SetMarkerColor(kRed);
  xsValue_MCnlo->SetLineColor(kRed);
  xsValue_MCnlo->SetLineWidth(1);
  xsValue_MCnlo->SetLineStyle(1);
  xsValue_MCnlo->SetMarkerStyle(24);
  xsValue_MCnlo->SetMarkerSize(1.0);

  hmcError->SetLineWidth(0);
  hmcError->SetMarkerSize (      0); 
  hmcError->SetFillColor  (kOrange);
 



  //-- Plot Data

  xsValue->SetMarkerStyle(kFullCircle);
      
  if (differential == 0) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dp_{T,max}^{l}}");
  if (differential == 1) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dp_{T}(ll)}");
  if (differential == 2) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{dm_{#font[12]{ll}}}");
  if (differential == 3) AxisFonts (xsValue->GetYaxis(), "#frac{1}{#sigma} #frac{d#sigma}{d#Delta#phi_{ll}}");

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
  DrawLegend (0.65, 0.80, xsValue_Powheg,   "Powheg", "PL");
  DrawLegend (0.65, 0.75, hmadError,   "", "F");
  DrawLegend (0.65, 0.75, xsValue_Madgraph,   "Madgraph", "PL");  
  DrawLegend (0.65, 0.70, hmcError,   "", "F");
  DrawLegend (0.65, 0.70, xsValue_MCnlo,   "MCNLO", "LP");

  canvas->GetFrame()->DrawClone();



  // Draw text 
  //----------------------------------------------------------------------------
  TLatex * CMSLabel = new TLatex (0.18, 0.96, "#bf{CMS}");
  CMSLabel->SetNDC ();
  CMSLabel->SetTextAlign (10);
  CMSLabel->SetTextFont (42);
  CMSLabel->SetTextSize (_tsize);
  CMSLabel->Draw ("same") ;


  TLatex * _lumiLabel = new TLatex (0.95, 0.96, "19.4fb#lower[0.3]{^{-1}} (8 TeV)");
  _lumiLabel->SetNDC ();
  _lumiLabel->SetTextAlign (30);
  _lumiLabel->SetTextFont (42);
  _lumiLabel->SetTextSize (_tsize);
  _lumiLabel->Draw ("same") ;


  // Draw also ratio
  //----------------------------------------------------------------------------
  if (drawRatio) {

     pad2->cd();
	
    TH1F* ratio_pow       = xsValue_Powheg->Clone("ratio");
    TH1F* ratio_mad       = xsValue_Madgraph->Clone("ratio");
    TH1F* ratio_mcnlo     = xsValue_MCnlo->Clone("ratio");
    TH1F* hratio_pow       = xsValue_Powheg->Clone("ratio");
    TH1F* hratio_mad       = xsValue_Madgraph->Clone("ratio");
    TH1F* hratio_mcnlo     = xsValue_MCnlo->Clone("ratio");
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
    ratioErr  ->Draw("e2");
    
   
    
    ratio_mad      ->SetLineColor(kPink-9);
    ratio_mad      ->SetMarkerSize(1.0);
    ratio_mad      ->SetLineWidth(1);
    ratio_mad      ->SetMarkerStyle(21);
    hratio_mad     ->SetLineWidth(0);
    hratio_mad     ->SetMarkerSize (      0);  
    hratio_mad     ->SetFillColor  (kPink+1);
    hratio_mad     ->SetFillStyle  (1001);
    hratio_mad     ->Draw("e2,same");
    ratio_mad      ->Draw("e1p,same");
    
    
    ratio_mcnlo     ->SetLineColor(kRed);
    ratio_mcnlo     ->SetMarkerSize(1.0);
    ratio_mcnlo      ->SetLineWidth(1);
    ratio_mcnlo     ->SetMarkerStyle(24);
    hratio_mcnlo    ->SetLineWidth(0);
    hratio_mcnlo    ->SetMarkerSize (      0);  
    hratio_mcnlo    ->SetFillColor  (kOrange);
    hratio_mcnlo     ->SetFillStyle  (1001);
    hratio_mcnlo    ->Draw("e2,same");
    ratio_mcnlo     ->Draw("ep,same");

    ratio_pow      ->SetLineColor(kAzure-3);
    ratio_pow      ->SetMarkerSize(1.2);
    ratio_pow      ->SetLineWidth(1);
    ratio_pow      ->SetMarkerStyle(22);
    hratio_pow     ->SetLineWidth(0);
    hratio_pow     ->SetMarkerSize (      0);  
    hratio_pow     ->SetFillColor  (kAzure-9);
    hratio_pow     ->SetFillStyle  (1001);
    hratio_pow     ->Draw("e2,same");
    ratio_pow      ->Draw("e1p,same");
    
  
    
    ratioErr->SetFillColor  (kGray+2);
    ratioErr->SetFillStyle  (   3345);
    ratioErr->SetLineColor  (kGray+2);
    ratioErr->SetMarkerColor(kGray+2);
    ratioErr->SetMarkerSize (      0);

    ratioErr->Draw("sameaxis");

    ratioErr->GetYaxis()->SetRangeUser(0.4, 1.6);

   
    AxisFontsRatio (ratioErr->GetYaxis(), "y", "MC/data");
    AxisFontsRatio (ratioErr->GetXaxis(), "x", XTitle);

  }
      

}




//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts (TAxis* axis,
		TString title)
{
  axis->SetLabelFont ( 42);
  axis->SetLabelOffset (0.015);
  axis->SetLabelSize (0.040);
  axis->SetNdivisions ( 505);
  axis->SetTitleFont ( 42);
  axis->SetTitleOffset ( 1.9);
  axis->SetTitleSize (0.040);
  axis->SetTitle (title);
}

//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFontsRatio (TAxis* axis, TString which = "x", 
		TString title)
{
  axis->SetLabelFont ( 42);
  axis->SetLabelOffset (0.015);
  axis->SetLabelSize (0.090);
  axis->SetNdivisions ( 505);
  axis->SetTitleFont ( 42);
  axis->SetTitleOffset ( 1.5);
  if (which == "y" )  axis->SetTitleOffset ( 0.9);
  axis->SetTitleSize (0.090);
  axis->SetTitle (title);
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
