#include "TFile.h"
#include "TH1F.h"
#include "TInterpreter.h"

#include "../DY.C"
#include "../Top.C"

#include "TSVDUnfold.h"
#include "TRandom.h"
#include "TH1D.h"
#include <iostream>
#include "TVectorD.h"
#include "TMatrixD.h"


using std::cout;
using std::endl;



gSystem->Load("libRooUnfold");

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"



const Double_t nlo8tev       = 1;//57.25;  // [pb]
const Double_t nlo8tevPlus   = 4.2 * nlo8tev / 1e2;
const Double_t nlo8tevMinus  = 2.2 * nlo8tev / 1e2;

const Double_t CMS8tev       = 1;//56.87;//63.53;  // [pb]
const Double_t CMS8tev_err   = 5.55;   // [pb]

const Double_t CMSFidutial8tev       = 57.25;  // [pb]
const Double_t BR_WW_to_lnln = (3 * 0.108) * (3 * 0.108);


const Double_t ggWW_xs = nlo8tev * BR_WW_to_lnln * 0.03;
const Double_t qqWW_xs = nlo8tev * BR_WW_to_lnln * 0.97;

const Double_t NTotalggWW =  109987;
const Double_t NTotalqqWW = 1933235;


const Double_t NFiducialggWW =   67977;  // pt > 20, |eta| < 2.5
const Double_t NFiducialqqWW = 1050807;  // pt > 20, |eta| < 2.5


const UInt_t Nsyst = 11;



const UInt_t nchannels = 7;

TString channel [nchannels] = {"All", "SF", "OF", "MuMu",   "EE", "EMu",  "MuE" };

enum {All, SF, OF, MuMu, EE, EMu, MuE};


// Initialize the systematics 

const UInt_t Nsyst =2;

Double_t Syst_data [Nsyst] = {1.3, 2.9}; // in % systematics on the efficiency from { PDF 1.3, Luminosity 2.9}


// Wjets --> normalization 28%
// Vg --> scale 30% 
// VV --> 4% UEPS + 5% QCD
// H125 --> 8% UEPS + 20% QCD
// TTbar nom --> 10 %
// DY error --> 30 %

Double_t Syst_back [6] = {0.28, 0.30, 0.06, 0.22, 0.10, 0.30};

//------------------------------------------------------------------------------
// XS
//------------------------------------------------------------------------------
void XSDiffUnfold(Double_t  luminosity = 19365,
		  Int_t     njet = 0,
		  TString   channel = channel[OF],
		  TString   directory = "rootfiles.2014Selection.4L.NewMC.NewLumi",
		  //TString   directory = "rootfiles.Systematics.FromDaniel/nominals",
		  TString   systematics = "141030_Systematics_v1.0",
		  Bool_t    useDataDriven = false,
		  Int_t     printLevel = 0,
		  Bool_t    fiducial = true, 
		  Int_t     differential = 3,
		  Bool_t    drawTheXS = true,
		  Bool_t    drawRatio = true,
		  Int_t     verbose = 0,
		  Bool_t    useOvf = 1
	   )

{
 

  /// cancell - out some systematics if we compute the normalized distributions. 

  if ( fiducial ) { 
    Syst_data [0] = 0.0;  Syst_data [1] = 0.0;
  }


  TH1::SetDefaultSumw2();
 
  gInterpreter->ExecuteMacro("../draw/ChargeRatioStyle.C");
  gStyle->SetPalette(1);

 TString dyChannel = "SF";

  if (channel.Contains("EE"))   dyChannel = "EE";
  if (channel.Contains("MuMu")) dyChannel = "MuMu";

  Double_t dyScaleFactor       = -999;
  Double_t dyScaleFactorForTop = 1;  Double_t topScaleFactor      = -999;

 
  Double_t NTop_incl [3];

  Double_t NDY_incl [3];




  // Needed to always have the SF DY scale factor in the top estimation
  //----------------------------------------------------------------------------
 
  if ( useDataDriven ) { 

    Top(NTop_incl[0],
	NTop_incl[1],
	NTop_incl[2],
	topScaleFactor,
	dyScaleFactorForTop,
	njet,
	"OF",  // "All" --> channel (individual topScaleFactor)
	directory,
	useDataDriven,
	printLevel);
  }
  

 

  //----------------------------------------
  // differential == 0 --> leading lepton pt 
  //----------------------------------------


  printf("\n");
  printf("   nlo8tev xs = %.2f pb", nlo8tev );
  printf("   CMS8tev xs = %.2f pb", CMS8tev );
  printf("   CMSFidutial8tev xs = %.2f pb", CMSFidutial8tev );

  
 
  // Input files
  //----------------------------------------------------------------------------
  TString path = Form("../%s/%djet/%s/", directory.Data(), njet, channel.Data());
  
  TFile* inputWW       = new TFile(path + "WW_pow_nnll_smear.root");  
  TFile* inputggWW     = new TFile(path + "ggWWto2L.root");  
  TFile* inputqqWW     = new TFile(path + "WWTo2L2Nu.root");  
  TFile* inputqqWW_pow = new TFile(path + "WWTo2L2Nu_pow_nnll_smear.root");  
  TFile* inputqqWW_mc  = new TFile(path + "WWTo2L2Nu_mcnlo.root");  
  TFile* inputTT       = new TFile(path + "TTbar.root");  
  TFile* inputTW       = new TFile(path + "TW.root");  
  TFile* inputWj       = new TFile(path + "WJetsFakes_Total.root");
  TFile* inputWZ       = new TFile(path + "WZ.root");  
  TFile* inputZZ       = new TFile(path + "ZZ.root");  
  TFile* inputDY       = new TFile(path + "DY.root");  
  TFile* inputDYtautau = new TFile(path + "DYtautau.root");  
  TFile* inputWg       = new TFile(path + "Wgamma.root");  
  TFile* inputH125     = new TFile(path + "HWW125.root");  
  TFile* inputZgamma   = new TFile(path + "Zgamma.root");  
  TFile* inputData     = new TFile(path + "DataRun2012_Total.root");  

  // TFile* inputWW_GEN_pow     = new TFile("../rootfiles.Systematics.FromDaniel/jetEnergyScale_up/0jet/WWGEN/WW_GEN_0jet_pow_full.root");

  //TFile* inputWW_GEN_pow     = new TFile("files/WW_GEN_0jet_pow_gg_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root"); //for doing unfolding
  

  //TFile* inputWW_GEN_pow     = new TFile("files/final/WW_GEN_0jet_gg_pow_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root");


  //TFile* inputWW_GEN_mad     = new TFile("files/WW_GEN_0jet_mad_gg_full_JetGenVeto_Eff_NNLOXsec_NewLumi.root"); //for doing unfolding
  //TFile* inputWW_GEN_mcnlo   = new TFile("files/WW_GEN_0jet_mcnlo_gg_full_JetGenVeto_Eff_NNLOXsec_NewLumi.root"); //for doing unfolding
  //TFile* inputWW_GEN_mad     = new TFile("files/WW_GEN_0jet_mad_gg_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root"); //for doing unfolding
  //TFile* inputWW_GEN_mcnlo   = new TFile("files/WW_GEN_0jet_mcnlo_gg_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root"); //for doing unfolding



  TFile* inputWW_GEN_pow     = new TFile("files/final/WW_GEN_0jet_gg_pow_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root");
  TFile* inputWW_GEN_mad     = new TFile("files/final/WW_GEN_0jet_gg_mad_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root");
  TFile* inputWW_GEN_mcnlo   = new TFile("files/final/WW_GEN_0jet_gg_mcnlo_full_NNLL_JetGenVeto_Eff_NNLOXsec_NewLumi.root");



  //----------------------------------------------------------------------------
  //
  // Estimate WW differential cross-section                   
  //
  //
  //----------------------------------------------------------------------------

 
  // [0] pt Lepton 1

  TString distribution = ""; 
  TString genDistribution = ""; 
  TString recoDistribution = ""; 
  TString response = "";
  int kreg = 1; 

  TFile* syst;
  TFile* syst_pow;
  TFile* syst_mc;

  if ( differential == 0 ) {
    distribution = "hPtLepton1WWLevel_Diff";
    response = "hPtLepton1WWLevel_RECO_hPtLepton1WWLevel_GEN";
    genDistribution = "hPtLepton1_GEN";
    kreg = 3;
    syst = new TFile(systematics+"/syst_PtLepton1.root"); 
    syst_pow = new TFile("141210_MCUncertainties_v0.9/mcuncert_pow_PtLepton1.root");
    syst_mc = new TFile("141210_MCUncertainties_v0.9/mcuncert_mcnlo_PtLepton1.root");
  }

  if ( differential == 1 ) {
    distribution = "hDileptonWWLevel_Diff";
    response = "hDileptonWWLevel_RECO_hDileptonWWLevel_GEN";
    genDistribution = "hDilepton_GEN";
    kreg = 3;
    syst = new TFile(systematics+"/syst_Dilepton.root"); 
    syst_pow = new TFile("141210_MCUncertainties_v0.9/mcuncert_pow_Dilepton.root");
    syst_mc = new TFile("141210_MCUncertainties_v0.9/mcuncert_mcnlo_Dilepton.root");
  }
  
  if ( differential == 2 ) {
    distribution = "hMinvWWLevel_Diff";
    response = "hmllWWLevel_RECO_hmllWWLevel_GEN";
    genDistribution = "hmll_GEN";
    kreg = 4;  
    syst = new TFile(systematics+"/syst_mll.root");
    syst_pow = new TFile("141210_MCUncertainties_v0.9/mcuncert_pow_mll.root");
    syst_mc = new TFile("141210_MCUncertainties_v0.9/mcuncert_mcnlo_mll.root");
}

  if ( differential == 3 ) {
    distribution = "hDeltaPhiWWLevel_Diff";
    response = "hdphiWWLevel_RECO_hdphiWWLevel_GEN";
    genDistribution = "hdphi_GEN";
    kreg = 3; 
    syst = new TFile(systematics+"/syst_dphi.root");
    syst_pow = new TFile("141210_MCUncertainties_v0.9/mcuncert_pow_dphi.root");
    syst_mc = new TFile("141210_MCUncertainties_v0.9/mcuncert_mcnlo_dphi.root");
  }

  if ( differential == 6 ) {
    distribution = "hWTopTagging";
    response = "hInclusiveWWLevel_RECO_hInclusiveWWLevel_GEN;";
    genDistribution = "hInclusive_GEN";
    kreg = 3; 
    syst = new TFile(systematics+"/syst_Inclusive.root");
    syst_pow = new TFile("141210_MCUncertainties_v0.9/mcuncert_pow_Inclusive.root");
    syst_mc = new TFile("141210_MCUncertainties_v0.9/mcuncert_mcnlo_Inclusive.root");
    //syst = new TFile(systematics+"/syst_Inclusive_leptonEfficiency.root");
    //syst = new TFile(systematics+"/syst_Inclusive_JER.root");
    //syst = new TFile(systematics+"/syst_Inclusive_jetEnergyScale.root");
    //syst = new TFile(systematics+"/syst_Inclusive_metResolution.root");
    //syst = new TFile(systematics+"/syst_Inclusive_electroScale.root");
    //syst = new TFile(systematics+"/syst_Inclusive_muonScale.root");
    //syst = new TFile(systematics+"/syst_Inclusive_electronResolution.root");

}
 

  TKey *key = inputWW->FindKey(distribution);
  if (key ==0){
    cout << "  " << endl;
    cout << "!!Histogram does not exist!!" << endl;
    cout << "  " << endl;
  }

  TH1F* hNWW, *hNggWW, *hNqqWW_mad, *hNqqWW_pow, *hNqqWW_mcnlo, *hNqqWW_mc, *hNTT, *hNTW, *hNZgamma;
  TH1F* hNWj, *hNWZ, *hNZZ, *hNDY, *hNDYtautau;
  TH1F* hNWg, *hNH125, *hNData;  


  TH1F* hNDYOF, * hNDYSF;

  Int_t NBins;

  
  
  // SIGNAL
  hNWW       = (TH1F*) inputWW      ->Get(distribution);
  NBins = hNWW->GetNbinsX();

  hNggWW         = (TH1F*) inputggWW             ->Get(distribution);
  hNqqWW_mad     = (TH1F*) inputWW_GEN_mad       ->Get(genDistribution);
  hNqqWW_pow     = (TH1F*) inputWW_GEN_pow       ->Get(genDistribution);
  hNqqWW_mcnlo   = (TH1F*) inputWW_GEN_mcnlo     ->Get(genDistribution);
  hNqqWW_mc      = (TH1F*) inputqqWW_mc          ->Get(distribution);

 
  // BACKGROUNDS
  hNTT       = (TH1F*) inputTT      ->Get(distribution);
  hNTW       = (TH1F*) inputTW      ->Get(distribution);
  hNWj       = (TH1F*) inputWj      ->Get(distribution);
  hNWZ       = (TH1F*) inputWZ      ->Get(distribution);
  hNZZ       = (TH1F*) inputZZ      ->Get(distribution);
  hNDY       = (TH1F*) inputDY      ->Get(distribution);
  hNDYtautau = (TH1F*) inputDYtautau->Get(distribution);
  hNWg       = (TH1F*) inputWg      ->Get(distribution);
  hNH125     = (TH1F*) inputH125    ->Get(distribution);
  hNZgamma   = (TH1F*) inputZgamma  ->Get(distribution);
  hNData     = (TH1F*) inputData    ->Get(distribution);
  
 
 
  // Define data - background distributions
  //----------------------------------------------------------------------------
  
  
  TH1F *WWdata = (TH1F*) hNData->Clone();

  WWdata->Add(hNTT,-1.0);
  WWdata->Add(hNTW,-1.0);
  WWdata->Add(hNWj,-1.0);
  WWdata->Add(hNWZ,-1.0);
  WWdata->Add(hNZZ,-1.0);
  WWdata->Add(hNDY,-1.0);
  WWdata->Add(hNDYtautau,-1.0);
  WWdata->Add(hNH125,-1.0);
  WWdata->Add(hNZgamma ,-1.0);
  WWdata->Add(hNWg,-1.0);
  

  WWdata->SetName(distribution);



  
  // DoUnfold of data- background
  //----------------------------------------------------------------------------

  const RooUnfoldResponse *responseFinal =  inputWW_GEN_pow->Get(response);


  /// To be included also the underflow bin
  if (useOvf ) responseFinal->UseOverflow();    
 

 //initialize SVD unfolding
 RooUnfoldSvd   unfold(responseFinal, WWdata,1);   
 unfold.SetRegParm(kreg);  

  //initialize Bayesian unfolding
  //RooUnfoldBayes  unfold(responseFinal, WWdata,10);

  //initialize ML unfolding
  //RooUnfoldInvert unfold(responseFinal, WWdata);
 
 
  // ----> Choose error treatment::  
  //       0: Errors are the square root of the bin content
  //       1: Errors from the diagonals of the covariance matrix given by the unfolding 
  //          (variance values) 
  //       2: Errors from the covariance matrix given by the unfolding
  //       3: Errors from the covariance matrix from the variation of the results in toy MC tests
  //==============================================================================

  //setup error type
  unfold.SetNToys(1000);
  RooUnfold::ErrorTreatment *doerror = 2;//RooUnfold::kCovariance;


  if (verbose > 0 ) { 

    RooUnfoldParms *parms= new RooUnfoldParms(&unfold,(RooUnfold::ErrorTreatment)doerror, hNqqWW_pow);
  
    parms->SetMinParm(1);
    parms->SetMaxParm(NBins+1);
    //    parms->SetStepSizeParm(1);

    TProfile *hParmChi2= parms->GetChi2();
    TProfile *hParmErr= parms->GetRMSError();  // Mean values of errors in each bin. 
    TProfile *hParmRes= parms->GetMeanResiduals();
    TH1      *hParmRms= parms->GetRMSResiduals();


   

    // ---->  create matrix to store covariance matrix
    //==============================================================================

    if (useOvf) NBins = NBins+2;

    TMatrixD m_covMat (NBins,NBins);
    m_covMat = unfold.Ereco((RooUnfold::ErrorTreatment)doerror); 

    
    TMatrixD m_errMat(NBins,NBins);

    for (Int_t i=0; i<NBins; i++) { 
      for (Int_t j=0; j<NBins; j++) { 
	
	m_errMat(i,j)= m_covMat(i,j)>=0 ? sqrt(m_covMat(i,j)) : -sqrt(-m_covMat(i,j));
      }
    }

    unfold.PrintTable (cout, hNqqWW_pow, (RooUnfold::ErrorTreatment)doerror);

    PrintMatrix(m_errMat,"","covariance matrix",10);

    TH2D *hCorr = CorrelationHist ( m_covMat,    
				    "corr", "Unfolded correlation matrix",
				    responseFinal->Hresponse()->GetYaxis()->GetXmin(),
				    responseFinal->Hresponse()->GetYaxis()->GetXmax(),
				    useOvf);
    

    //if (useOvf) NBins = NBins-2;

    
    // ---->   Returns d vector (for choosing appropriate regularisation)
    //==============================================================================
  
    TCanvas * param = new TCanvas("param", "param", 550, 550);
    param->cd();

    TSVDUnfold *myTSVD = (TSVDUnfold*) unfold.Impl();

    TH1D *svVector = myTSVD->GetSV();
    TH1D *dVector = myTSVD->GetD();
     
    TH1D *dzVector = (TH1D*) dVector->Clone();

    double tau = svVector->GetBinContent(kreg+1);

    for (int i=1; i < NBins; i++) {

    double Si = svVector->GetBinContent(i);  
    double scale = (Si*Si)/((Si*Si)+(tau*tau));
    double di = dVector->GetBinContent(i);
  
    dzVector->SetBinContent(i, di*scale);
    
  }


    dVector->GetYaxis()->SetTitle("d_{i}");
    dVector->GetXaxis()->SetTitle("i");
    

    if (differential == 0) dVector->SetTitle("SVD unfolding d_{i} for p_{T}^{max}");
    if (differential == 1) dVector->SetTitle("SVD unfolding d_{i} for dp_{T}(ll)");
    if (differential == 2) dVector->SetTitle("SVD unfolding d_{i} for m_{#font[12]{ll}}");
    if (differential == 3) dVector->SetTitle("SVD unfolding d_{i} for #Delta#phi_{ll}");

    dVector->SetMarkerStyle(4);
    dVector->DrawCopy("hist");  
    
    dzVector->SetLineStyle(2);
    dzVector->SetLineWidth(2);
    dzVector->SetLineColor(kRed);
    dzVector->Draw("histSAME");  
    
    TLine().DrawLine(dVector->GetBinLowEdge(1), 1.0, dVector->GetBinLowEdge(NBins+1), 1.0);  // draw a line at y=0;


    //For getting the response matrix

    TCanvas * resp = new TCanvas("resp", "response", 550, 550);
    resp->cd();

    TH2 *resM = responseFinal->Hresponse(); 

    resM->Draw("COLZ");



  } // end verbose



  // ---->  get unfolded histogram 
  //==============================================================================
  TH1F* h_dataReco_unfolded = (TH1F*) unfold.Hreco((RooUnfold::ErrorTreatment)doerror); 



  if (verbose > 0 ) { 


    //For plotting error bars from unfolding.
    
    TCanvas * errors = new TCanvas("errors", "errors", 550, 550);
    errors->cd();

    TH1F* h1 = (TH1F*) unfold.Hreco((RooUnfold::ErrorTreatment)1); 
    TH1F* h2 = (TH1F*) unfold.Hreco((RooUnfold::ErrorTreatment)2);
    TH1F* h3 = (TH1F*) unfold.Hreco((RooUnfold::ErrorTreatment)3);
    
   
    h1 ->SetFillColor(kRed);
    h1->SetFillStyle(3004);
    h1->Draw("E2");
    h2->Draw("E1same");
    h3 ->SetFillColor(kGreen);
    h3->SetFillStyle(3006);
    h3->Draw("E2same");



    TCanvas * residual = new TCanvas("residual", "residual", 550, 550);
    residual->cd();

    // ---->  Calculate pulls and residuals
    //==============================================================================
  
    TH1F *hRes = (TH1F*) h_dataReco_unfolded->Clone("res");
    hRes  ->Reset();
    hRes  ->SetTitle ("Residuals");

    TH1F *hRes_up1s   = (TH1F*) hNqqWW_pow->Clone("res");
    TH1F *hRes_up2s   = (TH1F*) hNqqWW_pow->Clone("res");
    TH1F *hRes_down1s = (TH1F*) hNqqWW_pow->Clone("res");
    TH1F *hRes_down2s = (TH1F*) hNqqWW_pow->Clone("res");

    TH1F *hPulls= (TH1F*) h_dataReco_unfolded->Clone("pulls");
    hPulls->Reset();
    hPulls->SetTitle ("Pulls");

    TH1F *h_mcTruth= hNqqWW_pow->Clone("res");

    for (Int_t i= 1; i<=NBins; i++) {
      if ((h_dataReco_unfolded->GetBinContent(i)!=0.0 || (doerror && h_dataReco_unfolded->GetBinError(i)>0.0)) &&
	  ( h_mcTruth->GetBinContent(i)!=0.0 || (doerror &&  h_mcTruth->GetBinError(i)>0.0))) {
      
	Double_t res= h_dataReco_unfolded->GetBinContent(i) -  h_mcTruth->GetBinContent(i);
	Double_t err= h_dataReco_unfolded->GetBinError  (i);
      
	Double_t Err = h_mcTruth->GetBinError  (i);
	Double_t Bias = Err; 


	hRes->SetBinContent (i, res);
	hRes->SetBinError   (i, err);

	hRes_up1s->SetBinContent(i, Bias);
      	hRes_up2s->SetBinContent(i, 2*Bias);
	hRes_down1s->SetBinContent(i, -Bias);
      	hRes_down2s->SetBinContent(i, -2*Bias);


	if (err>0.0) {
	  hPulls->SetBinContent (i, res/err);
	  hPulls->SetBinError   (i, 1.0);
	}
      }
    }
    
     hRes->Draw();
    
    hRes_up2s->Draw("histsame");
    hRes_down2s->Draw("histsame");
    hRes_up2s->SetFillColor(kYellow);
    hRes_up2s->SetFillStyle(3004);
    hRes_up2s->SetLineColor(kYellow);
    hRes_down2s->SetFillColor(kYellow);
    hRes_down2s->SetFillStyle(3004);
    hRes_down2s->SetLineColor(kYellow);


    hRes_up1s->Draw("histsame");
    hRes_down1s->Draw("histsame");
    hRes_up1s->SetFillColor(kGreen);
    hRes_up1s->SetFillStyle(3002);
    hRes_up1s->SetLineColor(kGreen);
    hRes_down1s->SetFillColor(kGreen);
    hRes_down1s->SetFillStyle(3002);
    hRes_down1s->SetFillStyle(3002);
    hRes_down1s->SetLineColor(kGreen);
  
    hRes->SetMarkerStyle(kFullDotLarge);
    hRes->Draw("same");

    TLine().DrawLine(hRes->GetBinLowEdge(1), 0.0, hRes->GetBinLowEdge(NBins+1), 0.0);  // draw a line at y=0;
    
}


  // ---->  get yields 
  //==============================================================================

  const unsigned int bin = NBins; 

  Double_t NData [bin][4];
  Double_t NqqWW_pow [bin][4];
  Double_t NqqWW_mad [bin][4];
  Double_t NqqWW_mcnlo [bin][4];

  Double_t NWj [bin][4];
  Double_t NWg [bin][4];
  Double_t NWZ [bin][4];
  Double_t NZZ [bin][4];
  Double_t NH125 [bin][4];
  Double_t NTT [bin][4];
  Double_t NTW [bin][4];
  Double_t NDY [bin][4];


  TH1F* systUp;  TH1F* systDown;
  TH1F* systUp_pow;  TH1F* systDown_pow;
  TH1F* systUp_mc;  TH1F* systDown_mc;

  float backgError [bin]; 

  systUp    = (TH1F*) syst   ->Get("up");
  systDown  = (TH1F*) syst   ->Get("down");

  systUp_pow    = (TH1F*) syst_pow   ->Get("up");
  systDown_pow  = (TH1F*) syst_pow   ->Get("down");

  systUp_mc    = (TH1F*) syst_mc   ->Get("up");
  systDown_mc  = (TH1F*) syst_mc   ->Get("down");


  for (int ib=0; ib < bin; ib++) { 

    NData[ib][0] = h_dataReco_unfolded    ->GetBinContent(ib+1); NData[ib][1] = h_dataReco_unfolded    ->GetBinWidth(ib+1);
    NData[ib][2] = h_dataReco_unfolded    ->GetBinError(ib+1);   NData[ib][3] = 0.0; 

    NqqWW_pow [ib][0] = hNqqWW_pow ->GetBinContent(ib+1); NqqWW_pow [ib][1] = hNqqWW_pow ->GetBinWidth(ib+1);
    NqqWW_pow [ib][2] = hNqqWW_pow ->GetBinError(ib+1); NqqWW_pow [ib][3] = 0.0;

    NqqWW_mad [ib][0] = hNqqWW_mad ->GetBinContent(ib+1); NqqWW_mad [ib][1] = hNqqWW_mad ->GetBinWidth(ib+1);
    NqqWW_mad [ib][2] = hNqqWW_mad ->GetBinError(ib+1); NqqWW_mad [ib][3] = 0.0;

    NqqWW_mcnlo [ib][0] = hNqqWW_mcnlo ->GetBinContent(ib+1); NqqWW_mcnlo [ib][1] = hNqqWW_mcnlo ->GetBinWidth(ib+1);
    NqqWW_mcnlo [ib][2] = hNqqWW_mcnlo ->GetBinError(ib+1); NqqWW_mcnlo [ib][3] = 0.0;


    // backgrounds
    NWj [ib][0] = hNWj->GetBinContent(ib+1); NWj [ib][1] = hNWj ->GetBinWidth(ib+1);
    NWj [ib][2] =   hNWj ->GetBinError(ib+1); NWj [ib][3] = NWj [ib][0] *Syst_back[0];

    NWg [ib][0] = hNWg->GetBinContent(ib+1); NWg [ib][1] = hNWg ->GetBinWidth(ib+1);
    NWg [ib][2] =   hNWg ->GetBinError(ib+1); NWg [ib][3] = NWg [ib][0] *Syst_back[1];

    NWZ [ib][0] = hNWZ->GetBinContent(ib+1); NWZ [ib][1] = hNWZ ->GetBinWidth(ib+1);
    NWZ [ib][2] =   hNWZ ->GetBinError(ib+1); NWZ [ib][3] = NWZ [ib][0] *Syst_back[2];
  
    NZZ [ib][0] = hNZZ->GetBinContent(ib+1); NZZ [ib][1] = hNZZ ->GetBinWidth(ib+1);
    NZZ [ib][2] =   hNZZ ->GetBinError(ib+1); NZZ [ib][3] = NZZ [ib][0] *Syst_back[2];

    NH125 [ib][0] = hNH125->GetBinContent(ib+1); NH125 [ib][1] = hNH125 ->GetBinWidth(ib+1);
    NH125 [ib][2] =   hNH125 ->GetBinError(ib+1); NH125 [ib][3] = NH125 [ib][0] *Syst_back[3];

    NTT [ib][0] = hNTT->GetBinContent(ib+1); NTT [ib][1] = hNTT ->GetBinWidth(ib+1);
    NTT [ib][2] =   hNTT ->GetBinError(ib+1); NTT [ib][3] = NTT [ib][0] *Syst_back[4];

    NTW [ib][0] = hNTW->GetBinContent(ib+1); NTW [ib][1] = hNTW ->GetBinWidth(ib+1);
    NTW [ib][2] =   hNTW ->GetBinError(ib+1); NTW [ib][3] = NTW [ib][0] *Syst_back[4];

    NDY [ib][0] = hNDY->GetBinContent(ib+1)+ hNDYtautau->GetBinContent(ib+1); NDY [ib][1] = hNDY ->GetBinWidth(ib+1);
    NDY [ib][2] =   hNDY ->GetBinError(ib+1); NDY [ib][3] = NDY [ib][0] *Syst_back[5];


    backgError [ib] =  sqrt(NWj [ib][3]*NWj [ib][3] + NWg [ib][3]*NWg [ib][3] + NWZ [ib][3]*NWZ [ib][3] 
			    + NZZ [ib][3]*NZZ [ib][3] +  NH125 [ib][3]* NH125 [ib][3] +  NTW [ib][3] *  NTW [ib][3] 
			    + NTT [ib][3] * NTT [ib][3] +  NDY [ib][3] * NDY [ib][3] ) / (WWdata->GetBinContent(ib+1));

    

  // Add relative systematic uncertainties
  //----------------------------------------------------------------------------
 
    for (UInt_t iSyst=0; iSyst<Nsyst; iSyst++) {
      
      NData[ib][3]+=(Syst_data [iSyst]*Syst_data [iSyst]);
      
    }


    /// ----------------- ON DATA 
    float nom = h_dataReco_unfolded    ->GetBinContent(ib+1);
    float up = systUp->GetBinContent(ib+1);
    float down = systDown->GetBinContent(ib+1);

    if ( nom <= 0) continue;

    cout << "ERROR!:  " << nom << "  " << up << "  " << down << endl;

    float errUp = up / nom;
    float errDown = down / nom;

    float maxErr = errUp; 
    
    if (errUp > errDown) {
      maxErr = errUp; 
    } else {
      maxErr = errDown; 
    }

    float total = sqrt((NData[ib][3]/10000) + (maxErr*maxErr) + backgError [ib]*backgError [ib]); 

    //float total = sqrt((NData[ib][3]/10000) +  backgError [ib]*backgError [ib]); 

    NData[ib][3]= total; //relative error 

    /// ----------------- 


    /// ----------------- ON POWHEG

    float nom_pow =  hNqqWW_pow   ->GetBinContent(ib+1);
    float up_pow = systUp_pow->GetBinContent(ib+1);
    float down_pow = systDown_pow->GetBinContent(ib+1);

    if ( nom_pow  <= 0) continue;

    float errUp_pow = up_pow / nom_pow;
    float errDown_pow = down_pow / nom_pow;

    float maxErr_pow = errUp_pow; 
    
    if (errUp_pow > errDown_pow) {
      maxErr_pow = errUp_pow; 
    } else {
      maxErr_pow = errDown_pow; 
    }

    NqqWW_pow [ib][3]= maxErr_pow; //relative error 

    /// ----------------- 


 /// ----------------- ON MCNLO

    float nom_mc =  hNqqWW_mcnlo   ->GetBinContent(ib+1);
    float up_mc = systUp_mc->GetBinContent(ib+1);
    float down_mc = systDown_mc->GetBinContent(ib+1);

    if ( nom_mc  <= 0) continue;

    float errUp_mc = up_mc / nom_mc;
    float errDown_mc = down_mc / nom_mc;

    float maxErr_mc = errUp_mc; 
    
    if (errUp_mc > errDown_mc) {
      maxErr_mc = errUp_mc; 
    } else {
      maxErr_mc = errDown_mc; 
    }

    NqqWW_mcnlo [ib][3]= maxErr_mc; //relative error 

   
    /// ----------------- 


  }



  // ---->  compute differential WW Xsec on data 
  //==============================================================================
  

  Double_t xsUnfold [bin];

  Double_t xsUnfold_stat [bin];
  
  TH1F *xsValue = (TH1F*) hNWW->Clone("xsValue");

  Double_t xsUnfold_fid = 5.24; 
  Double_t xsUnfold_fid_stat = 0.11; 
  Double_t xsUnfold_fid_syst = 0.21;  //removing everything that can cancell 
  Double_t xsUnfold_fid_err = 0.29; 



  for (int ib=0; ib < bin; ib++) { 
    //cout<< NData[ib][0]<< endl;
    xsUnfold [ib] = NData[ib][0] / (luminosity * NData[ib][1] *BR_WW_to_lnln );  
    xsUnfold_stat [ib] = NData[ib][2] / (luminosity  * NData[ib][1] *BR_WW_to_lnln );
 

    if (fiducial) { 
      xsValue->SetBinContent(ib+1, (xsUnfold[ib]/xsUnfold_fid));
      xsValue->SetBinError(ib+1,  (xsUnfold_stat [ib]/xsUnfold_fid));
    } else {
      xsValue->SetBinContent(ib+1, xsUnfold[ib]);
      xsValue->SetBinError(ib+1, xsUnfold_stat [ib]);
    }

  }



  // ---->  compute differential WW Xsec for MADGRAPH MC samples
  //==============================================================================
  
  Double_t xsMadgraph[bin];

  Double_t xsMadgraph_stat[bin];

  TH1F *xsValue_Madgraph = (TH1F*) hNWW->Clone("xsValue_mad");

  Double_t xsMadgraph_fid = 5.12; // 4.95 (5.12 for NNLL); 
  Double_t xsMadgraph_fid_err = 0.01;  // only statistical 


  for (int ib=0; ib < bin; ib++) { 
    
    xsMadgraph[ib] = NqqWW_mad [ib][0] / (luminosity  * NqqWW_mad[ib][1] *BR_WW_to_lnln);
    xsMadgraph_stat[ib] = NqqWW_mad [ib][2] / (luminosity * NqqWW_mad[ib][1] *BR_WW_to_lnln);

    /*    double xsMadgraph_stat_rel =  xsMadgraph_stat[ib] / xsMadgraph[ib]; 
    double total_rel_error = sqrt ( xsMadgraph_stat_rel*xsMadgraph_stat_rel + (xsMadgraph_fid_err * xsMadgraph_fid_err/(xsMadgraph_fid*xsMadgraph_fid)));
    double total_err = total_rel_error * xsMadgraph [ib]/xsMadgraph_fid;
    */

    if ( fiducial) {
      xsValue_Madgraph->SetBinContent(ib+1, (xsMadgraph [ib]/xsMadgraph_fid));
      xsValue_Madgraph->SetBinError(ib+1, (xsMadgraph_stat [ib]/xsMadgraph_fid));
    } else {
      xsValue_Madgraph->SetBinContent(ib+1, xsMadgraph [ib]);
      xsValue_Madgraph->SetBinError(ib+1, xsMadgraph_stat [ib]);
    }


  }

 // ---->  compute differential WW Xsec for POWHEG MC samples
  //==============================================================================
  
  Double_t xsPowheg[bin];

  Double_t xsPowheg_stat[bin], xsPowheg_tot[bin];

  TH1F *xsValue_Powheg = (TH1F*) hNWW->Clone("xsValue_pow");

  Double_t xsPowheg_fid = 5.29;  
  Double_t xsPowheg_fid_err = 0.03; //  statistical + systematics



  for (int ib=0; ib < bin; ib++) { 
    
    xsPowheg[ib] = NqqWW_pow [ib][0] / (luminosity  * NqqWW_pow[ib][1] *BR_WW_to_lnln);
    xsPowheg_stat[ib] = NqqWW_pow [ib][2] / (luminosity * NqqWW_pow[ib][1] *BR_WW_to_lnln);

    xsPowheg_tot[ib] = sqrt(NqqWW_pow [ib][3]*NqqWW_pow [ib][3]*xsPowheg[ib]*xsPowheg[ib]+xsPowheg_stat[ib]*xsPowheg_stat[ib]);

    /*    double xsPowheg_stat_rel =  xsPowheg_stat[ib] / xsPowheg[ib]; 
    double total_rel_error = sqrt ( xsPowheg_stat_rel*xsPowheg_stat_rel + (xsPowheg_fid_err * xsPowheg_fid_err/(xsPowheg_fid*xsPowheg_fid)));
    double total_err = total_rel_error * xsPowheg [ib]/xsPowheg_fid;
    */

    if (fiducial) {
      xsValue_Powheg->SetBinContent(ib+1, (xsPowheg [ib]/xsPowheg_fid));
      xsValue_Powheg->SetBinError(ib+1, (xsPowheg_tot [ib]/xsPowheg_fid));
    } else {
      xsValue_Powheg->SetBinContent(ib+1, xsPowheg [ib]);
      xsValue_Powheg->SetBinError(ib+1, xsPowheg_tot [ib]);
    }

  }

 // ---->  compute differential WW Xsec for MCNLO MC samples
  //==============================================================================
  
  Double_t xsMCnlo[bin];

  Double_t xsMCnlo_stat[bin];
  Double_t xsMCnlo_tot[bin];

  TH1F *xsValue_MCnlo = (TH1F*) hNWW->Clone("xsValue_mcnlo");


  Double_t xsMCnlo_fid = 5.25; // 5.20 (5.25 for NNLL) ; 
  Double_t xsMCnlo_fid_err = 0.03; // 0.02 (0.03 for NNLL) ; // only statistical 



  for (int ib=0; ib < bin; ib++) { 
    
    xsMCnlo[ib] = NqqWW_mcnlo [ib][0] / (luminosity  * NqqWW_mcnlo[ib][1] *BR_WW_to_lnln);
    xsMCnlo_stat[ib] = NqqWW_mcnlo [ib][2] / (luminosity  * NqqWW_mcnlo[ib][1] *BR_WW_to_lnln);
    xsMCnlo_tot[ib] = sqrt(NqqWW_mcnlo [ib][3]*NqqWW_mcnlo [ib][3]*xsMCnlo[ib]*xsMCnlo[ib]+xsMCnlo_stat[ib]*xsMCnlo_stat[ib]);

    /* double xsMCnlo_stat_rel =  xsMCnlo_stat[ib] / xsMCnlo[ib]; 
    double total_rel_error = sqrt ( xsMCnlo_stat_rel*xsMCnlo_stat_rel + (xsMCnlo_fid_err * xsMCnlo_fid_err/(xsMCnlo_fid*xsMCnlo_fid)));
    double total_err = total_rel_error * xsMCnlo [ib]/xsMCnlo_fid;
    */

    if (fiducial ) {
      xsValue_MCnlo->SetBinContent(ib+1, (xsMCnlo [ib]/xsMCnlo_fid));
      xsValue_MCnlo->SetBinError(ib+1, (xsMCnlo_tot [ib]/xsMCnlo_fid));
    } else {
      xsValue_MCnlo->SetBinContent(ib+1, xsMCnlo [ib]);
      xsValue_MCnlo->SetBinError(ib+1, xsMCnlo_stat [ib]);
    }  

}
  
  if( drawTheXS) { 

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

      //-- Plot Powheg
      xsValue_Powheg->SetLineColor(kGreen+2);
      xsValue_Powheg->SetMarkerColor(kGreen+2);
      xsValue_Powheg->SetLineWidth(2);

      //-- Plot Madgraph
      xsValue_Madgraph->SetLineColor(kRed);
      xsValue_Madgraph->SetMarkerColor(kRed);
      xsValue_Madgraph->SetLineWidth(2);
      xsValue_Madgraph->SetLineStyle(2);

   
      //-- Plot Madgraph
      xsValue_MCnlo->SetLineColor(kBlue);
      xsValue_MCnlo->SetMarkerColor(kBlue);
      xsValue_MCnlo->SetLineWidth(2);
      xsValue_MCnlo->SetLineStyle(3);


      //-- Plot Data

      xsValue->SetMarkerStyle(kFullCircle);
      
      TString title = Form("CMS preliminary (%s)", channel.Data() );

      //xsValue->SetTitle(title);
       xsValue->SetTitleSize(0.040);
       // #frac{d#sigma}
       if (differential == 0) xsValue->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}^{max}}");//("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{max}}");
      if (differential == 1) xsValue->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}(ll)}");
      if (differential == 2) xsValue->GetYaxis()->SetTitle("#frac{d#sigma}{dm_{#font[12]{ll}}}");
      if (differential == 3) xsValue->GetYaxis()->SetTitle("#frac{d#sigma}{d#Delta#phi_{ll}}");
      
      if (differential == 0 && fiducial) xsValue->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{max}}");//("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}^{max}}");
      if (differential == 1 && fiducial) xsValue->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dp_{T}(ll)}");
      if (differential == 2 && fiducial) xsValue->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{#font[12]{ll}}}");
      if (differential == 3 && fiducial) xsValue->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{d#Delta#phi_{ll}}");
      

      xsValue->GetYaxis()->SetTitleOffset(1.6);

      if (differential == 0) xsValue->GetXaxis()->SetTitle("p_{T}^{max}");
      if (differential == 1) xsValue->GetXaxis()->SetTitle("p_{T}(ll)");
      if (differential == 2) xsValue->GetXaxis()->SetTitle("m_{#font[12]{ll}}");
      if (differential == 3) xsValue->GetXaxis()->SetTitle("#Delta#phi_{ll}");
     

      // Include systematic errors 

      TH1F* systHisto = xsValue->Clone("syst");
     
      systHisto->SetFillColor  (kGray+2);
      systHisto->SetFillStyle  (   3345);
      systHisto->SetLineColor  (kGray+2);
      systHisto->SetMarkerColor(kGray+2);
      systHisto->SetMarkerSize (      0);
      
      for (int i=1; i<= NBins; i++) {

	Double_t binValue = 0;
	Double_t binError = 0;

	float value = xsValue->GetBinContent(i);
	float error = xsValue->GetBinError(i)/value;  // --> only relative statistical error 
	
       	float totalError = sqrt(NData[i-1][3]*NData[i-1][3] + error*error);

	systHisto->SetBinContent(i, value);
	systHisto->SetBinError(i, totalError*value);

	if (differential == 3 ) {

	  float errorbin1 = sqrt(totalError*totalError+0.041*0.041);
	  float errorbin2 = sqrt(totalError*totalError+0.078*0.078);
	
	  systHisto->SetBinError(NBins-1, errorbin1*value);
	  systHisto->SetBinError(NBins, errorbin2*value);
	}

	if ( differential == 6 && i == 2 ) {

	  cout << "statistical error: " << error*value << "; syst error: " <<  NData[1][3]*value << endl;

	}    
	
      }

     
      xsValue->GetXaxis()->SetTitleOffset(1.6);
      //      xsValue->Draw();
      


      xsValue->Draw("p");
      xsValue_Powheg->Draw("hist,same");
      xsValue_Madgraph->Draw("histe2,same");
      xsValue_MCnlo->Draw("histe2,same");
      systHisto->Draw("e2, same");
      xsValue->Draw("p,same");
      

      if ( differential == 6) {
	cout << "[DATA]   FIDUTIAL CROSS-SECTION:  " << xsValue->GetBinContent(2) << " +- " << systHisto->GetBinError(2) << endl;
	cout << "[POWHEG] FIDUTIAL CROSS-SECTION:  " << xsValue_Powheg->GetBinContent(2) << " +- " << xsValue_Powheg->GetBinError(2) << endl;
	cout << "[MADGRAPH] FIDUTIAL CROSS-SECTION:  " << xsValue_Madgraph->GetBinContent(2) << " +- " << xsValue_Madgraph->GetBinError(2) << endl;
	cout << "[MCNLO] FIDUTIAL CROSS-SECTION:  " << xsValue_MCnlo->GetBinContent(2) << " +- " << xsValue_MCnlo->GetBinError(2) << endl;
      }


      // Legend
      //----------------------------------------------------------------------------
      TLegend* legend = new TLegend(0.85, 0.65, 0.65, 0.89);
      
  
      SetLegendStyle(legend);

      legend->AddEntry(xsValue,   "Data", "P");
      legend->AddEntry(xsValue_Powheg,   "Powheg", "L");
      legend->AddEntry(xsValue_Madgraph,   "Madgraph", "L");      
      legend->AddEntry(xsValue_MCnlo,   "MCNLO", "L");
 


      // Put everything together
      //----------------------------------------------------------------------------
      legend->Draw("same");


      // Draw also ratio
      //----------------------------------------------------------------------------
      if (drawRatio) { 
	
	pad2->cd();
 
	TH1F* ratio_pow       = xsValue_Powheg->Clone("ratio");
	TH1F* ratio_mad       = xsValue_Madgraph->Clone("ratio");
	TH1F* ratio_mcnlo     = xsValue_MCnlo->Clone("ratio");
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

	  Double_t ratioValue_pow           = (powValue > 0) ? dataValue/powValue : 0.0;
	  Double_t ratioError_pow           = (powValue > 0) ? dataError/powValue : 0.0;

	  Double_t ratioValue_mad           = (madValue > 0) ? dataValue/madValue : 0.0;
	  Double_t ratioError_mad           = (madValue > 0) ? dataError/madValue : 0.0;

	  Double_t ratioValue_mcnlo         = (mcnloValue > 0) ? dataValue/mcnloValue : 0.0;
	  Double_t ratioError_mcnlo         = (mcnloValue > 0) ? dataError/mcnloValue : 0.0;

	  Double_t uncertaintyError         = (dataValue > 0) ? dataError/dataValue : 0.0;
	  
	  //dataError/dataValue 
	  ratio_pow->SetBinContent(ibin, ratioValue_pow);
	  ratio_pow->SetBinError  (ibin, ratioError_pow);

	  ratio_mad->SetBinContent(ibin, ratioValue_mad);
	  ratio_mad->SetBinError  (ibin, ratioError_mad);

	  ratio_mcnlo->SetBinContent(ibin, ratioValue_mcnlo);
	  ratio_mcnlo->SetBinError  (ibin, ratioError_mcnlo);

	  ratioErr->SetBinContent(ibin, 1.0);
	  ratioErr->SetBinError  (ibin, uncertaintyError);
	}

	ratioErr->SetTitle("");
	ratioErr  ->Draw("e2");

	ratio_pow      ->Draw("ep,same");
	ratio_pow      ->SetLineColor(kGreen+2);
	ratio_pow      ->SetMarkerSize(1.0);
	ratio_pow      ->SetLineWidth(2);
	ratio_pow      ->SetMarkerStyle(20);

	ratio_mad      ->Draw("ep,same");
	ratio_mad      ->SetLineColor(kRed);
	ratio_mad      ->SetMarkerSize(1.0);
	ratio_mad      ->SetLineWidth(2);
	ratio_mad      ->SetMarkerStyle(20);

	ratio_mcnlo      ->Draw("ep,same");
	ratio_mcnlo     ->SetLineColor(kBlue);
	ratio_mcnlo     ->SetMarkerSize(1.0);
	ratio_mcnlo      ->SetLineWidth(2);
	ratio_mcnlo     ->SetMarkerStyle(20);

	ratioErr->SetLineColor(kBlue);
	ratioErr->SetFillColor  (kGray+2);
	ratioErr->SetFillStyle  (   3345);
	ratioErr->SetLineColor  (kGray+2);
	ratioErr->SetMarkerColor(kGray+2);
	ratioErr->SetMarkerSize (      0);

	ratioErr->GetYaxis()->SetRangeUser(0.4, 1.6);

	Pad2TAxis(ratioErr, xsValue->GetXaxis()->GetTitle(), "unfolded / generated");


      }

      // And save it
      //----------------------------------------------------------------------------

    

      TString addedDistribution;

      if (differential == 0) addedDistribution = "pt1";
      if (differential == 1) addedDistribution = "ptll";
      if (differential == 2) addedDistribution = "mll";
      if (differential == 3) addedDistribution = "dphill";


      TString file = Form("pdf/wwxsec_%djet_%s_%s.pdf", njet, channel.Data(), addedDistribution.Data());      
      canvas->SaveAs(file);

      if (differential == 0 ) {
	TFile* output = new TFile("XSLeadingPt_AN.root", "update");
	output->cd();
	xsValue->Write("xsValue", TObject::kWriteDelete );
	xsValue_Powheg->Write("xsValue_Powheg", TObject::kWriteDelete );
	xsValue_Madgraph->Write("xsValue_Madgraph", TObject::kWriteDelete );
	xsValue_MCnlo->Write("xsValue_MCnlo", TObject::kWriteDelete );
	systHisto->Write("systHisto", TObject::kWriteDelete );
	output->Close();
      }

      if (differential == 1 ) {
	TFile* output = new TFile("XSDileptonPt_AN.root", "update");
	output->cd();
	xsValue->Write("xsValue", TObject::kWriteDelete );
	xsValue_Powheg->Write("xsValue_Powheg", TObject::kWriteDelete );
	xsValue_Madgraph->Write("xsValue_Madgraph", TObject::kWriteDelete );
	xsValue_MCnlo->Write("xsValue_MCnlo", TObject::kWriteDelete );
	systHisto->Write("systHisto", TObject::kWriteDelete );
	output->Close();
      }
      if (differential == 2 ) {
	TFile* output = new TFile("XSMll_AN.root", "update");
	output->cd();
	xsValue->Write("xsValue", TObject::kWriteDelete );
	xsValue_Powheg->Write("xsValue_Powheg", TObject::kWriteDelete );
	xsValue_Madgraph->Write("xsValue_Madgraph", TObject::kWriteDelete );
	xsValue_MCnlo->Write("xsValue_MCnlo", TObject::kWriteDelete );
	systHisto->Write("systHisto", TObject::kWriteDelete );
	output->Close();
      }
      if (differential == 3 ) {
	TFile* output = new TFile("XSDeltaPhi_AN.root", "update");
	output->cd();
	xsValue->Write("xsValue", TObject::kWriteDelete );
	xsValue_Powheg->Write("xsValue_Powheg", TObject::kWriteDelete );
	xsValue_Madgraph->Write("xsValue_Madgraph", TObject::kWriteDelete );
	xsValue_MCnlo->Write("xsValue_MCnlo", TObject::kWriteDelete );
	systHisto->Write("systHisto", TObject::kWriteDelete );
	output->Close();
      }
  }
  
 

}



//==============================================================================
//Tool for computing and printing correlation matrix 
//==============================================================================

TH2D* CorrelationHist (const TMatrixD& cov,
		       const char* name, const char* title,
		       Double_t lo, Double_t hi, Bool_t _useOvf)
{

  // correlation coef == cov(xy)/(sigmay*sigmay);

  Int_t nb= cov.GetNrows();

  TH2D* h; 

  h = new TH2D (name, title, nb, lo, hi, nb, lo, hi);

  if (_useOvf ) h = new TH2D (name, title, nb-2, lo, hi, nb-2, lo, hi);
  
  h->SetAxisRange (-1.0, 1.0, "Z");
  
  if (_useOvf) nb = nb-2;

  for(int i=0; i < nb; i++) {
  
    for(int j=0; j < nb; j++) {
      
      if (_useOvf ) {
      
	Double_t Viijj= cov(i+1,i+1)*cov(j+1,j+1);
	h->SetBinContent (i+1, j+1, cov(i+1,j+1)/sqrt(Viijj));
      } else {
	Double_t Viijj= cov(i,i)*cov(j,j);
	h->SetBinContent (i+1, j+1, cov(i,j)/sqrt(Viijj));
      }
    }
  }

  return h;
}




//==============================================================================
//Tool for printing covariance matrix 
//==============================================================================

void PrintMatrix(const TMatrixD& m, const char* format,
                                       const char* name, Int_t cols_per_sheet)
{

  // Print the matrix as a table of elements.
   // Based on TMatrixTBase<>::Print, but allowing user to specify name and cols_per_sheet (also option -> format).
   // By default the format "%11.4g" is used to print one element.
   // One can specify an alternative format with eg
   //  format ="%6.2f  "

   if (!m.IsValid()) {
     m.Error("PrintMatrix","%s is invalid",name);
     return;
   }

   const Int_t ncols  = m.GetNcols();
   const Int_t nrows  = m.GetNrows();
   const Int_t collwb = m.GetColLwb();
   const Int_t rowlwb = m.GetRowLwb();

   if (!(format && format[0])) format= "%11.4g ";
   char topbar[1000];
   snprintf(topbar,1000,format,123.456789);
   Int_t nch = strlen(topbar)+1;
   if (nch > 18) nch = 18;
   char ftopbar[20];
   for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
   Int_t nk = 1 + Int_t(log10(ncols));
   snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
   Int_t nch2 = strlen(ftopbar);
   for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
   ftopbar[nch] = '|';
   ftopbar[nch+1] = 0;

   printf("\n%dx%d %s is as follows",nrows,ncols,name);

   if (cols_per_sheet <= 0) {
     cols_per_sheet = 5;
     if (nch <= 8) cols_per_sheet =10;
   }
   nk = 5+nch*(cols_per_sheet<ncols ? cols_per_sheet : ncols);
   for (Int_t i = 0; i < nk; i++) topbar[i] = '-';
   topbar[nk] = 0;
   for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
      printf("\n\n     |");
      for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
         printf(ftopbar,j+collwb-1);
      printf("\n%s\n",topbar);
      if (m.GetNoElements() <= 0) continue;
      for (Int_t i = 1; i <= nrows; i++) {
         printf("%4d |",i+rowlwb-1);
         for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++)
            printf(format,m(i+rowlwb-1,j+collwb-1));
         printf("\n");
      }
   }
   printf("\n");

}






//------------------------------------------------------------------------------
// SetLegendStyle
//------------------------------------------------------------------------------


void SetLegendStyle(TLegend* legend) {

  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (   0.040);
}



//------------------------------------------------------------------------------
// Pad2TAxis
//------------------------------------------------------------------------------
void Pad2TAxis(TH1* hist, TString xtitle, TString ytitle)
{
  TAxis* xaxis = (TAxis*)hist->GetXaxis();
  TAxis* yaxis = (TAxis*)hist->GetYaxis();

  xaxis->SetLabelFont  (    42);
  xaxis->SetLabelOffset( 0.025);
  xaxis->SetLabelSize  (   0.08);
  xaxis->SetNdivisions (   505);
  xaxis->SetTitle      (xtitle);
  xaxis->SetTitleFont  (    42);
  xaxis->SetTitleOffset(  1.35);
  xaxis->SetTitleSize  (  0.1);

  yaxis->CenterTitle   (      );
  yaxis->SetLabelFont  (    42);
  yaxis->SetLabelOffset(  0.02);
  yaxis->SetLabelSize  (   0.08);
  yaxis->SetNdivisions (   505);
  yaxis->SetTitle      (ytitle);
  yaxis->SetTitleFont  (    42);
  yaxis->SetTitleOffset(  0.75);
  yaxis->SetTitleSize  (  0.08);
}
