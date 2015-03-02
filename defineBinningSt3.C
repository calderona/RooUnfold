
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"
#include <iomanip>
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom.h"

//gSystem->Load("libRooUnfold");

#include "src/RooUnfoldResponse.h"



//------------------------------------------------------------------------------
// defineBinning
//------------------------------------------------------------------------------


void defineBinning(Int_t   jetChannel = 0 , TString theSample = "WW", Bool_t jetGenVeto = 0 ) {

  gSystem->Load("libRooUnfold");

  TH1::SetDefaultSumw2();

  //TString path = Form("rootfiles/%djet/%s/", jetChannel, flavorChannel.Data());
 
  //gSystem->mkdir(path, kTRUE);
 
 
  //----------------------------------------------------------------------------
  // Input files
  //----------------------------------------------------------------------------

  TString filesPath;

  //filesPath = "/gpfs/csic_projects/tier3data/LatinosSkims/ReducedTrees/DiferentialXSection/";
  filesPath = "/gpfs/csic_projects/cms/calderon/WWGEN/";
  

  TChain* tree = new TChain("latino", "latino"); 


  //-*** Final for data 
  //tree->Add("/gpfs/csic_projects/cms/calderon/WWGEN/finalFilesForData/latino006_nll_ewk.root"); // POW 
  tree->Add("/gpfs/csic_projects/cms/calderon/WWGEN/finalFilesForData/latino_001_GGWWJets_jetSmear_METxy.root"); // GG
  //tree->Add("/gpfs/csic_projects/cms/calderon/WWGEN/finalFilesForData/latino000_nll_ewk.root"); // MAD
  //tree->Add("/gpfs/csic_projects/cms/calderon/WWGEN/finalFilesForData/latino002_nll_ewk.root"); // MCNLO



  //----------------------------------------------------------------------------
  // Define functions
  //----------------------------------------------------------------------------
  Float_t smear (Float_t xt); 

   //----------------------------------------------------------------------------
  // Output files
  //----------------------------------------------------------------------------
  

  TString path = Form("_gg_status3.root");

  TFile* output = new TFile( theSample+path, "recreate");


  // Defining binning
  //----------------------------------------------------------------------------

  //Double_t pt1bins[6] = {25,50,100,150,200,400};
  
  Double_t RECOpt1bins[11] = {10,20,40,60,80,100,125,150,175,200,210};
  Double_t GENpt1bins[7] = {10,20,50,100,150,200,210};

  const Int_t pt1Nbin = 9;
  const Int_t ptllNbin = 8; 
  const Int_t mllNbin = 9;
  const Int_t dphiNbin = 13;
  const Int_t jetEtNbin = 10;

  Double_t pt1bins[pt1Nbin] = {20,40,60,80,100,125,150,175,200}; 
  Double_t ptllbins[ptllNbin] = {30,40,50,60,70,85,120,150};
  Double_t mllbins[mllNbin] = {20,40,60,80,100,125,150,175,200};
  Double_t dphibins[dphiNbin] = {0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3};//{0,0.5,1,1.5,2,2.5,3};
  Double_t jetEtbins[jetEtNbin] = {30,40,50,60,70,80,90,100,110,120}; 


  // Pt, Dilepton, DeltaPhi, Mll

  // GEN level ( phase space)  differential histograms 
  //----------------------------------------------------------------------------

  TH1F* hPtLepton1_GEN  = new TH1F("hPtLepton1_GEN",       "", pt1Nbin-1, pt1bins);
  TH1F* hPtLepton1NoTaus_GEN  = new TH1F("hPtLepton1NoTaus_GEN",       "", pt1Nbin-1, pt1bins);
  TH1F* hPtLepton1S3NoTaus_GEN  = new TH1F("hPtLepton1S3NoTaus_GEN",       "", pt1Nbin-1, pt1bins);

  TH1F* hPtLepton1_GEN_C  = new TH1F("hPtLepton1_GEN_C",       "", 180, 20, 200);
  TH1F* hPtLepton1NoTaus_GEN_C  = new TH1F("hPtLepton1NoTaus_GEN_C",       "", 180, 20, 200);
  TH1F* hPtLepton1S3NoTaus_GEN_C  = new TH1F("hPtLepton1S3NoTaus_GEN_C",       "", 180, 20, 200);

  TH1F* hPtLepton1_RECO  = new TH1F("hPtLepton1_RECO",       "", pt1Nbin-1, pt1bins);
  TH2F* hPtLepton1_RECO_GEN =  new TH2F("hPtLepton1_RECO_GEN", "", pt1Nbin-1, pt1bins, pt1Nbin-1, pt1bins);

  TH1F* hDilepton_GEN  = new TH1F("hDilepton_GEN",       "", ptllNbin-1, ptllbins);
  TH1F* hDileptonNoTaus_GEN  = new TH1F("hDileptonNoTaus_GEN",       "", ptllNbin-1, ptllbins);
  TH1F* hDileptonS3NoTaus_GEN  = new TH1F("hDileptonS3NoTaus_GEN",       "", ptllNbin-1, ptllbins);

  TH1F* hDilepton_RECO  = new TH1F("hDilepton_RECO",       "", ptllNbin-1, ptllbins);
  TH2F* hDilepton_RECO_GEN =  new TH2F("hDilepton_RECO_GEN", "", ptllNbin-1, ptllbins, ptllNbin-1, ptllbins);

  TH1F* hmll_GEN  = new TH1F("hmll_GEN",       "", mllNbin-1, mllbins);
  TH1F* hmllNoTaus_GEN  = new TH1F("hmllNoTaus_GEN",       "", mllNbin-1, mllbins);
  TH1F* hmllS3NoTaus_GEN  = new TH1F("hmllS3NoTaus_GEN",       "", mllNbin-1, mllbins);

  TH1F* hmll_RECO  = new TH1F("hmll_RECO",       "",mllNbin-1,  mllbins);
  TH2F* hmll_RECO_GEN =  new TH2F("hmll_RECO_GEN", "", mllNbin-1, mllbins, mllNbin-1, mllbins);

  TH1F* hdphi_GEN  = new TH1F("hdphi_GEN",       "", dphiNbin-1, dphibins);
  TH1F* hdphiNoTaus_GEN  = new TH1F("hdphiNoTaus_GEN",       "", dphiNbin-1, dphibins);
  TH1F* hdphiS3NoTaus_GEN  = new TH1F("hdphiS3NoTaus_GEN",       "", dphiNbin-1, dphibins);


  TH1F* hdphi_RECO  = new TH1F("hdphi_RECO",       "", dphiNbin-1,  dphibins);
  TH1F* hdphiNoTaus_RECO  = new TH1F("hdphiNoTaus_RECO",       "", dphiNbin-1,  dphibins); 

  TH2F* hdphi_RECO_GEN =  new TH2F("hdphi_RECO_GEN", "", dphiNbin-1, dphibins, dphiNbin-1, dphibins);

  TH1F* hjetEt_GEN  = new TH1F("hjetEt_GEN",       "", jetEtNbin-1, jetEtbins);
  TH1F* hjetEt_RECO  = new TH1F("hjetEt_RECO",       "", jetEtNbin-1,  jetEtbins);
  TH2F* hjetEt_RECO_GEN =  new TH2F("hjetEt_RECO_GEN", "", jetEtNbin-1, jetEtbins, jetEtNbin-1, jetEtbins);

  TH1F* hInclusive_GEN  = new TH1F("hInclusive_GEN",       "", 3,0,3);

  TH1F* hjetEt_Nolepton_GEN  = new TH1F("hjetEt_Nolepton_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_lepton_GEN  = new TH1F("hjetEt_lepton_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_leptonID_GEN  = new TH1F("hjetEt_leptonID_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_leptonEta_GEN  = new TH1F("hjetEt_leptonEta_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_leptonEtaPt1_GEN  = new TH1F("hjetEt_leptonEtaPt1_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_leptonEtaPt1Pt2_GEN  = new TH1F("hjetEt_leptonEtaPt1Pt2_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_hadronicWW_GEN  = new TH1F("hjetEt_hadronicWW_GEN",       "", 100, 5, 100);

  TH1F* hjetEt_Nofid_GEN  = new TH1F("hjetEt_Nofid_GEN",       "", 100, 5, 100);
  TH1F* hjetEt_fid_GEN  = new TH1F("hjetEt_fid_GEN",       "", 100, 5, 100);

  TH2F* hjetEt_vs_leptonpt1_GEN = new TH2F("hjetEt_vs_leptonpt1_GEN", "", 100,5,100,100,5,100);
  TH2F* hjetEt_vs_leptonpt2_GEN = new TH2F("hjetEt_vs_leptonpt2_GEN", "", 100,5,100,100,5,100);

  RooUnfoldResponse responsePtLepton1GEN(hPtLepton1_RECO, hPtLepton1_GEN);
  RooUnfoldResponse responseMllGEN(hmll_RECO, hmll_GEN );
  RooUnfoldResponse responseJetEtGEN(hjetEt_RECO, hjetEt_GEN );


  // WW level differential histograms 
  //----------------------------------------------------------------------------

  TH1F* hPtLepton1WWLevel_RECO  = new TH1F("hPtLepton1WWLevel_RECO",       "", pt1Nbin-1, pt1bins);
  TH1F* hPtLepton1WWLevel_GEN  = new TH1F("hPtLepton1WWLevel_GEN",       "", pt1Nbin-1, pt1bins);
  TH2F* hPtLepton1WWLevel_RECO_GEN =  new TH2F("hPtLepton1WWLevel_RECO_GEN", "", pt1Nbin-1, pt1bins, pt1Nbin-1, pt1bins);

 
  TH1F* hDileptonWWLevel_RECO  = new TH1F("hDileptonWWLevel_RECO",       "", ptllNbin-1, ptllbins);
  TH1F* hDileptonWWLevel_GEN  = new TH1F("hDileptonWWLevel_GEN",       "", ptllNbin-1, ptllbins);
  TH2F* hDileptonWWLevel_RECO_GEN =  new TH2F("hDileptonWWLevel_RECO_GEN", "", ptllNbin-1, ptllbins, ptllNbin-1, ptllbins);

  TH1F* hmllWWLevel_RECO  = new TH1F("hmllWWLevel_RECO",       "", mllNbin-1, mllbins);
  TH1F* hmllWWLevel_GEN  = new TH1F("hmllWWLevel_GEN",       "", mllNbin-1, mllbins);
  TH2F* hmllWWLevel_RECO_GEN =  new TH2F("hmllWWLevel_RECO_GEN", "", mllNbin-1, mllbins,mllNbin-1, mllbins);

  TH1F* hdphiWWLevel_RECO  = new TH1F("hdphiWWLevel_RECO",       "", dphiNbin-1, dphibins);
  TH1F* hdphiWWLevel_GEN  = new TH1F("hdphiWWLevel_GEN",       "", dphiNbin-1, dphibins);
  TH2F* hdphiWWLevel_RECO_GEN =  new TH2F("hdphiWWLevel_RECO_GEN", "", dphiNbin-1, dphibins,dphiNbin-1, dphibins);

  TH1F* hjetEtWWLevel_GEN  = new TH1F("hjetEtWWLevel_GEN",       "", jetEtNbin-1, jetEtbins);
  TH1F* hjetEtWWLevel_RECO  = new TH1F("hjetEtWWLevel_RECO",       "", 1,  jetEtbins);
  TH2F* hjetEtWWLevel_RECO_GEN =  new TH2F("hjetEtWWLevel_RECO_GEN", "", jetEtNbin-1, jetEtbins, jetEtNbin-1, jetEtbins);
  
  TH1D* hPtLepton1WWLevel_nonselected = new TH1D("hPtLepton1WWLevel_nonselected", "hPtLepton1WWLevel_nonselected",pt1Nbin-1, pt1bins);

  TH1F* hdeltaR = new TH1F("hdeltaR", "hdeltaR", 50,0,0.1);

  TH1F* hgenJetEt = new TH1F("hgenJetEt", "hgenJetEt", 100,0,100);
  TH1F* hrecoJetEt = new TH1F("hrecoJetEt", "hrecoJetEt", 100,0,100);
  TH2F* hJetEt_Gen_Reco = new TH2F("hJetEt_Gen_Reco", "hJetEt_Gen_Reco", 100,0,100,100,0,100);

  TH1F* hInclusiveWWLevel_GEN  = new TH1F("hInclusiveWWLevel_GEN",       "",3,0,3);
  TH1F* hInclusiveWWLevel_RECO  = new TH1F("hInclusiveWWLevel_RECO",       "",3,0,3);


  // WW level  define response matrix
  
  RooUnfoldResponse responsePtLepton1(hPtLepton1WWLevel_RECO, hPtLepton1WWLevel_GEN);

  RooUnfoldResponse responseDilepton(hDileptonWWLevel_RECO, hDileptonWWLevel_GEN);

  RooUnfoldResponse responseMll(hmllWWLevel_RECO, hmllWWLevel_GEN );

  RooUnfoldResponse responseDphi(hdphiWWLevel_RECO , hdphiWWLevel_GEN);

  RooUnfoldResponse responseJetEt(hjetEtWWLevel_RECO , hjetEtWWLevel_GEN);

  RooUnfoldResponse responseInclusive(hInclusiveWWLevel_RECO, hInclusiveWWLevel_GEN );

  // Declaration of leaf types
  //----------------------------------------------------------------------------

  Float_t baseW;        tree->SetBranchAddress("baseW"       , &baseW);
  Float_t channel;      tree->SetBranchAddress("channel"     , &channel);
  Float_t chmet;        tree->SetBranchAddress("chmet"       , &chmet);
  Float_t dataset;      tree->SetBranchAddress("dataset"     , &dataset);
  Float_t dphill;       tree->SetBranchAddress("dphill"      , &dphill);
  Float_t dphilljet;    tree->SetBranchAddress("dphilljet"   , &dphilljet);
  Float_t dphilljetjet; tree->SetBranchAddress("dphilljetjet", &dphilljetjet);
  Float_t drll;         tree->SetBranchAddress("drll"        , &drll);
  Float_t effW;         tree->SetBranchAddress("effW"        , &effW);
  Float_t jetphi1;      tree->SetBranchAddress("jetphi1"     , &jetphi1);
  Float_t jetphi2;      tree->SetBranchAddress("jetphi2"     , &jetphi2);
  Float_t jetphi3;      tree->SetBranchAddress("jetphi3"     , &jetphi3);
  Float_t jeteta1;      tree->SetBranchAddress("jeteta1"     , &jeteta1);
  Float_t jeteta2;      tree->SetBranchAddress("jeteta2"     , &jeteta2);
  Float_t jeteta3;      tree->SetBranchAddress("jeteta3"     , &jeteta3);
  Float_t jetpt1;       tree->SetBranchAddress("jetpt1"      , &jetpt1);
  Float_t jetpt2;       tree->SetBranchAddress("jetpt2"      , &jetpt2);
  Float_t jetpt3;       tree->SetBranchAddress("jetpt3"      , &jetpt3);
  Float_t jettche1;     tree->SetBranchAddress("jettche1"    , &jettche1);
  Float_t jettche2;     tree->SetBranchAddress("jettche2"    , &jettche2);
  Float_t mctruth;      tree->SetBranchAddress("mctruth"     , &mctruth);
  Float_t mll;          tree->SetBranchAddress("mll"         , &mll);
  Float_t mpmet;        tree->SetBranchAddress("mpmet"       , &mpmet); 
  Float_t mth;          tree->SetBranchAddress("mth"         , &mth);
  Float_t nbjet;        tree->SetBranchAddress("nbjet"       , &nbjet);
  Float_t nbjettche;    tree->SetBranchAddress("nbjettche"   , &nbjettche);
  Float_t nextra;       tree->SetBranchAddress("nextra"      , &nextra);
  Float_t njet;         tree->SetBranchAddress("njet"        , &njet);
  Float_t nvtx;         tree->SetBranchAddress("nvtx"        , &nvtx);
  Float_t pchmet;       tree->SetBranchAddress("pchmet"      , &pchmet);
  Float_t pfmet;        tree->SetBranchAddress("pfmet"       , &pfmet);
  Float_t ppfmet;       tree->SetBranchAddress("ppfmet"      , &ppfmet);
  Float_t isomva1;      tree->SetBranchAddress("isomva1"     , &isomva1);
  Float_t isomva2;      tree->SetBranchAddress("isomva2"     , &isomva2);
  Float_t pt1;          tree->SetBranchAddress("pt1"         , &pt1);
  Float_t pt2;          tree->SetBranchAddress("pt2"         , &pt2);
  Float_t pt3;          tree->SetBranchAddress("pt3"         , &pt3);
  Float_t phi1;         tree->SetBranchAddress("phi1"         , &phi1);
  Float_t dymva1;       tree->SetBranchAddress("dymva1"         , &dymva1);
  Float_t phi2;         tree->SetBranchAddress("phi2"         , &phi2);
  Float_t eta1;         tree->SetBranchAddress("eta1"         , &eta1);
  Float_t eta2;         tree->SetBranchAddress("eta2"         , &eta2);
  Float_t ch1;          tree->SetBranchAddress("ch1"         , &ch1);
  Float_t ch2;          tree->SetBranchAddress("ch2"         , &ch2);
  Float_t ptll;         tree->SetBranchAddress("ptll"        , &ptll);
  Float_t softtche;     tree->SetBranchAddress("softtche"    , &softtche);
  Float_t trigger;      tree->SetBranchAddress("trigger"     , &trigger);
  Float_t triggW;       tree->SetBranchAddress("triggW"      , &triggW);
  Int_t   bveto;        tree->SetBranchAddress("bveto"       , &bveto);
  Int_t   bveto_ip;     tree->SetBranchAddress("bveto_ip"    , &bveto_ip);
  Int_t   bveto_mu;     tree->SetBranchAddress("bveto_mu"    , &bveto_mu);
  Int_t   bveto_nj30;   tree->SetBranchAddress("bveto_nj30"  , &bveto_nj30);
  Int_t   dphiveto;     tree->SetBranchAddress("dphiveto"    , &dphiveto);
  Int_t   sameflav;     tree->SetBranchAddress("sameflav"    , &sameflav);
  Int_t   zveto;        tree->SetBranchAddress("zveto"       , &zveto);
  UInt_t  event;        tree->SetBranchAddress("event"       , &event);
  UInt_t  lumi;         tree->SetBranchAddress("lumi"        , &lumi);
  UInt_t  run;          tree->SetBranchAddress("run"         , &run);
  Float_t puW;          tree->SetBranchAddress("puW"         , &puW);
  

  // Apply NNLL resummation 
  Float_t nllW = 1; tree->SetBranchAddress("nllW", &nllW);

// GEN info... 


//Define Status1 leptons 

  Float_t lepGenpt1, lepGenpt2, lepGenpt3;
  tree->SetBranchAddress("genVV_lepton1_pt", &lepGenpt1);
  tree->SetBranchAddress("genVV_lepton2_pt", &lepGenpt2);
  tree->SetBranchAddress("genVV_lepton3_pt", &lepGenpt3);

  Float_t lepGeneta1, lepGeneta2, lepGeneta3;
  tree->SetBranchAddress("genVV_lepton1_eta", &lepGeneta1);
  tree->SetBranchAddress("genVV_lepton2_eta", &lepGeneta2);
  tree->SetBranchAddress("genVV_lepton3_eta", &lepGeneta3);

  Float_t lepGenphi1, lepGenphi2, lepGenphi3;
  tree->SetBranchAddress("genVV_lepton1_phi", &lepGenphi1);
  tree->SetBranchAddress("genVV_lepton2_phi", &lepGenphi2);
  tree->SetBranchAddress("genVV_lepton3_phi", &lepGenphi3);

  Float_t lepGenM1, lepGenM2, lepGenM3;
  tree->SetBranchAddress("genVV_lepton1_oVpid", &lepGenM1); 
  tree->SetBranchAddress("genVV_lepton2_oVpid", &lepGenM2); 
  tree->SetBranchAddress("genVV_lepton3_oVpid", &lepGenM3); 

  Float_t lepGenimTau1, lepGenimTau2, lepGenimTau3;
  tree->SetBranchAddress("genVV_lepton1_imTau", &lepGenimTau1); 
  tree->SetBranchAddress("genVV_lepton2_imTau", &lepGenimTau2); 
  tree->SetBranchAddress("genVV_lepton3_imTau", &lepGenimTau3); 

  Float_t lepGenpid1, lepGenpid2, lepGenpid3;
  tree->SetBranchAddress("genVV_lepton1_pid", &lepGenpid1);
  tree->SetBranchAddress("genVV_lepton2_pid", &lepGenpid2);
  tree->SetBranchAddress("genVV_lepton3_pid", &lepGenpid3);

  Float_t jetGen1_pt, jetGen2_pt, jetGen3_pt, jetGen4_pt, jetGen5_pt;
  tree->SetBranchAddress("genVV_jet1_pt", &jetGen1_pt);
  tree->SetBranchAddress("genVV_jet2_pt", &jetGen2_pt);
  tree->SetBranchAddress("genVV_jet3_pt", &jetGen3_pt);
  tree->SetBranchAddress("genVV_jet4_pt", &jetGen4_pt);  
  tree->SetBranchAddress("genVV_jet5_pt", &jetGen5_pt);

  Float_t jetGen1_eta, jetGen2_eta, jetGen3_eta, jetGen4_eta, jetGen5_eta;
  tree->SetBranchAddress("genVV_jet1_eta", &jetGen1_eta);
  tree->SetBranchAddress("genVV_jet2_eta", &jetGen2_eta);
  tree->SetBranchAddress("genVV_jet3_eta", &jetGen3_eta);
  tree->SetBranchAddress("genVV_jet4_eta", &jetGen4_eta);  
  tree->SetBranchAddress("genVV_jet5_eta", &jetGen5_eta);

  Float_t jetGen1_phi, jetGen2_phi, jetGen3_phi, jetGen4_phi, jetGen5_phi;
  tree->SetBranchAddress("genVV_jet1_phi", &jetGen1_phi);
  tree->SetBranchAddress("genVV_jet2_phi", &jetGen2_phi);
  tree->SetBranchAddress("genVV_jet3_phi", &jetGen3_phi);
  tree->SetBranchAddress("genVV_jet4_phi", &jetGen4_phi);  
  tree->SetBranchAddress("genVV_jet5_phi", &jetGen5_phi);


 

 // Set the channel
  //----------------------------------------------------------------------------
  Float_t SelectedChannel = -999;

  /*  if      (flavorChannel == "MuMu") SelectedChannel =  0;
  else if (flavorChannel == "EE"  ) SelectedChannel =  1;
  else if (flavorChannel == "EMu" ) SelectedChannel =  2;
  else if (flavorChannel == "MuE" ) SelectedChannel =  3;
  else if (flavorChannel == "All" ) SelectedChannel = -1;
  */

  int kk = 0;

 //----------------------------------------------------------------------------
  // Loop
  //----------------------------------------------------------------------------
  

  // Float_t Nentries = 179545; 
  //Float_t Nentries = 187525;
  Float_t Nentries = 10745; // with powheg 
  //Float_t Nentries = 9792; // with madgraph 

  for (int ievent=0; ievent<tree->GetEntriesFast(); ievent++) {
  //for (int ievent=0; ievent<Nentries; ievent++) {
  //for (int ievent=Nentries; ievent<tree->GetEntriesFast(); ievent++) {
   
    tree->GetEntry(ievent);

    //Double_t mybaseW =  5984.0/1933235;//5812.3/1933235; // madgraph (1933232)
    //Double_t mybaseW = 5984.0/999864;//  5812.3/999864; // powheg (999860)
    Double_t mybaseW = 182.852 /109986; // GGWW 
    //Double_t mybaseW = 5984.0/539594; //5812.3/539594; // mcnlo



    Float_t luminosity = 19.365;

    Double_t efficiencyW =  effW * triggW * puW;
    // Double_t totalW = effW * triggW * baseW * efficiencyW * luminosity;

    Double_t totalW = puW * mybaseW * luminosity;//efficiencyW;

    Double_t totalWGen = nllW *  mybaseW * luminosity ; //Removing PU from gen level
    
    Double_t totalWReco =  effW * triggW * puW * nllW *  mybaseW * luminosity;// * effW * triggW ;//puW *  mybaseW * luminosity;//puW * effW * triggW * mybaseW * luminosity;

  
    // The GEN selection begins here
    //--------------------------------------------------------------------------
    
    /// ---> 1) Need status 1 leptons to define the same fiducial region
    /// ---> 2) Count how many GEN leptons we have in each bin, applying the fidual region cuts
    /// ---> 3) Apply also, OF, jetbin and opposite-charged cuts.
    
    
    bool genEvent = false; 

    bool genEventNoTaus  = false;

   
    
    Float_t dileptonGenPt;
    Float_t mllGen;
    Float_t dphiGen;

    TLorentzVector leptonGen1p4;
    TLorentzVector leptonGen2p4;
    leptonGen1p4.SetPtEtaPhiM(lepGenpt1, lepGeneta1, lepGenphi1, 0.0);
    leptonGen2p4.SetPtEtaPhiM(lepGenpt2, lepGeneta2, lepGenphi2, 0.0);   

  
    Float_t Genpt1S = smear(lepGenpt1);    
    
    dileptonGenPt = (leptonGen1p4+leptonGen2p4).Pt();
    mllGen = (leptonGen1p4+leptonGen2p4).M();
    dphiGen = fabs(leptonGen1p4.DeltaPhi(leptonGen2p4));




    if ( fabs(lepGenpid1) <= 20  && 
	 fabs(lepGenpid2) <= 20  && 
    	 fabs(lepGenM1)== 24 &&  
	 fabs(lepGenM2)== 24 &&
	 fabs(lepGenpid1) != 15 &&
	 fabs(lepGenpid1) != 15 &&
	 lepGenpt1 > 20 &&  
	 lepGenpt2 > 20 &&
	 fabs(lepGenpid1) != fabs(lepGenpid2) &&  
	 (lepGenpid1*lepGenpid2) < 0 && 
	 ((fabs(lepGenpid1) == 13 && fabs(lepGeneta1) < 2.4) || 
	  (fabs(lepGenpid1) == 11 && fabs(lepGeneta1) < 2.5)) &&  
	 ((fabs(lepGenpid2) == 13 && fabs(lepGeneta2) < 2.4) || 
	  (fabs(lepGenpid2) == 11 && fabs(lepGeneta2) < 2.5)) ) 


  { 
     
      hjetEt_lepton_GEN->Fill(jetGen1_pt, totalWGen); 


	   // If jet veto at GEN level
	   //--------------------------------------------------------------------------

	   Int_t nGenJets = 0, nGenJet1 = 0, nGenJet2 = 0, nGenJet3 = 0, nGenJet4 = 0, nGenJet5 = 0; 
  
	   if ( jetGen1_pt>=30 ) nGenJet1++;
	   if ( jetGen2_pt>=30 ) nGenJet2++;
	   if ( jetGen3_pt>=30 ) nGenJet3++;
	   if ( jetGen4_pt>=30 ) nGenJet4++;
	   if ( jetGen5_pt>=30 ) nGenJet5++;
	   
	   nGenJets = nGenJet1 + nGenJet2 + nGenJet3 + nGenJet4 + nGenJet5; 
    
	   if (jetGenVeto && nGenJets < 1 )  genEvent = true;

	   
    }




    TLorentzVector lepton1p4;
    TLorentzVector lepton2p4;
    lepton1p4.SetPtEtaPhiM(pt1, eta1, phi1, 0.0);
    lepton2p4.SetPtEtaPhiM(pt2, eta2, phi2, 0.0);
   

    Float_t dileptonPt = (lepton1p4+lepton2p4).Pt();
    Float_t mll = (lepton1p4+lepton2p4).M();
    Float_t deltaphill = fabs(lepton1p4.DeltaPhi(lepton2p4));



    if  (genEvent ) {

      if ( lepGenpt1 >= 200) hPtLepton1_GEN_C->Fill(199, totalWGen);

      hPtLepton1_GEN_C->Fill(lepGenpt1, totalWGen);
      hPtLepton1_GEN->Fill(lepGenpt1, totalWGen);//*baseW*luminosity*0.00300652); // leading pt ---> which pt should I store here? 
    
      hDilepton_GEN->Fill(dileptonGenPt,totalWGen); // ptll 
    
      hmll_GEN->Fill(mllGen,totalWGen); // mll
      
      hdphi_GEN->Fill(dphiGen,totalWGen); // deltaPhi

      hjetEt_GEN->Fill(jetGen1_pt, totalWGen); 

      hInclusive_GEN->Fill(1, totalWGen);
    
    
    // The RECO selection begins here. The RECO leptons are supposed to pass the ID+ISO selection already? 
    //----------------------------------------------------------------------------------------------------
  

 
    hPtLepton1_RECO_GEN->Fill(pt1, lepGenpt1, totalW);
    hPtLepton1_RECO->Fill(pt1, totalWReco);

    hDilepton_RECO->Fill(dileptonPt, totalW);
    hDilepton_RECO_GEN->Fill(dileptonPt, dileptonGenPt, totalW);

    hmll_RECO->Fill(mll, totalW);
    hmll_RECO_GEN->Fill(mll, mllGen , totalW);

    hdphi_RECO->Fill(dphill, totalW);
    hdphi_RECO_GEN->Fill(dphill,dphiGen , totalW);
 
    hjetEt_RECO->Fill(jetpt1, totalW);
    hjetEt_RECO_GEN->Fill(jetpt1, jetGen1_pt, totalW);
    }

    // ---->>> Fill Response matrix to do only unfolding on resolution
    //         without any selection efficiency propagated, but fiducial region as reco one.

    responsePtLepton1GEN.Fill(pt1, lepGenpt1, totalWReco);
    responseMllGEN.Fill(mll, mllGen, totalWReco);

 

    /// ---> 3) Going to apply the selection analysis cuts on RECO objects

    Int_t dphiv = (njet <= 1 || (njet > 1 && dphilljetjet < 165.*TMath::DegToRad()));
    
    Float_t metvar = (njet <= 1) ? mpmet : pfmet;
    
    Float_t jetbin = njet;
    
    Float_t dyMVA = ( !sameflav || ( (njet!=0 || dymva1>0.88) && (njet!=1 || dymva1>0.84) && ( njet==0 || njet==1 || (pfmet > 45.0)) ) );


    Bool_t isMatched = true; 
    Bool_t isMatchedGEN = true; 
    Bool_t isMatchedRECO = true; 
    Float_t deltaR = 9999.9;
   
    hgenJetEt->Fill(jetGen1_pt );
    
    

    if(  pt1 > 20 && pt2 > 20  && !sameflav && ch1*ch2 < 0 &&
	 trigger == 1                        &&
	 nextra == 0                         && 
	 pfmet > 20                          &&
	 mll > 12                            &&
	 (zveto==1 || !sameflav)             &&
	 (mpmet > 20  && dyMVA)              &&
	 (dphiv || !sameflav)                &&
	 bveto_mu                            &&
	 ptll>30 && (!sameflav || ptll>45 )  &&
	 jetbin == jetChannel                && 
	 (bveto_ip==1 &&  nbjettche==0)   
	 ){
     
	//hgenJetEt->Fill(jetGen1_pt );
	hrecoJetEt->Fill(jetpt1);
	hJetEt_Gen_Reco->Fill(jetGen1_pt, jetpt1);

	//	cout << jetGen1_pt << "  "  <<jetpt1  << endl;

	// if (true) {	     




      // DEFINE MATCHING RECO - GEN
      //--------------------------------------------------------------------------
    
	//if ( lepton1p4.DeltaR(leptonGen1p4) >= 0.15 )isMatchedGEN = false;//&& lepton1p4.DeltaR(leptonGen2p4) >= 0.15 )  isMatchedGEN = false;

	if ( lepton1p4.DeltaR(leptonGen1p4) >= 0.15 && lepton1p4.DeltaR(leptonGen2p4) >= 0.15 )  isMatchedGEN = false;
	
	if ( leptonGen1p4.DeltaR(lepton1p4) >= 0.15 && leptonGen1p4.DeltaR(lepton2p4) >= 0.15 )  isMatchedRECO = false; 
   


	Float_t pt1S = smear(pt1);

	//Float_t pt1S = linearW(pt1, 0.006);
	
	
	//cout << totalWReco  << endl;

	
	hPtLepton1WWLevel_RECO->Fill(pt1,    totalWReco); 
	hPtLepton1WWLevel_GEN->Fill(lepGenpt1,totalW);
	hPtLepton1WWLevel_RECO_GEN->Fill(pt1,lepGenpt1, totalW);
	
	hDileptonWWLevel_RECO->Fill(dileptonPt,        totalWReco); 
	hDileptonWWLevel_GEN->Fill(dileptonGenPt,totalW);
	hDileptonWWLevel_RECO_GEN->Fill(dileptonPt,dileptonGenPt , totalW);
	
	hmllWWLevel_RECO->Fill(mll,        totalWReco); 
	hmllWWLevel_GEN->Fill(mllGen,totalW);
	hmllWWLevel_RECO_GEN->Fill(mll , mllGen, totalW);
	
	hdphiWWLevel_RECO->Fill(dphill,        totalWReco); 
	hdphiWWLevel_GEN->Fill(dphiGen,totalW);
	hdphiWWLevel_RECO_GEN->Fill(dphill,dphiGen, totalW);
	
	hjetEtWWLevel_RECO->Fill(jetpt1, totalWReco);
	hjetEtWWLevel_GEN->Fill( jetGen1_pt, totalWReco);
	hjetEtWWLevel_RECO_GEN->Fill(jetpt1, jetGen1_pt, totalW);

	hInclusiveWWLevel_GEN->Fill(1, totalWGen);
	hInclusiveWWLevel_RECO->Fill(1,   totalWReco);

	
	//---- Fill response matrix ( we have always with our selection 2 gen leptons and 2 reco leptons)
	if  (genEvent) {

	  responsePtLepton1.Fill(pt1, lepGenpt1, totalWReco);
	  responseDilepton.Fill(dileptonPt, dileptonGenPt, totalWReco);
	  responseMll.Fill(mll, mllGen, totalWReco);
	  responseDphi.Fill(dphill, dphiGen, totalWReco);
	  responseInclusive.Fill(1, 1, totalWReco);

	  responsePtLepton1.Miss(lepGenpt1, (1-efficiencyW)*totalWGen );
	  responseDilepton.Miss(dileptonGenPt, (1-efficiencyW)*totalWGen);
	  responseMll.Miss(mllGen, (1-efficiencyW)*totalWGen);
	  responseDphi.Miss(dphiGen, (1-efficiencyW)*totalWGen);
	  responseInclusive.Miss(1, (1-efficiencyW)*totalWGen);
	  
	  hjetEt_fid_GEN->Fill(jetGen1_pt, totalWGen); 

	} else {
	  
	  responsePtLepton1.Fake(pt1, totalWReco);
	  responseDilepton.Fake(dileptonPt, totalWReco);
	  responseMll.Fake(mll,  totalWReco);
	  responseDphi.Fake(dphill, totalWReco);
	  responseInclusive.Fake(1, totalWReco);

	  hjetEt_Nofid_GEN->Fill(jetGen1_pt, totalWGen); 
	}
   
    } else {	

      if  (genEvent) {
	hPtLepton1WWLevel_nonselected->Fill(lepGenpt1, totalW);
	
	//---- Fill response matrix with efficiency 
	
	responsePtLepton1.Miss(lepGenpt1, totalWGen);
	responseDilepton.Miss(dileptonGenPt, totalWGen);
	responseMll.Miss(mllGen, totalWGen);
	responseDphi.Miss(dphiGen, totalWGen);

	responseInclusive.Miss(1, totalWGen);
   
      }

    }
	   





    /*
	//---- Fill response matrix 	
	
	if (isMatchedGEN ) {
	  
	  responsePtLepton1.Fill(pt1, lepGenpt1, totalWReco);
	  responseDilepton.Fill(dileptonPt, dileptonGenPt, totalWReco);
	  responseMll.Fill(mll, mllGen, totalWReco);
	  responseDphi.Fill(dphill, dphiGen, totalWReco);
	  
	  //      }	else if (!isMatchedGEN ) { 


	} else { 
	   
	  responsePtLepton1.Fake(pt1, totalWReco);
	  responseDilepton.Fake(dileptonPt, totalWReco);
	  responseMll.Fake(mll,  totalWReco);
	  responseDphi.Fake(dphill, totalWReco);
	  
	  if (!isMatchedRECO) responsePtLepton1.Miss(lepGenpt1, totalWGen);
	  
	  //cout << "Not found matched GEN !!! " << endl;
	    
	  //	responsePtLepton1.Miss(lepGenpt1, totalWGen);
	        
	} 

	 } else {	
 
	hPtLepton1WWLevel_nonselected->Fill(lepGenpt1, totalW);
	
	//---- Fill response matrix with efficiency 
	
	responsePtLepton1.Miss(lepGenpt1, totalWGen);
	responseDilepton.Miss(dileptonGenPt, totalWGen);
	responseMll.Miss(mllGen, totalWGen);
	responseDphi.Miss(dphiGen, totalWGen);
	 }
	 
*/

    
    }
 
  // Save the histograms
  //----------------------------------------------------------------------------
  output->cd();
  responsePtLepton1GEN.Write();
  responseMllGEN.Write();
  responseJetEtGEN.Write();
  responsePtLepton1.Write();
  responseDilepton.Write();
  responseMll.Write();
  responseDphi.Write();
  responseInclusive.Write();
  output->Write("", TObject::kOverwrite);
  output->Close();



  // Define binning 

  


  
}





//==============================================================================
// Gaussian smearing, systematic translation, and variable inefficiency
//==============================================================================

Float_t smear (Float_t xt)
{

  Float_t cutdummy= -99999.0;

  Float_t xeff= 0.3 + (1.0-0.3)/20*(xt+10.0);  // efficiency
  Float_t x= gRandom->Rndm();
  //  if (x>xeff) return cutdummy;
  Float_t xsmear= gRandom->Gaus(25,10);     // bias and smear
  return xt+xsmear;
}


//==============================================================================
// Linear Weight
//==============================================================================

Double_t linearW (Double_t xt, Double_t slope)
{

  Double_t weight = 1;

  Double_t inter = 1+ (xt-50)*slope; 

  if ( inter*inter > 0.1) { 

    weight = inter*inter ;
  } else {
    weight = 0.1;
  }
  
     
  return xt*weight;
}






