#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TLegend.h"
#include "THStack.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <tuple>


void AnalyzeOutput(){
  Int_t NumOfEvt=100000;
  Int_t bw_ini=3;
  Int_t bw_fin=10;
  Double_t bw[bw_fin-bw_ini+1];
  Double_t Acceptance_em[bw_fin-bw_ini+1];
  Double_t Acceptance_ep[bw_fin-bw_ini+1];
  Double_t RecoRate_em[bw_fin-bw_ini+1];
  Double_t RecoRate_ep[bw_fin-bw_ini+1];

  memset (Acceptance_em, 0, sizeof (Acceptance_em));
  memset (Acceptance_ep, 0, sizeof (Acceptance_ep));
  memset (RecoRate_em, 0, sizeof (RecoRate_em));
  memset (RecoRate_ep, 0, sizeof (RecoRate_ep));


  ///////////////////////////////////////////////////////////////////////////////
  //
  //
  // e- Part
  //
  //
  ///////////////////////////////////////////////////////////////////////////////

  for (Int_t i_bw=bw_ini; i_bw<bw_fin+1; i_bw++){
    
    bw[i_bw-bw_ini] = i_bw;
    std::string dir_em = "./FindedEvents/";
    std::string fileName_em = "bw"+std::to_string(i_bw)+"_find_trig_em104_onlyPrimary.root";
    TFile *f_em = TFile::Open(TString(dir_em+fileName_em));
    TTree *t_em = (TTree*)f_em->Get("trdata");

    Int_t nCALCDCHit;
    Int_t nRecoHit;
    Int_t RecoMaxWireLayerId;
    Bool_t Reco_ifCL3;
    Bool_t ifSingleTurn;
    Bool_t ifMultiTurn;
    Double_t fittedR_even;
    Double_t fittedR_odd;
    Double_t drEvenToOdd;

    t_em->SetBranchAddress("nCALCDCHit", &nCALCDCHit);
    t_em->SetBranchAddress("nRecoHit", &nRecoHit);
    t_em->SetBranchAddress("RecoMaxWireLayerId", &RecoMaxWireLayerId);
    t_em->SetBranchAddress("Reco_ifCL3", &Reco_ifCL3);
    t_em->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
    t_em->SetBranchAddress("ifMultiTurn", &ifMultiTurn);
    t_em->SetBranchAddress("fittedR_even", &fittedR_even);
    t_em->SetBranchAddress("fittedR_odd", &fittedR_odd);
    t_em->SetBranchAddress("drEvenToOdd", &drEvenToOdd);

    Int_t TotalHit=0;
    Int_t Total_RecoHit=0;

    for (Int_t i_evt=0; i_evt<t_em->GetEntries(); i_evt++){      
      t_em->GetEntry(i_evt);
      
      // Acceptance

      if (Reco_ifCL3==1 && RecoMaxWireLayerId>=4 && nRecoHit>=30 && drEvenToOdd<14 && fittedR_even<40 && fittedR_odd<40){
	Acceptance_em[i_bw-bw_ini]++;
      }

      // Reco Rate

      TotalHit+=nCALCDCHit;
      Total_RecoHit+=nRecoHit;
    }

    Acceptance_em[i_bw-bw_ini]=Acceptance_em[i_bw-bw_ini]/NumOfEvt*100;
    RecoRate_em[i_bw-bw_ini]=Double_t(Total_RecoHit)/(TotalHit)*100.0;

    f_em->cd();
    f_em->Close();
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  //
  // e+ Part
  //
  //
  ///////////////////////////////////////////////////////////////////////////////
  
  for (Int_t i_bw=bw_ini; i_bw<bw_fin+1; i_bw++){
    
    bw[i_bw-bw_ini] = i_bw;
    std::string dir_ep = "./FindedEvents/";
    std::string fileName_ep = "bw"+std::to_string(i_bw)+"_find_trig_ep92_onlyPrimary.root";
    TFile *f_ep = TFile::Open(TString(dir_ep+fileName_ep));
    TTree *t_ep = (TTree*)f_ep->Get("trdata");

    Int_t nCALCDCHit;
    Int_t nRecoHit;
    Int_t RecoMaxWireLayerId;
    Bool_t Reco_ifCL3;
    Bool_t ifSingleTurn;
    Bool_t ifMultiTurn;
    Double_t fittedR_even;
    Double_t fittedR_odd;
    Double_t drEvenToOdd;

    t_ep->SetBranchAddress("nCALCDCHit", &nCALCDCHit);
    t_ep->SetBranchAddress("nRecoHit", &nRecoHit);
    t_ep->SetBranchAddress("RecoMaxWireLayerId", &RecoMaxWireLayerId);
    t_ep->SetBranchAddress("Reco_ifCL3", &Reco_ifCL3);
    t_ep->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
    t_ep->SetBranchAddress("ifMultiTurn", &ifMultiTurn);
    t_ep->SetBranchAddress("fittedR_even", &fittedR_even);
    t_ep->SetBranchAddress("fittedR_odd", &fittedR_odd);
    t_ep->SetBranchAddress("drEvenToOdd", &drEvenToOdd);

    Int_t TotalHit=0;
    Int_t Total_RecoHit=0;

    for (Int_t i_evt=0; i_evt<t_ep->GetEntries(); i_evt++){      
      t_ep->GetEntry(i_evt);
      
      // Acceptance

      if (Reco_ifCL3==1 && RecoMaxWireLayerId>=4 && nRecoHit>=30 && drEvenToOdd<14 && fittedR_even<45 && fittedR_odd<45){
	Acceptance_ep[i_bw-bw_ini]++;
      }

      // Reco Rate

      TotalHit+=nCALCDCHit;
      Total_RecoHit+=nRecoHit;
    }

    Acceptance_ep[i_bw-bw_ini]=Acceptance_ep[i_bw-bw_ini]/NumOfEvt*100;
    RecoRate_ep[i_bw-bw_ini]=Double_t(Total_RecoHit)/(TotalHit)*100.0;

    std::cout <<      Acceptance_em[i_bw-bw_ini] << "   " << RecoRate_em[i_bw-bw_ini] << "   " <<  Acceptance_ep[i_bw-bw_ini] << "   " << RecoRate_ep[i_bw-bw_ini] << std::endl;


    f_ep->cd();
    f_ep->Close();
  }
  

  TCanvas *c_div = new TCanvas("canvas", "canvas", 1000, 500);
  c_div->Divide(2,1);
  TPad *pad1_1 = new TPad("pad1","",0,0,1,1);
  TPad *pad1_2 = new TPad("pad2","",0,0,1,1);
  pad1_1->SetFillStyle(4000);
  pad1_1->SetFrameFillStyle(0);
  pad1_2->SetFillStyle(4000);
  pad1_2->SetFrameFillStyle(0);

  TPad *pad2_1 = new TPad("pad1","",0,0,1,1);
  TPad *pad2_2 = new TPad("pad2","",0,0,1,1);
  pad2_2->SetFillStyle(4000);
  pad2_2->SetFrameFillStyle(0);

  TMultiGraph *mg1 = new TMultiGraph(); 
  TGraph *grAcceptance_em= new TGraph(bw_fin-bw_ini+1, bw, Acceptance_em);
  TGraph *grRecoRate_em= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_em);
  TGraph *grAcceptance_ep= new TGraph(bw_fin-bw_ini+1, bw, Acceptance_ep);
  TGraph *grRecoRate_ep= new TGraph(bw_fin-bw_ini+1, bw, RecoRate_ep);

  c_div->cd(1);
  pad1_1->Draw();
  pad1_1->cd();
  grAcceptance_em->SetTitle("Acceptance");
  grAcceptance_em->SetMarkerStyle(20);
  grAcceptance_em->SetMarkerColor(4);
  grAcceptance_em->GetYaxis()->SetRangeUser(9,20);
  grAcceptance_em->GetYaxis()->SetAxisColor(4);
  grAcceptance_em->GetYaxis()->SetTitle("Acceptance of e- (%)");
  grAcceptance_em->GetXaxis()->SetTitle("Band Width of Circle");
  grAcceptance_em->Draw("ALP");
  
  pad1_2->Draw();
  pad1_2->cd();
  grAcceptance_ep->SetMarkerStyle(20);
  grAcceptance_ep->SetMarkerColor(2);
  grAcceptance_ep->GetYaxis()->SetRangeUser(2,3.5);
  grAcceptance_ep->GetYaxis()->SetAxisColor(2);
  grAcceptance_ep->GetYaxis()->SetTitle("Acceptance of e+ (%)");
  grAcceptance_ep->Draw("LPAY+");

  ////////////////////////////////////////////////////////////////////

  c_div->cd(2);
  pad2_1->Draw();
  pad2_1->cd();
  grRecoRate_em->SetTitle("Signal Hits Recognition Rate");
  grRecoRate_em->SetMarkerStyle(20);
  grRecoRate_em->SetMarkerColor(38);
  grRecoRate_em->GetYaxis()->SetRangeUser(93,99);
  grRecoRate_em->GetYaxis()->SetTitle("Recognition Rate (%)");
  grRecoRate_em->GetXaxis()->SetTitle("Band Width of Circle");
  grRecoRate_em->Draw("ALP");
  
  pad2_2->Draw();
  pad2_2->cd();
  grRecoRate_ep->SetTitle("Signal Hits Recognition Rate");
  grRecoRate_ep->SetMarkerStyle(20);
  grRecoRate_ep->SetMarkerColor(46);
  grRecoRate_ep->GetYaxis()->SetRangeUser(93,99);
  grRecoRate_ep->Draw("LPA");

  //////////////////////////////////
  

}

