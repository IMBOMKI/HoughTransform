#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TTree.h"
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
#include <cmath>
#include <tuple>

Float_t ConfTransX(Float_t x, Float_t y);
Float_t ConfTransY(Float_t x, Float_t y);

void HoughTransform(){
  
  TCanvas *c1 = new TCanvas("c1", "Even", 600, 600);

  TCanvas *c2 = new TCanvas("c2", "Odd", 600, 600);

  TCanvas *c3 = new TCanvas("c3", "Truth", 600, 600);
  

  std::string dir = "./";
  std::string fileName = "trig_ep92.root";

  TFile *f = TFile::Open(TString(dir+fileName));
  TTree *t = (TTree*)f->Get("trdata");

  std::string outputdir = "../TrackFinding/";
  std::string outputfileName = "trig_"+fileName;
  
  TFile *f_out = new TFile(TString(outputdir+outputfileName), "recreate");
  TTree *t_out = t->CloneTree(0);

  Int_t nCDCHit;
  Float_t CDCHitX[100000];
  Float_t CDCHitY[100000];
  Float_t CDCHitZ[100000];
  Float_t CDCHitT[100000];
  Float_t CDCEDep[100000];
  Int_t nCALCDCHit;
  Float_t CDCDriftDist[100000];
  Int_t CDCCharge[100000];
  Float_t WireEnd0X[100000];
  Float_t WireEnd0Y[100000];  
  Float_t WireEnd0Z[100000];
  Float_t WireEnd1X[100000];
  Float_t WireEnd1Y[100000];  
  Float_t WireEnd1Z[100000];
  Int_t WireLayerId[100000];
  Int_t WireMaxLayerId;

  Float_t genTrE;
  Float_t genTrPx;
  Float_t genTrPy;
  Float_t genTrPz;
  Float_t genTrT;

  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[1000];
  Int_t TrigCRK[1000];
  Float_t TrigSTLTime[1000];
  Float_t TrigCRKTime[1000];
  Int_t MaxCDCLayerId;
  Int_t nCDCHitPrimary;
  Bool_t ifSingleTurn;
  Bool_t ifMultiTurn;

  t->SetBranchAddress("nCDCHit", &nCDCHit);
  t->SetBranchAddress("CDCHitX", CDCHitX);  
  t->SetBranchAddress("CDCHitY", CDCHitY);  
  t->SetBranchAddress("CDCHitZ", CDCHitZ); 
  t->SetBranchAddress("CDCEDep", CDCEDep);   
  t->SetBranchAddress("nCALCDCHit", &nCALCDCHit);
  t->SetBranchAddress("CDCDriftDist", CDCDriftDist);
  t->SetBranchAddress("CDCCharge", CDCCharge);
  t->SetBranchAddress("WireEnd0X", WireEnd0X);
  t->SetBranchAddress("WireEnd0Y", WireEnd0Y);
  t->SetBranchAddress("WireEnd0Z", WireEnd0Z);
  t->SetBranchAddress("WireEnd1X", WireEnd1X);
  t->SetBranchAddress("WireEnd1Y", WireEnd1Y);
  t->SetBranchAddress("WireEnd1Z", WireEnd1Z);
  t->SetBranchAddress("WireLayerId", WireLayerId);
  t->SetBranchAddress("WireMaxLayerId", &WireMaxLayerId);

  t->SetBranchAddress("genTrE", &genTrE);
  t->SetBranchAddress("genTrPx", &genTrPx);
  t->SetBranchAddress("genTrPy", &genTrPy);
  t->SetBranchAddress("genTrPz", &genTrPz);
  t->SetBranchAddress("genTrT", &genTrT);

  t->SetBranchAddress("nTrigSTL", &nTrigSTL);
  t->SetBranchAddress("nTrigCRK", &nTrigCRK);    
  t->SetBranchAddress("TrigSTL", TrigSTL);
  t->SetBranchAddress("TrigCRK", TrigCRK);
  t->SetBranchAddress("TrigSTLTime", TrigSTLTime);
  t->SetBranchAddress("TrigCRKTime", TrigCRKTime);
  t->SetBranchAddress("MaxCDCLayerId", &MaxCDCLayerId);
  t->SetBranchAddress("ifSingleTurn", &ifSingleTurn);
  t->SetBranchAddress("ifMultiTurn", &ifMultiTurn);
  t->SetBranchAddress("nCDCHitPrimary", &nCDCHitPrimary);




  TH1F *fitpT_hist  = new TH1F("fitpT","fitpT",100,60,110) ;
  TH1F *truthpT_hist = new TH1F("truthpT","truthpT",100,60,110) ;
  TH1F *diff_hist = new TH1F("diff","truthpT-fitpT",300,-100,50) ;
    


  for (Int_t i=0; i < t->GetEntries(); i++){
    t->GetEntry(i);

    if (nCALCDCHit>=30 && WireMaxLayerId>=4){
    Float_t EvenConfX[100000];
    Float_t EvenConfY[100000];
    Float_t OddConfX[100000];
    Float_t OddConfY[100000];
    Int_t EvenHit=0;
    Int_t OddHit=0;

    for (Int_t nHits=0; nHits<nCALCDCHit ; nHits++){

      if (WireLayerId[nHits]%2==0){
	EvenConfX[EvenHit] = ConfTransX(WireEnd0X[EvenHit],WireEnd0Y[EvenHit]);
	EvenConfY[EvenHit] = ConfTransY(WireEnd0X[EvenHit],WireEnd0Y[EvenHit]);
	EvenHit++;
      }
      else if (WireLayerId[nHits]%2==1){
	OddConfX[OddHit] = ConfTransX(WireEnd0X[OddHit],WireEnd0Y[OddHit]);
	OddConfY[OddHit] = ConfTransY(WireEnd0X[OddHit],WireEnd0Y[OddHit]);
	OddHit++;
      }
      
    }


    TGraph *grEven = new TGraph(EvenHit, EvenConfX, EvenConfY);
    TGraph *grOdd = new TGraph(OddHit, OddConfX, OddConfY);
    TGraph *grTruth = new TGraph(nCALCDCHit, WireEnd0X, WireEnd0Y);


    TF1 *fitEven; TF1 *fitOdd;
    Float_t offsetEven; Float_t offsetOdd;
    Float_t slopeEven; Float_t slopeOdd;

    ///////////////// Conformal Fitting /////////////////////

    /**** Even Layer ****/
    c1->cd();
    grEven->SetTitle("Even");
    grEven->SetMarkerStyle(20);
 
    grEven->Fit("pol1");
    fitEven = grEven->GetFunction("pol1");
    offsetEven = fitEven->GetParameter(0);
    slopeEven  = fitEven->GetParameter(1);
    //grEven->Draw("AP");

    /**** Odd Layer ****/
    c2->cd();
    grOdd->SetTitle("Odd");
    grOdd->SetMarkerStyle(20);

    grOdd->Fit("pol1");
    fitOdd = grOdd->GetFunction("pol1");
    offsetOdd = fitOdd->GetParameter(0);
    slopeOdd  = fitOdd->GetParameter(1);
    //grOdd->Draw("AP");

    Float_t slopeAvg=(slopeEven+slopeOdd)/2;
    Float_t offsetAvg=(offsetEven+offsetOdd)/2;
    Float_t fitrad=sqrt(pow(slopeAvg/(-2*offsetAvg),2)+pow(1/(2*offsetAvg),2));
    Float_t fitpT = 10*fitrad/3.3356;
    Float_t truthpT = sqrt(pow(genTrPy,2)+pow(genTrPz,2));
    
    fitpT_hist->Fill(fitpT);
    truthpT_hist->Fill(truthpT);
    diff_hist->Fill(truthpT-fitpT);

    ////////////////////////////////////////////////////////





    ///////////////// MC Truth Drawing /////////////////////

    c3->cd();
    c3->DrawFrame(-90,-90,90,90);

    Float_t WireRad[18] = {53.00,54.60,56.20,57.80,59.40,61.00,62.60,64.20,65.80,67.40,69.00,70.60,72.20,73.80,75.40,77.00,78.60,80.20};
    for (Int_t i=0; i<18; i++){
      TEllipse *WireCircle = new TEllipse(0,0,WireRad[i],WireRad[i]);
      WireCircle->SetFillColor(0);
      WireCircle->SetFillStyle(4000);
      WireCircle->SetLineColor(40);
      //WireCircle->Draw();
    }

    Float_t DiskRad=10.00;
    TEllipse *Disk = new TEllipse(0,0,DiskRad,DiskRad);
    Disk->SetFillColor(40);
    Disk->SetLineColor(40);
    //Disk->Draw();  

    grTruth->SetTitle("Truth");
    grTruth->SetMarkerStyle(20);
    grTruth->SetMarkerSize(1);
    //grTruth->Draw("P");
 
    std::cout << "Even Layer NDF: " << EvenHit << std::endl;
    std::cout << " Odd Layer NDF: " << OddHit << std::endl;
    std::cout << "  MC Truth NDF: " << nCALCDCHit << std::endl;
    std::cout << i <<"-th Event"<< std::endl;

    }
  }

  c1->cd();
  fitpT_hist->Draw();
  c2->cd();
  truthpT_hist->Draw();
  c3->cd();
  diff_hist->Draw();

  std::cout << "Finish!" << std::endl;
  f->Close();
  
}

Float_t ConfTransX(Float_t x, Float_t y){
  return x/(pow(x,2)+pow(y,2));
} 

Float_t ConfTransY(Float_t x, Float_t y){
  return y/(pow(x,2)+pow(y,2));
} 





