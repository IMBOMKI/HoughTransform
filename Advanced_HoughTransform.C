#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TLegend.h"
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
#include <cmath>
#include <tuple>

Float_t ConfTransX(Float_t x, Float_t y);
Float_t ConfTransY(Float_t x, Float_t y);
Float_t HoughTrans(Float_t x, Float_t y, Float_t theta);
Bool_t ifInsideDisk(Float_t x, Float_t y);
std::vector<std::pair<Float_t,Float_t> > makeOrigins(Float_t rad, std::pair<Float_t,Float_t> ref);
Int_t findMaxpoint(std::vector<Int_t> vec);

void Advanced_HoughTransform(){
  std::string dir = "./";
  std::string fileName = "trig_ep92.root";
  TFile *f = TFile::Open(TString(dir+fileName));
  TTree *t = (TTree*)f->Get("trdata");
  std::string outputdir = "../TrackFinding/";
  std::string outputfileName = "fit_"+fileName;
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

  TCanvas *c_div = new TCanvas ("c_div", "c divided", 600,600);
  c_div->Divide(2,4);
  TH1F *fitpT_hist  = new TH1F("fitpT","fitpT",100,60,110) ;
  TH1F *truthpT_hist = new TH1F("truthpT","truthpT",100,60,110) ;
  TH1F *diff_hist = new TH1F("diff","truthpT-fitpT",300,-100,50) ;
  TH1F *maxvote_even_hist = new TH1F("maxvote_even_hist","max_even_vote_hist", 100, 0, 100);
  TH1F *maxvote_odd_hist = new TH1F("maxvote_odd_hist","max_odd_vote_hist", 100, 0, 100);
  TH1F *rad_even_hist = new TH1F("rad_even_hist","rad_even_hist", 300, 0, 300);
  TH1F *rad_odd_hist = new TH1F("rad_odd_hist","rad_odd_hist", 300, 0, 300);

  TCanvas *c_hits = new TCanvas("c_hits", "c_hits", 600,600);

  for (Int_t i_evt=0; i_evt < 30; i_evt++){

    t->GetEntry(i_evt);    

    Float_t ConfX[100000];
    Float_t ConfY[100000];

    Int_t nBins=250;
    Float_t nPt=1000.0;
    Float_t rhomax=0.02;
    Float_t rhomin=-0.02;

    Float_t rad_uncertainty=5;
    std::pair <Float_t,Float_t> ref = std::make_pair(0,0);
    //std::pair <Float_t,Float_t> ref_odd = std::make_pair(0,0);
    //std::pair <Float_t,Float_t> ref_even = std::make_pair(0,0);

    Int_t rho_even_index; Int_t deg_even_index;
    Int_t rho_odd_index; Int_t deg_odd_index;    
    Int_t do_not_use;


    if (nCALCDCHit>=30 && WireMaxLayerId>=4 && (ifSingleTurn==1 || ifMultiTurn==1)){ 

      for (Int_t it2=0; it2<3; it2++){

	std::vector<std::pair<Float_t,Float_t> > origins = makeOrigins(rad_uncertainty,ref);     
	//std::vector<std::pair<Float_t,Float_t> > origins_odd = makeOrigins(rad_uncertainty,ref_odd);     
	//std::vector<std::pair<Float_t,Float_t> > origins_even = makeOrigins(rad_uncertainty,ref_even);     
	std::vector<Int_t> vote_sum;
	//std::vector<Int_t> maxvote_even;
	//std::vector<Int_t> maxvote_odd;
	std::vector<std::pair<Int_t,Int_t> > even_index_vec; 
	std::vector<std::pair<Int_t,Int_t> > odd_index_vec;   

	
	for (Int_t it=0; it<origins.size(); it++){
	  
	  Float_t oriX=origins[it].first;
	  Float_t oriY=origins[it].second;

	  if  (ifInsideDisk(oriX,oriY)==0){
	    vote_sum.push_back(0);
	    //maxvote_even.push_back(0);
	    //maxvote_odd.push_back(0);
	  }

	  else if (ifInsideDisk(oriX,oriY)==1){
	    
	    /*---------------------------------------
	      |                                        |
	      |      Conformal Transformation          |
	      |                                        |
	      -----------------------------------------*/
	    
	  for (Int_t nHits=0; nHits<nCALCDCHit ; nHits++){   
	    ConfX[nHits] = ConfTransX(WireEnd0X[nHits]-oriX,WireEnd0Y[nHits]-oriY);
	    ConfY[nHits] = ConfTransY(WireEnd0X[nHits]-oriX,WireEnd0Y[nHits]-oriY);
	  }
	    
	    /*-------------------------------------------
	      |                                         |
	      |      Hough Transformation & Voting      |
	      |                                         |
	      ------------------------------------------*/
	    
	  Bool_t if_already_vote[nBins][nBins][2];
	  TNtuple * hough = new TNtuple("hough","hough", "rho:theta:hitId:wireId:is_even");
	  TH2F *vote_even = new TH2F("vote_even", "vote_even", nBins, 0, nBins, nBins, 0, nBins);
	  TH2F *vote_odd = new TH2F("vote_odd", "vote_odd", nBins, 0, nBins, nBins, 0, nBins);

	  for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	    memset(if_already_vote,0,sizeof(if_already_vote));
	    for (Int_t i_Pt=0 ; i_Pt<nPt; i_Pt++){
	      
	      Float_t deg = i_Pt/nPt*180.0;
	      Float_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);
	      Bool_t is_even = WireLayerId[i_hit]%2;	  
	      Int_t deg_index = (deg)*nBins/180;
	      Int_t rho_index = (rho+rhomax)*nBins/(rhomax-rhomin);
	      
	      if (is_even==1 && if_already_vote[deg_index][rho_index][0]==0){
		vote_even->Fill(deg_index,rho_index);
		if_already_vote[deg_index][rho_index][0]=1;
	      }
	      
	      else if (is_even==0 && if_already_vote[deg_index][rho_index][1]==0){
		vote_odd->Fill(deg_index,rho_index);
		if_already_vote[deg_index][rho_index][1]=1;
	      }
	    }
	  }
	  
	  Int_t sum = vote_even->GetMaximum()+vote_odd->GetMinimum();
	  vote_sum.push_back(sum);
	  //maxvote_even.push_back(vote_even->GetMaximum());
	  //maxvote_odd.push_back(vote_odd->GetMaximum());

	  vote_even->GetMaximumBin(deg_even_index, rho_even_index, do_not_use);
	  vote_odd->GetMaximumBin(deg_odd_index, rho_odd_index,do_not_use);

	  even_index_vec.push_back(std::make_pair(deg_even_index,rho_even_index));
	  odd_index_vec.push_back(std::make_pair(deg_odd_index,rho_odd_index));	  

	  delete vote_even;
	  delete vote_odd;
	  }
	}
	
	Int_t maxIndex=findMaxpoint(vote_sum);
	//Int_t maxeven_Index=findMaxpoint(maxvote_even);
	//Int_t maxodd_Index=findMaxpoint(maxvote_odd);

	ref.first = origins[maxIndex].first; ref.second = origins[maxIndex].second;
	rad_uncertainty = rad_uncertainty/2.0;
	deg_even_index = even_index_vec[maxIndex].first; rho_even_index = even_index_vec[maxIndex].second;
	deg_odd_index = odd_index_vec[maxIndex].first; rho_odd_index = odd_index_vec[maxIndex].second;

	//deg_even_index = even_index_vec[maxeven_Index].first; 
	//rho_even_index = even_index_vec[maxeven_Index].second;
	//deg_odd_index = odd_index_vec[maxodd_Index].first; 
	//rho_odd_index = odd_index_vec[maxodd_Index].second;
      }
      
      //std::cout << ref.first << "   " << ref.second << std::endl;
      
      //////////////////////////////////////////////////////////////////////
      /*
      Int_t rho_even_index; Int_t deg_even_index;
      Int_t rho_odd_index; Int_t deg_odd_index;    
      Int_t do_not_use;
      vote_even->GetMaximumBin(rho_even_index, deg_even_index, do_not_use);
      vote_odd->GetMaximumBin(rho_odd_index, deg_odd_index,do_not_use);
      */
      Float_t deg_even_eval = 180.0/nBins*(2*deg_even_index+1)/2;
      Float_t rho_even_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_even_index+1)/2;
      Float_t cX_even = TMath::Cos(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      Float_t cY_even = TMath::Sin(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      
      Float_t deg_odd_eval = 180.0/nBins*(2*deg_odd_index+1)/2;
      Float_t rho_odd_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_odd_index+1)/2;
      Float_t cX_odd = TMath::Cos(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
      Float_t cY_odd = TMath::Sin(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);

      Float_t rad_even = sqrt(pow(cX_even,2)+pow(cY_even,2));
      Float_t rad_odd = sqrt(pow(cX_odd,2)+pow(cY_odd,2));
      rad_even_hist->Fill(rad_even);
      rad_odd_hist->Fill(rad_odd);

      std::cout << rho_even_index<< std::endl;


      /*
      Int_t maxvote_even =vote_even->GetMaximum();
      Int_t maxvote_odd =vote_odd->GetMaximum();    
      maxvote_even_hist->Fill(maxvote_even);
      maxvote_odd_hist->Fill(maxvote_odd);  
            
      c_div->cd(5);
      vote_even->Draw("colz");
      c_div->cd(6);
      vote_odd->Draw("colz");
      */

      /*-----------------------------------------
	|                                         |
	|   Wire Hit Drawing with Hough Circle    |
	|                                         |
	------------------------------------------*/
      
      c_hits->cd();
      c_hits->DrawFrame(-90,-90,90,90);
      
      // Wire Layers
      
      Float_t WireRad[18] = {53.00,54.60,56.20,57.80,59.40,61.00,62.60,64.20,65.80,67.40,69.00,70.60,72.20,73.80,75.40,77.00,78.60,80.20};
      
      for (Int_t i=0; i<18; i++){
	TEllipse *WireCircle = new TEllipse(0,0,WireRad[i],WireRad[i]);
	WireCircle->SetFillColor(0);
	WireCircle->SetFillStyle(4000);
	WireCircle->SetLineColor(40);
	WireCircle->Draw();
      }
      
      // Stopping target
      
      Float_t DiskRad=10.00;
      TEllipse *Disk = new TEllipse(0,0,DiskRad,DiskRad);
      Disk->SetFillColor(40);
      Disk->SetLineColor(40);
      Disk->Draw();  
      
      // Hough Circles (Even, Odd)
      
      TEllipse *circle_even = new TEllipse(cX_even+ref.first,cY_even+ref.second,rad_even,rad_even);
      TEllipse *circle_odd = new TEllipse(cX_odd+ref.first,cY_odd+ref.second,rad_odd,rad_odd);

      circle_even->SetFillColor(0);
      circle_even->SetFillStyle(4000);
      circle_even->SetLineColor(2);
      circle_even->Draw();
      
      circle_odd->SetFillColor(0);
      circle_odd->SetFillStyle(4000);
      circle_odd->SetLineColor(4);      
      circle_odd->Draw();
      

      // Hits
      
      TGraph *grTruth = new TGraph(nCALCDCHit, WireEnd0X, WireEnd0Y);
      grTruth->SetTitle("Truth");
      grTruth->SetMarkerStyle(20);
      grTruth->SetMarkerSize(1);
      grTruth->Draw("P");
      
      std::cout << "  MC Truth NDF: " << nCALCDCHit << std::endl;
      std::cout << i_evt <<"-th Event"<< std::endl;
      
      //if (rad_even>33 && rad_odd>33){break;}
    }
  }
  
  c_div->cd(1);
  maxvote_even_hist->SetLineColor(4);
  maxvote_even_hist->SetFillColor(4);
  maxvote_even_hist->SetFillStyle(3004);
  maxvote_even_hist->Draw();
  maxvote_odd_hist->SetLineColor(2);
  maxvote_odd_hist->SetFillColor(2);
  maxvote_odd_hist->SetFillStyle(3005);
  maxvote_odd_hist->Draw("same");

  c_div->cd(7);
  rad_even_hist->SetLineColor(4);
  rad_even_hist->SetFillColor(4);
  rad_even_hist->SetFillStyle(3004);
  rad_even_hist->Draw();
  rad_odd_hist->SetLineColor(2);
  rad_odd_hist->SetFillColor(2);
  rad_odd_hist->SetFillStyle(3005);
  rad_odd_hist->Draw("same");


  std::cout << "Finish!" << std::endl;
  f->Close();
  
}

Float_t ConfTransX(Float_t x, Float_t y){
  return x/(pow(x,2)+pow(y,2));
} 

Float_t ConfTransY(Float_t x, Float_t y){
  return y/(pow(x,2)+pow(y,2));
} 

Float_t HoughTrans(Float_t x, Float_t y, Float_t theta){
  return x*TMath::Cos(theta*TMath::Pi()/180)+y*TMath::Sin(theta*TMath::Pi()/180);
}

Bool_t ifInsideDisk(Float_t x, Float_t y){
  if (pow(x,2)+pow(y,2)<100){ return 1;}
  else {return 0;}
};

std::vector<std::pair<Float_t,Float_t> > makeOrigins(Float_t rad, std::pair<Float_t,Float_t> ref){
  std::vector<std::pair<Float_t,Float_t> > origins;
  origins.push_back(std::make_pair(ref.first+rad, ref.second+rad));
  origins.push_back(std::make_pair(ref.first-rad, ref.second+rad));
  origins.push_back(std::make_pair(ref.first-rad, ref.second-rad));
  origins.push_back(std::make_pair(ref.first+rad, ref.second-rad));
  return origins;
}

Int_t findMaxpoint(std::vector<Int_t> vec){
  Int_t maxpoint;
  Int_t val=0;
  for (int i=0; i<vec.size(); i++){
    if (vec[i]>val){
      val=vec[i];
      maxpoint=i;
    }
  }
  return maxpoint;
}
