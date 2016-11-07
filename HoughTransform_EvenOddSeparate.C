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
#include <iterator>
#include <cmath>
#include <tuple>

Float_t ConfTransX(Float_t x, Float_t y);
Float_t ConfTransY(Float_t x, Float_t y);
Float_t HoughTrans(Float_t x, Float_t y, Float_t theta);
Bool_t ifInsideDisk(Float_t x, Float_t y);
std::vector<std::pair<Float_t,Float_t> > makeOrigins(Float_t rad_uncertainty, std::pair<Float_t,Float_t> ref);
Int_t findMaxpoint(std::vector<Int_t> vec);
Bool_t ifInsideBand(Float_t x, Float_t y, Float_t cX, Float_t cY, Float_t outerR, Float_t interR);

void HoughTransform_EvenOddSeparate(){

  std::string dir = "../";
  std::string fileName = "trig_ep92_onlyPrimary.root";
  TFile *f = TFile::Open(TString(dir+fileName));
  TTree *t = (TTree*)f->Get("trdata");
  std::string outputdir = "../";
  std::string outputfileName = "fit_"+fileName;
  TFile *f_out = new TFile(TString(outputdir+outputfileName),"recreate");
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
  
  Int_t nRecoHit;
  Float_t RecoWireEnd0X[50000];
  Float_t RecoWireEnd0Y[50000];
  Float_t RecoWireEnd1X[50000];
  Float_t RecoWireEnd1Y[50000];
  Float_t RecoCDCDriftDist[50000];
  Float_t RecoWireLayerId[50000];

  t_out->Branch("nRecoHit", &nRecoHit, "nRecoHit/I");
  t_out->Branch("RecoWireEnd0X", RecoWireEnd0X, "RecoWireEnd0X[nRecoHit]/F");
  t_out->Branch("RecoWireEnd0Y", RecoWireEnd0Y, "RecoWireEnd0Y[nRecoHit]/F");
  t_out->Branch("RecoWireEnd1X", RecoWireEnd1X, "RecoWireEnd1X[nRecoHit]/F");
  t_out->Branch("RecoWireEnd1Y", RecoWireEnd1Y, "RecoWireEnd1Y[nRecoHit]/F");
  t_out->Branch("RecoCDCDriftDist", RecoCDCDriftDist, "RecoCDCDriftDist[nRecoHit]/F");
  t_out->Branch("RecoWireLayerId", RecoWireLayerId, "RecoWireLayerId[nRecoHit]/F");

  TCanvas *c_div = new TCanvas ("c_div", "c divided", 600,900);
  c_div->Divide(2,4);

  //TH1F *maxvote_even_hist = new TH1F("maxvote_even_hist","max_even_vote_hist", 100, 0, 100);
  //TH1F *maxvote_odd_hist = new TH1F("maxvote_odd_hist","max_odd_vote_hist", 100, 0, 100);
  TH1F *fitpT_hist  = new TH1F("fitpT","fitpT",100,60,160) ;
  TH1F *truthpT_hist = new TH1F("truthpT","truthpT",100,60,160) ;
  TH1F *diff_hist = new TH1F("diff","fitpT-truthpT",200,-100,100) ;
  TH1F *rad_even_hist = new TH1F("rad_even_hist","rad_even_hist", 100, 0, 50);
  TH1F *rad_odd_hist = new TH1F("rad_odd_hist","rad_odd_hist", 100, 0, 50);
  TH2F *ref_even_dist = new TH2F("ref_even_dist", "ref_even_dist", 30, -10, 10, 30, -10, 10);
  TH2F *ref_odd_dist = new TH2F("ref_odd_dist", "ref_odd_dist", 30, -10, 10, 30, -10, 10);
 
  TCanvas *c_hits = new TCanvas("c_hits", "c_hits", 1000,1000);

  Int_t niter=3;
  Int_t nBins=100;
  Float_t nPt=1000.0;
  Float_t rhomax=0.02;
  Float_t rhomin=-0.02;

  //////// Bandwidth Variables ////////

  Int_t TotalHits_Single=0;
  Int_t TotalHits_Multi=0;
  Int_t RecoHits_Single=0;
  Int_t RecoHits_Multi=0;
  
  Float_t bandwidth=3;
  Int_t nIter_band = 10;
  Float_t RecoRate_Single[10];
  Float_t RecoRate_Multi[10];

  /////////////////////////////////////
  /*
  for (Int_t i_band=0; i_band<nIter_band+1; i_band++){
    bandwidth=i_band+1;
  */
  for (Int_t i_evt=0; i_evt < 5; i_evt++){

    nRecoHit=0;
    memset (RecoWireEnd0X, 0, sizeof(RecoWireEnd0X)); 
    memset (RecoWireEnd0Y, 0, sizeof(RecoWireEnd0Y)); 
    memset (RecoWireEnd1X, 0, sizeof(RecoWireEnd1X)); 
    memset (RecoWireEnd1Y, 0, sizeof(RecoWireEnd1Y)); 
    memset (RecoCDCDriftDist, 0, sizeof(RecoCDCDriftDist));
    memset (RecoWireLayerId, 0, sizeof(RecoWireLayerId));
    
    t->GetEntry(i_evt);    
    
    std::vector<std::pair<Float_t,Float_t> > WireEnd0;
    std::vector<Int_t> LayerId;
    Int_t nEvenhits=0;
    Int_t nOddhits=0;    
    Float_t WireEnd0X_even[10000];
    Float_t WireEnd0Y_even[10000];  
    Float_t WireEnd0X_odd[10000];
    Float_t WireEnd0Y_odd[10000];

    Float_t ConfX[100000];
    Float_t ConfY[100000];

    Int_t deg_index;    
    Int_t rho_index;
    Int_t do_not_use;

    Int_t deg_even_index; 
    Int_t rho_even_index;
    std::pair <Float_t,Float_t> ref_even;
    Int_t deg_odd_index;
    Int_t rho_odd_index;
    std::pair <Float_t,Float_t> ref_odd;

    /////////////////////////////////////////////////////////////////////
    // Remove MultiPosition Hit & Sort Even, Odd Hit for Visulaization //
    /////////////////////////////////////////////////////////////////////

    for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
      if ((WireLayerId[i_hit]+1)%2==0){
	WireEnd0X_even[nEvenhits]=WireEnd0X[i_hit];
	WireEnd0Y_even[nEvenhits]=WireEnd0Y[i_hit];
	nEvenhits++;
      }
      else if((WireLayerId[i_hit]+1)%2==1){
	WireEnd0X_odd[nOddhits]=WireEnd0X[i_hit];
	WireEnd0Y_odd[nOddhits]=WireEnd0Y[i_hit];
	nOddhits++;
      }
    }

    ////////////////////////
    //    Main Analysis   //
    ////////////////////////
       
    if (nCALCDCHit>=30 && WireMaxLayerId>=4 && (ifMultiTurn==1 || ifSingleTurn==1)){ 

      for (Int_t is_even=0; is_even<2; is_even++){

	Float_t rad_uncertainty=5;
	std::pair <Float_t,Float_t> ref = std::make_pair(0,0);

	for (Int_t it2=0; it2<niter; it2++){
	  
	  std::vector<std::pair<Float_t,Float_t> > origins = makeOrigins(rad_uncertainty,ref);     
	  std::vector<Int_t> vote_max;
	  std::vector<std::pair<Int_t,Int_t> > index_vec; 
	  
	  for (Int_t it=0; it<origins.size(); it++){
	    
	    Float_t oriX=origins[it].first; 
	    Float_t oriY=origins[it].second;
	    
	    if  (ifInsideDisk(oriX,oriY)==0){
	      vote_max.push_back(0);
	    }
	    
            else if (ifInsideDisk(oriX,oriY)==1){
	      
	      /*---------------------------------------
		|                                        |
		|      Conformal Transformation          |
		|                                        |
		-----------------------------------------*/
	      	      	      
	      for (Int_t i_hit=0; i_hit<nCALCDCHit ; i_hit++){   
		ConfX[i_hit] = ConfTransX(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
		ConfY[i_hit] = ConfTransY(WireEnd0X[i_hit]-oriX,WireEnd0Y[i_hit]-oriY);
	      }
	      
	      /*-------------------------------------------
		|                                         |
		|      Hough Transformation & Voting      |
		|                                         |
		------------------------------------------*/
	      
	      Bool_t if_already_vote[nBins][nBins][2];
	      //TNtuple * hough = new TNtuple("hough","hough", "rho:theta:hitId:wireId:is_even");
	      TH2F *vote = new TH2F("vote_even", "vote_even", nBins, 0, nBins, nBins, 0, nBins);

	      
	      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
		memset(if_already_vote,0,sizeof(if_already_vote));
		for (Int_t i_Pt=0 ; i_Pt<nPt; i_Pt++){
		  
		  Float_t deg = i_Pt/nPt*180.0;
		  Float_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);	  
		  Int_t tmp_deg_index = (deg)*nBins/180;
		  Int_t tmp_rho_index = (rho+rhomax)*nBins/(rhomax-rhomin);
	      
		  if (is_even==WireLayerId[i_hit]%2 && if_already_vote[tmp_deg_index][tmp_rho_index][0]==0){
		    if_already_vote[tmp_deg_index][tmp_rho_index][0]=1;
		    vote->Fill(tmp_deg_index,tmp_rho_index);
		  }
		}
	      }
   
	      vote_max.push_back(vote->GetMaximum());	      
	      vote->GetMaximumBin(deg_index, rho_index, do_not_use);
	      index_vec.push_back(std::make_pair(deg_index,rho_index));
	      delete vote;
	    }
	  }

	  
	  Int_t maxIndex=findMaxpoint(vote_max);
	  ref.first = origins[maxIndex].first; 
	  ref.second = origins[maxIndex].second;

	  if (is_even==1){	  
	    ref_even.first = ref.first; 
	    ref_even.second = ref.second;
	    deg_even_index = index_vec[maxIndex].first; 
	    rho_even_index = index_vec[maxIndex].second;

	    if (it2==niter-1){
	      ref_even_dist->Fill(ref_even.first,ref_even.second);}

	  }

	  else if (is_even==0){
	    ref_odd.first = ref.first; 
	    ref_odd.second = ref.second;
	    deg_odd_index = index_vec[maxIndex].first; 
	    rho_odd_index = index_vec[maxIndex].second;

	    if (it2==niter-1){
	      ref_odd_dist->Fill(ref_odd.first,ref_odd.second);}

	  }
	  
	  rad_uncertainty = rad_uncertainty/2.0;
	
	}
      }
      
      Float_t deg_even_eval = 180.0/nBins*(2*deg_even_index+1)/2;
      Float_t rho_even_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_even_index-1)/2;   
      Float_t cX_even = TMath::Cos(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      Float_t cY_even = TMath::Sin(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      
      Float_t deg_odd_eval = 180.0/nBins*(2*deg_odd_index+1)/2;
      Float_t rho_odd_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_odd_index-1)/2;
      Float_t cX_odd = TMath::Cos(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
      Float_t cY_odd = TMath::Sin(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);

      Float_t rad_even = sqrt(pow(cX_even,2)+pow(cY_even,2));
      Float_t rad_odd = sqrt(pow(cX_odd,2)+pow(cY_odd,2));
      rad_even_hist->Fill(rad_even);
      rad_odd_hist->Fill(rad_odd);


      Float_t fitpT = (rad_even/0.3356+rad_odd/0.3356)/2;
      Float_t truthpT = sqrt(pow(genTrPy,2)+pow(genTrPz,2));
      fitpT_hist->Fill(fitpT);
      truthpT_hist->Fill(truthpT);
      diff_hist->Fill(fitpT-truthpT);

      
      /*-------------------------------------------------
	|                                               |
	|    Count Signal Hit Number in Circle Band     |
	|                                               |
	------------------------------------------------*/
      
      Int_t nRecoHit_even=0;
      Int_t nRecoHit_odd=0;
      Float_t RecoWireEnd0X_even[10000];
      Float_t RecoWireEnd0Y_even[10000];
      Float_t RecoWireEnd0X_odd[10000];
      Float_t RecoWireEnd0Y_odd[10000];
      Float_t outerR_even = rad_even + bandwidth/2.0;
      Float_t interR_even = rad_even - bandwidth/2.0;
      Float_t outerR_odd = rad_odd + bandwidth/2.0;
      Float_t interR_odd = rad_odd - bandwidth/2.0;


      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	if ((WireLayerId[i_hit]+1)%2 == 0){ // even
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_even+ref_even.first,cY_even+ref_even.second,outerR_even,interR_even)==1){
	    RecoWireEnd0X_even[nRecoHit_even]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_even[nRecoHit_even]=WireEnd0Y[i_hit];
	    nRecoHit_even++;	   

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[nRecoHit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[nRecoHit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[nRecoHit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[nRecoHit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[nRecoHit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[nRecoHit];
	    nRecoHit++;
	  }
	}
	else if ((WireLayerId[i_hit]+1)%2 == 1){ // odd
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_odd+ref_odd.first,cY_odd+ref_odd.second,outerR_odd,interR_odd)==1){
	    RecoWireEnd0X_odd[nRecoHit_odd]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_odd[nRecoHit_odd]=WireEnd0Y[i_hit];
	    nRecoHit_odd++;

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[nRecoHit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[nRecoHit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[nRecoHit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[nRecoHit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[nRecoHit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[nRecoHit];
	    nRecoHit++;
	  }
	}		
      } 


      if (ifSingleTurn==1){
	TotalHits_Single+=nCALCDCHit;
	RecoHits_Single+=(nRecoHit_even+nRecoHit_odd);
      }

      else if (ifMultiTurn==1){
	TotalHits_Multi+=nCALCDCHit;
	RecoHits_Multi+=(nRecoHit_even+nRecoHit_odd);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////
      //                                                                                               //
      //                                         PLOTTING SECTION                                      //  
      //                                                                                               //
      ///////////////////////////////////////////////////////////////////////////////////////////////////

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
      
      TEllipse *circle_even = new TEllipse(cX_even+ref_even.first,cY_even+ref_even.second,rad_even,rad_even);
      TEllipse *circle_odd = new TEllipse(cX_odd+ref_odd.first,cY_odd+ref_odd.second,rad_odd,rad_odd);

      circle_even->SetFillColor(0);
      circle_even->SetFillStyle(4000);
      circle_even->SetLineColor(4);
      circle_even->SetLineWidth(1);
      circle_even->Draw();
      
      circle_odd->SetFillColor(0);
      circle_odd->SetFillStyle(4000);
      circle_odd->SetLineColor(2);      
      circle_odd->SetLineWidth(1);
      circle_odd->Draw();
      
      // Hits
      
      TGraph *grEvenhits = new TGraph(nEvenhits, WireEnd0X_even, WireEnd0Y_even);
      grEvenhits->SetTitle("Evenhits");
      grEvenhits->SetMarkerStyle(20);
      grEvenhits->SetMarkerSize(1);
      grEvenhits->SetMarkerColor(4);  //even - blue
      grEvenhits->Draw("P");
      
      TGraph *grOddhits = new TGraph(nOddhits, WireEnd0X_odd, WireEnd0Y_odd);
      grOddhits->SetTitle("Oddhits");
      grOddhits->SetMarkerStyle(20);
      grOddhits->SetMarkerSize(1);
      grOddhits->SetMarkerColor(2);  //odd - red
      grOddhits->Draw("P");
      
      // Recognized Hits

      TGraph *grRecoEvenhits = new TGraph(nRecoHit_even, RecoWireEnd0X_even, RecoWireEnd0Y_even);
      grRecoEvenhits->SetTitle("RecoEvenhits");  
      grRecoEvenhits->SetMarkerStyle(20);
      grRecoEvenhits->SetMarkerSize(1);
      grRecoEvenhits->SetMarkerColor(38); // Reco_even - right blue
      grRecoEvenhits->Draw("P");
      
      TGraph *grRecoOddhits = new TGraph(nRecoHit_odd, RecoWireEnd0X_odd, RecoWireEnd0Y_odd);
      grRecoOddhits->SetTitle("RecoOddhits");
      grRecoOddhits->SetMarkerStyle(20);
      grRecoOddhits->SetMarkerSize(1);
      grRecoOddhits->SetMarkerColor(46); // Reco_odd - right red
      grRecoOddhits->Draw("P");      

      std::cout << "  MC Even hit: " << nEvenhits << std::endl;      
      std::cout << "  MC Odd hit: " << nOddhits << std::endl;      
      std::cout << "  MC Hit NDF: " << nCALCDCHit << std::endl;
      std::cout << i_evt <<"-th Event"<< std::endl;      

      t_out->Fill();
    }
  }
  
  c_div->cd(1);
  diff_hist->SetLineColor(4);
  diff_hist->SetFillColor(4);
  diff_hist->SetFillStyle(3004);
  diff_hist->Draw();

  c_div->cd(2);
  rad_even_hist->SetLineColor(4);
  rad_even_hist->SetFillColor(4);
  rad_even_hist->SetFillStyle(3004);
  rad_even_hist->Draw();
  rad_odd_hist->SetLineColor(2);
  rad_odd_hist->SetFillColor(2);
  rad_odd_hist->SetFillStyle(3005);
  rad_odd_hist->Draw("same");

  c_div->cd(3);
  truthpT_hist->SetLineColor(1);
  truthpT_hist->SetFillColor(4);
  truthpT_hist->SetFillStyle(3004);
  truthpT_hist->Draw();

  c_div->cd(4);
  fitpT_hist->SetLineColor(1);
  fitpT_hist->SetFillColor(4);
  fitpT_hist->SetFillStyle(3004);
  fitpT_hist->Draw();

  c_div->cd(5);
  ref_even_dist->Draw("colz");
  c_div->cd(6);
  ref_odd_dist->Draw("colz");
  
  /*
  std::cout << RecoHits_Single << std::endl;
  std::cout << TotalHits_Single << std::endl;
  std::cout << RecoHits_Multi << std::endl;
  std::cout << TotalHits_Multi << std::endl;
  std::cout << Float_t(RecoHits_Single)/TotalHits_Single*100.0 << "   " << Float_t(RecoHits_Multi)/TotalHits_Multi*100.0 << std::endl; 
  */

  std::cout << "Finish!" << std::endl;
  
  f_out->cd();
  t_out->Write();
  //f_out->Close();
  f->cd();
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

std::vector<std::pair<Float_t,Float_t> > makeOrigins(Float_t rad_uncertainty, std::pair<Float_t,Float_t> ref){
  std::vector<std::pair<Float_t,Float_t> > origins;
  origins.push_back(std::make_pair(ref.first, ref.second));  
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first, ref.second+rad_uncertainty));
  origins.push_back(std::make_pair(ref.first, ref.second-rad_uncertainty));
  origins.push_back(std::make_pair(ref.first+rad_uncertainty, ref.second));
  origins.push_back(std::make_pair(ref.first-rad_uncertainty, ref.second));
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

Bool_t ifInsideBand(Float_t x, Float_t y, Float_t cX, Float_t cY, Float_t outerR, Float_t interR){
  if ((pow(x-cX,2)+pow(y-cY,2))<pow(outerR,2) && (pow(x-cX,2)+pow(y-cY,2))>pow(interR,2)){
    return 1;
  }
  return 0;
}
