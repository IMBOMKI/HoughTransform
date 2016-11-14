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

Double_t ConfTransX(Double_t x, Double_t y);
Double_t ConfTransY(Double_t x, Double_t y);
Double_t HoughTrans(Double_t x, Double_t y, Double_t theta);
Bool_t ifInsideDisk(Double_t x, Double_t y);
std::vector<std::pair<Double_t,Double_t> > makeOrigins(Double_t rad_uncertainty, std::pair<Double_t,Double_t> ref);
Int_t findMaxpoint(std::vector<Int_t> vec);
Bool_t ifInsideBand(Double_t x, Double_t y, Double_t cX, Double_t cY, Double_t outerR, Double_t interR);
Bool_t ifInsideVec(Int_t element, std::vector<Int_t> vec);
/*
Double_t twoPtDistance(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
Bool_t ifInsideClusterSet(std::tuple< Double_t, Double_t, Int_t> tmpXYL,std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet);
void FindMatchedCluster(std::tuple< Double_t, Double_t, Int_t> tmpXYL,std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet);
Bool_t ifLayerIsEmptyInCluster(Int_t layer,std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet);
*/
bool sortFunction (int i,int j) { return (i<j); }
void MakeCluster(std::vector<Int_t> &WireIds, std::vector<std::vector< Int_t> > &ClusterSet, Int_t WireNumberInLayer);
bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size);


void HoughTransform_BandIter_em(int bw){

  std::string dir = "../";
  std::string fileName = "trig_em104_onlyPrimary.root";
  TFile *f = TFile::Open(TString(dir+fileName));
  TTree *t = (TTree*)f->Get("trdata");
  std::string outputdir = "FindedEvents/";
  std::string outputfileName = "bw"+std::to_string(bw)+"_"+"find_"+fileName;
  TFile *f_out = new TFile(TString(outputdir+outputfileName),"recreate");
  TTree *t_out = t->CloneTree(0);

  Int_t nCDCHit;
  Double_t CDCHitX[30000];
  Double_t CDCHitY[30000];
  Double_t CDCHitZ[30000];
  Double_t CDCHitT[30000];
  Double_t CDCEDep[30000];
  Int_t nCALCDCHit;
  Double_t CDCDriftDist[30000];
  Int_t CDCCharge[30000];
  Double_t WireEnd0X[30000];
  Double_t WireEnd0Y[30000];  
  Double_t WireEnd0Z[30000];
  Double_t WireEnd1X[30000];
  Double_t WireEnd1Y[30000];  
  Double_t WireEnd1Z[30000];
  Int_t WireLayerId[30000];
  Int_t WireMaxLayerId;
  Int_t WireId[30000];


  Double_t genTrE;
  Double_t genTrPx;
  Double_t genTrPy;
  Double_t genTrPz;
  Double_t genTrT;
  Int_t nTrigSTL;
  Int_t nTrigCRK;
  Int_t TrigSTL[1000];
  Int_t TrigCRK[1000];
  Double_t TrigSTLTime[1000];
  Double_t TrigCRKTime[1000];
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
  t->SetBranchAddress("WireId", WireId);

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
  Double_t RecoWireEnd0X[30000];
  Double_t RecoWireEnd0Y[30000];
  Double_t RecoWireEnd1X[30000];
  Double_t RecoWireEnd1Y[30000];
  Double_t RecoCDCDriftDist[30000];
  Int_t RecoWireLayerId[30000];
  Int_t RecoWireId[30000];
  Int_t RecoPID;
  Int_t RecoMaxWireLayerId;
  Bool_t Reco_ifCL3;
  Bool_t Reco_ifSingleTurn;
  Bool_t Reco_ifMultiTurn;
  Double_t fittedR_even;
  Double_t fittedR_odd;
  Double_t drEvenToOdd;
  Double_t TruthZ1;

  t_out->Branch("nRecoHit", &nRecoHit, "nRecoHit/I");
  t_out->Branch("RecoWireEnd0X", RecoWireEnd0X, "RecoWireEnd0X[nRecoHit]/D");
  t_out->Branch("RecoWireEnd0Y", RecoWireEnd0Y, "RecoWireEnd0Y[nRecoHit]/D");
  t_out->Branch("RecoWireEnd1X", RecoWireEnd1X, "RecoWireEnd1X[nRecoHit]/D");
  t_out->Branch("RecoWireEnd1Y", RecoWireEnd1Y, "RecoWireEnd1Y[nRecoHit]/D");
  t_out->Branch("RecoCDCDriftDist", RecoCDCDriftDist, "RecoCDCDriftDist[nRecoHit]/D");
  t_out->Branch("RecoWireLayerId", RecoWireLayerId, "RecoWireLayerId[nRecoHit]/I");
  t_out->Branch("RecoWireId", RecoWireId, "RecoWireId[nRecoHit]/I");
  t_out->Branch("RecoPID", &RecoPID, "RecoPID/I");
  t_out->Branch("TruthZ1", &TruthZ1, "TruthZ1/D");
  t_out->Branch("drEvenToOdd", &drEvenToOdd, "drEvenToOdd/D");
  t_out->Branch("TruthZ1", &TruthZ1, "TruthZ1/D");
  t_out->Branch("RecoMaxWireLayerId", &RecoMaxWireLayerId, "RecoMaxWireLayerId/I");
  t_out->Branch("Reco_ifCL3", &Reco_ifCL3, "Reco_ifCL3/O");
  t_out->Branch("Reco_ifSingleTurn", &Reco_ifSingleTurn, "Reco_ifSingleTurn/O");
  t_out->Branch("Reco_ifMultiTurn", &Reco_ifMultiTurn, "Reco_ifMultiTurn/O");  
  t_out->Branch("fittedR_even", &fittedR_even, "fittedR_even/D");  
  t_out->Branch("fittedR_odd", &fittedR_odd, "fittedR_odd/D");

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

  Int_t NumOfLayers=18;
  Int_t NumOfWiresPerLayer[18]={198,204,210,216,222,228,234,240,246,252,258,264,270,276,282,288,294,300};

  Int_t niter=3;
  Int_t nBins=100;
  Double_t nPt=1000.0;
  Double_t rhomax=0.02;
  Double_t rhomin=-0.02;

  //////// Bandwidth Variables ////////

  Int_t TotalHits_Single=0;
  Int_t TotalHits_Multi=0;
  Int_t RecoHits_Single=0;
  Int_t RecoHits_Multi=0;
  
  Double_t bandwidth=bw;
  Int_t nIter_band = 10;
  Double_t RecoRate_Single[10];
  Double_t RecoRate_Multi[10];

  /////////////////////////////////////

  for (Int_t i_evt=0; i_evt < t->GetEntries(); i_evt++){
    t->GetEntry(i_evt);    

    std::cout << "-------------------------" << std::endl;
    std::cout << i_evt <<"-th Event"<< std::endl;           
    std::cout << "MC Hit NDF: " << nCALCDCHit << std::endl;

    ///////////////////
    // Pre Track Cut //
    ///////////////////
    
    Bool_t ifThereisOddhit=0;
    Bool_t ifThereisEvenhit=0;
    if (nCALCDCHit<2) continue; // Require At least more than 2 hits
    for (Int_t i_layer=0; i_layer<nCALCDCHit; i_layer++){  
      if ((WireLayerId[i_layer]+1)%2==0) {ifThereisEvenhit=1;}
      else if ((WireLayerId[i_layer]+1)%2==1) {ifThereisOddhit=1;}
    }

    if (ifThereisEvenhit==0 || ifThereisOddhit==0) continue; // Require At least one hit at odd and even layer

    ///////////////////
    // Nullification //
    ///////////////////

    nRecoHit=0;
    memset (RecoWireEnd0X, 0, sizeof(RecoWireEnd0X)); 
    memset (RecoWireEnd0Y, 0, sizeof(RecoWireEnd0Y)); 
    memset (RecoWireEnd1X, 0, sizeof(RecoWireEnd1X)); 
    memset (RecoWireEnd1Y, 0, sizeof(RecoWireEnd1Y)); 
    memset (RecoCDCDriftDist, 0, sizeof(RecoCDCDriftDist));
    memset (RecoWireLayerId, 0, sizeof(RecoWireLayerId));
    RecoPID=0;
    TruthZ1=0;
    RecoMaxWireLayerId=0;
    Reco_ifCL3=0;
    Reco_ifSingleTurn=0;
    Reco_ifMultiTurn=0;

    
    std::vector<std::pair<Double_t,Double_t> > WireEnd0;
    std::vector<Int_t> LayerId;
    Int_t nEvenhits=0;
    Int_t nOddhits=0;    
    Double_t WireEnd0X_even[10000];
    Double_t WireEnd0Y_even[10000];  
    Double_t WireEnd0X_odd[10000];
    Double_t WireEnd0Y_odd[10000];

    Double_t ConfX[100000];
    Double_t ConfY[100000];

    Int_t deg_index;    
    Int_t rho_index;
    Int_t do_not_use;

    Int_t deg_even_index; 
    Int_t rho_even_index;
    std::pair <Double_t,Double_t> ref_even;
    Int_t deg_odd_index;
    Int_t rho_odd_index;
    std::pair <Double_t,Double_t> ref_odd;

    //////////////////////////////////////////////////////////
    // Sort Even, Odd Hit for Visulaization  &&  Clustering //
    //////////////////////////////////////////////////////////

    std::vector<Int_t> WireIdsPerLayer[18];
    std::vector< std::vector <Int_t> > tmpClusterSet;
    std::vector< std::vector <Int_t> > ClusterSet;

    for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){

      // Sort Even, Odd Hit
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
      
      // Sort Layer of Hits
      Int_t layerid=WireLayerId[i_hit];
      if(std::find(WireIdsPerLayer[layerid].begin(), WireIdsPerLayer[layerid].end(), WireId[i_hit]) != WireIdsPerLayer[layerid].end()){
	continue;
      }
      WireIdsPerLayer[layerid].push_back(WireId[i_hit]);    	  
      
    }
    
    // Sort one more time from smallest to biggest 
    // + Make Clusters and push into ClusterSet
    
    for (Int_t i_layer=0; i_layer<NumOfLayers; i_layer++){
      std::sort (WireIdsPerLayer[i_layer].begin(), WireIdsPerLayer[i_layer].end(), sortFunction);
      for (Int_t i_ele=0; i_ele<WireIdsPerLayer[i_layer].size(); i_ele++){
      }
      MakeCluster(WireIdsPerLayer[i_layer], tmpClusterSet, NumOfWiresPerLayer[i_layer]);
    }

    // Reduce size 1 Cluster
    for (Int_t i_clu=0; i_clu<tmpClusterSet.size(); i_clu++){
      if (tmpClusterSet[i_clu].size()>1){
	ClusterSet.push_back(tmpClusterSet[i_clu]);
      }
    }

   
    ////////////////////////
    //    Main Analysis   //
    ////////////////////////
       
      for (Int_t is_even=0; is_even<2; is_even++){
	  
	Double_t rad_uncertainty=5;
	std::pair <Double_t,Double_t> ref = std::make_pair(0,0);

	for (Int_t it2=0; it2<niter; it2++){
	  
	  std::vector<std::pair<Double_t,Double_t> > origins = makeOrigins(rad_uncertainty,ref);     
	  std::vector<Int_t> vote_max;
	  std::vector<std::pair<Int_t,Int_t> > index_vec; 
	  
	  for (Int_t it=0; it<origins.size(); it++){
	    
	    Double_t oriX=origins[it].first; 
	    Double_t oriY=origins[it].second;
	    
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
		  
		  Double_t deg = i_Pt/nPt*180.0;
		  Double_t rho =HoughTrans(ConfX[i_hit],ConfY[i_hit],deg);	  
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
      
      Double_t deg_even_eval = 180.0/nBins*(2*deg_even_index+1)/2;
      Double_t rho_even_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_even_index-1)/2;   
      Double_t cX_even = TMath::Cos(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      Double_t cY_even = TMath::Sin(deg_even_eval*TMath::Pi()/180.0)/(2*rho_even_eval);
      
      Double_t deg_odd_eval = 180.0/nBins*(2*deg_odd_index+1)/2;
      Double_t rho_odd_eval = rhomin+(rhomax-rhomin)/nBins*(2*rho_odd_index-1)/2;
      Double_t cX_odd = TMath::Cos(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);
      Double_t cY_odd = TMath::Sin(deg_odd_eval*TMath::Pi()/180.0)/(2*rho_odd_eval);

      Double_t rad_even = sqrt(pow(cX_even,2)+pow(cY_even,2));
      Double_t rad_odd = sqrt(pow(cX_odd,2)+pow(cY_odd,2));
      rad_even_hist->Fill(rad_even);
      rad_odd_hist->Fill(rad_odd);
      fittedR_even=rad_even;
      fittedR_odd=rad_odd;


      Double_t fitpT = (rad_even/0.3356+rad_odd/0.3356)/2;
      Double_t truthpT = sqrt(pow(genTrPy,2)+pow(genTrPz,2));
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
      Double_t RecoWireEnd0X_even[10000];
      Double_t RecoWireEnd0Y_even[10000];
      Double_t RecoWireEnd0X_odd[10000];
      Double_t RecoWireEnd0Y_odd[10000];
      Double_t outerR_even = rad_even + bandwidth/2.0;
      Double_t interR_even = rad_even - bandwidth/2.0;
      Double_t outerR_odd = rad_odd + bandwidth/2.0;
      Double_t interR_odd = rad_odd - bandwidth/2.0;
      
      RecoMaxWireLayerId=0;
      
      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	if ((WireLayerId[i_hit]+1)%2 == 0){ // even
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_even+ref_even.first,cY_even+ref_even.second,outerR_even,interR_even)==1){
	    RecoWireEnd0X_even[nRecoHit_even]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_even[nRecoHit_even]=WireEnd0Y[i_hit];
	    nRecoHit_even++;	   

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];

	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	  }
	}
	else if ((WireLayerId[i_hit]+1)%2 == 1){ // odd
	  if (ifInsideBand(WireEnd0X[i_hit],WireEnd0Y[i_hit],cX_odd+ref_odd.first,cY_odd+ref_odd.second,outerR_odd,interR_odd)==1){
	    RecoWireEnd0X_odd[nRecoHit_odd]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_odd[nRecoHit_odd]=WireEnd0Y[i_hit];
	    nRecoHit_odd++;

	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];

	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	  }
	}		
      }

      /////////////////////////////////////////////////////
      if (ifSingleTurn==1){
	TotalHits_Single+=nCALCDCHit;
	RecoHits_Single+=(nRecoHit_even+nRecoHit_odd);
      }

      else if (ifMultiTurn==1){
	TotalHits_Multi+=nCALCDCHit;
	RecoHits_Multi+=(nRecoHit_even+nRecoHit_odd);
      }
      ////////////////////////////////////////////////////

      /*-----------------------------------------------------------
	|                                                         |
	|           Count Non-Recognized Hits in Cluster          |
	|                                                         |
	----------------------------------------------------------*/

      std::vector <Int_t> WireIdsToBeReco;

      for (Int_t i_clu=0; i_clu<ClusterSet.size(); i_clu++){
	Bool_t isItRecoCluster=0;

	for(Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){	  
	  if (ifInsideArray(ClusterSet[i_clu][i_ele], RecoWireId, nRecoHit)==1){
	    isItRecoCluster=1;
	    break;
	  }
	}
	
	if (isItRecoCluster==1){
	  for(Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){	  
	    if (ifInsideArray(ClusterSet[i_clu][i_ele], RecoWireId, nRecoHit)==0){
	      WireIdsToBeReco.push_back(ClusterSet[i_clu][i_ele]);
	    }
	  }
	}
      }
  
      std::cout << "Num Of Ids to be Recognized: " << WireIdsToBeReco.size() << std::endl;

      for (Int_t i_hit=0; i_hit<nCALCDCHit; i_hit++){
	if (ifInsideVec(WireId[i_hit],WireIdsToBeReco)==1){

	  if ((WireLayerId[i_hit]+1)%2 == 0){ // even
	    RecoWireEnd0X_even[nRecoHit_even]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_even[nRecoHit_even]=WireEnd0Y[i_hit];
	    nRecoHit_even++;	   
	    
	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];
	    
	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;
	    
	  }
	  else if ((WireLayerId[i_hit]+1)%2 == 1){ // odd
	    RecoWireEnd0X_odd[nRecoHit_odd]=WireEnd0X[i_hit];
	    RecoWireEnd0Y_odd[nRecoHit_odd]=WireEnd0Y[i_hit];
	    nRecoHit_odd++;
	    
	    RecoWireEnd0X[nRecoHit]=WireEnd0X[i_hit];
	    RecoWireEnd0Y[nRecoHit]=WireEnd0Y[i_hit];
	    RecoWireEnd1X[nRecoHit]=WireEnd1X[i_hit];
	    RecoWireEnd1Y[nRecoHit]=WireEnd1Y[i_hit];
	    RecoCDCDriftDist[nRecoHit]=CDCDriftDist[i_hit];
	    RecoWireLayerId[nRecoHit]=WireLayerId[i_hit];
	    RecoWireId[nRecoHit]=WireId[i_hit];
	    
	    if (RecoWireLayerId[nRecoHit]>RecoMaxWireLayerId){
	      RecoMaxWireLayerId=RecoWireLayerId[nRecoHit];
	    }
	    nRecoHit++;	    
	  }			  
	} 
      }
      

      /*--------------------
	|                  |                                      
	|    CL3 Check     |
	|                  |
	-------------------*/

      Double_t abs_cX_even=cX_even+ref_even.first;
      Double_t abs_cY_even=cY_even+ref_even.second;
      Double_t abs_cX_odd=cX_odd+ref_odd.first;
      Double_t abs_cY_odd=cY_odd+ref_odd.second;      
      std::vector<Int_t> domain1_layer;
      std::vector<Int_t> domain2_layer;

      Double_t slope_even = cY_even/cX_even;
      Double_t slope_odd  = cY_odd/cX_odd;

      for (Int_t i_reco=0; i_reco<nRecoHit; i_reco++){
	if ((RecoWireLayerId[i_reco]+1)%2 == 0){ //even

	  Double_t LHS= RecoWireEnd0Y[i_reco];
	  Double_t RHS= slope_even*(RecoWireEnd0X[i_reco]-abs_cX_even)+abs_cY_even;
	  if (LHS>RHS){
	    domain1_layer.push_back(RecoWireLayerId[i_reco]);
	  }
	  else if (LHS<RHS){
	    domain2_layer.push_back(RecoWireLayerId[i_reco]);
	  }
	}

	else if ((RecoWireLayerId[i_reco]+1)%2 == 1){ //odd
	  Double_t LHS= RecoWireEnd0Y[i_reco];
	  Double_t RHS= slope_odd*(RecoWireEnd0X[i_reco]-abs_cX_odd)+abs_cY_odd;
	  if (LHS>RHS){
	    domain1_layer.push_back(RecoWireLayerId[i_reco]);
	  }
	  else if (LHS<RHS){
	    domain2_layer.push_back(RecoWireLayerId[i_reco]);
	  }
	}
      }

      if (ifInsideVec(0,domain1_layer)==1 && ifInsideVec(1,domain1_layer)==1 && ifInsideVec(2,domain1_layer)==1){
	if (ifInsideVec(0,domain2_layer)==1 && ifInsideVec(1,domain2_layer)==1 && ifInsideVec(2,domain2_layer)==1){
	  Reco_ifCL3=1;
	}
      }


      /*--------------------------------
	|                              |
	|    Charge Identification     | // Need to Correct
	|                              |
	-------------------------------*/

      /*
      Double_t crossZ=abs_cX_odd*abs_cY_even-abs_cY_odd*abs_cX_even;
      Int_t charge;

      if (crossZ>=0){ // e-
	RecoPID=11; 
	charge=-1;
      } 
      else if (crossZ<0){ // e+
	RecoPID=-11; 
	charge=1;
      } 
      */

      drEvenToOdd = sqrt(pow(abs_cX_even-abs_cX_odd,2)+pow(abs_cY_even-abs_cY_odd,2));
      TruthZ1=CDCHitZ[0];

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
      /*      
      c_hits->cd();
      c_hits->DrawFrame(-90,-90,90,90);
      
      // Wire Layers
      
      Double_t WireRad[18] = {53.00,54.60,56.20,57.80,59.40,61.00,62.60,64.20,65.80,67.40,69.00,70.60,72.20,73.80,75.40,77.00,78.60,80.20};
      
      for (Int_t i=0; i<18; i++){
	TEllipse *WireCircle = new TEllipse(0,0,WireRad[i],WireRad[i]);
	WireCircle->SetFillColor(0);
	WireCircle->SetFillStyle(4000);
	WireCircle->SetLineColor(40);
	WireCircle->Draw();
      }
      
      // Stopping target
      
      Double_t DiskRad=10.00;
      TEllipse *Disk = new TEllipse(0,0,DiskRad,DiskRad);
      Disk->SetFillColor(40);
      Disk->SetLineColor(40);
      Disk->Draw();  
      
      // Hough Circles (Even, Odd)
      
      TEllipse *circle_even = new TEllipse(abs_cX_even,abs_cY_even,rad_even,rad_even);
      TEllipse *circle_odd = new TEllipse(abs_cX_odd,abs_cY_odd,rad_odd,rad_odd);
      TEllipse *center_even = new TEllipse(abs_cX_even,abs_cY_even,0.1,0.1);
      TEllipse *center_odd = new TEllipse(abs_cX_odd,abs_cY_odd,0.1,0.1);

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

      center_even->SetFillColor(0);
      center_even->SetFillStyle(4000);
      center_even->SetLineColor(4);
      center_even->SetLineWidth(1);
      center_even->Draw();
      
      center_odd->SetFillColor(0);
      center_odd->SetFillStyle(4000);
      center_odd->SetLineColor(2);      
      center_odd->SetLineWidth(1);
      center_odd->Draw();
      
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

      // Center crossing line

      Int_t nPoints = 8;
      Double_t xLine_even[nPoints];
      Double_t yLine_even[nPoints];
      Double_t xLine_odd[nPoints];
      Double_t yLine_odd[nPoints];

      for (Int_t i_line=0; i_line<nPoints; i_line++) {	
	xLine_even[i_line] = i_line*10;
	yLine_even[i_line] = slope_even*(xLine_even[i_line]-abs_cX_even)+abs_cY_even;
      }
      for (Int_t i_line=0; i_line<nPoints; i_line++) {
	xLine_odd[i_line] = i_line*10;
	yLine_odd[i_line] = slope_odd*(xLine_odd[i_line]-abs_cX_odd)+abs_cY_odd;
      }

      TGraph *crossingLine_even = new TGraph(nPoints, xLine_even, yLine_even);
      crossingLine_even->SetLineColor(4);
      crossingLine_even->Draw("L");

      TGraph *crossingLine_odd = new TGraph(nPoints, xLine_odd, yLine_odd);
      crossingLine_odd->SetLineColor(2);
      crossingLine_odd->Draw("L");
      */
      t_out->Fill();

  }
  /*
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
  
  */

  std::cout << "Finish!" << std::endl;
  
  f_out->cd();
  f_out->Write();
  f_out->Close();
  f->cd();
  f->Close();

}

Double_t ConfTransX(Double_t x, Double_t y){
  return x/(pow(x,2)+pow(y,2));
} 

Double_t ConfTransY(Double_t x, Double_t y){
  return y/(pow(x,2)+pow(y,2));
} 

Double_t HoughTrans(Double_t x, Double_t y, Double_t theta){
  return x*TMath::Cos(theta*TMath::Pi()/180)+y*TMath::Sin(theta*TMath::Pi()/180);
}

Bool_t ifInsideDisk(Double_t x, Double_t y){
  if (pow(x,2)+pow(y,2)<100){ return 1;}
  else {return 0;}
};

std::vector<std::pair<Double_t,Double_t> > makeOrigins(Double_t rad_uncertainty, std::pair<Double_t,Double_t> ref){
  std::vector<std::pair<Double_t,Double_t> > origins;
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

Bool_t ifInsideBand(Double_t x, Double_t y, Double_t cX, Double_t cY, Double_t outerR, Double_t interR){
  if ((pow(x-cX,2)+pow(y-cY,2))<pow(outerR,2) && (pow(x-cX,2)+pow(y-cY,2))>pow(interR,2)){
    return 1;
  }
  return 0;
}

Bool_t ifInsideVec(Int_t element, std::vector<Int_t> vec){
  for (Int_t i=0; i<vec.size(); i++){
    if (vec[i]==element){
      return 1;
    }
  }
  return 0;
}
/*
Double_t twoPtDistance(Double_t x1, Double_t y1, Double_t x2, Double_t y2){
  return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}

Bool_t ifInsideClusterSet(std::tuple< Double_t, Double_t, Int_t> tmpXYL,std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet){
  Bool_t isInside=0;
  if (ClusterSet.empty()==1){return isInside;}

  for (Int_t i_clu=0; i_clu<ClusterSet.size(); i_clu++){
    for (Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){
      if (tmpXYL == ClusterSet[i_clu][i_ele]){
	isInside=1;
	break;
      }
    }
  } 

  return isInside;  
}

void FindMatchedCluster(std::tuple< Double_t, Double_t, Int_t> tmpXYL, std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet){

  if (!ClusterSet.empty()){
      for (Int_t i_clu=0; i_clu<ClusterSet.size(); i_clu++){
	for (Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){
	  if (std::get<2>(tmpXYL)==std::get<2>(ClusterSet[i_clu][i_ele])){
	    if (twoPtDistance(std::get<0>(tmpXYL),std::get<1>(tmpXYL),std::get<0>(ClusterSet[i_clu][i_ele]),std::get<1>(ClusterSet[i_clu][i_ele]))<2){
	      ClusterSet[i_clu].push_back(tmpXYL);	     	      
	    }
	  }
	}
      }   
    }
}

Bool_t ifLayerIsEmptyInCluster(Int_t layer,std::vector < std::vector <std::tuple <Double_t, Double_t, Int_t> > > ClusterSet){
  Bool_t isEmpty=1;
  if (ClusterSet.empty()==1) { return isEmpty;}
  else if (!ClusterSet.empty()){
    for (Int_t i_clu=0; i_clu<ClusterSet.size(); i_clu++){
      for (Int_t i_ele=0; i_ele<ClusterSet[i_clu].size(); i_ele++){
	if (layer==std::get<2>(ClusterSet[i_clu][i_ele])){
	  isEmpty=0;
	  return isEmpty;
	}
      }   
    }
  }
  return isEmpty;
}
*/
void MakeCluster(std::vector<Int_t> &WireIds, std::vector<std::vector< Int_t> > &ClusterSet, Int_t WireNumberInLayer)
{
  if (!WireIds.empty()){
    std::vector<Int_t> tmpCluster;
    std::vector<Int_t> firstCluster;
    tmpCluster.push_back(WireIds[0]);
    firstCluster.push_back(WireIds[0]);

    for(Int_t i_id = 1; i_id < WireIds.size(); ++i_id){
      if(WireIds[i_id] == WireIds[i_id-1] + 1){
	tmpCluster.push_back(WireIds[i_id]);
	
	if (i_id==WireIds.size()-1){ // Last element 
	  // 1. Check if it satisfies periodic boundary condition with first id of layer
	  if (WireIds[i_id]-WireIds[0]+1==WireNumberInLayer){
	    firstCluster.insert(firstCluster.end(), tmpCluster.begin(),tmpCluster.end());	    
	    ClusterSet.push_back(firstCluster);	
	    //	    std::cout << "Found Periodic Boundary!" << std::endl;
	  }
	  // 2. Otherwise, just push it.
	  else{
	    ClusterSet.push_back(firstCluster);
	    ClusterSet.push_back(tmpCluster);
	  }
  
	}

      }
      else if (WireIds[i_id] != WireIds[i_id-1] + 1){
	if (std::find(tmpCluster.begin(), tmpCluster.end(), WireIds[0]) != tmpCluster.end()){
	  firstCluster=tmpCluster;
	}
	else {
	  ClusterSet.push_back(tmpCluster);
	}
	tmpCluster.clear();
	tmpCluster.push_back(WireIds[i_id]);
      }
    }

  }
}

bool ifInsideArray(Int_t wireid, Int_t *reco_ids, Int_t arr_size){
  for (Int_t i=0; i<arr_size; i++){
    if (reco_ids[i]==wireid) return 1;
  }
  return 0;
}
