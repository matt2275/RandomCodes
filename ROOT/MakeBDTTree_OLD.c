/// \file
/// \ingroup tutorial_tree
/// \notebook -nodraw
/// Example of Root macro to copy a subset of a Tree to a new Tree, selecting entries.
///
/// Only selected entries are copied to the new Tree.
/// The input file has been generated by the program in `$ROOTSYS/test/Event`
/// with `Event 1000 1 99 1`
///
/// \macro_code
///
/// \author Rene Brun

#include "TMath.h"
using namespace std;


//int isGoodEvent(vector<float> Muon_Pt, vector<float> Muon_Eta, vector<float> Muon_Phi, vector<float> Ele_Pt, vector<float> Ele_Eta, vector<float> Ele_Phi, vector<float> Trk_Pt, vector<float> Trk_Eta, vector<float> Trk_Phi){

int isGoodEvent(vector<float> Muon[3], vector<float> Ele[3], vector<float> Trk[3]){
  int EventType = 0;
  string lepton_type;
  int lepton_index;
  int matched_track_index = -1;
  int unmatched_track_index = -1;
  int nMuons = Muon[0].size(); 
  int nTracks = Trk[0].size();
  int nElectrons = Ele[0].size();
  float Max_Muon_Pt = 0.0;
  float Max_Electron_Pt = 0.0;
  int Max_Muon_index = -1;
  int Max_Electron_index = -1;
  int nUnmatchedTracks = 0;
  for(int i = 0; i < nMuons; i++){
	if(Muon[0].at(i) > Max_Muon_Pt){
	  Max_Muon_Pt = Muon[0].at(i);
	  Max_Muon_index = i;
	}
  }

  for(int i = 0; i < nElectrons; i++){
	if(Ele[0].at(i) > Max_Electron_Pt){
	  Max_Electron_Pt = Ele[0].at(i);
	  Max_Electron_index = i;
	}
  }
  
  if(Max_Electron_Pt > Max_Muon_Pt){
	lepton_index = Max_Electron_index;
	lepton_type = "electron";
  }
 
	
  if(Max_Electron_Pt < Max_Muon_Pt){
	lepton_index = Max_Muon_index;
	lepton_type = "muon";
  }

    if(lepton_type == "muon"){
	for(int i = 0; i < nTracks; i++){
	  float diff_phi = fabs(Muon[1].at(lepton_index) - Trk[1].at(i));
 	  float diff_eta = fabs(Muon[2].at(lepton_index) - Trk[2].at(i));
      if( diff_eta < .01 and diff_phi < .01){
		matched_track_index = i;
	  }
      if(!( diff_eta < .01 and diff_phi < .01)){
		unmatched_track_index = i;
		nUnmatchedTracks ++;
	  }
	}
	}

    if(lepton_type == "electron"){
	for(int i = 0; i < nTracks; i++){
	  float diff_phi = fabs(Ele[1].at(lepton_index) - Trk[1].at(i));
 	  float diff_eta = fabs(Ele[2].at(lepton_index) - Trk[2].at(i));
      if( diff_eta < .01 and diff_phi < .01){
		matched_track_index = i;
	  }
      if(!( diff_eta < .01 and diff_phi < .01)){
		unmatched_track_index = i;
		nUnmatchedTracks ++;
	  }
	}
	}

    if(lepton_type == "muon" && unmatched_track_index != -1 && matched_track_index != -1 && nUnmatchedTracks == 1){
		EventType = 1;
	}	
	
    if(lepton_type == "electron" && unmatched_track_index != -1 && matched_track_index != -1 && nUnmatchedTracks == 1){
		EventType = 2;
	}
	return EventType;
}


void MakeBDTTree(string inputfile, string outfile)
{
	TFile *inFile = TFile::Open(inputfile.c_str());
	TTree *eventTree     = (TTree*)inFile->Get("ggHiNtuplizer/EventTree");

	int nEle = 0;
	vector<float> *elePt = nullptr;
	vector<float> *eleEta = nullptr;
	vector<float> *elePhi = nullptr;
	int nMu = 0;
	vector<float> *muPt = nullptr;
	vector<float> *muEta = nullptr;
	vector<float> *muPhi = nullptr;	
	int nTrk = 0;
	vector<float> *trkPt = nullptr;
	vector<float> *trkEta = nullptr;
	vector<float> *trkPhi = nullptr; 
	vector<float> *trkvx = nullptr;
	vector<float> *trkvy = nullptr;
	vector<float> *trkvz = nullptr;
	vector<float> *trkP = nullptr; 
	vector<float> *trkd0 = nullptr;
	vector<float> *trkdxy = nullptr;
	vector<float> *trkdz = nullptr;
	float trkDist = 0 ;
	float trkDeltaPhi = 0;

	float trkPt_0 = 0 ;
	float trkPt_1 = 0;
	float trkP_0 = 0 ;
	float trkP_1 = 0;
	float trkEta_0 = 0 ;
	float trkEta_1 = 0;
	float trkPhi_0 = 0;
	float trkPhi_1 = 0 ;
	float trkvx_0 = 0;
	float trkvx_1 = 0 ;
	float trkvy_0 = 0;
	float trkvy_1 = 0 ;
	float trkvz_0 = 0;
	float trkvz_1 = 0 ;
	float trkd0_0 = 0;
	float trkd0_1 = 0 ;
	float trkdxy_0 = 0;
	float trkdxy_1 = 0 ;
	float trkdz_0 = 0;
	float trkdz_1 = 0 ;

	eventTree->SetBranchStatus("*",0);
	eventTree->SetBranchStatus("nMC",1);
	eventTree->SetBranchStatus("mcPID",1);
	eventTree->SetBranchStatus("mcVtx_x",1);
	eventTree->SetBranchStatus("mcVtx_y",1);
	eventTree->SetBranchStatus("mcVtx_z",1);
	eventTree->SetBranchStatus("mcPt",1);
	eventTree->SetBranchStatus("mcEta",1);
	eventTree->SetBranchStatus("mcPhi",1);
	eventTree->SetBranchStatus("mcE",1);
	eventTree->SetBranchStatus("mcEt",1);
	eventTree->SetBranchStatus("mcMass",1);
	eventTree->SetBranchStatus("nEle",1);
	eventTree->SetBranchStatus("eleCharge",1);
	eventTree->SetBranchStatus("elePt",1);
	eventTree->SetBranchStatus("eleEta",1);
	eventTree->SetBranchStatus("elePhi",1);
	eventTree->SetBranchStatus("nMu",1);
	eventTree->SetBranchStatus("muPt",1);
	eventTree->SetBranchStatus("muEta",1);
	eventTree->SetBranchStatus("muPhi",1);
	eventTree->SetBranchStatus("muCharge",1);
	eventTree->SetBranchStatus("nTrk",1);
	eventTree->SetBranchStatus("trkPt",1);
	eventTree->SetBranchStatus("trkP",1);
	eventTree->SetBranchStatus("trkEta",1);
	eventTree->SetBranchStatus("trkPhi",1);
	eventTree->SetBranchStatus("trkcharge",1);
	eventTree->SetBranchStatus("trkvx",1);
	eventTree->SetBranchStatus("trkvy",1);
	eventTree->SetBranchStatus("trkvz",1);
	eventTree->SetBranchStatus("trknormchi2",1);
	eventTree->SetBranchStatus("trkchi2",1);
	eventTree->SetBranchStatus("trkd0",1);
	eventTree->SetBranchStatus("trkdxy",1);
	eventTree->SetBranchStatus("trkdz",1);
	
	eventTree->SetBranchAddress("nEle", &nEle);
	eventTree->SetBranchAddress("elePt", &elePt);
	eventTree->SetBranchAddress("eleEta", &eleEta);
	eventTree->SetBranchAddress("elePhi", &elePhi);
	eventTree->SetBranchAddress("nMu", &nMu);
	eventTree->SetBranchAddress("muPt", &muPt);
	eventTree->SetBranchAddress("muEta", &muEta);
	eventTree->SetBranchAddress("muPhi", &muPhi);	
	eventTree->SetBranchAddress("nTrk", &nTrk);
	eventTree->SetBranchAddress("trkPt", &trkPt);
	eventTree->SetBranchAddress("trkEta", &trkEta);
	eventTree->SetBranchAddress("trkPhi", &trkPhi);
	eventTree->SetBranchAddress("trkvx", &trkvx);
	eventTree->SetBranchAddress("trkvy", &trkvy);
	eventTree->SetBranchAddress("trkvz", &trkvz);
	eventTree->SetBranchAddress("trkd0", &trkd0);
	eventTree->SetBranchAddress("trkdxy", &trkdxy);
	eventTree->SetBranchAddress("trkdz", &trkdz);
	eventTree->SetBranchAddress("trkP", &trkP);


   int nentries = eventTree->GetEntries();

   // Create a new file + a clone of old tree in new file
   TFile newfile(outfile.c_str(), "recreate");
   auto newtree = eventTree->CloneTree(0);
   int index = 0;

	newtree->Branch("trkPt_0", &trkPt_0, "trkPt_0/F");
	newtree->Branch("trkPt_1", &trkPt_1, "trkPt_1/F");
	newtree->Branch("trkP_0", &trkP_0, "trkP_0/F");
	newtree->Branch("trkP_1", &trkP_1, "trkP_1/F");
	newtree->Branch("trkEta_0", &trkEta_0, "trkEta_0/F");
	newtree->Branch("trkEta_1", &trkEta_1, "trkEta_1/F");
	newtree->Branch("trkPhi_0", &trkPhi_0, "trkPhi_0/F");
	newtree->Branch("trkPhi_1", &trkPhi_1, "trkPhi_1/F");
	newtree->Branch("trkvx_0", &trkvx_0, "trkvx_0/F");
	newtree->Branch("trkvx_1", &trkvx_1, "trkvx_1/F");
	newtree->Branch("trkvy_0", &trkvy_0, "trkvy_0/F");
	newtree->Branch("trkvy_1", &trkvy_1, "trkvy_1/F");
	newtree->Branch("trkvz_0", &trkvz_0, "trkvz_0/F");
	newtree->Branch("trkvz_1", &trkvz_1, "trkvz_1/F");
	newtree->Branch("trkd0_0", &trkd0_0, "trkd0_0/F");
	newtree->Branch("trkd0_1", &trkd0_1, "trkd0_1/F");
	newtree->Branch("trkdxy_0", &trkdxy_0, "trkdxy_0/F");
	newtree->Branch("trkdxy_1", &trkdxy_1, "trkdxy_1/F");
	newtree->Branch("trkdz_0", &trkdz_0, "trkdz_0/F");	
	newtree->Branch("trkdz_1", &trkdz_1, "trkdz_1/F");
	newtree->Branch("trkDist", &trkDist, "trkDist/F");	
	newtree->Branch("trkDeltaPhi", &trkDeltaPhi, "trkdDeltaPhi/F");

   for (int i=0;i<nentries;i++) {
      eventTree->GetEntry(i);
	  vector<float> Muon[3] = {*muPt, *muEta, *muPhi};
	  vector<float> Ele[3] = {*elePt, *eleEta, *elePhi};
	  vector<float> Trk[3] = {*trkPt, *trkEta, *trkPhi};
      if (isGoodEvent(Muon, Ele, Trk) != 0 and nTrk == 2){
		 trkPt_0 = trkPt->at(0);
		 trkPt_1 = trkPt->at(1);
		 trkP_0 = trkP->at(0);
		 trkP_1 = trkP->at(1);
		 trkEta_0 = trkEta->at(0);
		 trkEta_1 = trkEta->at(1);
		 trkPhi_0 = trkPhi->at(0);
		 trkPhi_1 = trkPhi->at(1);
		 trkvx_0 = trkvx->at(0);
		 trkvx_1 = trkvx->at(1);
		 trkvy_0 = trkvy->at(0);
		 trkvy_1 = trkvy->at(1);
		 trkvz_0 = trkvz->at(0);
		 trkvz_1 = trkvz->at(1);
		 trkd0_0 = trkd0->at(0);
		 trkd0_1 = trkd0->at(1);
		 trkdxy_0 = trkdxy->at(0);
		 trkdxy_1 = trkdxy->at(1);
		 trkdz_0 = trkdz->at(0);
		 trkdz_1 = trkdz->at(1);
		 trkDeltaPhi = fabs(trkPhi->at(0) - trkPhi->at(1));
		 trkDist = fabs(trkvx->at(0) - trkvx->at(1)) + fabs(trkvy->at(0) - trkvy->at(1));
         newtree->Fill();
		 index ++;
	  }
	  if(index > 10000){
		cout << "breaking" << endl;
		break;
	  }
	nEle = 0;
	elePt->clear();
	eleEta->clear();
	elePhi->clear();
	nMu = 0;
	muPt->clear();
	muEta->clear();
	muPhi->clear();	
	nTrk = 0;
	trkPt->clear();
	trkEta->clear();
	trkPhi->clear(); 
	trkvx->clear();
	trkvy->clear();
	trkvz->clear();
	trkP->clear(); 
	trkd0->clear();
	trkdxy->clear();
	trkdz->clear();
	trkDist = 0 ;
	trkDeltaPhi = 0;
	trkPt_0 = 0 ;
	trkPt_1 = 0;
	trkP_0 = 0 ;
	trkP_1 = 0;
	trkEta_0 = 0 ;
	trkEta_1 = 0;
	trkPhi_0 = 0;
	trkPhi_1 = 0 ;
	trkvx_0 = 0;
	trkvx_1 = 0 ;
	trkvy_0 = 0;
	trkvy_1 = 0 ;
	trkvz_0 = 0;
	trkvz_1 = 0 ;
	trkd0_0 = 0;
	trkd0_1 = 0 ;
	trkdxy_0 = 0;
	trkdxy_1 = 0 ;
	trkdz_0 = 0;
	trkdz_1 = 0 ;
	  
   }

   newtree->Print();
   newfile.Write();
}
