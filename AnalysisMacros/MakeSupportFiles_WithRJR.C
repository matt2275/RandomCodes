#include <TSystem.h>
#include <TROOT.h>
#include <TH1D.h>
#include <TTree.h>
#include <TChain.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPave.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <iostream>
#include<fstream>
#include <TClonesArray.h>
#include <cmath>//includes the C/C++ math library
#include <string>//includes C/C++ way of making strings
#include "TRandom.h"//includes random number generator as part of the program
#include "TString.h"//includes Root way of making strings
#include "TTree.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <vector>
using namespace std;


// structure for getting photon info around a track 
struct photonInfo {
    TLorentzVector pho_sum;
    int sum_nPho = 0;
};
 

// const bool DoZeroInvisMass = false;
const bool DoZeroVisMass = false;
const int Max_Events = 1000000;
const string prong_type = "3Prong";


float GetAcoplanarity(float Phi_1, float Phi_2)
{
  float phi1 = Phi_1;
  float phi2 = Phi_2;
  
  // Make sure that angles are in range [0, 2π)
  while(phi1 < 0)               phi1 += 2*TMath::Pi();
  while(phi1 >= 2*TMath::Pi())  phi1 -= 2*TMath::Pi();
  
  while(phi2 < 0)               phi2 += 2*TMath::Pi();
  while(phi2 >= 2*TMath::Pi())  phi2 -= 2*TMath::Pi();
  
  float deltaPhi = fabs(phi2-phi1);
  
  // Make sure that Δφ is in range [0, π]
  if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi() - deltaPhi;
  
  // Calcualte acoplanarity
  float aco = 1 - deltaPhi/TMath::Pi();
  
  return aco;
}

float GetDeltaPhi(float Phi_1, float Phi_2)
{
  float phi1 = Phi_1;
  float phi2 = Phi_2;
  
  return TMath::Pi() - fabs(fabs(phi1 - phi2) - TMath::Pi());
}


float GetDeltaR(float Phi_1, float Phi_2, float Eta_1, float Eta_2)
{
  return sqrt(pow(Eta_1-Eta_2, 2) + pow(Phi_1-Phi_2, 2));
}

// returns 4 vector info and num of photons of photons within dR of max_DR of a track
photonInfo GetPhotonInfo(int nPho, vector<float> *phoEt, vector<float> *phoEta, vector<float> *phoPhi, TLorentzVector Track, float max_DR){
  photonInfo s;
  TLorentzVector pho_sum;
  TLorentzVector tmp_pho;
  int pho_num = 0;
  for(int ipho =0 ; ipho < nPho ; ipho ++){
     tmp_pho.SetPtEtaPhiM(phoEt->at(ipho),phoEta->at(ipho),phoPhi->at(ipho),0);
     if(Track.DeltaR(tmp_pho) > max_DR) continue;
     pho_sum = pho_sum + tmp_pho;
     pho_num++;     
  }
  
  
  s.pho_sum = pho_sum;
  s.sum_nPho = pho_num;
  
  return s;
   
   
}

// returns TChain of all files that match string. Typically string ends in *.root 
TChain *CreateChainMC(string inFile, const char* treeName)//creates a TChain in which a tree name is fed into the argument
{
TChain* chain = new TChain(treeName);//declares a TChain called chain where a tree name is fed into the argument


chain->Add(inFile.c_str());//adds all data files in the directory in quotes to the chain, * is recursive command which will including every file whose name includes "HiForestAOD.root"
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0000/*.root");
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0001/*.root");
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0002/*.root");
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0003/*.root");
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0004/*.root");
//chain->Add("/eos/cms/store/group/phys_diffraction/lbyl_2018/HIForward_Reco/ntuples/ntuples_data/HIForward/ntuples_data_lbl/200617_140125/0005/*.root");
return chain;//returns the tree that is fed into the TChain arguments above from all of the chained data files
}

// returns TChain of Data files. strings are specified since more than one input folder
TChain *CreateChainData(const char* treeName)//creates a TChain in which a tree name is fed into the argument
{

 TChain* chain = new TChain(treeName);//declares a TChain called chain where a tree name is fed into the argument


chain->Add("/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/HIForward/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_Run2018/230325_0346280000/output_full_data_*.root");
chain->Add("/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/HIForward/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_Run2018/230325_034628/0001/output_full_data_*.root");
return chain;//returns the tree that is fed into the TChain arguments above from all of the chained data files
}


// main function returns Entries where the TnP or 3 prong system has the highest ranked pair of 
// different metrics highest mass, lowest acoplanarity, etc
// also adds addition variable like Recursive Jigsaw Reconstruction of 2 Particles and summed photon info

void MakeSupportFiles_WithRJR(){


vector< tuple<string, string>> DataTuple = {
// data type, file, color,  marker type, scale	
// {"Data","Histogram_ReduceTree_Starlight_Data.root" , 1 ,20, 1},
// {"starlight_MuMu","Histogram_ReduceTree_Starlight_MuMu.root" ,876 ,21, 1},
{"Data","/eos/user/m/mnickel/TauTauDecays_3Prong/output_data_3prong.root"},
{"DiTau_SuperChic","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggTauTau_TuneCP5_5p02TeV_SuperChic_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauSuperChic/230325_034441/0000/output_full_TauTau_mc_*.root"},
{"DiTau_Madgraph","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggTauTau_TuneCP5_5p02TeV_amcatnlo_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauMadgraph/230325_034351/0000/output_full_TauTau_mc_*.root"},
// {"DiTau_BSM1","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggTauTau_atau_m0p01_dtau_0p0_TuneCP5_5p02TeV_amcatnlo_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauBSM1/230201_220503/0000/output_full_TauTau_mc_*.root"},
// {"DiTau_BSM2","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggTauTau_atau_0p005_dtau_0p0_TuneCP5_5p02TeV_amcatnlo_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauBSM2/230201_220534/0000/output_full_TauTau_mc_*.root"},
{"BBbar","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggBBbar_4f_TuneCP5_5p02TeV_MG5_aMCatNLO_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_BBbar/230325_034527/0000/output_full_Other_mc_*.root"},
{"CCbar","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/ggCCbar_TuneCP5_5p02TeV_MG5_aMCatNLO_pythia8/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_CCbar/230325_034553/0000/output_full_Other_mc_*.root"},
// {"DiMuon_Gamma_OLD","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoMuMu_5p02TeV_gammaUPCEDFF-pLHE-v1/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_MuMuGamma/230201_220048/0000/output_full_MuMu_mc_*.root"},
{"DiMuon_Gamma","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/gammaUPCmumuFSR/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_MuMuGamma/230325_032902/0000/output_full_MuMu_mc_*.root"},
{"DiMuon_noFSR","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoMuMu_5p02TeV_STARlight/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_MuMu/230325_032740/0000/output_full_MuMu_mc_*.root"},
{"DiElectron_SuperChic","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/QEDGammaGamma_5p02TeV_SuperChic/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_EESuperChic/230324_214152/0000/output_full_EE_mc_*.root"},
{"DiElectron_Starlight","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/QEDGammaGamma_5p02TeV_STARlight/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_EEStarlight/230325_032420/0000/output_full_EE_mc_*.root"},
{"DiMuon_Arash","/eos/user/m/mnickel/CrabOut_MultiTag/DiMuon_Arash/output_full_MuMu_mc_*.root"},
// {"DiTau_GammaUPC_ChFF_kTSmearing","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoTauTau_5p02TeV_gammaUPCChFFkTSmearing-pLHE-v1/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauGammaUPC_ChFF_kTSmearing/230315_170943/0000/output_full_TauTau_mc_*.root"},
// {"DiTau_GammaUPC_ChFF","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoTauTau_5p02TeV_gammaUPCChFF-pLHE-v3/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauGammaUPC_ChFF/230315_170722/0000/output_full_TauTau_mc_*.root"},
// {"DiTau_GammaUPC_EDFF_kTSmearing","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoTauTau_5p02TeV_gammaUPCEDFFkTSmearing-pLHE-v1/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauGammaUPC_EDFF_kTSmearing/230315_163311/0000/output_full_TauTau_mc_*.root"},
// {"DiTau_GammaUPC_EDFF","/eos/user/m/mnickel/CrabOut_MultiTag/TnP_ntuples/muon/HIUPC/Run2018/AOD/GammaGammatoTauTau_5p02TeV_gammaUPCEDFF-pLHE-v1/crab_TnP_ntuplizer_muon_HIUPC_Run2018_AOD_TauTauGammaUPC_EDFF/230315_170546/0000/output_full_TauTau_mc_*.root"},


};


   map< string , string > DataTypes; 
   
   for(auto &[datatype, filename] : DataTuple){
      DataTypes[datatype] = filename;
      
   }
   
  vector<string> input_MCTypes = {"Data","DiTau_SuperChic","DiTau_Madgraph","DiMuon_Gamma","DiElectron_SuperChic","DiElectron_Starlight", "DiMuon_noFSR","BBbar","CCbar","DiMuon_Arash"};
  // vector<string> input_MCTypes = {"DiMuon_Gamma"};
  // vector<string> input_MCTypes = {"DiMuon_Arash"};
  // vector<string> input_MCTypes = {"Data"};
  if(prong_type == "3Prong"){
  for(string input_MCType : input_MCTypes){
  // string input_MCType = "Data";
string inputfile = DataTypes[input_MCType];
string outfile_begin = "/eos/user/m/mnickel/MultiTagData/3Prong/";
cout << inputfile <<" "<< outfile_begin <<endl;

TChain *t5;
if (input_MCType != "Data") t5 = CreateChainMC(inputfile, "muon/Events_3prong");
if (input_MCType == "Data") t5 = CreateChainData("muon/Events_3prong");
t5->SetBranchStatus("*",1);

int evt = t5->GetEntries();
   bool tag_isElectron;
   bool tag_isMuon;

   float tag_pt;
   float tag_eta;
   float tag_phi;
   float tau_3prong_total_pt;   
   float tau_3prong_total_eta;   
   float tau_3prong_total_phi;   
   float tau_3prong_total_M;   
   float tau_3prong_DR_1_2;
   float tau_3prong_DR_1_3;
   float tau_3prong_DR_2_3;
   float tau_3prong_DR_tag_1;
   float tau_3prong_DR_tag_2;
   float tau_3prong_DR_tag_3;
   float tau_3prong_total_vtx_prob;
   int pair_rank_vtx_prob_3prong;
   int pair_rank_dPhi_muons_3prong;
   int pair_rank_Mass_Mmumu_3prong;
   int iprobe_3prong;
   

  t5->SetBranchAddress("tag_isElectron", &tag_isElectron );
  t5->SetBranchAddress("tag_isMuon", &tag_isMuon );

  t5->SetBranchAddress("tag_pt", &tag_pt );
  t5->SetBranchAddress("tag_eta", &tag_eta );
  t5->SetBranchAddress("tag_phi", &tag_phi );
  t5->SetBranchAddress("tau_3prong_total_pt", &tau_3prong_total_pt );  
  t5->SetBranchAddress("tau_3prong_total_eta", &tau_3prong_total_eta );  
  t5->SetBranchAddress("tau_3prong_total_phi", &tau_3prong_total_phi );  
  t5->SetBranchAddress("tau_3prong_total_M", &tau_3prong_total_M );  
  t5->SetBranchAddress("tau_3prong_DR_1_2", &tau_3prong_DR_1_2 );
  t5->SetBranchAddress("tau_3prong_DR_1_3", &tau_3prong_DR_1_3 );
  t5->SetBranchAddress("tau_3prong_DR_2_3", &tau_3prong_DR_2_3 );
  t5->SetBranchAddress("tau_3prong_DR_tag_1", &tau_3prong_DR_tag_1 );
  t5->SetBranchAddress("tau_3prong_DR_tag_2", &tau_3prong_DR_tag_2 );
  t5->SetBranchAddress("tau_3prong_DR_tag_3", &tau_3prong_DR_tag_3 );
  t5->SetBranchAddress("tau_3prong_total_vtx_prob", &tau_3prong_total_vtx_prob );
  t5->SetBranchAddress("pair_rank_vtx_prob_3prong", &pair_rank_vtx_prob_3prong );    
  t5->SetBranchAddress("pair_rank_dPhi_muons_3prong", &pair_rank_dPhi_muons_3prong );    
  t5->SetBranchAddress("pair_rank_Mass_Mmumu_3prong", &pair_rank_Mass_Mmumu_3prong );    
  t5->SetBranchAddress("iprobe_3prong", &iprobe_3prong);  


   int nPho = 0;
   vector<float> *phoE = nullptr;
   vector<float> *phoEt = nullptr;
   vector<float> *phoEta = nullptr;
   vector<float> *phoPhi = nullptr;
   
  t5->SetBranchAddress("nPho", &nPho);  
  t5->SetBranchAddress("phoE", &phoE);  
  t5->SetBranchAddress("phoEt", &phoEt);  
  t5->SetBranchAddress("phoEta", &phoEta);  
  t5->SetBranchAddress("phoPhi", &phoPhi);  

string Invis_Mass_string = "";
string Vis_Mass_string = "";
// if(!DoZeroInvisMass) Invis_Mass_string = "HasInvisMass/";
// if(DoZeroInvisMass) Invis_Mass_string = "NoInvisMass/";
if(!DoZeroVisMass) Vis_Mass_string = "HasVisMass/";
if(DoZeroVisMass) Vis_Mass_string = "NoVisMass/";
string outfile;
outfile = outfile_begin + Vis_Mass_string +"MinimalCuts_"+input_MCType + ".root";
TFile f1(outfile.c_str(),"recreate");
f1.cd();
TTree *outTree = t5->CloneTree(0);

TTree *MaxEventTree = new TTree("MaxEventTree", "Event Tree with Max Event Info");

bool ExceedMaxEvents = false;

  MaxEventTree->Branch("TotalEvents", &evt);
  MaxEventTree->Branch("ExceedMaxEvents", &ExceedMaxEvents);

 bool hasHighestPt = false;
 bool hasBestVertex = false;
 bool hasBestDeltaPhi = false;
 bool hasHighestMass = false;

float Delta_Phi = 0;
float Aco = 0;
float Rapidity = 0;

  float Reco_Vis_dphi = 0;
  float RJR_tau1_Reco_pt = 0;
  float RJR_tau1_Reco_eta = 0;
  float RJR_tau1_Reco_phi = 0;
  float RJR_tau1_Reco_M = 0;
  float RJR_tau1_Invis_pt = 0;
  float RJR_tau1_Invis_eta = 0;
  float RJR_tau1_Invis_phi = 0;
  float RJR_tau1_Invis_M = 0;


  float RJR_tau2_Reco_pt = 0;
  float RJR_tau2_Reco_eta = 0;
  float RJR_tau2_Reco_phi = 0;
  float RJR_tau2_Reco_M = 0;
  float RJR_tau2_Invis_pt = 0;
  float RJR_tau2_Invis_eta = 0;
  float RJR_tau2_Invis_phi = 0;
  float RJR_tau2_Invis_M = 0;
  
  float RJR_Reco_Total_pt = 0;
  float RJR_Reco_Total_eta = 0;
  float RJR_Reco_Total_phi = 0;
  float RJR_Reco_Total_M = 0;
  float RJR_Reco_Total_Pz = 0;

outTree->Branch("Delta_Phi", &Delta_Phi);
outTree->Branch("Aco", &Aco);
outTree->Branch("Rapidity", &Rapidity);

  outTree->Branch("Reco_Vis_dphi", &Reco_Vis_dphi );
  outTree->Branch("RJR_tau1_Reco_pt", &RJR_tau1_Reco_pt );
  outTree->Branch("RJR_tau1_Reco_eta", &RJR_tau1_Reco_eta );
  outTree->Branch("RJR_tau1_Reco_phi", &RJR_tau1_Reco_phi );
  outTree->Branch("RJR_tau1_Reco_M", &RJR_tau1_Reco_M );
  outTree->Branch("RJR_tau1_Invis_pt", &RJR_tau1_Invis_pt );
  outTree->Branch("RJR_tau1_Invis_eta", &RJR_tau1_Invis_eta );
  outTree->Branch("RJR_tau1_Invis_phi", &RJR_tau1_Invis_phi );
  outTree->Branch("RJR_tau1_Invis_M", &RJR_tau1_Invis_M );
  
  outTree->Branch("RJR_tau2_Reco_pt", &RJR_tau2_Reco_pt );
  outTree->Branch("RJR_tau2_Reco_eta", &RJR_tau2_Reco_eta );
  outTree->Branch("RJR_tau2_Reco_phi", &RJR_tau2_Reco_phi );
  outTree->Branch("RJR_tau2_Reco_M", &RJR_tau2_Reco_M );
  outTree->Branch("RJR_tau2_Invis_pt", &RJR_tau2_Invis_pt );
  outTree->Branch("RJR_tau2_Invis_eta", &RJR_tau2_Invis_eta );
  outTree->Branch("RJR_tau2_Invis_phi", &RJR_tau2_Invis_phi );
  outTree->Branch("RJR_tau2_Invis_M", &RJR_tau2_Invis_M );
  
  outTree->Branch("RJR_Reco_Total_pt", &RJR_Reco_Total_pt );
  outTree->Branch("RJR_Reco_Total_eta", &RJR_Reco_Total_eta );
  outTree->Branch("RJR_Reco_Total_phi", &RJR_Reco_Total_phi );
  outTree->Branch("RJR_Reco_Total_M", &RJR_Reco_Total_M );
  
  outTree->Branch("hasHighestPt", &hasHighestPt);
  outTree->Branch("hasBestVertex", &hasBestVertex);
  outTree->Branch("hasBestDeltaPhi", &hasBestDeltaPhi);
  outTree->Branch("hasHighestMass", &hasHighestMass);


  int tag_nPho;
  float tag_phoPt;
  float tag_phoEta;
  float tag_phoPhi;
  float tag_phoM;
  
  int probe_nPho;
  float probe_phoPt;
  float probe_phoEta;
  float probe_phoPhi;
  float probe_phoM;
  
  outTree->Branch("tag_nPho", &tag_nPho);
  outTree->Branch("tag_phoPt", &tag_phoPt);
  outTree->Branch("tag_phoEta", &tag_phoEta);
  outTree->Branch("tag_phoPhi", &tag_phoPhi);
  outTree->Branch("tag_phoM", &tag_phoM);

  outTree->Branch("probe_nPho", &probe_nPho);
  outTree->Branch("probe_phoPt", &probe_phoPt);
  outTree->Branch("probe_phoEta", &probe_phoEta);
  outTree->Branch("probe_phoPhi", &probe_phoPhi);
  outTree->Branch("probe_phoM", &probe_phoM);


cout << " Total Events " << evt << endl;
for(int i=0; i<t5->GetEntries(); i++){
  if(i > Max_Events && input_MCType != "Data"){
     ExceedMaxEvents = true;
     break;
     }
  hasHighestPt = false;
  hasBestVertex = false;
  hasBestDeltaPhi = false;
  hasHighestMass = false;
  if(i%10000 == 0) cout << " Event " << i << endl;
  // if(i >= 5000000) break;
  // bool hasValidEvent = false;
  t5->GetEntry(i);
  // if(tau_3prong_DR_1_2 > 2) continue;
  // if(tau_3prong_DR_1_3 > 2) continue;
  // if(tau_3prong_DR_2_3 > 2) continue;
  // if(tau_3prong_DR_tag_1 < 1) continue;
  // if(tau_3prong_DR_tag_2 < 1) continue;
  // if(tau_3prong_DR_tag_3 < 1) continue;
  Delta_Phi = GetDeltaPhi(tag_phi,tau_3prong_total_phi);
  Aco = GetAcoplanarity(tag_phi,tau_3prong_total_phi);
  
  
  float pi = 3.14159265359;  
  float tag_M = 0;
  float probe_M = 0;
  TLorentzVector gentau1_visible;
  TLorentzVector gentau2_visible;
  TLorentzVector Reco_Vis;
  if(!DoZeroVisMass){
  if(tag_isElectron) tag_M = 0.510998E-3;
  if(tag_isMuon) tag_M = 105.6583755E-3;
  
  probe_M = tau_3prong_total_M;
  }


  // Initialize RJR Vectors
  gentau1_visible.SetPtEtaPhiM(tag_pt,tag_eta,tag_phi,tag_M);
  gentau2_visible.SetPtEtaPhiM(tau_3prong_total_pt,tau_3prong_total_eta,tau_3prong_total_phi,probe_M);
  Reco_Vis = gentau1_visible + gentau2_visible;
  Reco_Vis_dphi = gentau1_visible.DeltaPhi(gentau2_visible);
  
  
  // Get Photon Info
  photonInfo tag_phoInfo, probe_phoInfo;

  tag_phoInfo = GetPhotonInfo(nPho, phoEt,phoEta,phoPhi,gentau1_visible, .4);
  probe_phoInfo = GetPhotonInfo(nPho, phoEt,phoEta,phoPhi,gentau2_visible, .4);

  tag_nPho = tag_phoInfo.sum_nPho;
  tag_phoPt = tag_phoInfo.pho_sum.Pt();
  tag_phoEta = tag_phoInfo.pho_sum.Eta();
  tag_phoPhi = tag_phoInfo.pho_sum.Phi();
  tag_phoM = tag_phoInfo.pho_sum.M();

  probe_nPho = probe_phoInfo.sum_nPho;
  probe_phoPt = probe_phoInfo.pho_sum.Pt();
  probe_phoEta = probe_phoInfo.pho_sum.Eta();
  probe_phoPhi = probe_phoInfo.pho_sum.Phi();
  probe_phoM = probe_phoInfo.pho_sum.M();
  
  // Resume RJR things
  
  
  float Reco_Vis_M = Reco_Vis.M();
  float Reco_Vis_pt = Reco_Vis.Pt();
  float Reco_Vis_phi = Reco_Vis.Phi();
  float Reco_Vis_Pz = Reco_Vis.Pz();

  

  float RJR_Invis_pt = Reco_Vis_pt;
  float RJR_Invis_phi = Reco_Vis_phi+ pi ;
  float RJR_Invis_M = 0;
  float RJR_Invis_Pz = 0;
  RJR_Invis_M = Reco_Vis_M;

  RJR_Invis_Pz = Reco_Vis_Pz * sqrt(pow(RJR_Invis_pt,2) + pow(RJR_Invis_M,2))/sqrt(pow(Reco_Vis_pt,2) + pow(Reco_Vis_M,2));
  // float RJR_Invis_Pz = Reco_Invis_Pz;  
   // cout << " RJR Invis Pz " << RJR_Invis_Pz << " REco_Vis_Pz " << Reco_Vis_Pz << " REco Invis Pz " << Reco_Invis_Pz <<  endl;
  TLorentzVector RJR_Reco_Invis;
  
  // cout << " RECO_VIS_Eta " << Reco_Vis_eta << " TEST ETA " << asinh(Reco_Vis_Pz/ Reco_Vis_pt)<< endl;
  float RJR_Invis_eta = asinh(RJR_Invis_Pz/ RJR_Invis_pt);
  RJR_Reco_Invis.SetPtEtaPhiM(RJR_Invis_pt,RJR_Invis_eta ,RJR_Invis_phi,RJR_Invis_M);
  TLorentzVector RJR_Boost = RJR_Reco_Invis + Reco_Vis;
  gentau1_visible.Boost(-RJR_Boost.BoostVector());
  gentau2_visible.Boost(-RJR_Boost.BoostVector());
  TLorentzVector RJR_Reco_Vis = gentau1_visible + gentau2_visible;
 float RJR_C = .5*(1 + sqrt(pow(RJR_Reco_Vis.E(),2) - pow(Reco_Vis_M,2) + pow(RJR_Invis_M,2) )/ RJR_Reco_Vis.E());
 // cout << "RJR_C " << RJR_C << endl;
  TLorentzVector RJR_Invis_1;
  TLorentzVector RJR_Invis_2;
  RJR_Invis_1.SetVect((RJR_C -1)*gentau1_visible.Vect() - RJR_C*gentau2_visible.Vect());
  RJR_Invis_2.SetVect((RJR_C -1)*gentau2_visible.Vect() - RJR_C*gentau1_visible.Vect());
  
  // if(DoZeroInvisMass){
  if(false){
  RJR_Invis_1.SetE(RJR_Invis_1.P());
  RJR_Invis_2.SetE(RJR_Invis_2.P());
  }
  
  else{
  RJR_Invis_1.SetE((RJR_C -1)*gentau1_visible.E() + RJR_C*gentau2_visible.E());
  RJR_Invis_2.SetE((RJR_C -1)*gentau2_visible.E() + RJR_C*gentau1_visible.E());
  }
  
  // cout << "RJR_BOOST "<< RJR_Boost.BoostVector().Z() << " GEN_VIS BOOST " << gen_sum.BoostVector().Z() << endl;
  gentau1_visible.Boost(RJR_Boost.BoostVector());
  gentau2_visible.Boost(RJR_Boost.BoostVector());
  RJR_Invis_1.Boost(RJR_Boost.BoostVector());
  RJR_Invis_2.Boost(RJR_Boost.BoostVector());

   TLorentzVector RJR_Reco_tau1 = RJR_Invis_1 + gentau1_visible;
   TLorentzVector RJR_Reco_tau2 = RJR_Invis_2 + gentau2_visible;
   TLorentzVector RJR_Reco_Total = RJR_Reco_tau1 + RJR_Reco_tau2;
   
  RJR_tau1_Reco_pt = RJR_Reco_tau1.Pt();
  RJR_tau1_Reco_eta = RJR_Reco_tau1.Eta();
  RJR_tau1_Reco_phi = RJR_Reco_tau1.Phi();
  RJR_tau1_Reco_M = RJR_Reco_tau1.M();
  RJR_tau1_Invis_pt = RJR_Invis_1.Pt();
  RJR_tau1_Invis_eta =  RJR_Invis_1.Eta();
  RJR_tau1_Invis_phi =  RJR_Invis_1.Phi();
  RJR_tau1_Invis_M = RJR_Invis_1.M();


  RJR_tau2_Reco_pt = RJR_Reco_tau2.Pt();
  RJR_tau2_Reco_eta = RJR_Reco_tau2.Eta();
  RJR_tau2_Reco_phi = RJR_Reco_tau2.Phi();
  RJR_tau2_Reco_M = RJR_Reco_tau2.M();
  RJR_tau2_Invis_pt = RJR_Invis_2.Pt();
  RJR_tau2_Invis_eta =  RJR_Invis_2.Eta();
  RJR_tau2_Invis_phi =  RJR_Invis_2.Phi();
  RJR_tau2_Invis_M = RJR_Invis_2.M();
  
  RJR_Reco_Total_pt = RJR_Reco_Total.Pt();
  RJR_Reco_Total_eta = RJR_Reco_Total.Eta();
  RJR_Reco_Total_phi = RJR_Reco_Total.Phi();
  RJR_Reco_Total_M =RJR_Reco_Total.M();
  RJR_Reco_Total_Pz = RJR_Reco_Total.Pz();
  
  Rapidity = Reco_Vis.Rapidity();
  // if(tau_3prong_total_vtx_prob < 0.005) continue;
  // if(pair_rank_vtx_prob_3prong != 0) continue;
  if(iprobe_3prong == 1) hasHighestPt = true;
  if(pair_rank_vtx_prob_3prong == 0) hasBestVertex = true;
  if(pair_rank_dPhi_muons_3prong == 0) hasBestDeltaPhi = true;
  if(pair_rank_Mass_Mmumu_3prong == 0) hasHighestMass = true;
  
  if(hasHighestPt || hasHighestMass || hasBestDeltaPhi || hasBestVertex) outTree->Fill(); 
 
	}
   // t5->ResetBranchAddresses();
 MaxEventTree->Fill();

 f1.cd();
 outTree->Write();
 MaxEventTree->Write();
 f1.Close();  

  }
  }
  if(prong_type == "1Prong"){
 for(string input_MCType : input_MCTypes){
   
     // string input_MCType = "Data";
string inputfile = DataTypes[input_MCType];
string outfile_begin = "/eos/user/m/mnickel/MultiTagData/1Prong/";
cout << inputfile <<" "<< outfile_begin <<endl;
TChain *t4;
if (input_MCType != "Data") t4 = CreateChainMC(inputfile, "muon/Events");
if (input_MCType == "Data") t4 = CreateChainData("muon/Events");

int evt = t4->GetEntries();


   
   float tag_pt;
   float tag_eta;
   float tag_phi;
   float probe_pt;
   float probe_eta;
   float probe_phi;
   
   bool tag_isElectron;
   bool tag_isMuon;
   bool probe_isElectron;
   bool probe_isMuon;
   bool probe_isPion;

   int pair_rank_vtx_prob;
   int pair_rank_dPhi_muons;
   int pair_rank_Mass_Mmumu;

  t4->SetBranchAddress("tag_pt", &tag_pt );
  t4->SetBranchAddress("tag_eta", &tag_eta );
  t4->SetBranchAddress("tag_phi", &tag_phi );
  t4->SetBranchAddress("probe_pt", &probe_pt );
  t4->SetBranchAddress("probe_eta", &probe_eta );
  t4->SetBranchAddress("probe_phi", &probe_phi );

  t4->SetBranchAddress("tag_isElectron", &tag_isElectron );
  t4->SetBranchAddress("tag_isMuon", &tag_isMuon );
  t4->SetBranchAddress("probe_isElectron", &probe_isElectron );
  t4->SetBranchAddress("probe_isMuon", &probe_isMuon );
  t4->SetBranchAddress("probe_isPion", &probe_isPion );

  t4->SetBranchAddress("pair_rank_vtx_prob", &pair_rank_vtx_prob );    
  t4->SetBranchAddress("pair_rank_dPhi_muons", &pair_rank_dPhi_muons );    
  t4->SetBranchAddress("pair_rank_Mass_Mmumu", &pair_rank_Mass_Mmumu );    

   int nPho = 0;
   vector<float> *phoE = nullptr;
   vector<float> *phoEt = nullptr;
   vector<float> *phoEta = nullptr;
   vector<float> *phoPhi = nullptr;
   
  t4->SetBranchAddress("nPho", &nPho);  
  t4->SetBranchAddress("phoE", &phoE);  
  t4->SetBranchAddress("phoEt", &phoEt);  
  t4->SetBranchAddress("phoEta", &phoEta);  
  t4->SetBranchAddress("phoPhi", &phoPhi);  


t4->SetBranchStatus("*", 1 );

string Invis_Mass_string = "";
string Vis_Mass_string = "";
// if(!DoZeroInvisMass) Invis_Mass_string = "HasInvisMass/";
// if(DoZeroInvisMass) Invis_Mass_string = "NoInvisMass/";
if(!DoZeroVisMass) Vis_Mass_string = "HasVisMass/";
if(DoZeroVisMass) Vis_Mass_string = "NoVisMass/";
string outfile;
outfile = outfile_begin + Vis_Mass_string +"MinimalCuts_"+input_MCType + ".root";
TFile f2(outfile.c_str(),"recreate");
f2.cd();
TTree *outTree = t4->CloneTree(0);

TTree *MaxEventTree = new TTree("MaxEventTree", "Event Tree with Max Event Info");

bool ExceedMaxEvents = false;

  MaxEventTree->Branch("TotalEvents", &evt);
  MaxEventTree->Branch("ExceedMaxEvents", &ExceedMaxEvents);

float Delta_Phi = 0;
float Aco = 0;
float Rapidity = 0;


  float Reco_Vis_dphi = 0;
  float RJR_tau1_Reco_pt = 0;
  float RJR_tau1_Reco_eta = 0;
  float RJR_tau1_Reco_phi = 0;
  float RJR_tau1_Reco_M = 0;
  float RJR_tau1_Invis_pt = 0;
  float RJR_tau1_Invis_eta = 0;
  float RJR_tau1_Invis_phi = 0;
  float RJR_tau1_Invis_M = 0;


  float RJR_tau2_Reco_pt = 0;
  float RJR_tau2_Reco_eta = 0;
  float RJR_tau2_Reco_phi = 0;
  float RJR_tau2_Reco_M = 0;
  float RJR_tau2_Invis_pt = 0;
  float RJR_tau2_Invis_eta = 0;
  float RJR_tau2_Invis_phi = 0;
  float RJR_tau2_Invis_M = 0;
  
  float RJR_Reco_Total_pt = 0;
  float RJR_Reco_Total_eta = 0;
  float RJR_Reco_Total_phi = 0;
  float RJR_Reco_Total_M = 0;
  float RJR_Reco_Total_Pz = 0;

bool hasHighestPt = false;
bool hasBestVertex = false;
bool hasBestDeltaPhi = false;
bool hasHighestMass = false;

outTree->Branch("Delta_Phi", &Delta_Phi);
outTree->Branch("Aco", &Aco);
outTree->Branch("Rapidity", &Rapidity);

  outTree->Branch("Reco_Vis_dphi", &Reco_Vis_dphi );
  outTree->Branch("RJR_tau1_Reco_pt", &RJR_tau1_Reco_pt );
  outTree->Branch("RJR_tau1_Reco_eta", &RJR_tau1_Reco_eta );
  outTree->Branch("RJR_tau1_Reco_phi", &RJR_tau1_Reco_phi );
  outTree->Branch("RJR_tau1_Reco_M", &RJR_tau1_Reco_M );
  outTree->Branch("RJR_tau1_Invis_pt", &RJR_tau1_Invis_pt );
  outTree->Branch("RJR_tau1_Invis_eta", &RJR_tau1_Invis_eta );
  outTree->Branch("RJR_tau1_Invis_phi", &RJR_tau1_Invis_phi );
  outTree->Branch("RJR_tau1_Invis_M", &RJR_tau1_Invis_M );
  
  outTree->Branch("RJR_tau2_Reco_pt", &RJR_tau2_Reco_pt );
  outTree->Branch("RJR_tau2_Reco_eta", &RJR_tau2_Reco_eta );
  outTree->Branch("RJR_tau2_Reco_phi", &RJR_tau2_Reco_phi );
  outTree->Branch("RJR_tau2_Reco_M", &RJR_tau2_Reco_M );
  outTree->Branch("RJR_tau2_Invis_pt", &RJR_tau2_Invis_pt );
  outTree->Branch("RJR_tau2_Invis_eta", &RJR_tau2_Invis_eta );
  outTree->Branch("RJR_tau2_Invis_phi", &RJR_tau2_Invis_phi );
  outTree->Branch("RJR_tau2_Invis_M", &RJR_tau2_Invis_M );
  
  outTree->Branch("RJR_Reco_Total_pt", &RJR_Reco_Total_pt );
  outTree->Branch("RJR_Reco_Total_eta", &RJR_Reco_Total_eta );
  outTree->Branch("RJR_Reco_Total_phi", &RJR_Reco_Total_phi );
  outTree->Branch("RJR_Reco_Total_M", &RJR_Reco_Total_M );
  
  outTree->Branch("hasBestVertex", &hasBestVertex);
  outTree->Branch("hasBestDeltaPhi", &hasBestDeltaPhi);
  outTree->Branch("hasHighestMass", &hasHighestMass);

  int tag_nPho;
  float tag_phoPt;
  float tag_phoEta;
  float tag_phoPhi;
  float tag_phoM;
  
  int probe_nPho;
  float probe_phoPt;
  float probe_phoEta;
  float probe_phoPhi;
  float probe_phoM;
  
  outTree->Branch("tag_nPho", &tag_nPho);
  outTree->Branch("tag_phoPt", &tag_phoPt);
  outTree->Branch("tag_phoEta", &tag_phoEta);
  outTree->Branch("tag_phoPhi", &tag_phoPhi);
  outTree->Branch("tag_phoM", &tag_phoM);

  outTree->Branch("probe_nPho", &probe_nPho);
  outTree->Branch("probe_phoPt", &probe_phoPt);
  outTree->Branch("probe_phoEta", &probe_phoEta);
  outTree->Branch("probe_phoPhi", &probe_phoPhi);
  outTree->Branch("probe_phoM", &probe_phoM);


cout << " Total Events " << evt << endl;
for(int i=0; i<t4->GetEntries(); i++){
  hasBestVertex = false;
  hasBestDeltaPhi = false;
  hasHighestMass = false;
  if(i > Max_Events && input_MCType != "Data"){
     ExceedMaxEvents = true;
     break;
  }
  if(i%10000 == 0) cout << " Event " << i << endl;
  // if(i >= 5000000) break;
  // bool hasValidEvent = false;
  t4->GetEntry(i);
  Delta_Phi = GetDeltaPhi(tag_phi,probe_phi);
  Aco = GetAcoplanarity(tag_phi,probe_phi);

 
  float pi = 3.14159265359;  
  float tag_M = 0;
  float probe_M = 0;
  TLorentzVector gentau1_visible;
  TLorentzVector gentau2_visible;
  TLorentzVector Reco_Vis;
  if(!DoZeroVisMass){
  if(tag_isElectron) tag_M = 0.510998E-3;
  if(tag_isMuon) tag_M = 105.6583755E-3;
  
  if(probe_isElectron) probe_M = 0.510998E-3;
  if(probe_isMuon) probe_M = 105.6583755E-3;
  if(probe_isPion) probe_M = 139.57039E-3;
  }
  
  // Initialize RJR Vectors
  
  gentau1_visible.SetPtEtaPhiM(tag_pt,tag_eta,tag_phi,tag_M);
  gentau2_visible.SetPtEtaPhiM(probe_pt,probe_eta,probe_phi,probe_M);
  Reco_Vis = gentau1_visible + gentau2_visible;
  Reco_Vis_dphi = gentau1_visible.DeltaPhi(gentau2_visible);
  

  // Get Photon Info
  
  photonInfo tag_phoInfo, probe_phoInfo;
  
  tag_phoInfo = GetPhotonInfo(nPho, phoEt,phoEta,phoPhi,gentau1_visible, .4);
  probe_phoInfo = GetPhotonInfo(nPho, phoEt,phoEta,phoPhi,gentau2_visible, .4);
  
  tag_nPho = tag_phoInfo.sum_nPho;
  tag_phoPt = tag_phoInfo.pho_sum.Pt();
  tag_phoEta = tag_phoInfo.pho_sum.Eta();
  tag_phoPhi = tag_phoInfo.pho_sum.Phi();
  tag_phoM = tag_phoInfo.pho_sum.M();
  
  probe_nPho = probe_phoInfo.sum_nPho;
  probe_phoPt = probe_phoInfo.pho_sum.Pt();
  probe_phoEta = probe_phoInfo.pho_sum.Eta();
  probe_phoPhi = probe_phoInfo.pho_sum.Phi();
  probe_phoM = probe_phoInfo.pho_sum.M();
  
  // Resume RJR things
  
  
  float Reco_Vis_M = Reco_Vis.M();
  float Reco_Vis_pt = Reco_Vis.Pt();
  float Reco_Vis_phi = Reco_Vis.Phi();
  float Reco_Vis_Pz = Reco_Vis.Pz();

  

  float RJR_Invis_pt = Reco_Vis_pt;
  float RJR_Invis_phi = Reco_Vis_phi+ pi ;
  float RJR_Invis_M = 0;
  float RJR_Invis_Pz = 0;
  RJR_Invis_M = Reco_Vis_M;

  RJR_Invis_Pz = Reco_Vis_Pz * sqrt(pow(RJR_Invis_pt,2) + pow(RJR_Invis_M,2))/sqrt(pow(Reco_Vis_pt,2) + pow(Reco_Vis_M,2));
  // float RJR_Invis_Pz = Reco_Invis_Pz;  
   // cout << " RJR Invis Pz " << RJR_Invis_Pz << " REco_Vis_Pz " << Reco_Vis_Pz << " REco Invis Pz " << Reco_Invis_Pz <<  endl;
  TLorentzVector RJR_Reco_Invis;
  
  // cout << " RECO_VIS_Eta " << Reco_Vis_eta << " TEST ETA " << asinh(Reco_Vis_Pz/ Reco_Vis_pt)<< endl;
  float RJR_Invis_eta = asinh(RJR_Invis_Pz/ RJR_Invis_pt);
  RJR_Reco_Invis.SetPtEtaPhiM(RJR_Invis_pt,RJR_Invis_eta ,RJR_Invis_phi,RJR_Invis_M);
  TLorentzVector RJR_Boost = RJR_Reco_Invis + Reco_Vis;
  gentau1_visible.Boost(-RJR_Boost.BoostVector());
  gentau2_visible.Boost(-RJR_Boost.BoostVector());
  TLorentzVector RJR_Reco_Vis = gentau1_visible + gentau2_visible;
 float RJR_C = .5*(1 + sqrt(pow(RJR_Reco_Vis.E(),2) - pow(Reco_Vis_M,2) + pow(RJR_Invis_M,2) )/ RJR_Reco_Vis.E());
 // cout << "RJR_C " << RJR_C << endl;
  TLorentzVector RJR_Invis_1;
  TLorentzVector RJR_Invis_2;
  RJR_Invis_1.SetVect((RJR_C -1)*gentau1_visible.Vect() - RJR_C*gentau2_visible.Vect());
  RJR_Invis_2.SetVect((RJR_C -1)*gentau2_visible.Vect() - RJR_C*gentau1_visible.Vect());
  
  // if(DoZeroInvisMass){
  if(false){
  RJR_Invis_1.SetE(RJR_Invis_1.P());
  RJR_Invis_2.SetE(RJR_Invis_2.P());
  }
  
  else{
  RJR_Invis_1.SetE((RJR_C -1)*gentau1_visible.E() + RJR_C*gentau2_visible.E());
  RJR_Invis_2.SetE((RJR_C -1)*gentau2_visible.E() + RJR_C*gentau1_visible.E());
  }
  
  // cout << "RJR_BOOST "<< RJR_Boost.BoostVector().Z() << " GEN_VIS BOOST " << gen_sum.BoostVector().Z() << endl;
  gentau1_visible.Boost(RJR_Boost.BoostVector());
  gentau2_visible.Boost(RJR_Boost.BoostVector());
  RJR_Invis_1.Boost(RJR_Boost.BoostVector());
  RJR_Invis_2.Boost(RJR_Boost.BoostVector());

   TLorentzVector RJR_Reco_tau1 = RJR_Invis_1 + gentau1_visible;
   TLorentzVector RJR_Reco_tau2 = RJR_Invis_2 + gentau2_visible;
   TLorentzVector RJR_Reco_Total = RJR_Reco_tau1 + RJR_Reco_tau2;
   
  RJR_tau1_Reco_pt = RJR_Reco_tau1.Pt();
  RJR_tau1_Reco_eta = RJR_Reco_tau1.Eta();
  RJR_tau1_Reco_phi = RJR_Reco_tau1.Phi();
  RJR_tau1_Reco_M = RJR_Reco_tau1.M();
  RJR_tau1_Invis_pt = RJR_Invis_1.Pt();
  RJR_tau1_Invis_eta =  RJR_Invis_1.Eta();
  RJR_tau1_Invis_phi =  RJR_Invis_1.Phi();
  RJR_tau1_Invis_M = RJR_Invis_1.M();


  RJR_tau2_Reco_pt = RJR_Reco_tau2.Pt();
  RJR_tau2_Reco_eta = RJR_Reco_tau2.Eta();
  RJR_tau2_Reco_phi = RJR_Reco_tau2.Phi();
  RJR_tau2_Reco_M = RJR_Reco_tau2.M();
  RJR_tau2_Invis_pt = RJR_Invis_2.Pt();
  RJR_tau2_Invis_eta =  RJR_Invis_2.Eta();
  RJR_tau2_Invis_phi =  RJR_Invis_2.Phi();
  RJR_tau2_Invis_M = RJR_Invis_2.M();
  
  RJR_Reco_Total_pt = RJR_Reco_Total.Pt();
  RJR_Reco_Total_eta = RJR_Reco_Total.Eta();
  RJR_Reco_Total_phi = RJR_Reco_Total.Phi();
  RJR_Reco_Total_M =RJR_Reco_Total.M();
  RJR_Reco_Total_Pz = RJR_Reco_Total.Pz();
  
  
  Rapidity = Reco_Vis.Rapidity();

  if(pair_rank_vtx_prob == 0) hasBestVertex = true;
  if(pair_rank_dPhi_muons == 0) hasBestDeltaPhi = true;
  if(pair_rank_Mass_Mmumu == 0) hasHighestMass = true;
  
    if(hasHighestMass || hasBestDeltaPhi || hasBestVertex) outTree->Fill(); 
  
 
	}
   
   
   // t4->ResetBranchAddresses();
  MaxEventTree->Fill();
 f2.cd();
 outTree->Write();
 MaxEventTree->Write();
 f2.Close();
}
  }
}
