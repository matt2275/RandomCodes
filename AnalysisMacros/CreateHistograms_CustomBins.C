
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <typeinfo>
#include <string>
#include <algorithm>
#include<boost/algorithm/string.hpp>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <dirent.h>
#include <errno.h>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

// this has the custom bins information
#include "variables.h"

// need this for RNode to work to add new variables
//  . /cvmfs/sft.cern.ch/lcg/releases/LCG_96/ROOT/6.18.00/x86_64-centos7-gcc8-opt/bin/thisroot.sh


 
using namespace std;
using namespace ROOT;

// Name of Folder for Where histograms will be put
string RJR_Method = "HasVisMass";

// variables that are modified via this config file
string prong_type = "";
string save_folder = "";
string infile_name_beginning = "";

// Variables to save histograms
bool do1DHists = true;
bool do2DHists = true;


// maps for different cut selection
 map<string,string> CutTypesMap;
 map<string,string> SelectionTypesMap;
 vector<string> CutTypesVec;



// functions used to add variables in RDataframe. i.e. seeing the sum of dR of trks around a lead trk 

auto CalcPt = [](ROOT::VecOps::RVec< double > &px, ROOT::VecOps::RVec< double > &py, ROOT::VecOps::RVec< double > &E) {
   ROOT::VecOps::RVec< double > v;
   for (auto i=0U;i < px.size(); ++i) {
      if (E[i] > 100) {
         v.emplace_back(sqrt(px[i]*px[i] + py[i]*py[i]));
      }
   }
   return v;
};

bool GetBool(ROOT::VecOps::RVec< int > &output_vec, ROOT::VecOps::RVec< bool > &check_vec) {
   bool v = false;
   for (int i=0;i < check_vec.size(); i++) {
      if (check_vec[i]) {
         v = output_vec[i];
      }
   }
   return v;
}

int GetInt(ROOT::VecOps::RVec< int > &output_vec, ROOT::VecOps::RVec< bool > &check_vec) {
   int v = -9999;
   for (int i=0;i < check_vec.size(); i++) {
      if (check_vec[i]) {
         v = output_vec[i];
      }
   }
   return v;
}

float GetFloat(ROOT::VecOps::RVec< float > &output_vec, ROOT::VecOps::RVec< bool > &check_vec) {
   float v = -9999;
   for (int i=0;i < check_vec.size(); i++) {
      if (check_vec[i]) {
         v = output_vec[i];
      }
   }
   return v;
}

float GetMinMaxEta(ROOT::VecOps::RVec< float > &eta_vec, ROOT::VecOps::RVec< float> &pt_vec, float pt_cut, bool MaxEta) {
   float min = 9999;
   float max = -9999;
   for (int i=0;i < pt_vec.size(); i++) {
      float eta = eta_vec[i];
      float pt = pt_vec[i];
      if (pt < pt_cut) continue;
      if (eta > max) max = eta;
      if (eta < min) min = eta;

   }
   if(max > 2.5 ) max = 2.5;
   if(min < -2.5 ) min = -2.5;
   if(MaxEta) return max;
   else return min;
}

float GetSumDR(ROOT::VecOps::RVec< float > &eta_vec, ROOT::VecOps::RVec< float> &phi_vec, float tag_eta, float tag_phi, float probe_eta, float probe_phi, string sum_type = "all") {
   float sum = 0;
   float tag_sum = 0;
   float probe_sum = 0;
   for (int i=0;i < eta_vec.size(); i++) {
      float eta = eta_vec[i];
      float phi = phi_vec[i];
      float tag_DR =  pow(tag_eta - eta,2) + pow(tag_phi - phi,2);
      float probe_DR = pow(probe_eta - eta,2) + pow(probe_phi - phi,2);
      if(tag_DR < probe_DR) tag_sum = tag_sum + tag_DR;
      else if(probe_DR < tag_DR) probe_sum = probe_sum + probe_DR;

   }
   sum = tag_sum + probe_sum;
   
   if(sum_type == "probe") return sqrt(probe_sum);
   else if(sum_type == "tag") return sqrt(tag_sum);
   else return sqrt(sum);
}



// Creates new varibles to plot and use for selections
ROOT::RDF::RNode AddVariables(ROOT::RDF::RNode node, string TreeName)
{
   if(TreeName == "Events"){
   return node.Define("abs_tag_eta","abs(tag_eta)").Define("abs_probe_eta", "abs(probe_eta)")
   .Define("probe_muon_softID","GetBool(muIDSoft,mu_isProbe)").Define("tag_muon_softID","GetBool(muIDSoft,mu_isTag)")
   .Define("probe_electron_dEtaAtVtx","GetFloat(eledEtaAtVtx,ele_isProbe)").Define("tag_electron_dEtaAtVtx","GetFloat(eledEtaAtVtx,ele_isTag)")
   .Define("probe_electron_HoverE","GetFloat(eleHoverE,ele_isProbe)").Define("tag_electron_HoverE","GetFloat(eleHoverE,ele_isTag)")
   .Define("probe_electron_MissHits","GetInt(eleMissHits,ele_isProbe)").Define("tag_electron_MissHits","GetInt(eleMissHits,ele_isTag)")
   .Define("Rel_Diff_HF","(maxHFp - maxHFm)/(maxHFp + maxHFm)").Define("Rel_Diff_ZDC","(ZDC_P_Total_Energy - ZDC_M_Total_Energy)/(ZDC_P_Total_Energy + ZDC_M_Total_Energy)")
   .Define("PosRapGap","2.5 -  GetMinMaxEta(trkEta,trkPt,.2,true)").Define("NegRapGap","2.5 +  GetMinMaxEta(trkEta,trkPt,.2,false)")
   .Define("SumDeltaR", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,probe_eta,probe_phi,\"all\")")
   .Define("SumDeltaR_Tag", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,probe_eta,probe_phi,\"tag\")")
   .Define("SumDeltaR_Probe", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,probe_eta,probe_phi,\"probe\")");
      return node;
   }
   else if(TreeName == "Events_3prong"){
   return node.Define("abs_tag_eta","abs(tag_eta)")
   .Define("tag_electron_dEtaAtVtx","GetFloat(eledEtaAtVtx,ele_isTag)")
   .Define("tag_muon_softID","GetBool(muIDSoft,mu_isTag)").Define("tag_electron_HoverE","GetFloat(eleHoverE,ele_isTag)").Define("tag_electron_MissHits","GetInt(eleMissHits,ele_isTag)")
   .Define("Rel_Diff_HF","(maxHFp - maxHFm)/(maxHFp + maxHFm)").Define("Rel_Diff_ZDC","(ZDC_P_Total_Energy - ZDC_M_Total_Energy)/(ZDC_P_Total_Energy + ZDC_M_Total_Energy)")
   .Define("PosRapGap","2.5 -  GetMinMaxEta(trkEta,trkPt,.2,true)").Define("NegRapGap","2.5 +  GetMinMaxEta(trkEta,trkPt,.2,false)")
   .Define("SumDeltaR", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,tau_3prong_total_eta,tau_3prong_total_phi,\"all\")")
   .Define("SumDeltaR_Tag", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,tau_3prong_total_eta,tau_3prong_total_phi,\"tag\")")
   .Define("SumDeltaR_Probe", "GetSumDR(trkEta,trkPhi,tag_eta,tag_phi,tau_3prong_total_eta,tau_3prong_total_phi,\"probe\")");
      return node;
   }
   else{
      return node;
   }
}


// for cleaning strings provided in Selection and CutType files
string eraseSubStr(std::string cut_file, std::string  toErase)
{
    // Search for the substring in string
    string tmp_cut_string = "VOID";
    ifstream file(cut_file);
    string full_cut_string;
    while(getline(file, full_cut_string)){
    size_t pos = full_cut_string.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        full_cut_string.erase(pos, toErase.length());
        tmp_cut_string = full_cut_string;
    }
}
return tmp_cut_string;
}


// Removes spaces and \r \n from string
string CleanString(string dirty_string){
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\n'), dirty_string.cend()); 
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\r'), dirty_string.cend());
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), ' '), dirty_string.cend());
  return dirty_string;
}


// Maps the Specific Cuts that can be applied
// included in this is a universal cut that is applied to all selections
void MakeCutTypeMap (std::string cut_type_file)
{
    // Search for the substring in string
    string toErase = ":";
    string selection_type = "VOID";
    string selection_info = "VOID";
    ifstream file(cut_type_file);
    string full_cut_string;
    while(getline(file, full_cut_string)){
    size_t pos = full_cut_string.find(toErase);
    if (pos != std::string::npos)
    {
        // If found then erase it from string
        selection_type = full_cut_string.substr( 0, pos);
        selection_info = full_cut_string.substr(pos +1 , full_cut_string.length());
        selection_type.erase(std::remove(selection_type.begin(), selection_type.end(), '\n'), selection_type.cend());
        selection_info.erase(std::remove(selection_info.begin(), selection_info.end(), '\n'), selection_info.cend());
        selection_type.erase(std::remove(selection_type.begin(), selection_type.end(), '\r'), selection_type.cend());
        selection_info.erase(std::remove(selection_info.begin(), selection_info.end(), '\r'), selection_info.cend());
        CutTypesMap[selection_type] = selection_info;
        CutTypesVec.push_back(selection_type);
    }
}
}

//Maps the Selections to be applied. Takes information from Selection file
void MakeSelectionMap(std::string Selection_List, std::string Selection, string folder_loc)
{
    string Full_Selection = "";
    vector<string> CutTypes_contained{};
    vector<char>  valid_chars{' ', ')','(','!','&','|','\n','\r'};
    vector<char> valid_selection_chars{':', ' '};
    // Search for the substring in string
    ifstream file(Selection_List);
    string full_cut_string;
    while(getline(file, full_cut_string)){
    int selection_num = 0;
    bool valid_start = false;
    bool valid_end = false;
    size_t pos_selection = full_cut_string.find(Selection);
    if(pos_selection  != std::string::npos){
    if(pos_selection == 0) valid_start = true;
    if(pos_selection != 0)
    {
    if(std::find(valid_selection_chars.begin(), valid_selection_chars.end(), full_cut_string[pos_selection -1]) != valid_selection_chars.end()) valid_start = true;
    }
    if(pos_selection + Selection.length() < full_cut_string.length()){
    if(std::find(valid_selection_chars.begin(), valid_selection_chars.end(), full_cut_string[pos_selection + Selection.length()]) != valid_selection_chars.end()) valid_end = true; 
    }
    if( !(valid_start && valid_end)) continue;
    size_t pos_colon = full_cut_string.find(':');
    if(pos_colon  == std::string::npos) continue;
    Full_Selection = full_cut_string.substr(pos_colon + 1, full_cut_string.size());
    for(string selection_type : CutTypesVec){
    int str_length = selection_type.length();
    // size_t pos = Full_Selection.find(selection_type,pos_selection+Selection.length());
    size_t pos = Full_Selection.find(selection_type);
    if (pos != std::string::npos)
    {
       valid_start = false;
       valid_end = false;
       if(pos == 0) valid_start = true;
       if(pos != 0){
       if(std::find(valid_chars.begin(), valid_chars.end(), Full_Selection[pos -1]) != valid_chars.end()) valid_start = true;
       }
       if(pos + str_length >= Full_Selection.length())valid_end = true;
       if(pos + str_length < Full_Selection.length()){
        if(std::find(valid_chars.begin(), valid_chars.end(), Full_Selection[pos + str_length]) != valid_chars.end())valid_end = true;         
       }
       
       
       if( !(valid_start && valid_end)) continue;
       CutTypes_contained.push_back(selection_type);
       string replacement_string = "(" + CutTypesMap[selection_type] + ")";
       // Full_Selection.replace(pos,selection_type.length(),replacement_string);
       boost::algorithm::replace_all(Full_Selection, selection_type , replacement_string);
       selection_num ++;  
    }

    }
    
    // Writes the Selection file description into the location where histograms are written
    SelectionTypesMap[Selection] = Full_Selection;
    ofstream myfile;
    myfile.open(folder_loc + "/" + Selection +"_selection_info.txt");
    myfile << Selection + "\n";
    myfile << "Unprocessed Selection is " + full_cut_string + "\n";
    myfile << "has Selections : \n";
      myfile << "Universal : " + CutTypesMap["Universal"] + "\n";    
    for(string Cut : CutTypes_contained){
      myfile << Cut + " : " + CutTypesMap[Cut] + "\n";
    }
    myfile << " Full_Selection : \n";
    myfile << Full_Selection + "\n";
    myfile.close();
    return;
    }
}
}

// Main Function input variables are config file, and the list of selections used and datatypes to be analyzed
void CreateHistograms_CustomBins(std::string config_file_name, std::string list_select = "list_selections.txt", std::string list_data = "list_datatypes.txt"){
    string file_line;
    vector<string> config_inputs = {};
    vector<string> datatypes = {};
    vector<string> selections = {};
    
    list_select = CleanString(list_select);
    list_data = CleanString(list_data);
    config_file_name = CleanString(config_file_name);


    // sets up config
    ifstream file_config(config_file_name);
    while(getline(file_config, file_line)){
     config_inputs.push_back(CleanString(file_line));   
    }
    // minimum number of config variables to run CreateHistogram Macro, Note that other Macros may require more 
    // variables but the variables will always be in same line number to allow for a universal config file
    if(config_inputs.size() < 7){ cout << " Bad config file " << endl; return;}

   prong_type = config_inputs.at(2);
   save_folder = config_inputs.at(1)+ prong_type+ "/";
   infile_name_beginning = config_inputs.at(0)+ prong_type+ "/";
   string paraTree_filename = save_folder + config_inputs.at(3);
   string paraTree_filename_2D = save_folder + config_inputs.at(4);
   string cuts_types_file = save_folder + config_inputs.at(5);
   string selection_file = save_folder + config_inputs.at(6);
   
   
   
   // Gets list of selections and datatypes. assumes the selection and datatypes are
   // in the save_folder location if absolute path isn't given. save_folder in config file 
   
    if( list_select.at(0) != '/' ) list_select = save_folder + list_select;
    if( list_data.at(0) != '/') list_data = save_folder + list_data;
    ifstream file_list_select(list_select);
    while(getline(file_list_select, file_line)){
      selections.push_back(CleanString(file_line));   
    }
   
    ifstream file_list_data(list_data);
    while(getline(file_list_data, file_line)){
      datatypes.push_back(CleanString(file_line));   
    }


// ParameterTree is the histogram parameters for 1D histograms. CSV file provided in config file

TTree *ParameterTree = new TTree("ParameterTree", "tree to fill hist parameters");
ParameterTree->ReadFile(paraTree_filename.c_str(),"KeepHists/I:BranchName/C:HistName/C:HistTitle/C:CustomBinsName/C:nBins/I:min/D:max/D");
int KeepHists;
Char_t HistName[1000]; 
Char_t HistTitle[1000]; 
Char_t BranchName[1000]; 
Char_t CustomBinsName[1000]; 
int nBins;
double min;
double max;
ParameterTree->SetBranchAddress("KeepHists", &KeepHists);
ParameterTree->SetBranchAddress("BranchName", &BranchName);
ParameterTree->SetBranchAddress("HistName", &HistName);
ParameterTree->SetBranchAddress("HistTitle", &HistTitle);
ParameterTree->SetBranchAddress("CustomBinsName", &CustomBinsName);
ParameterTree->SetBranchAddress("nBins", &nBins);
ParameterTree->SetBranchAddress("min", &min);
ParameterTree->SetBranchAddress("max", &max); 


// 2d hists
// ParameterTree_2D is the histogram parameters for 2D histograms. CSV file provided in config file

TTree *ParameterTree_2D = new TTree("ParameterTree", "tree to fill hist parameters");
ParameterTree_2D->ReadFile(paraTree_filename_2D.c_str(),"KeepHists/I:BranchName_X/C:BranchName_Y/C:HistName/C:HistTitle/C:CustomBinsName_X/C:nBins_X/I:min_X/D:max_X/D:CustomBinsName_Y/C:nBins_Y/I:min_Y/D:max_Y/D");
int KeepHists_2D;
Char_t HistName_2D[1000]; 
Char_t HistTitle_2D[1000]; 
Char_t BranchName_X_2D[1000]; 
Char_t BranchName_Y_2D[1000];
Char_t CustomBinsName_X_2D[1000];  
int nBins_X_2D;
double min_X_2D;
double max_X_2D;
Char_t CustomBinsName_Y_2D[1000]; 
int nBins_Y_2D;
double min_Y_2D;
double max_Y_2D;
ParameterTree_2D->SetBranchAddress("KeepHists", &KeepHists_2D);
ParameterTree_2D->SetBranchAddress("BranchName_X", &BranchName_X_2D);
ParameterTree_2D->SetBranchAddress("BranchName_Y", &BranchName_Y_2D);
ParameterTree_2D->SetBranchAddress("HistName", &HistName_2D);
ParameterTree_2D->SetBranchAddress("HistTitle", &HistTitle_2D);
ParameterTree_2D->SetBranchAddress("CustomBinsName_X", &CustomBinsName_X_2D);
ParameterTree_2D->SetBranchAddress("nBins_X", &nBins_X_2D);
ParameterTree_2D->SetBranchAddress("min_X", &min_X_2D);
ParameterTree_2D->SetBranchAddress("max_X", &max_X_2D); 
ParameterTree_2D->SetBranchAddress("CustomBinsName_Y", &CustomBinsName_Y_2D);
ParameterTree_2D->SetBranchAddress("nBins_Y", &nBins_Y_2D);
ParameterTree_2D->SetBranchAddress("min_Y", &min_Y_2D);
ParameterTree_2D->SetBranchAddress("max_Y", &max_Y_2D); 



MakeCutTypeMap(cuts_types_file);


// Makes the selection map and creates folder for each of selections.
   for( string selection : selections){
      string folder_loc = save_folder +  RJR_Method +"/" +selection;
      DIR* dir = opendir(folder_loc.c_str());
      if (dir) {
       cout << "Folder Exists " << endl;
       closedir(dir);
      } else if (ENOENT == errno) {
      int status;     
      string mkdir_string = "mkdir -p " + folder_loc;
      status = system(mkdir_string.c_str()); // Creating a directory
      if (status == -1)
         cerr << "Error : " << strerror(errno) << endl;
      else
         cout << "Directories are created: " << folder_loc << endl;
      }
      MakeSelectionMap(selection_file, selection, folder_loc); 
   }

// loops over datatype
for( string datatype : datatypes){
      // data input location. Assumes input file has format MinimalCuts_ + dataname + .root
      // this should be changed to accomidate your naming convention. 
      string infile_name = infile_name_beginning + RJR_Method +"/MinimalCuts_"+ datatype + ".root";
      TFile *_file0 = new TFile(infile_name.c_str(), "READ");

      // This is where the tree name is determined. Also should be changed to your root files
      string TreeName;
      if(prong_type == "1Prong") TreeName = "Events";
      // if(prong_type == "1Prong") TreeName = "GenEfficiency";
      if(prong_type == "3Prong") TreeName = "Events_3prong";
     
      // root file is loaded into RDataFrame     
      ROOT::RDataFrame d(TreeName, _file0); 
      if(CutTypesMap["Universal"] == ""){
         cout << "Null universal selection " << endl;
         continue;}  
      // universal selection is applied         
      auto d1 = d.Filter(CutTypesMap["Universal"], "df with selection 1");
      
      // DataFrame adds new variables here. 
         
      auto d2 = AddVariables(d1, TreeName);
      
      // loops over selectinos and created histograms
      for( string selection : selections){
         cout << "Data Type : " << datatype << endl;
         cout << " Selection : " << selection << endl;
         string folder_loc = save_folder +  RJR_Method +"/" +selection;
         if(SelectionTypesMap[selection] == ""){
         cout << "Null selection :  " << selection << endl;
         continue;}  
         // final selections are applied 
         auto d3 = d2.Filter(SelectionTypesMap[selection], "df with selection");
         
         // output histogram name is created. Can be changed but down the line macros should be changed as well.
         
         string outfile_name = folder_loc+"/MinCuts_" + datatype+".root";
         TFile *outfile = new TFile(outfile_name.c_str(), "recreate"); 
         
         // creates 1D histograms
         if(do1DHists){         
         for(int i = 0; i < ParameterTree->GetEntries(); i++){
            ParameterTree->GetEntry(i);
            if(KeepHists==0) continue;
            // cout << HistName << endl;
            string custom_bin_name = CustomBinsName;
            if(custom_bin_name == "NA"){
            auto myHist1 = d3.Histo1D({HistName, HistTitle, nBins, min, max}, BranchName);
            myHist1->Write();
            }
            if(custom_bin_name != "NA"){
            int size = constants::mapvectors[custom_bin_name].size();
            double testarray[size];
            for(int j =0; j < size; j++) testarray[j] = constants::mapvectors[custom_bin_name].at(j);
            auto myHist1 = d3.Histo1D({HistName, HistTitle, size -1, testarray}, BranchName);
            myHist1->Write();
            }
         if(i == ParameterTree->GetEntries()-1) cout << "Final 1D hist done " << endl;
        }
         }
         
         // Creates 2D histograms
         if(do2DHists){
         for(int i = 0; i < ParameterTree_2D->GetEntries(); i++){
            ParameterTree_2D->GetEntry(i);
            if(KeepHists_2D==0) continue;
            // cout << HistName_2D << endl;
            string custom_bin_name_X = CustomBinsName_X_2D;
            string custom_bin_name_Y = CustomBinsName_Y_2D;
            string hist2d_name = string("Hist2D_") +HistName_2D;
            string profile_name = string("Profile_") +HistName_2D;
            if(custom_bin_name_X == "NA"){
            if(custom_bin_name_Y == "NA"){
            auto myHist2D = d3.Histo2D({hist2d_name.c_str(),HistTitle_2D, nBins_X_2D, min_X_2D, max_X_2D, nBins_Y_2D, min_Y_2D, max_Y_2D}, BranchName_X_2D,BranchName_Y_2D);
            myHist2D->Write();
            }
            if(custom_bin_name_Y != "NA"){
            int size_Y = constants::mapvectors[custom_bin_name_Y].size();
            double testarray_Y[size_Y];
            for(int j =0; j < size_Y; j++) testarray_Y[j] = constants::mapvectors[custom_bin_name_Y].at(j);
            auto myHist2D = d3.Histo2D({hist2d_name.c_str(), HistTitle_2D, nBins_X_2D, min_X_2D, max_X_2D, size_Y -1 , testarray_Y}, BranchName_X_2D,BranchName_Y_2D);
            myHist2D->Write();
            }
            auto myProfile = d3.Profile1D({profile_name.c_str(), HistTitle_2D, nBins_X_2D, min_X_2D, max_X_2D}, BranchName_X_2D,BranchName_Y_2D);
            myProfile->Write();
            }
            if(custom_bin_name_X != "NA"){
            int size_X = constants::mapvectors[custom_bin_name_X].size();
            double testarray_X[size_X];
            for(int j =0; j < size_X; j++) testarray_X[j] = constants::mapvectors[custom_bin_name_X].at(j);
            if(custom_bin_name_Y == "NA"){
            auto myHist2D = d3.Histo2D({hist2d_name.c_str(), HistTitle_2D, size_X -1 , testarray_X, nBins_Y_2D, min_Y_2D, max_Y_2D}, BranchName_X_2D,BranchName_Y_2D);
            myHist2D->Write();
            }
            if(custom_bin_name_Y != "NA"){
            int size_Y = constants::mapvectors[custom_bin_name_Y].size();
            double testarray_Y[size_Y];
            for(int j =0; j < size_Y; j++) testarray_Y[j] = constants::mapvectors[custom_bin_name_Y].at(j);
            auto myHist2D = d3.Histo2D({hist2d_name.c_str(), HistTitle_2D, size_X -1 , testarray_X, size_Y -1 , testarray_Y}, BranchName_X_2D,BranchName_Y_2D);
            myHist2D->Write();
            }
            auto myProfile = d3.Profile1D({profile_name.c_str(), HistTitle_2D, size_X -1 , testarray_X}, BranchName_X_2D,BranchName_Y_2D);
            myProfile->Write();
            }
         if(i == ParameterTree_2D->GetEntries()-1) cout << "Final 2D hist done " << endl;
        } 
         }
         
         // closes output file after all histograms are created
        outfile->Close();
      }
      
      // closes datafile after all selections are applied to data. 
      _file0->Close();
}



}