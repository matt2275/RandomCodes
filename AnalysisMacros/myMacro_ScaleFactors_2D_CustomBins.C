#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

#include "TCanvas.h"
#include "TStyle.h"
#include "TGaxis.h"
#include "TRandom.h"

#include <iomanip>
#include <sstream>

#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <dirent.h>
#include <errno.h>

 using namespace std;
 
  float T_margin = 0.08;
  float B_margin = 0.12; 
  float L_margin = 0.14;
  float R_margin = 0.04;

// vector<string> selections ={};
vector<string> RJR_Methods = {};

string prong_type = "";
string save_folder = "";

 vector<string> DataTypes = {};
 
vector<string>  selections = {};
 
string DENOM = ""; // denominator of efficiency calculation. selections is nominator

string SF_NOM_datatype = ""; // numerator of scale type. usually data


string outName_start = "";

 
bool HasRatioPlot = false;
string normalization = "XS";
float epsi = 0.000000000001;
float Lumi = .432553397338;
int num_bins = 100;
float min_bin = 0;
float max_bin = 40;
float hist_max = 5000;
int num_bins_X = 100;
float min_bin_X = 0;
float max_bin_X = 40;
float hist_max_X = 5000;
int num_bins_Y = 100;
float min_bin_Y = 0;
float max_bin_Y = 40;
float hist_max_Y = 5000;
TString x_title = "m_{lep^{+}lep^{-}} (GeV)";
TString y_title = "Events / 0.4 GeV";
TString hist_name = "pair_mass";
TString outName = "pair_mass_hist";
string DoLogScale = "no";
bool DoLegend = true;
int label_pos = 0;
float x_legend_max = .92;
float x_legend_width = .30;
float y_legend_max = .3;
float y_legend_width = .18;


// "Data","DiTau_SuperChic","DiTau_Madgraph","BBbar","CCbar","DiMuon_Gamma","DiElectron_SuperChic","DiElectron_Starlight","DiMuon_noFSR"
 // if(prong_type == "3Prong" ) selections = {"El3Prong","Mu3Prong"};
 // if(prong_type == "1Prong" ) selections = {"MuMu_0n0n","MuMu_1n0n","MuMu_2n0n","MuMu_0n1n","MuMu_0n2n","MuMu_1n","MuMu_2n"};
 // selections = {"MuMu_0n0n","MuMu_1n0n","MuMu_2n0n","MuMu_0n1n","MuMu_0n2n","MuMu_1n","MuMu_2n"};
 // vector<string> DataTypes = {"Data","DiElectron_SuperChic","DiElectron_Starlight"};
 // vector<string> DataTypes = {"Data","DiTau_SuperChic","BBbar","CCbar","DiElectron_SuperChic","DiMuon_Gamma","DiMuon_noFSR"};
 // vector<string> MCTypes = {"DiTau_SuperChic","BBbar","CCbar","DiElectron_SuperChic","DiMuon_Gamma","DiMuon_noFSR"};
 // vector<string> MCTypes = {"DiMuon_noFSR","DiMuon_Arash"};
map< string, tuple<bool, TString, int, int, float>> histParams;
 map<string,TString> FancyNames;


string Global_dataType;
string Global_selection;
// choose colors with https://root.cern/doc/master/classTColor.html
// choose markers with https://root.cern.ch/doc/master/classTAttMarker.html


void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  // pad1->Draw();  
  pad1->cd();  pad1->SetLeftMargin(L_margin);  pad1->SetBottomMargin(B_margin); pad1->SetRightMargin(R_margin); pad1->SetTopMargin(T_margin);
  //pad1->SetGrid(); 
  //pad1->SetGrid(10,10); 
}

string CleanString(string dirty_string){
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\n'), dirty_string.cend()); 
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\r'), dirty_string.cend());
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), ' '), dirty_string.cend());
  return dirty_string;
}

TGraphErrors* Hists_Diff(TH1* hist1, TH1* hist2){
	   int nbins = hist1->GetNbinsX();
      cout << " Nbins ratio diff " << nbins << endl;
       double X_val[nbins];
       double X_err[nbins];
       double Y_val[nbins];
       double Y_err[nbins];     	   	   
       double val1;
       double val2;
       double sum;
       for (Int_t i=0; i<=nbins ;i++) {
          val1 = hist1->GetBinContent(i+1);
		  val2 = hist2->GetBinContent(i+1);	
		  sum = val1 + val2;
		  X_val[i] = hist1->GetXaxis()->GetBinCenter(i+1);
		  X_err[i] = 0;
		  Y_val[i] = (val1 - val2) / (sum + epsi);
		  Y_err[i] = 4 * val1 * val2/ (sum*sum*sum + epsi);	
	   }

	 TGraphErrors *gr = new TGraphErrors(nbins, X_val, Y_val,X_err,Y_err);
	 return gr;
}


TGraphErrors* Hists_Ratio(TH1* hist1, TH1* hist2){
	   int nbins = hist1->GetNbinsX();
       double X_val[nbins];
       double X_err[nbins];
       double Y_val[nbins];
       double Y_err[nbins];     	   	   
       double val1;
       double val2;
       double sum;
       for (Int_t i=0; i<=nbins -1 ;i++) {
          val1 = hist1->GetBinContent(i + 1);
		  val2 = hist2->GetBinContent(i + 1);	
		  sum = val1 + val2;
		  X_val[i] = hist1->GetXaxis()->GetBinCenter(i);
		  X_err[i] = 0;
		  Y_val[i] = (val1 /( val2 + epsi));
		  Y_err[i] = (val1 /( val2 + epsi))*(1/(val1+epsi) + 1/(val2+epsi),.5);	
	   }

	 TGraphErrors *gr = new TGraphErrors(nbins, X_val, Y_val,X_err,Y_err);
	 return gr;
}


// TH2* GetScaleFactor_2D(TH2*  hist1, TH2* hist2){
      // int nbins_x = hist1->GetNbinsX();
      // int nbins_y = hist1->GetNbinsY();
       // double X_val[nbins];
       // double X_err[nbins];
       // double Y_val[nbins];
       // double Y_err[nbins];     	   	   
       // double val1;
       // double err1;
       // double rel_err1;
       // double val2;
       // double err2;
       // double rel_err2;
       // double sum;
       // for (Int_t i=0; i<=nbins_x +1 ;i++) {
         // for (Int_t j=0; j<=nbins_y +1 ;j++) {          
        // val1 = hist1->GetBinContent(i,j);
		  // val2 = hist2->GetBinContent(i,j);
        // // if(error_type == "down"){
        // // err1 = hist1->GetEfficiencyErrorLow(i,j);
		  // // err2 = hist2->GetEfficiencyErrorLow(i,j);           
        // // }
        // // else if(error_type == "up"){
        // // err1 = hist1->GetEfficiencyErrorUp(i,j);
		  // // err2 = hist2->GetEfficiencyErrorUp(i,j);              
        // // }
        // // else{
        // // err1 = 0;
		  // // err2 = 0;   
        // // }
        // // rel_err1 = pow(err1/val1,2);
        // // rel_err2 = pow(err2/val2,2);
		  // sum = val1 + val2;
		  // // X_val[i] = copy->GetXaxis()->GetBinCenter(i);
		  // // X_err[i] = copy->GetXaxis()->GetBinWidth(i)/2;
		  // hist1->SetBinContent(i,j,1,(val1 /( val2 + epsi)));
		  // // copy->SetBinError(i,j,(val1 /( val2 + epsi))*pow(rel_err1 + rel_err2, .5));	
	     // }
       // }

	 // // TGraphErrors *gr = new TGraphErrors(nbins, X_val, Y_val,X_err,Y_err);
	 // return hist1;
// }

TString replace_first(
    TString char_array,
    std::string const& toReplace,
    std::string const& replaceWith
) {
    std::string s(char_array);
    // cout << s << endl;
    std::size_t pos = s.find(toReplace);
    if (pos == std::string::npos) return char_array;
    s.replace(pos, toReplace.length(), replaceWith);
    // cout << s <<endl;
    return s.c_str();
}
void fillMaps(string prong_type, string RJR_Method, string selection){
   string Hists_Folder = save_folder;
// data type, isMC ,file, color,  marker type, scale	(Lumi*Cross section nb/ Number of Events Generated)   
histParams["Data"] = {false ,Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_Data.root" , 1 ,20, 1};
histParams["DiTau_SuperChic"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiTau_SuperChic.root" ,860 ,21, Lumi*8.83E6/9.97/1184600};
histParams["DiTau_Madgraph"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiTau_Madgraph.root" ,860 ,21, Lumi*570.0E3/950000};
histParams["BBbar"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_BBbar.root" ,397 ,66, Lumi*1.5E3/20000000*1};
histParams["CCbar"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_CCbar.root" ,419 ,44, Lumi*300.0E3/19976000*1};
histParams["DiElectron_SuperChic"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiElectron_SuperChic.root" ,910 ,33, Lumi*8.83E6/67810000/2};
histParams["DiElectron_Starlight"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiElectron_Starlight.root" ,614 ,32, Lumi*7.92E6/66750000/2};
histParams["DiMuon_Gamma"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Gamma.root" ,890 ,21, Lumi*140E3/2390000};
histParams["DiMuon_noFSR"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_noFSR.root" ,880 ,22, Lumi*189.08E3/30000000*17.97 };
histParams["DiMuon_Starlight"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Starlight.root" ,880 ,22, Lumi*189.08E3/30000000*17.97 };
histParams["DiMuon_Private"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Private.root" ,880 ,22, Lumi*189.08E3/30000000*17.97 };
histParams["DiMuon_Arash"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Arash.root" ,871 ,23, Lumi*8.83E6/250000/2.5 };
 // DiMuon_FSR 17.97 comes from only taking 1M of 18 M events passed primary selections
 // DiTau_SuperChic X-sec comes from DiElectron_SuperChic / (2*M_tau / 2( mass cut ee))^4 , divided by 2 for fudge factor
 // DiElectron_Starlight X-sec Divided by 9.11 for Single EG3 Prescale, times 3 for fudge factor
 // DiElectron_SuperChic X-sec Divided by 9.11 for Single EG3 Prescale,times 3 for fudge factor

// data type, Fancy Name (latex version)

FancyNames["Data"] = "Data";
FancyNames["DiTau_SuperChic"] = "#gamma #gamma #rightarrow #tau^{+}#tau^{-} (MC SC)";
FancyNames["DiTau_Madgraph"] = "#gamma #gamma #rightarrow #tau^{+}#tau^{-} (MC MG)";
FancyNames["BBbar"] = "B #bar{B} (MC)";
FancyNames["CCbar"] = "C #bar{C} (MC)";
FancyNames["DiElectron_SuperChic"] = "#gamma #gamma #rightarrow e^{+}e^{-} (MC SC)";
FancyNames["DiElectron_Starlight"] = "#gamma #gamma #rightarrow e^{+}e^{-} (MC SL)";
FancyNames["DiMuon_Gamma"] = "#gamma #gamma #rightarrow #mu^{+}#mu^{-} + #gamma (MC SL)";
FancyNames["DiMuon_noFSR"] = "#gamma #gamma #rightarrow #mu^{+}#mu^{-} (MC)";
FancyNames["DiMuon_Starlight"] = "#gamma #gamma #rightarrow #mu^{+}#mu^{-} (MC)";
FancyNames["DiMuon_Private"] = "#gamma #gamma #rightarrow #mu^{+}#mu^{-} (MC)";
FancyNames["DiMuon_Arash"] = "#gamma #gamma #rightarrow #mu^{+}#mu^{-} (MC Arash)";
}


// fills map for getting histogram parameters using csv file Note there is a modified version in a different macro for plotting multiple selections in same plot
void fillMaps_CSV(string RJR_Method,string selection){
   // string Hists_Folder = save_folder;
   
string csvname = save_folder + "DataSamplesInfo.csv";
TTree *SampleInfoTree = new TTree("SampleInfoTree", "tree to sample info");
SampleInfoTree->ReadFile(csvname.c_str(),"SampleName/C:SampleIsMC/I:SampleFileName/C:FillColor/I:MarkerType/I:CrossSection/D:NumofEvents/D:FudgeFactor/D:FancyName/C");
Char_t SampleName[1000]; 
int SampleIsMC;
Char_t SampleFileName[1000]; 
int FillColor;
int MarkerType;
double CrossSection;
double NumofEvents;
double FudgeFactor;
Char_t FancyName[1000]; 
SampleInfoTree->SetBranchAddress("SampleName", &SampleName);
SampleInfoTree->SetBranchAddress("SampleIsMC", &SampleIsMC);
SampleInfoTree->SetBranchAddress("SampleFileName", &SampleFileName);
SampleInfoTree->SetBranchAddress("FillColor", &FillColor);
SampleInfoTree->SetBranchAddress("MarkerType", &MarkerType);
SampleInfoTree->SetBranchAddress("CrossSection", &CrossSection);
SampleInfoTree->SetBranchAddress("NumofEvents", &NumofEvents);
SampleInfoTree->SetBranchAddress("FudgeFactor", &FudgeFactor);
SampleInfoTree->SetBranchAddress("FancyName", &FancyName);

for(int i = 0; i < SampleInfoTree->GetEntries(); i++){
    SampleInfoTree->GetEntry(i);
    // string dataname = convertToString(SampleName);
    // cout << " SampleName " << SampleName << endl;
    // cout << " SampleIsMC " << (SampleIsMC==1) << endl;
    // cout << " SampleFileName " << SampleFileName << endl;
    // cout << " FillColor " << FillColor << endl;
    // cout << " MarkerType " << MarkerType << endl;
    // cout << " CrossSection " << CrossSection << endl;
    // cout << " NumofEvents " << NumofEvents << endl;
    // cout << " FudgeFactor " << FudgeFactor << endl;
    // cout << " FancyName " << FancyName << endl;
    // cout << save_folder +  RJR_Method +"/" +selection+ SampleFileName << endl;
// data type, isMC ,file, color,  marker type, scale	(Lumi*Cross section nb/ Number of Events Generated)   
if(SampleIsMC==1) histParams[convertToString(SampleName)] = {(SampleIsMC==1) ,save_folder +  RJR_Method +"/" +selection+ SampleFileName , FillColor ,MarkerType, Lumi*FudgeFactor*CrossSection/NumofEvents};
if(SampleIsMC!=1) histParams[convertToString(SampleName)] = {(SampleIsMC==1) ,save_folder +  RJR_Method +"/" +selection+ SampleFileName , FillColor ,MarkerType, FudgeFactor*CrossSection/NumofEvents};

// cout << dataname << "  :  " << histParams[dataname] << endl;
// data type, Fancy Name (latex version)

FancyNames[convertToString(SampleName)] = FancyName;

}
}


TCanvas* example_plot( int iPeriod, int iPos );
 
 // vector<string> selections = {"Muon_Trigger_All","Muon_Trigger_Pass","Muon_Trigger_Fail","Electron_EG3_Trigger_All","Electron_EG3_Trigger_Pass","Electron_EG3_Trigger_Fail", "Electron_EG5_Trigger_All","Electron_EG5_Trigger_Pass","Electron_EG5_Trigger_Fail"};
 
void myMacro_ScaleFactors_2D_CustomBins(std::string config_file_name, std::string list_select = "list_selections.txt", std::string list_data = "list_datatypes.txt"){
    string file_line;
    vector<string> config_inputs = {};
    // vector<string> datatypes = {};
    // vector<string> selections = {};
    
    list_select = CleanString(list_select);
    list_data = CleanString(list_data);
    config_file_name = CleanString(config_file_name);

    ifstream file_config(config_file_name);
    while(getline(file_config, file_line)){
     config_inputs.push_back(CleanString(file_line));   
    }
    if(config_inputs.size() < 10){ cout << " Bad config file " << endl; return;}

   prong_type = config_inputs.at(2);
   save_folder = config_inputs.at(1)+ prong_type+ "/";
   string infile_name_beginning = config_inputs.at(0)+ prong_type+ "/";
   string paraTree_filename = save_folder + config_inputs.at(3);
   string paraTree_filename_2D = save_folder + config_inputs.at(4);
   string cuts_types_file = save_folder + config_inputs.at(5);
   string selection_file = save_folder + config_inputs.at(6);
   outName_start = config_inputs.at(7);
   DENOM = config_inputs.at(8); // denominator of efficiency calculation. selections is nominator
   SF_NOM_datatype = config_inputs.at(9); // numerator of scale type. usually data

   
   // string RJR_Method = "HasVisMass";
   
    if( list_select.at(0) != '/' ) list_select = save_folder + list_select;
    if( list_data.at(0) != '/') list_data = save_folder + list_data;
    ifstream file_list_select(list_select);
    while(getline(file_list_select, file_line)){
      selections.push_back(CleanString(file_line));   
    }
   
    ifstream file_list_data(list_data);
    while(getline(file_list_data, file_line)){
      DataTypes.push_back(CleanString(file_line));   
    }
    
 if(normalization=="Unity") HasRatioPlot = false;
   
 if(prong_type == "3Prong" )RJR_Methods ={"HasVisMass"};
 if(prong_type == "1Prong" )RJR_Methods ={"HasVisMass"};
 
 

TTree *ParameterTree_2D = new TTree("ParameterTree", "tree to fill hist parameters");
ParameterTree_2D->ReadFile(paraTree_filename_2D.c_str(),"KeepHists/I:BranchName_X/C:BranchName_Y/C:HistName/C:HistTitle/C:CustomBinsName_X/C:nBins_X/I:min_X/D:max_X/D:CustomBinsName_Y/C:nBins_Y/I:min_Y/D:max_Y/D:xaxis/C:yaxis/C");
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
Char_t xaxis_2D[1000]; 
Char_t yaxis_2D[1000];
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
ParameterTree_2D->SetBranchAddress("xaxis", &xaxis_2D); 
ParameterTree_2D->SetBranchAddress("yaxis", &yaxis_2D); 


// (lep , #tau _{3prong})
  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  
for( string RJR_Method  :  RJR_Methods ){ 
for( string selection : selections){

// fillMaps(prong_type, RJR_Method, selection);
fillMaps_CSV(RJR_Method, selection);

  //  gROOT->LoadMacro("CMS_lumi.C");
std::stringstream stream;
stream << std::fixed << std::setprecision(2) << Lumi;
string cms_text = "#sqrt{s_{NN}} = 5.02 TeV ( " + stream.str() +  " nb^{-1} )";
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = cms_text;       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  // second parameter in example_plot is iPos, which drives the position of the CMS logo in the plot
  // iPos=11 : top-left, left-aligned
  // iPos=33 : top-right, right-aligned
  // iPos=22 : center, centered
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
   for(string dataType : DataTypes){
      Global_dataType = dataType;
      Global_selection = selection;
      
      
    string out_folder;
   TString eff_folder = "SF_" + selection + "_OVER_" + DENOM;
   out_folder = save_folder +  RJR_Method+ "/" +eff_folder + "/" + Global_dataType;
   
   
   DIR* dir = opendir(out_folder.c_str());
   if (dir) {
    cout << "Folder Exists " << endl;
    closedir(dir);
   } else if (ENOENT == errno) {
   int status;     
   string mkdir_string = "mkdir -p " + out_folder;
   status = system(mkdir_string.c_str()); // Creating a directory
   if (status == -1)
      cerr << "Error : " << strerror(errno) << endl;
   else
      cout << "Directories are created: " << out_folder << endl;
   }
   
for(int i = 0; i < ParameterTree_2D->GetEntries(); i++){
    ParameterTree_2D->GetEntry(i);
   // xaxis = replace_first(xaxis,"#comma", ",");
   // yaxis = replace_first(yaxis,"#comma", ",");
    num_bins_X = nBins_X_2D;
    min_bin_X = min_X_2D;
    max_bin_X = max_X_2D;
    num_bins_Y = nBins_Y_2D;
    min_bin_Y = min_Y_2D;
    max_bin_Y = max_Y_2D;
    x_title = replace_first(xaxis_2D,"#comma", ",");
    y_title = replace_first(yaxis_2D,"#comma", ",");
    hist_name = HistName_2D;

    // outName = out_folder + "/" + hist_name;

    DoLogScale = "both";

    outName = out_folder + "/" + Global_dataType + "_"+ hist_name;
      cout << Global_dataType << "  " << Global_selection << endl;
   example_plot( iPeriod, label_pos ); 
   }
        
}
}
}
}

TCanvas* example_plot( int iPeriod, int iPos )
{ 
  //  if( iPos==0 ) relPosX = 0.12;
  int W = 800;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 600; 
  int W_ref = 800; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.14*W_ref;
  float R = 0.04*W_ref;
  
  TString canvName = "FigExample_";
  canvName += W;
  canvName += "-";
  canvName += H;
  canvName += "_";  
  canvName += iPeriod;
  if( writeExtraText ) canvName += "-prelim";
  if( iPos%10==0 ) canvName += "-out";
  else if( iPos%10==1 ) canvName += "-left";
  else if( iPos%10==2 )  canvName += "-center";
  else if( iPos%10==3 )  canvName += "-right";

  TCanvas* canv = new TCanvas(canvName,canvName,0,0,W,H);
  canv->SetFillColor(0);
  // canv->SetBorderMode(0);
  // canv->SetFrameFillStyle(0);
  // canv->SetFrameBorderMode(0);
  // canv->SetLeftMargin( L/W );
  // canv->SetRightMargin( R/W );
  // canv->SetTopMargin( T/H );
  // canv->SetBottomMargin( B/H );
  // canv->SetTickx(0);
  // canv->SetTicky(0);
  
  // TCanvas* c1 = new TCanvas("c1","A Simple Graph with asymmetric error bars",1200,1200);
   TPad* pad1;
    // TPad* pad2;
  // c1->SetFillColor(0);
  //c1->SetGrid();

  pad1 = new TPad("pad1","The pad",0,0,1,1);
    pad1->Draw(); 
  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);
        pad1->SetLeftMargin(.1);
         pad1->SetRightMargin(.15);
         pad1->SetBottomMargin(.1);
         pad1->SetTopMargin(0.1);
  canv->cd();
  // pad1->Draw();
  // pad2->Draw();
  pad1->cd();
    TString hist_name_2D = "Hist2D_" + hist_name;
auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[SF_NOM_datatype]; 

   TFile file_NOM_Data(FileName,"READ");
    cout << FileName << " "  << hist_name << endl;

    TH2F *h_Data = static_cast<TH2F*>(file_NOM_Data.Get(hist_name_2D)->Clone());
    TString FileName_DENOM_Data = replace_first(FileName,Global_selection, DENOM);


    TFile file_DENOM_Data(FileName_DENOM_Data,"READ");
    TH2F *tmp_hist_DENOM_Data= static_cast<TH2F*>(file_DENOM_Data.Get(hist_name_2D)->Clone());
    h_Data->Divide(tmp_hist_DENOM_Data);
 
   auto &[isMC_MC, FileName_MC, fill_color_MC, marker_style_MC, scale_MC]  = histParams[Global_dataType];
   TFile file_NOM_MC(FileName_MC,"READ");
    cout << FileName_MC<< " "  << hist_name << endl;

    TH2F *h_MC = static_cast<TH2F*>(file_NOM_MC.Get(hist_name_2D)->Clone());
    TString FileName_DENOM_MC = replace_first(FileName_MC,Global_selection, DENOM);


    TFile file_DENOM_MC(FileName_DENOM_MC,"READ");
    TH2F *tmp_hist_DENOM_MC= static_cast<TH2F*>(file_DENOM_MC.Get(hist_name_2D)->Clone());
    h_MC->Divide(tmp_hist_DENOM_MC);
    h_Data->Divide(h_MC);
  // TH2* h = new TH2F("h","h",num_bins_X,min_bin_X,max_bin_X,num_bins_Y,min_bin_Y,max_bin_Y);
  // h->GetXaxis()->SetNdivisions(6,5,0);
  h_Data->GetXaxis()->SetTitle(x_title);  
  // h->GetYaxis()->SetNdivisions(6,5,0);
  h_Data->GetYaxis()->SetTitleOffset(1.25);
  // TString Eff_y_title = "Efficiency of " + x_title;
  h_Data->GetYaxis()->SetTitle(y_title);
  // h->GetXaxis()->SetTitleSize(.05);   
  // h->GetYaxis()->SetTitleSize(.05); 



// h->GetYaxis()->SetNdivisions(6,5,0);
// h->GetXaxis()->SetNdivisions(6,5,0);



  // h->SetMaximum(1.05);
  // h->SetMinimum(0);
  // if( iPos==1 ) h->SetMaximum( 10000 );
  h_Data->Draw("colz");
  



  
  

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos );
 // outName = save_folder + prong_type +"/" + "HasVisMass/Electron_EG3_Trigger_Efficiency/" + hist_name;
  canv->Update();
  // canv->RedrawAxis();
  canv->GetFrame()->Draw();
  if(DoLogScale == "yes"){
  pad1->SetLogz(1);
  canv->Print(outName+"_log.pdf",".pdf");
  canv->Print(outName+"_log.png",".png");
  }
  if(DoLogScale == "no"){
  pad1->SetLogz(0);
  canv->Print(outName+".pdf",".pdf");
  canv->Print(outName+".png",".png");
  }  
   if(DoLogScale == "both"){
  pad1->SetLogz(0);
  canv->Print(outName+".pdf",".pdf");
  canv->Print(outName+".png",".png");
  pad1->SetLogz(1);
  canv->Print(outName+"_log.pdf",".pdf");
  canv->Print(outName+"_log.png",".png");
  }
  
  // canv->Close();


  return canv;
}