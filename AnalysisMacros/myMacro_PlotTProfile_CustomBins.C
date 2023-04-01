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


#include "variables.h"

 using namespace std;
 
string save_folder = "";
string prong_type = "";

string normalization = "";

vector<string> DataTypes = {};
vector<string> selections ={};
  float T_margin = 0.08;
  float B_margin = 0.12; 
  float L_margin = 0.14;
  float R_margin = 0.04;


bool HasRatioPlot = false;
float epsi = 0.000000000001;
float Lumi = .432553397338;
int num_bins = 100;
float min_bin = 0;
float max_bin = 40;
float hist_max = 5000;
string custom_bin_name_X;
int num_bins_X = 100;
float min_bin_X = 0;
float max_bin_X = 40;
float hist_max_X = 5000;
string custom_bin_name_Y;
int num_bins_Y = 100;
float min_bin_Y = 0;
float max_bin_Y = 40;
float hist_max_Y = 5000;
TString x_title = "m_{lep^{+}lep^{-}} (GeV)";
TString y_title = "Events / 0.4 GeV";
TString hist_name = "pair_mass";
TString outName = "pair_mass_hist";
string DoLogScale = "no";
bool DoLegend = false;
int label_pos = 0;
float x_legend_max = .92;
float x_legend_width = .30;
float y_legend_max = .3;
float y_legend_width = .18;
vector<string> RJR_Methods = {};


map< string, map<string, tuple<bool, TString, int, int, float>>> histParams;
 map<string,map<string,TString>> FancyNames;


// choose colors with https://root.cern/doc/master/classTColor.html
// choose markers with https://root.cern.ch/doc/master/classTAttMarker.html

string CleanString(string dirty_string){
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\n'), dirty_string.cend()); 
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\r'), dirty_string.cend());
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), ' '), dirty_string.cend());
  return dirty_string;
}

string convertToString(char* a)
{
    string s(a);
 
    // we cannot use this technique again
    // to store something in s
    // because we use constructors
    // which are only called
    // when the string is declared.
 
    // Remove commented portion
    // to see for yourself
 
    /*
    char demo[] = "gfg";
    s(demo); // compilation error
    */
 
    return s;
}

void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  // pad1->Draw();  
  pad1->cd();  pad1->SetLeftMargin(L_margin);  pad1->SetBottomMargin(B_margin); pad1->SetRightMargin(R_margin); pad1->SetTopMargin(T_margin);
  //pad1->SetGrid(); 
  //pad1->SetGrid(10,10); 
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


TString replace_first(
    TString char_array,
    std::string const& toReplace,
    std::string const& replaceWith
) {
    std::string s(char_array);
    cout << s << endl;
    std::size_t pos = s.find(toReplace);
    if (pos == std::string::npos) return char_array;
    s.replace(pos, toReplace.length(), replaceWith);
    cout << s <<endl;
    return s.c_str();
}


// void fillMaps_CSV(string RJR_Method,string selection){
   // // string Hists_Folder = save_folder;
   
// string csvname = save_folder + "DataSamplesInfo.csv";
// TTree *SampleInfoTree = new TTree("SampleInfoTree", "tree to sample info");
// SampleInfoTree->ReadFile(csvname.c_str(),"SampleName/C:SampleIsMC/I:SampleFileName/C:FillColor/I:MarkerType/I:CrossSection/D:NumofEvents/D:FudgeFactor/D:FancyName/C");
// Char_t SampleName[1000]; 
// int SampleIsMC;
// Char_t SampleFileName[1000]; 
// int FillColor;
// int MarkerType;
// double CrossSection;
// double NumofEvents;
// double FudgeFactor;
// Char_t FancyName[1000]; 
// SampleInfoTree->SetBranchAddress("SampleName", &SampleName);
// SampleInfoTree->SetBranchAddress("SampleIsMC", &SampleIsMC);
// SampleInfoTree->SetBranchAddress("SampleFileName", &SampleFileName);
// SampleInfoTree->SetBranchAddress("FillColor", &FillColor);
// SampleInfoTree->SetBranchAddress("MarkerType", &MarkerType);
// SampleInfoTree->SetBranchAddress("CrossSection", &CrossSection);
// SampleInfoTree->SetBranchAddress("NumofEvents", &NumofEvents);
// SampleInfoTree->SetBranchAddress("FudgeFactor", &FudgeFactor);
// SampleInfoTree->SetBranchAddress("FancyName", &FancyName);

// for(int i = 0; i < SampleInfoTree->GetEntries(); i++){
    // SampleInfoTree->GetEntry(i);
    // // string dataname = convertToString(SampleName);
    // // cout << " SampleName " << SampleName << endl;
    // // cout << " SampleIsMC " << (SampleIsMC==1) << endl;
    // // cout << " SampleFileName " << SampleFileName << endl;
    // // cout << " FillColor " << FillColor << endl;
    // // cout << " MarkerType " << MarkerType << endl;
    // // cout << " CrossSection " << CrossSection << endl;
    // // cout << " NumofEvents " << NumofEvents << endl;
    // // cout << " FudgeFactor " << FudgeFactor << endl;
    // // cout << " FancyName " << FancyName << endl;
    // // cout << save_folder +  RJR_Method +"/" +selection+ SampleFileName << endl;
// // data type, isMC ,file, color,  marker type, scale	(Lumi*Cross section nb/ Number of Events Generated)   
// histParams[convertToString(SampleName)] = {(SampleIsMC==1) ,save_folder +  RJR_Method +"/" +selection+ SampleFileName , FillColor ,MarkerType, Lumi*FudgeFactor*CrossSection/NumofEvents};

// // cout << dataname << "  :  " << histParams[dataname] << endl;
// // data type, Fancy Name (latex version)

// FancyNames[convertToString(SampleName)] = FancyName;

// }
// }

void fillMaps_CSV_withSelections(string CSVFileName, string RJR_Method){
   // string Hists_Folder = save_folder;
   
string csvname = save_folder + CSVFileName;
TTree *SampleInfoTree = new TTree("SampleInfoTree", "tree to sample info");
SampleInfoTree->ReadFile(csvname.c_str(),"SampleName/C:SelectionName/C:SampleIsMC/I:SampleFileName/C:FillColor/I:MarkerType/I:CrossSection/D:NumofEvents/D:FudgeFactor/D:FancyName/C");
Char_t SampleName[1000]; 
Char_t SelectionName[1000]; 
int SampleIsMC;
Char_t SampleFileName[1000]; 
int FillColor;
int MarkerType;
double CrossSection;
double NumofEvents;
double FudgeFactor;
Char_t FancyName[1000]; 
SampleInfoTree->SetBranchAddress("SampleName", &SampleName);
SampleInfoTree->SetBranchAddress("SelectionName", &SelectionName);
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
if (SampleIsMC==1) histParams[convertToString(SampleName)][convertToString(SelectionName)] = {(SampleIsMC==1) ,save_folder +  RJR_Method +"/" +SelectionName+ SampleFileName , FillColor ,MarkerType, Lumi*FudgeFactor*CrossSection/NumofEvents};
if (SampleIsMC!=1) histParams[convertToString(SampleName)][convertToString(SelectionName)] = {(SampleIsMC==1) ,save_folder +  RJR_Method +"/" +SelectionName+ SampleFileName , FillColor ,MarkerType, FudgeFactor*CrossSection/NumofEvents};

// cout << dataname << "  :  " << histParams[dataname] << endl;
// data type, Fancy Name (latex version)

FancyNames[convertToString(SampleName)][convertToString(SelectionName)] = FancyName;

}
}



TCanvas* example_plot( int iPeriod, int iPos );
 
void myMacro_PlotTProfile_CustomBins(std::string config_file_name, std::string list_select = "list_selections.txt", std::string list_data = "list_datatypes.txt"){
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
    if(config_inputs.size() < 12){ cout << " Bad config file " << endl; return;}

   prong_type = config_inputs.at(2);
   save_folder = config_inputs.at(1)+ prong_type+ "/";
   string infile_name_beginning = config_inputs.at(0)+ prong_type+ "/";
   string paraTree_filename = save_folder + config_inputs.at(3);
   string paraTree_filename_2D = save_folder + config_inputs.at(4);
   string cuts_types_file = save_folder + config_inputs.at(5);
   string selection_file = save_folder + config_inputs.at(6);
   string outName_start = config_inputs.at(7);
   // DENOM = config_inputs.at(8); // denominator of efficiency calculation. selections is nominator
   // SF_NOM_datatype = config_inputs.at(9); // numerator of scale type. usually data
   normalization = config_inputs.at(10);
   string PlotFolderName = config_inputs.at(11);
   string CSVFileName = config_inputs.at(12);
   
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
ParameterTree_2D->ReadFile(paraTree_filename_2D.c_str(),"KeepHists/I:BranchName_X/C:BranchName_Y/C:HistName/C:HistTitle/C:CustomBinsName_X/C:nBins_X/I:min_X/D:max_X/D:CustomBinsName_Y/C:nBins_Y/I:min_Y/D:max_Y/D:xaxis/C:yaxis/C:label_pos/I:x_leg_max/D:x_leg_width/D:y_leg_max/D:y_leg_width/D");
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
int label_position; 
double x_leg_max; 
double x_leg_width;
double y_leg_max; 
double y_leg_width;  

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
ParameterTree_2D->SetBranchAddress("label_pos", &label_position); 
ParameterTree_2D->SetBranchAddress("x_leg_max", &x_leg_max); 
ParameterTree_2D->SetBranchAddress("x_leg_width", &x_leg_width); 
ParameterTree_2D->SetBranchAddress("y_leg_max", &y_leg_max); 
ParameterTree_2D->SetBranchAddress("y_leg_width", &y_leg_width);  

 
// TTree *ParameterTree = new TTree("ParameterTree", "tree to fill hist parameters");
// ParameterTree->ReadFile(paraTree_filename.c_str(),"KeepHists/I:BranchName/C:HistName/C:HistTitle/C:CustomBinsName/C:nBins/I:min/D:max/D:xaxis/C:yaxis/C:label_pos/I:x_leg_max/D:x_leg_width/D:y_leg_max/D:y_leg_width/D");
// int KeepHists;
// Char_t HistName[1000]; 
// Char_t HistTitle[1000]; 
// Char_t BranchName[1000]; 
// Char_t CustomBinsName[1000]; 
// int nBins;
// double bin_min;
// double bin_max;
// Char_t xaxis[1000]; 
// Char_t yaxis[1000];
// int label_position; 
// double x_leg_max; 
// double x_leg_width;
// double y_leg_max; 
// double y_leg_width;  
// ParameterTree->SetBranchAddress("KeepHists", &KeepHists);
// ParameterTree->SetBranchAddress("BranchName", &BranchName);
// ParameterTree->SetBranchAddress("HistName", &HistName);
// ParameterTree->SetBranchAddress("HistTitle", &HistTitle);
// ParameterTree->SetBranchAddress("CustomBinsName", &CustomBinsName);
// ParameterTree->SetBranchAddress("nBins", &nBins);
// ParameterTree->SetBranchAddress("min", &bin_min);
// ParameterTree->SetBranchAddress("max", &bin_max); 
// ParameterTree->SetBranchAddress("xaxis", &xaxis);
// ParameterTree->SetBranchAddress("yaxis", &yaxis); 
// ParameterTree->SetBranchAddress("label_pos", &label_position); 
// ParameterTree->SetBranchAddress("x_leg_max", &x_leg_max); 
// ParameterTree->SetBranchAddress("x_leg_width", &x_leg_width); 
// ParameterTree->SetBranchAddress("y_leg_max", &y_leg_max); 
// ParameterTree->SetBranchAddress("y_leg_width", &y_leg_width); 


// (lep , #tau _{3prong})
  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  
for( string RJR_Method  :  RJR_Methods ){ 

// fillMaps(prong_type, RJR_Method, selection);
fillMaps_CSV_withSelections(CSVFileName,RJR_Method);

  //  gROOT->LoadMacro("CMS_lumi.C");
std::stringstream stream;
stream << std::fixed << std::setprecision(2) << Lumi;
string cms_text = "PbPb #sqrt{s_{NN}} = 5.02 TeV ( " + stream.str() +  " nb^{-1} )";
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
  
   string out_folder = save_folder +  RJR_Method+"/" +PlotFolderName+ "/" + "ProfilePlot";
   
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
    outName = out_folder +"/" + outName_start + hist_name;
    DoLogScale = "no";
    DoLegend = false;
    label_pos = label_position;
    x_legend_max = x_leg_max;
    x_legend_width = x_leg_width;
    y_legend_max = y_leg_max;
    y_legend_width = y_leg_width;

    custom_bin_name_X = CustomBinsName_X_2D;
    custom_bin_name_Y = CustomBinsName_Y_2D;
    
    
    float unity_max = 0;
    // TH1F *tmp_hist;
    // if(normalization == "Unity"){
    // for(string DataType : DataTypes ){
    // for(string selection : selections){
    // auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[DataType][selection];
    // TFile file_(FileName,"READ");
    // tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
       // float entries = tmp_hist->Integral() + epsi;
       // tmp_hist->Scale(1/entries);
    // if (unity_max < tmp_hist->GetMaximum()) unity_max = tmp_hist->GetMaximum();
    // file_.Close();
    // }
    // }
    // }
    // if(normalization == "XS"){
    // // auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams["Data"];
    // // TFile file_(FileName,"READ");
    // // tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
    // // unity_max = tmp_hist->GetMaximum();
    // // file_.Close();
    // }
    hist_max = unity_max*1.2;
  example_plot( iPeriod, label_pos ); 
  
        }
  //  example_plot( iPeriod, 11 );  // left-aligned
  //  example_plot( iPeriod, 33 );  // right-aligned

  //  writeExtraText = false;       // remove Preliminary
  
  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)

  //  example_plot( iPeriod, 11 );  // default: left-aligned
  //  example_plot( iPeriod, 22 );  // centered
  //  example_plot( iPeriod, 33 );  // right-aligned  
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

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
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
    TPad* pad2;
  // c1->SetFillColor(0);
  //c1->SetGrid();
  if(HasRatioPlot){
  pad1 = new TPad("pad1","The pad",0,.3,1,1);
  pad2 = new TPad("pad2","The second pad",0,.0,1,.3);
    pad1->Draw();
  pad2->Draw();
    applyPadStyle(pad2);
  }
   if(!HasRatioPlot){
  pad1 = new TPad("pad1","The pad",0,0,1,1);
    pad1->Draw();
  }
  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // gStyle->SetEndErrorSize(4);
        // pad1->SetLeftMargin(.18);
         // pad1->SetRightMargin(.1);
         // pad1->SetBottomMargin(.2);
         // pad1->SetTopMargin(0.1);
  canv->cd();
  // pad1->Draw();
  // pad2->Draw();
  pad1->cd();
 
 if(HasRatioPlot){
    pad1->SetBottomMargin(.02);
    pad2->SetTopMargin(.02);
    pad2->SetBottomMargin(.2);
 }
 
 if(!HasRatioPlot){
    pad1->SetBottomMargin(.15);
    // pad1->SetTopMargin(.12);
 }


  TH1* h;
   if(custom_bin_name_X == "NA"){
   h = new TH1F("h","h",num_bins_X,min_bin_X,max_bin_X);
   }
   if(custom_bin_name_X != "NA"){
   int size = constants::mapvectors[custom_bin_name_X].size();
   double testarray[size];
   for(int j =0; j < size; j++) testarray[j] = constants::mapvectors[custom_bin_name_X].at(j);
   h = new TH1F("h","h", size -1, testarray);
   }
  // h->GetXaxis()->SetNdivisions(6,5,0);
  h->GetXaxis()->SetTitle(x_title);  
  // h->GetYaxis()->SetNdivisions(6,5,0);
  h->GetYaxis()->SetTitleOffset(1.25);
  // TString Unity_y_title = y_title + " (Normalized to Unity)";
  // if(normalization == "Unity")  h->GetYaxis()->SetTitle(Unity_y_title);  
  h->GetYaxis()->SetTitle(y_title);  

  h->GetXaxis()->SetTitleSize(.05);   
  h->GetYaxis()->SetTitleSize(.05); 



h->GetYaxis()->SetNdivisions(6,5,0);
h->GetXaxis()->SetNdivisions(6,5,0);



  // h->SetMaximum(hist_max);
  // // if( iPos==1 ) h->SetMaximum( 10000 );
  // h->Draw();
  


  float markerSize  = 1.0;

    float x1_l = x_legend_max;
    float y1_l = y_legend_max;
    float dx_l = x_legend_width;
    float dy_l = y_legend_width;
    float x0_l = x1_l-dx_l;
    float y0_l = y1_l-dy_l;

    auto legend = new TLegend(x0_l,y0_l,x1_l, y1_l );
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
  float tmp_hist_max = -9E10;
  float tmp_hist_min = 9E10;

  if(1==1){
    // Observed data
THStack *hs_stack = new THStack("hs_stack","Stacked Histogram"); 
    for(string dataType : DataTypes){
    for(string selection: selections){
    auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[dataType][selection];
    TFile file_(FileName,"READ");
    TString hist_name_2D = "Profile_" + hist_name;
    TProfile *tmp_hist = static_cast<TProfile*>(file_.Get(hist_name_2D)->Clone());
    tmp_hist->SetDirectory(0);
    tmp_hist->SetMarkerStyle(marker_style);
    tmp_hist->SetMarkerSize(markerSize);
    tmp_hist->SetLineColor(fill_color);
    if(tmp_hist->GetMaximum() > tmp_hist_max) tmp_hist_max = tmp_hist->GetMaximum();
    if(tmp_hist->GetMinimum() < tmp_hist_min) tmp_hist_min = tmp_hist->GetMinimum();
    // if(isMC) tmp_hist->SetFillColor(fill_color);
    tmp_hist->SetMarkerColor(fill_color);
    // tmp_hist->Scale(scale);
    legend->AddEntry(tmp_hist,FancyNames[dataType][selection]);
    // tmp_hist->Draw();
    hs_stack->Add(tmp_hist);

    file_.Close();
    }
  }
    pad1->cd();
    
    // hs_stack->GetXaxis()->SetTitle(x_title);  
    // hs_stack->GetYaxis()->SetTitleOffset(1.25); 
    // hs_stack->GetYaxis()->SetTitle(y_title);  

    // hs_stack->GetXaxis()->SetTitleSize(.05);   
    // hs_stack->GetYaxis()->SetTitleSize(.05); 



    // hs_stack->GetYaxis()->SetNdivisions(6,5,0);
    // hs_stack->GetXaxis()->SetNdivisions(6,5,0);
    h->SetMaximum(tmp_hist_max + .2*abs(tmp_hist_max));
    h->SetMinimum(tmp_hist_min - .2*abs(tmp_hist_min));
    h->Draw();
    hs_stack->Draw("NOSTACK same");
  }
  
  
  
  if(DoLegend) legend->Draw();
  
  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos );

  canv->Update();
  // canv->RedrawAxis();
  // canv->GetFrame()->Draw();
  if(DoLogScale == "yes"){
  pad1->SetLogy(1);
  canv->Print(outName+"_log.pdf",".pdf");
  canv->Print(outName+"_log.png",".png");
  }
  if(DoLogScale == "no"){
  pad1->SetLogy(0);
  canv->Print(outName+".pdf",".pdf");
  canv->Print(outName+".png",".png");
  }  
   if(DoLogScale == "both"){
  pad1->SetLogy(0);
  canv->Print(outName+".pdf",".pdf");
  canv->Print(outName+".png",".png");
  pad1->SetLogy(1);
  canv->Print(outName+"_log.pdf",".pdf");
  canv->Print(outName+"_log.png",".png");
  }
  
  // canv->Close();


  return canv;
}