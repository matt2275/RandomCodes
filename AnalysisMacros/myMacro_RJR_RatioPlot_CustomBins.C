// myMacro_RJR_RatioPlot_CustomBins draws histograms which were created 
// in the CreateHistogram_CustomBins Macro 
// This only does 1D histograms and will plot all datatypes in same plot
//  have the option to add a ratio plot below the main plot as well

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

// header file with custom histogram bin information
#include "variables.h"

 using namespace std;
 
 
 // variable initialization
string save_folder = "";
string prong_type = "";

string normalization = "";

vector<string> DataTypes = {};
vector<string> selections ={};

// margins of TCanvas
  float T_margin = 0.08;
  float B_margin = 0.12; 
  float L_margin = 0.14;
  float R_margin = 0.04;

// variable to plot ratio plot below main plot
bool HasRatioPlot = false;


float epsi = 0.000000000001;
float Lumi = .432553397338;
int num_bins = 100;
float min_bin = 0;
float max_bin = 40;
float hist_max = 5000;
TString x_title = "m_{lep^{+}lep^{-}} (GeV)";
TString y_title = "Events / 0.4 GeV";
TString hist_name = "pair_mass";
TString outName = "pair_mass_hist";
string DoLogScale = "no";
bool DoLegend = true;
int label_pos = 33;
string custom_bin_name;
float x_legend_max = .92;
float x_legend_width = .30;
float y_legend_max = .60;
float y_legend_width = .18;


map< string, tuple<bool, TString, int, int, float>> histParams;
 map<string,TString> FancyNames;


// choose colors with https://root.cern/doc/master/classTColor.html
// choose markers with https://root.cern.ch/doc/master/classTAttMarker.html


// cleans string of spaces and \r and \n

string CleanString(string dirty_string){
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\n'), dirty_string.cend()); 
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), '\r'), dirty_string.cend());
  dirty_string.erase(std::remove(dirty_string.begin(), dirty_string.end(), ' '), dirty_string.cend());
  return dirty_string;
}

// converts char array to string
string convertToString(char* a)
{
    string s(a);
 
    return s;
}


// applies margins and fill color to Tpad
void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  // pad1->Draw();  
  pad1->cd();  pad1->SetLeftMargin(L_margin);  pad1->SetBottomMargin(B_margin); pad1->SetRightMargin(R_margin); pad1->SetTopMargin(T_margin);
  //pad1->SetGrid(); 
  //pad1->SetGrid(10,10); 
}

// One of ratio plot methods, this one plots the Relative difference of MC vs Data (data - mc) /( data + mc) 
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

// One of ratio plot methods, this one plots the Ratio of MC vs Data (data /mc) 
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

// used to replace #comma with , in plot titles. Do to issues with csv using , as a delimiter
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


// old method for getting histogram file names and cross section. Now use fillMaps_CSV
void fillMaps(string prong_type, string RJR_Method, string selection){
   string Hists_Folder = save_folder;
// data type, isMC ,file, color,  marker type, scale	(Lumi*Cross section nb/ Number of Events Generated)   
histParams["Data"] = {false ,Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_Data.root" , 1 ,20, 1};
histParams["DiTau_SuperChic"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiTau_SuperChic.root" ,860 ,21, Lumi*8.83E6/9.97/1184600};
histParams["DiTau_Madgraph"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiTau_Madgraph.root" ,865 ,24, Lumi*570.0E3/950000};
histParams["BBbar"] = {true , Hists_Folder +  RJR_Method +"/" +selection+"/MinCuts_BBbar.root" ,397 ,66, Lumi*1.5E3/20000000*1};
histParams["CCbar"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_CCbar.root" ,419 ,44, Lumi*300.0E3/19976000*1};
histParams["DiElectron_SuperChic"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiElectron_SuperChic.root" ,910 ,33, Lumi*8.83E6/67810000/2};
histParams["DiElectron_Starlight"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiElectron_Starlight.root" ,614 ,32, Lumi*7.92E6/66750000/2};
histParams["DiMuon_Gamma"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Gamma.root" ,890 ,21, Lumi*140E3/2390000};
histParams["DiMuon_noFSR"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_noFSR.root" ,880 ,22, Lumi*189.08E3/30000000*17.97 };
histParams["DiMuon_Arash"] = {true , Hists_Folder  +  RJR_Method +"/" +selection+"/MinCuts_DiMuon_Arash.root" ,871 ,23, Lumi*8.83E6/250000/2.5 };
 // DiMuon_FSR 17.97 comes from only taking 1M of 18 M events passed primary selections
 // DiTau_SuperChic X-sec comes from DiElectron_SuperChic / (2*M_tau / 2( mass cut ee))^4 , divided by 2 for fudge factor
 // DiElectron_Starlight X-sec Divided by 9.11 for Single EG3 Prescale, times 3 for fudge factor
 // DiElectron_SuperChic X-sec Divided by 9.11 for Single EG3 Prescale,times 3 for fudge factor
 //DiTau_SuperChic X-sec comes from DiElectron_SuperChic / (2.5 (mass cut MuMu) / 2( mass cut ee))^4 , divided by 2 for fudge factor
 
 
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

// declaring example_plot, its defined later
TCanvas* example_plot( int iPeriod, int iPos );
 
 
void myMacro_RJR_RatioPlot_CustomBins(std::string config_file_name, std::string list_select = "list_selections.txt", std::string list_data = "list_datatypes.txt"){
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
    if(config_inputs.size() < 11){ cout << " Bad config file " << endl; return;}

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
   
   string RJR_Method = "HasVisMass";
   
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
   
 
// parameter tree for 1D histograms 
TTree *ParameterTree = new TTree("ParameterTree", "tree to fill hist parameters");
ParameterTree->ReadFile(paraTree_filename.c_str(),"KeepHists/I:BranchName/C:HistName/C:HistTitle/C:CustomBinsName/C:nBins/I:min/D:max/D:xaxis/C:yaxis/C:label_pos/I:x_leg_max/D:x_leg_width/D:y_leg_max/D:y_leg_width/D");
int KeepHists;
Char_t HistName[1000]; 
Char_t HistTitle[1000]; 
Char_t BranchName[1000]; 
Char_t CustomBinsName[1000]; 
int nBins;
double bin_min;
double bin_max;
Char_t xaxis[1000]; 
Char_t yaxis[1000];
int label_position; 
double x_leg_max; 
double x_leg_width;
double y_leg_max; 
double y_leg_width;  
ParameterTree->SetBranchAddress("KeepHists", &KeepHists);
ParameterTree->SetBranchAddress("BranchName", &BranchName);
ParameterTree->SetBranchAddress("HistName", &HistName);
ParameterTree->SetBranchAddress("HistTitle", &HistTitle);
ParameterTree->SetBranchAddress("CustomBinsName", &CustomBinsName);
ParameterTree->SetBranchAddress("nBins", &nBins);
ParameterTree->SetBranchAddress("min", &bin_min);
ParameterTree->SetBranchAddress("max", &bin_max); 
ParameterTree->SetBranchAddress("xaxis", &xaxis);
ParameterTree->SetBranchAddress("yaxis", &yaxis); 
ParameterTree->SetBranchAddress("label_pos", &label_position); 
ParameterTree->SetBranchAddress("x_leg_max", &x_leg_max); 
ParameterTree->SetBranchAddress("x_leg_width", &x_leg_width); 
ParameterTree->SetBranchAddress("y_leg_max", &y_leg_max); 
ParameterTree->SetBranchAddress("y_leg_width", &y_leg_width); 
// (lep , #tau _{3prong})
  //  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  
for( string selection : selections){

// fillMaps(prong_type, RJR_Method, selection);
fillMaps_CSV(RJR_Method, selection);

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
  
   string out_folder = save_folder +  RJR_Method+"/" +selection+ "/" + normalization;
   
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
  
for(int i = 0; i < ParameterTree->GetEntries(); i++){
    ParameterTree->GetEntry(i);
   // xaxis = replace_first(xaxis,"#comma", ",");
   // yaxis = replace_first(yaxis,"#comma", ",");
    num_bins = nBins;
    min_bin = bin_min;
    max_bin = bin_max;
    x_title = replace_first(xaxis,"#comma", ",");
    y_title = replace_first(yaxis,"#comma", ",");
    hist_name = HistName;
    outName = out_folder +"/" + outName_start + hist_name;
    DoLogScale = "both";
    DoLegend = true;
    label_pos = label_position;
    x_legend_max = x_leg_max;
    x_legend_width = x_leg_width;
    y_legend_max = y_leg_max;
    y_legend_width = y_leg_width;

    custom_bin_name = CustomBinsName;
    
    
    
    // finding max and min for plots 
    float unity_max = 0;
    TH1F *tmp_hist;
    if(normalization == "Unity"){
    for(string DataType : DataTypes ){
    auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[DataType];
    TFile file_(FileName,"READ");
    tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
       float entries = tmp_hist->Integral() + epsi;
       tmp_hist->Scale(1/entries);
    if (unity_max < tmp_hist->GetMaximum()) unity_max = tmp_hist->GetMaximum();
    file_.Close();
    }
    }
    if(normalization == "XS"){
    auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams["Data"];
    TFile file_(FileName,"READ");
    tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
    unity_max = tmp_hist->GetMaximum();
    file_.Close();
    }
    hist_max = unity_max*1.2;
    
  // function for creating plots
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

//creates canvas and draws histograms in CMS style
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

 // creates blank histogram with same bins as other histograms
  TH1* h;
   if(custom_bin_name == "NA"){
   h = new TH1F("h","h",num_bins,min_bin,max_bin);
   }
   if(custom_bin_name != "NA"){
   int size = constants::mapvectors[custom_bin_name].size();
   double testarray[size];
   for(int j =0; j < size; j++) testarray[j] = constants::mapvectors[custom_bin_name].at(j);
   h = new TH1F("h","h", size -1, testarray);
   }
  // h->GetXaxis()->SetNdivisions(6,5,0);
  h->GetXaxis()->SetTitle(x_title);  
  // h->GetYaxis()->SetNdivisions(6,5,0);
  h->GetYaxis()->SetTitleOffset(1.25);
  TString Unity_y_title = y_title + " (Normalized to Unity)";
  if(normalization == "Unity")  h->GetYaxis()->SetTitle(Unity_y_title);  
  if(normalization == "XS")  h->GetYaxis()->SetTitle(y_title);  

  h->GetXaxis()->SetTitleSize(.05);   
  h->GetYaxis()->SetTitleSize(.05); 



h->GetYaxis()->SetNdivisions(6,5,0);
h->GetXaxis()->SetNdivisions(6,5,0);



  h->SetMaximum(hist_max);
  // if( iPos==1 ) h->SetMaximum( 10000 );
  h->Draw();
  


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



// plot data and MC according to their given cross section, MC is stacked.
  if(normalization == "XS" ){
    // Observed data
THStack *mc_stack = new THStack("mc_stack","Stacked Histogram"); 
TH1F *data ;
    for(string dataType : DataTypes){
    auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[dataType];
    TFile file_(FileName,"READ");
    TH1F *tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
    tmp_hist->SetDirectory(0);
    tmp_hist->SetMarkerStyle(marker_style);
    tmp_hist->SetMarkerSize(markerSize);
    tmp_hist->SetLineColor(fill_color);
    if(isMC) tmp_hist->SetFillColor(fill_color);
    tmp_hist->SetMarkerColor(fill_color);
    tmp_hist->Scale(scale);
    legend->AddEntry(tmp_hist,FancyNames[dataType]);
    // tmp_hist->Draw();
    if(isMC)  mc_stack->Add(tmp_hist);
    if(!isMC)  data = tmp_hist;

    file_.Close();
    }
    if(HasRatioPlot){
       pad2->cd();
    TList *stackHists = mc_stack->GetHists();
 
    TH1* tmpHist = (TH1*)stackHists->At(0)->Clone();
    tmpHist->Reset();
    for (int i=0;i<stackHists->GetSize();++i) {
       tmpHist->Add((TH1*)stackHists->At(i));
    }
    auto ratio_graph = Hists_Diff( data, tmpHist);
   float nbins = data->GetNbinsX();
   float RG_max = TMath::MaxElement(nbins, ratio_graph->GetY());
   float RG_min = TMath::MinElement(nbins, ratio_graph->GetY());
	h->GetXaxis()->SetTitle("");
	h->GetXaxis()->SetLabelSize(0);
	ratio_graph->GetYaxis()->SetTitleSize(.1);
	ratio_graph->GetYaxis()->SetTitle("#frac{Data - MC}{Data + MC}");	
	// ratio_graph->GetYaxis()->SetTitle("#frac{Data}{MC}");	
//   ratio_graph->GetYaxis()->SetRangeUser(RG_min -fabs(RG_min)*.2, RG_max + fabs(RG_max)*.2);
   ratio_graph->GetYaxis()->SetRangeUser(-1.1, 1.1);
   ratio_graph->GetXaxis()->SetLimits(min_bin, max_bin);
	ratio_graph->GetXaxis()->SetTitleSize(.1);
	ratio_graph->GetXaxis()->SetTitle(x_title);
	ratio_graph->GetXaxis()->SetLabelSize(.1);	
	ratio_graph->GetXaxis()->SetLabelSize(.1);	
	ratio_graph->GetYaxis()->SetNdivisions(6,5,0);
	ratio_graph->GetXaxis()->SetNdivisions(6,5,0);
	ratio_graph->SetTitle("");
	ratio_graph->Draw("AP");
	// pad1->cd();
    }
    pad1->cd();
    mc_stack->Draw("histsame");
    data->Draw("esamex0");

  }
  
  // drawing plots where all datatypes are normalized to 1, there is no stacking of MC
  if(normalization == "Unity" ){
         pad1->cd();
    // Observed data
    for(string dataType : DataTypes){
    auto &[isMC, FileName, fill_color, marker_style, scale]  = histParams[dataType];
    TFile file_(FileName,"READ");
    TH1F *tmp_hist = static_cast<TH1F*>(file_.Get(hist_name)->Clone());
    float entries = tmp_hist->Integral() + epsi;
    tmp_hist->SetDirectory(0);
    tmp_hist->SetMarkerStyle(marker_style);
    tmp_hist->SetMarkerSize(markerSize);
    tmp_hist->SetLineColor(fill_color);
    // tmp_hist->SetFillColor(fill_color);
    tmp_hist->SetMarkerColor(fill_color);
    tmp_hist->Scale(1/entries);
    // tmp_hist->Draw();
    legend->AddEntry(tmp_hist,FancyNames[dataType]);
    tmp_hist->Draw("same");

    file_.Close();
    }
  }
  
  
  if(DoLegend) legend->Draw();


// printing plots in linear and log scale
  
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