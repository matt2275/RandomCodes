#!/bin/bash
# . /cvmfs/sft.cern.ch/lcg/releases/LCG_96/ROOT/6.18.00/x86_64-centos7-gcc8-opt/bin/thisroot.sh
for i in 1 2 3
# for i in 1
do
  root -x -q -b "myMacro_RJR_RatioPlot_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/Updated_Selections/1Prong/config.txt\", \"list_selections_$i.txt\", \"list_datatypes_plotting.txt\" )"
  
  # root -x -q -b "myMacro_RJR_RatioPlot_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/PionCheck/1Prong/config_plotting.txt\", \"SelectionFiles/list_selections_MuEl.txt\", \"SelectionFiles/list_datatypes_plotting.txt\" )"
  
  # root -x -q -b "myMacro_RJR_RatioPlot_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/GenTauCheck/1Prong/config.txt\", \"list_selections_$i.txt\", \"list_datatypes.txt\" )"
  
    
  # root -x -q -b "myMacro_RJR_RatioPlot_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/Efficiency_WithIDs/1Prong/Efficiency_Configs/config_Electron_EG5_$i.txt\", \"Efficiency_Configs/list_selections_Electron.txt\", \"Efficiency_Configs/list_datatypes_Electron.txt\" )"

  # root -x -q -b "myMacro_RJR_RatioPlot_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/Efficiency_WithIDs/1Prong/Efficiency_Configs/config_Muon_$i.txt\", \"Efficiency_Configs/list_selections_Muon.txt\", \"Efficiency_Configs/list_datatypes_Muon.txt\" )"
  
done