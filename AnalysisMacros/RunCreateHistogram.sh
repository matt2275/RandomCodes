#!/bin/bash
. /cvmfs/sft.cern.ch/lcg/releases/LCG_96/ROOT/6.18.00/x86_64-centos7-gcc8-opt/bin/thisroot.sh
for i in 1 2 3 4
# for i in 1
do
    for j in 1 2 3 
    do 
        # root -x -q -b "CreateHistograms_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/Efficiency_WithIDs/1Prong/config.txt\", \"list_selections_$j.txt\", \"list_datatypes_$i.txt\" )"
        
        # root -x -q -b "CreateHistograms_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/Updated_Selections/1Prong/config.txt\", \"list_selections_$j.txt\", \"list_datatypes_$i.txt\" )"
        
        root -x -q -b "CreateHistograms_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/TestSumDR/1Prong/config.txt\", \"list_selections_$j.txt\", \"list_datatypes_$i.txt\" )"
        
        # root -x -q -b "CreateHistograms_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/UpsilonJPsiStudy/1Prong/config_DiMuon.txt\", \"SelectionsJPsi/list_selections_DiMuon_$j.txt\", \"SelectionsJPsi/list_datatypes_$i.txt\" )"
        
        # root -x -q -b "CreateHistograms_CustomBins.C(\"/eos/user/m/mnickel/MultiTagData/Outputs/GenTauCheck/1Prong/config.txt\", \"list_selections_$j.txt\", \"list_datatypes_NEW_$i.txt\" )"
    done
done
