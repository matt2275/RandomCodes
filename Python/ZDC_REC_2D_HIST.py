import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import ROOT
import numpy as np
from array import array



f = ROOT.TFile.Open("/eos/user/m/mnickel/TauTau/TEST_ZDC.root")

h1f = ROOT.TH2F( 'h1f', 'ZDCNEnergy vs NEMtot', 2000, 0, 10000,2000, 0, 1000)
h2f = ROOT.TH2F( 'h2f', 'ZDCPEnergy vs PEMtot', 2000, 0, 10000, 2000, 0, 4000)
h3f = ROOT.TH2F( 'h3f', 'ZDCNEnergy vs NHADtot', 1000, 0, 10000,1000, 0, 10000)
h4f = ROOT.TH2F( 'h4f', 'ZDCPEnergy vs PHADtot', 1000, 0, 10000, 1000, 0, 10000)


for event in f.outTree :
	NEMtot = event.NEMtot
	NHADtot = event.NHADtot
	PEMtot = event.PEMtot
	PHADtot = event.PHADtot
	ZDCNEnergy = event.ZDCNEnergy
	ZDCPEnergy = event.ZDCPEnergy
	h1f.Fill(ZDCNEnergy,NEMtot)
	h2f.Fill(ZDCPEnergy,PEMtot)
	h3f.Fill(ZDCNEnergy,NHADtot)
	h4f.Fill(ZDCPEnergy,PHADtot)	
	
outfile = ROOT.TFile( 'test_zdc.root', 'recreate' )

h1f.Write()	
h2f.Write()	
h3f.Write()	
h4f.Write()		
outfile.Close()

	
	