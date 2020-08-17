import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
import ROOT
import numpy as np
from array import array


outfile = ROOT.TFile( 'test.root', 'recreate' )
t = ROOT.TTree( 'tautree', 'tree with simulated tau decays' )
 
maxn = 10
DiLeptonMass_OG = array( 'f', [ 0 ] )
DiLeptonMass_NEW = array( 'f', [ 0 ] )
DiLeptonRapidity_OG = array( 'f', [ 0 ] )
DiLeptonRapidity_NEW = array( 'f', [ 0 ] )
DiLeptonPt_OG = array( 'f', [ 0 ] )
DiLeptonPt_NEW = array( 'f', [ 0 ] )
DiLeptonAcoplanarity_OG = array( 'f', [ 0 ] )
DiLeptonAcoplanarity_NEW = array( 'f', [ 0 ] )

MuonPt = array( 'f', 2*[ 0 ] )
ElectronPt = array( 'f', 2*[ 0 ] )
MuonEta = array( 'f', 2*[ 0 ] )
ElectronEta = array( 'f', 2*[ 0 ] )
MuonPhi = array( 'f', 2*[ 0 ] )
ElectronPhi = array( 'f', 2*[ 0 ] )

t.Branch( 'DiLeptonMass_OG', DiLeptonMass_OG, 'DiLeptonMass_OG/F' )
t.Branch( 'DiLeptonMass_NEW', DiLeptonMass_NEW, 'DiLeptonMass_NEW/F' )
t.Branch( 'DiLeptonRapidity_OG', DiLeptonRapidity_OG, 'DiLeptonRapidity_OG/F' )
t.Branch( 'DiLeptonRapidity_NEW', DiLeptonRapidity_NEW, 'DiLeptonRapidity_NEW/F' )
t.Branch( 'DiLeptonPt_OG', DiLeptonPt_OG, 'DiLeptonPt_OG/F' )
t.Branch( 'DiLeptonPt_NEW', DiLeptonPt_NEW, 'DiLeptonPt_NEW/F' )
t.Branch( 'DiLeptonAcoplanarity_OG', DiLeptonAcoplanarity_OG, 'DiLeptonAcoplanarity_OG/F' )
t.Branch( 'DiLeptonAcoplanarity_NEW', DiLeptonAcoplanarity_NEW, 'DiLeptonAcoplanarity_NEW/F' )

t.Branch( 'LeptonPt_OG', LeptonPt_OG, 'LeptonPt_OG[2]/F' )
t.Branch( 'LeptonPt_NEW', LeptonPt_NEW, 'LeptonPt_NEW[2]/F' )
t.Branch( 'LeptonEta_OG', LeptonEta_OG, 'LeptonEta_OG[2]/F' )
t.Branch( 'LeptonEta_NEW', LeptonEta_NEW, 'LeptonEta_NEW[2]/F' )
t.Branch( 'LeptonPhi_OG', LeptonPhi_OG, 'LeptonPhi_OG[2]/F' )
t.Branch( 'LeptonPhi_NEW', LeptonPhi_NEW, 'LeptonPhi_NEW[2]/F' )



f = ROOT.TFile.Open("/eos/user/m/mnickel/TauTau/merge_DATA_PLOTS_mumu_tree-2.root")



for event in f.outTree :
	aVec = ROOT.TLorentzVector()
	bVec = ROOT.TLorentzVector()
	aPT = event.lepton_pt[0]
	bPT = event.lepton_pt[1]
	aeta = event.lepton_eta[0]
	beta = event.lepton_eta[1]
	aphi = event.lepton_phi[0]
	bphi = event.lepton_phi[1]
	aVec.SetPtEtaPhiM(aPT, aeta, aphi, Tau_mass)
	bVec.SetPtEtaPhiM(bPT, beta, bphi, Tau_mass)
	aVec = aVec + return_4Vec()
	bVec = bVec + return_4Vec()
	cVec = aVec + bVec
	DiLeptonMass_OG[0] = event.dilepton_mass
	DiLeptonRapidity_OG[0] = event.dilepton_rapidity
	DiLeptonPt_OG[0] = event.dilepton_pt
	DiLeptonAcoplanarity_OG[0] = event.dilepton_acoplanarity
	DiLeptonMass_NEW[0] = cVec.M()
	DiLeptonRapidity_NEW[0] = cVec.Rapidity()
	DiLeptonPt_NEW[0] = cVec.Pt()
	DiLeptonAcoplanarity_NEW[0] = return_aco(aVec.Phi(), bVec.Phi())
	LeptonPt_OG[0] = aPT
	LeptonPt_OG[1] = bPT
	LeptonEta_OG[0] = aeta
	LeptonEta_OG[1] = beta
	LeptonPhi_OG[0] = aphi
	LeptonPhi_OG[1] = bphi
	LeptonPt_NEW[0] = aVec.Pt()
	LeptonPt_NEW[1] = bVec.Pt()
	LeptonEta_NEW[0] = aVec.Eta()
	LeptonEta_NEW[1] = bVec.Eta()
	LeptonPhi_NEW[0] = aVec.Phi()
	LeptonPhi_NEW[1] = bVec.Phi()
	t.Fill()
	
outfile.Write()
outfile.Close()

	
	