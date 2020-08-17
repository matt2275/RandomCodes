import numpy as np
import ROOT
import numpy as np
from array import array


outfile = ROOT.TFile( 'TauTauDiLeptonMass.root', 'recreate' )
t = ROOT.TTree( 'tautree', 'tree with simulated tau decays' )
 
maxn = 10
DiLeptonMass = array( 'f', [ 0 ] )
DiLeptonRapidity = array( 'f', [ 0 ] )
DiLeptonPt = array( 'f', [ 0 ] )
DiLeptonAcoplanarity = array( 'f', [ 0 ] )
DiLeptonDeltaPhi = array( 'f', [ 0 ] )
DiLeptonCharge = array( 'f', [ 0 ] )
DiLeptonOpeningAngle = array( 'f', [ 0 ] )

t.Branch( 'DiLeptonMass', DiLeptonMass, 'DiLeptonMass/F' )
t.Branch( 'DiLeptonRapidity', DiLeptonRapidity, 'DiLeptonRapidity/F' )
t.Branch( 'DiLeptonPt', DiLeptonPt, 'DiLeptonPt/F' )
t.Branch( 'DiLeptonAcoplanarity', DiLeptonAcoplanarity, 'DiLeptonAcoplanarity/F' )
t.Branch( 'DiLeptonDeltaPhi', DiLeptonDeltaPhi, 'DiLeptonDeltaPhi/F' )
t.Branch( 'DiLeptonCharge', DiLeptonCharge, 'DiLeptonCharge/F' )
t.Branch( 'DiLeptonOpeningAngle', DiLeptonOpeningAngle, 'DiLeptonOpeningAngle/F' )

Tau_mass = 1.77686
Muon_mass = .1056583755
El_mass = .511E-3 

	
def return_aco(phi1, phi2):
	phi1 = phi1%(2*np.pi)
	phi2 = phi2%(2*np.pi)
	deltaphi = abs(phi1- phi2)
	if deltaphi > np.pi :
		deltaphi = 2 * np.pi - deltaphi
	aco = 1 - deltaphi/np.pi
	return(aco)
	

f = ROOT.TFile.Open("/eos/user/m/mnickel/TauTau/TauTau_MC.root ")



for event in f.EventTree :
	aVec = ROOT.TLorentzVector()
	bVec = ROOT.TLorentzVector()
	nEle = event.nEle
	nMu = event.nMu
	acharge = 0
	bcharge = 0
	aphi = 0
	bphi = 0
	saveTree = False
	if(nMu == 2 and nEle == 0):
		aPT = event.muPt[0]
		bPT = event.muPt[1]
		aeta = event.muEta[0]
		beta = event.muEta[1]
		aphi = event.muPhi[0]
		bphi = event.muPhi[1]
		acharge = event.muCharge[0]
		bcharge = event.muCharge[1]
		aVec.SetPtEtaPhiM(aPT, aeta, aphi, Muon_mass)
		bVec.SetPtEtaPhiM(bPT, beta, bphi, Muon_mass)
		saveTree = True
	if(nMu == 0 and nEle == 2):
		aPT = event.elePt[0]
		bPT = event.elePt[1]
		aeta = event.eleEta[0]
		beta = event.eleEta[1]
		aphi = event.elePhi[0]
		bphi = event.elePhi[1]
		acharge = event.eleCharge[0]
		bcharge = event.eleCharge[1]
		aVec.SetPtEtaPhiM(aPT, aeta, aphi, El_mass)
		bVec.SetPtEtaPhiM(bPT, beta, bphi, El_mass)
		saveTree = True
	if(nMu == 1 and nEle == 1):
		aPT = event.muPt[0]
		bPT = event.elePt[0]
		aeta = event.muEta[0]
		beta = event.eleEta[0]
		aphi = event.muPhi[0]
		bphi = event.elePhi[0]
		acharge = event.muCharge[0]
		bcharge = event.eleCharge[0]
		aVec.SetPtEtaPhiM(aPT, aeta, aphi, Muon_mass)
		bVec.SetPtEtaPhiM(bPT, beta, bphi, El_mass)		
		saveTree = True
	if(saveTree):
		DiLeptonDeltaPhi = abs(aphi - bphi)
		DiLeptonCharge[0] = abs(acharge + bcharge)
		cVec = aVec + bVec
		DiLeptonMass[0] = cVec.M()
		DiLeptonRapidity[0] = cVec.Rapidity()
		DiLeptonPt[0] = cVec.Pt()
		DiLeptonAcoplanarity[0] = return_aco(aphi, bphi)
		DiLeptonOpeningAngle[0] = aVec.Angle(bVec.Vect())
		t.Fill()
	
outfile.Write()
outfile.Close()

	
	