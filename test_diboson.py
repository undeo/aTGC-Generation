import ROOT
from DataFormats.FWLite import Events, Handle
import math
import time
import numpy as np

'''
Quarks		Leptons		Gauge and
				Higgs Bosons
d 1		e- 11		g (9) 21
u 2		nue 12		gamma 22
s 3		mu- 13		Z0 23
c 4		numu 14		W+ 24
b 5		tau- 15
t 6		nutau 16
		tau- 17
		nutau 18
'''



timestart = time.time()

events = Events ('lhe.root')
handle  = Handle ('LHEEventProduct')
label = ("source")

ROOT.gROOT.SetStyle('Plain') # white background
histolist = []

cuts = [0,0,0,0,0,0,0]
first = True
n_event = 0
n_used = 0
n_tot	= events.size()
n_base = 0

for event in events:
	n_event += 1
	if n_event%5000 == 0:
		print str(n_used) + ' / ' + str(n_event) + ' / ' + str(n_tot)

	event.getByLabel (label, handle)
	    
	lhe = handle.product()
	hepeup = lhe.hepeup()
	    
	if first:
	    first = False
	    for i in xrange(lhe.weights().size()):
		name = 'zww'+str(i)
		histolist.append(ROOT.TH1F(name, "WW Mass", 20, 1000, 3500))


	nw = 0
	nlep = 0
	keep_Event = True
	nowhad = True
	ww 	= ROOT.TLorentzVector(0.,0.,0.,0.)
	wlep 	= ROOT.TLorentzVector(0.,0.,0.,0.)
	whad 	= ROOT.TLorentzVector(0.,0.,0.,0.)
	whadNlist = []
	for i in xrange(hepeup.NUP):
		#w+-
		if abs(hepeup.IDUP[i])==24:
			nw 		= nw + 1
			pup 		= hepeup.PUP[i]
			wvector		= ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
			ww 		= ww + wvector
		#el+-
		if abs(hepeup.IDUP[i])==11:
			channel		= 'el'
			nlep		= nlep + 1
			pup	 	= hepeup.PUP[i]
			elvector	= ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
			if elvector.PseudoRapidity() > 2.5 or elvector.Pt() < 50:
				keep_Event = False
			wlep		= wlep + elvector
		#mu+-
		if abs(hepeup.IDUP[i])==13:
			channel		= 'mu'
			nlep		= nlep + 1
			pup	 	= hepeup.PUP[i]
			muvector	= ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
			if muvector.PseudoRapidity() > 2.1 or muvector.Pt() < 50:
				keep_Event = False
			wlep		= wlep + muvector
		#tau+-
		if abs(hepeup.IDUP[i])==15:
			keep_Event = False
		#n_el
		if abs(hepeup.IDUP[i])==12:
			pup 		= hepeup.PUP[i]
			nelvector	= ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
			if nelvector.Pt() < 80:
				keep_Event = False
			wlep 		= wlep + nelvector
		#n_mu
		if abs(hepeup.IDUP[i])==14:
			pup 		= hepeup.PUP[i]
			nmuvector	= ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
			if nmuvector.Pt() < 40:
				keep_Event = False
			wlep 		= wlep + nmuvector
		#w_had
		if nowhad and abs(hepeup.IDUP[i]) < 7:
			motherID = hepeup.MOTHUP[i].first
			if abs(hepeup.IDUP[motherID-1])==24:
				whadNlist.append(motherID)
				nowhad = False
				

	if len(whadNlist) == 1:
		pup	= hepeup.PUP[whadNlist[0]]
		whad	= whad + ROOT.TLorentzVector(pup[0],pup[1],pup[2],pup[3])
	
	if channel == 'el':
		delphi_met	= nelvector.DeltaPhi(whad)
	elif channel == 'mu':
		delphi_met	= nmuvector.DeltaPhi(whad)



	if nw != 2:
		print "didn't find 2 w's ???"
		keep_Event = False
	if nlep > 1:
		print "more than 1 lepton ???"
		keep_Event = False
	
##
	if keep_Event:	
		n_base =n_base +1
	
	if ww.M() <= 1 and keep_Event:
		keep_Event = False
		cuts[0] += 1
	if wlep.Pt() < 200. and keep_Event:
		keep_Event = False
		cuts[1] += 1
	if whad.Pt() < 200. and keep_Event:
		keep_Event = False
		cuts[2] += 1
	if abs(wlep.DeltaR(whad)) < math.pi/2. and keep_Event:
		keep_Event = False
		cuts[3] += 1
	if abs(delphi_met) < 2. and keep_Event:
		keep_Event = False
		cuts[4] += 1
	if abs(wlep.DeltaPhi(whad)) < 2. and keep_Event:
		keep_Event = False
		cuts[5] += 1
	if len(whadNlist) != 1:
		keep_Event = False
		cuts[6] += 1


	if keep_Event: 
		n_used += 1
		for i in xrange(lhe.weights().size()):
			weights = lhe.weights()
			histolist[i].Fill(ww.M(),weights[i].wgt)


'''if ww.M() >= 1000 and wlep.Pt() > 200. and whad.Pt() > 200. and abs(wlep.DeltaR(whad)) > math.pi/2. \
	and abs(delphi_met) > 2. and abs(wlep.DeltaPhi(whad)) > 2. and keep_Event: 
		for i in xrange(lhe.weights().size()):
			weights = lhe.weights()
			histolist[i].Fill(ww.M(),weights[i].wgt)
			n_used += 1'''



rfile = ROOT.TFile.Open("plots.root","RECREATE")
n = 0
c = ROOT.TCanvas()
c.SetLogy()
for hist in histolist:
    hist.Write()
    hist.SetLineColor(n)
    if n==0:
        first = False
        hist.Draw("h")
    else:
        hist.Draw("h,same")
        
c.Write()
rfile.Close()

timeend = time.time()
timeneeded = timeend-timestart

print 'used events: ' + str(n_used) + ' / ' + str(n_tot)
print 'used time : ' + str(round(timeneeded/60.,2)) + ' min'
print 'pass basline selection: ' +str(n_base)
print cuts
    
