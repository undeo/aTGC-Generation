from ROOT import *
from array import array

cwwws	= [-12,-6,-1,0,1,6,12]
cwwws1	= [-12,-6,0,6,12]
cwwws2	= [-1,0,1]
ccws	= [-20,-10,-2,0,2,10,20]
ccws1	= [-20,-10,0,10,20]
ccws2	= [-2,0,2]
cbs	= [-60,-30,-6,0,6,30,60]
cbs1	= [-60,-30,0,30,60]
cbs2	= [-6,0,6]
	
pvals	= [[0,0,0]]
pcount 	= 0
for cwww in cwwws1:
	for ccw in ccws1:
		for cb in cbs1:
			if cwww!=12 or ccw!=20 or cb!=60:
				if cwww!=0 or ccw!=0 or cb!=0:
					pvals_tmp	= []
					pvals_tmp.append(cwww)
					pvals_tmp.append(ccw)
					pvals_tmp.append(cb)
					pvals.append(pvals_tmp)
					pcount 		+= 1
for cwww in cwwws2:
	for ccw in ccws2:
		for cb in cbs2:
			if cwww!=0 or ccw!=0 or cb!=0:
				pvals_tmp	= []
				pvals_tmp.append(cwww)
				pvals_tmp.append(ccw)
				pvals_tmp.append(cb)
				pvals.append(pvals_tmp)
				pcount 		+= 1


fileIn		= TFile.Open('atgc_tree.root')
treeIn		= fileIn.Get('tree')

treeIn.GetEntry(1)
n_weights	= treeIn.weight.size()
#n_weights	= 5
n_entries	= treeIn.GetEntries()

print "creating " + str(n_weights) + " histograms"

histos 		= []
n_histos	= 0

for n in range(n_weights):
	if n_histos%5==0:
		print str(n_histos) + " / " + str(n_weights)


	histoname 	= "atgc" + str(n)
	hist 		= TH1F(histoname,histoname,20,900,3500)
	hist.GetXaxis().SetTitle("M_{WW}")
	hist.GetYaxis().SetTitle("N_{Events}")
	hist.SetTitle(str(pvals[n]))

	for i in range(n_entries):
		treeIn.GetEntry(i)
		weight_tmp	= treeIn.weight[n]
		hist.Fill(treeIn.MWW,weight_tmp)
	
	histos.append(hist)
	n_histos 	+= 1

print "done, " + str(n_histos) + " histograms created"

norm_SM		= histos[0].Integral()

bins_cwww	= array('d',[-15,-9,-3,-1,1,3,9,15])
bins_ccw	= array('d',[-25,-15,-5,-2,2,5,15,25])

hist_3d		= TH3F("hist_3d","hist_3d",25,-12.5,12.5,41,-20.5,20.5,121,-60.5,60.5)
hist_cwww	= TH1F("hist_cwww","hist_cwww",25,-12.5,12.5)
hist_ccw	= TH1F("hist_ccw","hist_ccw",41,-20.5,20.5)
hist_cb		= TH1F("hist_cb","hist_cb",121,-60.5,60.5)

hists_cwww_ccw	= []
for i in range(7):
	hist_cwww_ccw = TH2F("hist_cwww_ccw_cb%s"%i,"hist_cwww_ccw_cb%s"%i,25,-12.5,12.5,41,-20.5,20.5)
	hist_cwww_ccw.GetXaxis().Set(7,bins_cwww)
	hist_cwww_ccw.GetYaxis().Set(7,bins_ccw)
	hists_cwww_ccw.append(hist_cwww_ccw)


for i in range(n_weights):
	rel_yield	= histos[i].Integral()/norm_SM
	hist_3d.Fill(pvals[i][0],pvals[i][1],pvals[i][2],rel_yield)
	if pvals[i][1]==0 and pvals[i][2]==0:
		hist_cwww.Fill(pvals[i][0],rel_yield)
	if pvals[i][0]==0 and pvals[i][2]==0:
		hist_ccw.Fill(pvals[i][1],rel_yield)
	if pvals[i][0]==0 and pvals[i][1]==0:
		hist_cb.Fill(pvals[i][2],rel_yield)
	for j in range(7):
		if pvals[i][2]==cbs[j]:
			hists_cwww_ccw[j].Fill(pvals[i][0],pvals[i][2],rel_yield)

fileOut		= TFile.Open("atgc_histos.root","RECREATE")
#for i in range(n_histos):
#	histos[i].Write()
for i in range(7):
	hists_cwww_ccw[i].Write()
hist_3d.Write()
hist_cwww.Write()
hist_ccw.Write()
hist_cb.Write()
hist_cwww_ccw.Write()


fitfunc 	= TF3('parabel3D','[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y+[6]*z+[7]*z*z+[8]*x*z+[9]*y*z',-20,20,-30,30,-70,70)
hist_3d.Fit(fitfunc)

fitfunc.Write()
fileOut.Close()

