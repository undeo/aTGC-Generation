#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "math.h"
#include "iostream"
#include "ctime"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/EDProductGetter.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"



void read_lhe()
{

	clock_t begin = clock();

	gSystem->Load("libFWCoreFWLite.so"); 
	AutoLibraryLoader::enable();
	gSystem->Load("libDataFormatsFWLite.so");

	TFile file("root_files/lhe_WW.root");
	fwlite::Event ev(&file);

	TFile * outfile = TFile::Open("root_files/atgc_tree.root","RECREATE");
	TTree * tree 	= new TTree("tree","tree");
	double MWW_tree, DPhi_MET_tree, DPhi_Wlep_tree, whadpt_tree, wleppt_tree, deltaR_tree;
	int channel_tree;
	std::vector<double> atgc_weights_tree;
	tree->Branch("MWW",&MWW_tree);
	tree->Branch("weight",&atgc_weights_tree);
	tree->Branch("deltaPhi_WjetMet",&DPhi_MET_tree);
	tree->Branch("jet_pt",&whadpt_tree);
	tree->Branch("W_pt",&wleppt_tree);
	tree->Branch("deltaR_LeptonWJet",&deltaR_tree);
	tree->Branch("deltaPhi_WJetWlep",&DPhi_Wlep_tree);
	tree->Branch("channel",&channel_tree);

	int n_events = 0, n_used = 0, n_el = 0, n_mu = 0, corr_ev = 0;
	


	for( ev.toBegin(); ! ev.atEnd(); ++ev) 
	{

		bool keep_Event = true;
	 	fwlite::Handle<LHEEventProduct> lhe;
		//now can access data
		lhe.getByLabel(ev,"source");

		TLorentzVector ww, wlep, elvector, muvector, nelvector, nmuvector;	
		int nwhad = 0, nlep = 0, nw = 0, whadID = -1000, channel = -1;
		double delphi_met = 0, delphi_lep, delR = 0;

		if(lhe.product()->weights().size()!=150)
			{
				//std::cout<<"only "<<weights.size()<<" weights!"<<std::endl;
				corr_ev++;
				continue;
			}
		for(int i=0; i<150; i++)
		{
			if(lhe.product()->weights()[i].wgt>1000)
				corr_ev++;
				break;
		}

		n_events++;
		for(int i = 0; i < lhe->hepeup().NUP; i++)
		{
			//W+-
			if(abs(lhe->hepeup().IDUP[i]) == 24)
			{
				nw++;
				TLorentzVector wvector(	lhe->hepeup().PUP[i][0],
							lhe->hepeup().PUP[i][1],
							lhe->hepeup().PUP[i][2],
							lhe->hepeup().PUP[i][3]);
				ww += wvector;
			}
			//el+-
			if(abs(lhe->hepeup().IDUP[i])==11)
			{
				channel		= 1;	
				nlep++;
				TLorentzVector elvector_tmp(lhe->hepeup().PUP[i][0],
								lhe->hepeup().PUP[i][1],
								lhe->hepeup().PUP[i][2],
								lhe->hepeup().PUP[i][3]);
				elvector	+= elvector_tmp;
				if(abs(elvector.PseudoRapidity()) > 2.5 or elvector.Pt() < 50)
					keep_Event = false;
				wlep		+= elvector;
			}
			//mu+-
			if(abs(lhe->hepeup().IDUP[i])==13)
			{
				channel		= 2;
				nlep++;
				TLorentzVector muvector_tmp(lhe->hepeup().PUP[i][0],
								lhe->hepeup().PUP[i][1],
								lhe->hepeup().PUP[i][2],
								lhe->hepeup().PUP[i][3]);
				muvector	+= muvector_tmp;
				if(abs(muvector.PseudoRapidity()) > 2.1 or muvector.Pt() < 50)
					keep_Event = false;
				wlep		+= muvector;
			}
			//tau+-
			if(abs(lhe->hepeup().IDUP[i])==15)
			{
				nlep++;
				keep_Event = false;
			}
			
			//n_el
			if(abs(lhe->hepeup().IDUP[i])==12)
			{	
				TLorentzVector nelvector_tmp(lhe->hepeup().PUP[i][0],
								lhe->hepeup().PUP[i][1],
								lhe->hepeup().PUP[i][2],
								lhe->hepeup().PUP[i][3]);
				nelvector 	+= nelvector_tmp;
				if(nelvector.Pt() < 80)
					keep_Event = false;
				wlep 		+= nelvector;
			}
			//n_mu
			if(abs(lhe->hepeup().IDUP[i])==14)
			{
				TLorentzVector nmuvector_tmp(lhe->hepeup().PUP[i][0],
								lhe->hepeup().PUP[i][1],
								lhe->hepeup().PUP[i][2],
								lhe->hepeup().PUP[i][3]);
				nmuvector	+= nmuvector_tmp;
				if(nmuvector.Pt() < 40)
					keep_Event = false;
				wlep 		= wlep + nmuvector;
			}
			//w_had
			if(abs(lhe->hepeup().IDUP[i]) < 7)
			{
				int motherID = lhe->hepeup().MOTHUP[i].first - 1;
				if(abs(lhe->hepeup().IDUP[motherID])==24 and motherID!=whadID)
				{
						nwhad++;
						whadID = motherID;
				}
			}
		}
		
		TLorentzVector whad(lhe->hepeup().PUP[whadID][0],
					lhe->hepeup().PUP[whadID][1],
					lhe->hepeup().PUP[whadID][2],
					lhe->hepeup().PUP[whadID][3]);


		if(nw != 2 or nlep != 1 or nwhad != 1)
			std::cout << "something went wrong! (-> nw,nlep,nwhad)" << "(" << nw << " , " << nlep << " , " << nwhad << ")" << std::endl;	


		if(ww.M() < 600. or ww.M() > 3500.)
			keep_Event = false;

		if(keep_Event)
		{
			if(channel == 1)
			{
				delphi_met	= nelvector.DeltaPhi(whad);
				delR		= elvector.DeltaR(whad);
			}
			if(channel == 2)
			{
				delphi_met 	= nmuvector.DeltaPhi(whad);
				delR		= muvector.DeltaR(whad);

			}

			delphi_lep	= wlep.DeltaPhi(whad);

			if(wlep.Pt() > 200. and whad.Pt() > 200. and abs(whad.PseudoRapidity()) < 2.4 and delR > M_PI/2. and abs(delphi_lep) > 2. and abs(delphi_met) > 2.)		
			{
				
				n_used++;

				std::vector <double> weights;
				int n_weights = lhe.product()->weights().size();

				for(int i = 0; i < n_weights; i++)
				{
					double weight = lhe->weights()[i].wgt;
					weights.push_back(weight);
				}
				
				
				atgc_weights_tree	= weights;
				MWW_tree 		= ww.M();
				DPhi_MET_tree		= delphi_met;
				DPhi_Wlep_tree		= delphi_lep;
				whadpt_tree		= whad.Pt();
				wleppt_tree		= wlep.Pt();
				deltaR_tree		= delR;
				channel_tree		= channel;
				tree->Fill();
			}
		
		}

		if(n_events%50000==0)
			std::cout << n_used << " / " << n_events << " , corrupted events: " << corr_ev << std::endl;

		if(n_events==2000000)
			break;	

	}

	tree->Write();
	std::cout << n_used << " events used " << std::endl;
	outfile->Close();
	clock_t end = clock();
	double time_needed = double(end - begin) / CLOCKS_PER_SEC;
	std::cout << "total time nedded: " << time_needed/60 << " min" << std::endl;
	//exit(0);
}
