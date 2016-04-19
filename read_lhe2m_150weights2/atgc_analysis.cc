#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TF2.h>
#include <TH3.h>
#include <TF3.h>
#include <TAxis.h>
#include <TString.h>
#include <iostream>
#include <string>




void atgc_analysis()
{
/*
	double cwwws[]		= {-12.,-6.,0.,6.,12.};
	double cwwws2[]		= {-2.,0.,2.};
	double ccws[]		= {-20.,-10.,0.,10.,20.};
	double ccws2[]		= {-3.5,0.,3.5};
	double cbs[]		= {-60.,-30.,0.,30.,60.};
	double cbs2[]		= {-5.,0.,5.};
	double vals[150][3];
	int count 		= 1; 

	vals[0][0]	= 0.;
	vals[0][1]	= 0.;
	vals[0][2]	= 0.;

	for(int i = 0; i<5; i++)
		for(int j = 0; j<5; j++)
			for(int k = 0; k<5; k++)
				{
					if(cwwws[i]==12 and ccws[j]==20 and cbs[k]==60)
						continue;
					if(cwwws[i]!=0 or ccws[j]!=0 or cbs[k]!=0)
					{
						vals[count][0] = cwwws[i];
						vals[count][1] = ccws[j];
						vals[count][2] = cbs[k];
						count++;
					}
				}
	for(int i = 0; i<3; i++)
		for(int j = 0; j<3; j++)
			for(int k = 0; k<3; k++)
				if(cwwws2[i]!=0 or ccws2[j]!=0 or cbs2[k]!=0)
				{
					vals[count][0] = cwwws2[i];
					vals[count][1] = ccws2[j];
					vals[count][2] = cbs2[k];
					count++;
				}
*/
	TFile * fileIn 	= TFile::Open("atgc_tree.root");
	TTree * treeIn	= (TTree*)fileIn->Get("tree");

	double MWW;
	std::vector<double> * weights;

	std::cout<<1<<std::endl;	
	treeIn->SetBranchAddress("weight", &weights);
	std::cout<<1.5<<std::endl;
	treeIn->SetBranchAddress("MWW", &MWW);
	treeIn->GetEntry(1);
	std::cout<<2<<std::endl;


	std::cout<<weights->size()<<std::endl;

	int n_weights	= weights->size();
	std::cout<<"creating "<<n_weights<<" histograms"<<std::endl;
/*
	std::vector<TH1F*> histos;
	int n_histos	= 0;
	double weight_tmp = (*weights)[2];
	
	for(int n = 0; n<n_weights; n++)
	{	
		if(n%10==0)
			std::cout<<n_histos<<" / "<<n_weights<<std::endl;

		TString histoname;
		histoname 	= "atgc" + to_string(n);
		TH1F *hist	= new TH1F(histoname,histoname,20,900,3500);
		hist->Sumw2();
		hist->GetXaxis()->SetTitle("M_{WW}");
		hist->GetYaxis()->SetTitle("N_{Events}");
		hist->SetTitle(TString(to_string(vals[n][0])) + "," + TString(to_string(vals[n][1]))+ "," + TString(to_string(vals[n][2])));
		for(int i = 0; i<treeIn->GetEntries(); i++)
		{
			treeIn->GetEntry(i);
			double weight_tmp = (*weights)[n];
			hist->Fill(MWW,weight_tmp);
		}
		histos.push_back(hist);
		n_histos++;

	}
	std::cout<<"done, "<<n_histos<<" histograms created"<<std::endl;
	

	int count_cwww_ccw		= 0;
	int count_cwww_cb		= 0;
	int count_ccw_cb		= 0;

	double norm_SM	= histos[0]->Integral();
	std::cout<<norm_SM<<std::endl;

	Float_t bins_tmp1[8]		= {-15,-9,-3,-1,1,3,9,15};
	Float_t bins_tmp2[8]		= {-25,-15,-5,-2,2,5,15,25};
	
	TH3F * hist_3D			= new TH3F("hist_3D","hist_3D",25,-12.5,12.5,41,-20.5,20.5,121,-60.5,60.5);
	TH1F * hist_cwww		= new TH1F("hist_cwww","hist_cwww",25,-12.5,12.5);
	TH1F * hist_ccw			= new TH1F("hist_ccw","hist_ccw",41,-20.5,20.5);
	TH1F * hist_cb			= new TH1F("hist_cb","hist_cb",121,-60.5,60.5);
	TGraph2D * graph_cwww_ccw 	= new TGraph2D();
	TGraph2D * graph_cwww_cb	= new TGraph2D();
	TGraph2D * graph_ccw_cb		= new TGraph2D();
	TH2F * hist_cwww_ccw		= new TH2F("hist_cwww_ccw","hist_cwww_ccw",25,-15,15,41,-25,25);

	hist_cwww_ccw->GetXaxis()->Set(7,bins_tmp1);
	hist_cwww_ccw->GetYaxis()->Set(7,bins_tmp2);

	for(int i = 0; i<n_weights; i++)
	{
		double rel_yield	= histos[i]->Integral()/norm_SM;
		hist_3D->Fill(vals[i][0],vals[i][1],vals[i][2],rel_yield);
		if(vals[i][0]==0)
		{
			graph_ccw_cb->SetPoint(count_ccw_cb,vals[i][1],vals[i][2],rel_yield);
			count_ccw_cb++;
			if(vals[i][1]==0)
				hist_cb->Fill(vals[i][2],rel_yield);
			if(vals[i][2]==0)
				hist_ccw->Fill(vals[i][1],rel_yield);
		}
		if(vals[i][1]==0)
		{
			graph_cwww_cb->SetPoint(count_cwww_cb,vals[i][0],vals[i][2],rel_yield);
			count_cwww_cb++;
			if(vals[i][2]==0)
				hist_cwww->Fill(vals[i][0],rel_yield);
		}
		if(vals[i][2]==0)
		{
			graph_cwww_ccw->SetPoint(count_cwww_ccw,vals[i][0],vals[i][1],rel_yield);
			count_cwww_ccw++;
			hist_cwww_ccw->Fill(vals[i][0],vals[i][1],rel_yield);
		}
	}

	TF3 * fitfunc 		= new TF3("parabel3D","[0]+[1]*x+[2]*y+[3]*z+[4]*x*x+[5]*y*y+[6]*z*z+[7]*x*y+[8]*x*z+[9]*y*z",-20,20,-30,30,-70,70);
	hist_3D->Fit(fitfunc);


	TFile * fileOut = TFile::Open("atgc_histos2.root","RECREATE");
	//for(unsigned int i = 0; i<histos.size(); i++)
	//	histos[i]->Write();
	hist_3D->Write();
	hist_cwww->Write();
	hist_ccw->Write();
	hist_cb->Write();
	hist_cwww_ccw->Write();
	graph_cwww_ccw->Write();
	fitfunc->Write();
	fileOut->Close();
	
*/
	std::cout<<3<<std::endl;
}
