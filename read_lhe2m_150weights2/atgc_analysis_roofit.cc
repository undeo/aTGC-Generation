#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TBranch.h>
#include <TFile.h>
#include <TH3.h>
#include <TF1.h>
#include <TF3.h>
#include <TString.h>
#include <TMinuit.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooCategory.h>
#include <RooArgList.h>
#include <RooBinning.h>


double cwwws[]		= {-12.,-6.,-2.,0.,2.,6.,12.};
double cwwws1[]		= {-12.,-6.,0.,6.,12.};
double cwwws2[]		= {-2.,0.,2.};
double ccws[]		= {-20.,-10.,-3.5,0.,3.5,10.,20.};
double ccws1[]		= {-20.,-10.,0.,10.,20.};
double ccws2[]		= {-3.5,0.,3.5};
double cbs[]		= {-60.,-30.,-10.,0.,10.,30.,60.};
double cbs1[]		= {-60.,-30.,0.,30.,60.};
double cbs2[]		= {-10.,0.,10.};
double vals[150][3];
const unsigned int npar	= 10;



double_t rel_yields[150];
double normSM;
RooWorkspace w("w","w");



void chi2_minuit(int &nDim, double_t *gout, double_t &result, double_t par[], int flag)
{
	result = 0;
	for(int i=0;i<150;i++)
	{
		if(i==0)
			normSM = w.data("hist0")->sumEntries();
		double_t cwww		= vals[i][0];
		double_t ccw		= vals[i][1];
		double_t cb		= vals[i][2];
		double_t parabel	= par[0] + par[1]*cwww+par[2]*cwww*cwww + par[3]*ccw+par[4]*ccw*ccw + par[5]*cb+par[6]*cb*cb + par[7]*cwww*ccw + par[8]*cwww*cb + par[9]*ccw*cb;
		if(rel_yields[i]!=0)
			result	+= ((rel_yields[i]-parabel)*(rel_yields[i]-parabel) )/( (0.005*rel_yields[i])*(0.005*rel_yields[i]));
	}	
}





void atgc_analysis_roofit()
{
	
	vals[0][0]	= 0.;
	vals[0][1]	= 0.;
	vals[0][2]	= 0.;

	int count = 1; 
	for(int i = 0; i<5; i++)
		for(int j = 0; j<5; j++)
			for(int k = 0; k<5; k++)
			{
				if(cwwws1[i]==12 and ccws1[j]==20 and cbs1[k]==60)
					continue;
				if(cwwws1[i]!=0 or ccws1[j]!=0 or cbs1[k]!=0)
				{
					vals[count][0] = cwwws1[i];
					vals[count][1] = ccws1[j];
					vals[count][2] = cbs1[k];
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

	TFile * fileIn 	= TFile::Open("atgc_tree.root");

	TTreeReader reader("tree",fileIn);
	TTreeReaderValue<std::vector<double>> weights(reader,"weight");
	TTreeReaderValue<double> MWW_tree(reader,"MWW");

	RooRealVar cwww("cwww","cwww",0,-15,15);
	RooRealVar ccw("ccw","ccw",0,-25,25);
	RooRealVar cb("cb","cb",0,-70,70);
	RooRealVar MWW("MWW","MWW",2500,900,3500);
	RooRealVar weight("weight","weight",1,0,500);
	RooCategory cat_cwww("cat_cwww","cat_cwww");
	RooCategory cat_ccw("cat_ccw","cat_ccw");
	RooCategory cat_cb("cat_cb","cat_cb");
	w.import(weight);
	for(int i = 0; i<7; i++)
	{
		cat_cwww.defineType((to_string(cwwws[i])).c_str(),i);
		cat_ccw.defineType((to_string(ccws[i])).c_str(),i);
		cat_cb.defineType((to_string(cbs[i])).c_str(),i);
	}
	RooArgSet parset(cwww,ccw,cb,cat_cwww,cat_ccw,cat_cb,weight);

	RooDataSet data3D("data3D","data3D",parset,"weight");

	double normSM		= -1000;
	int tmp = 1;
	for(unsigned int i = 0; i<150; i++)
	{
		RooDataHist hist(("hist" + to_string(i)).c_str(),("hist" + to_string(i)).c_str(),RooArgSet(MWW));
		w.import(hist);
	}
	
	while(reader.Next())
	{
		MWW.setVal(*MWW_tree);
		for (unsigned int i = 0; i<150; i++)
			w.data(("hist"+to_string(i)).c_str())->add(RooArgSet(MWW),(*weights)[i]);
		tmp++;
		if(tmp%5000==0)
			std::cout<<tmp<<std::endl;
	}

	normSM	= w.data("hist0")->sumEntries();

	TH3F * hist4fit	= new TH3F("hist4fit","hist4fit",25,-12.5,12.5,41,-20.5,20.5,121,-60.5,60.5);


	for(int i = 0; i<150; i++)
	{
		double rel_yield	= (w.data(("hist"+to_string(i)).c_str())->sumEntries())/normSM;
		rel_yields[i]		= rel_yield;
		cat_cwww.setLabel((to_string(vals[i][0])).c_str());
		cat_ccw.setLabel((to_string(vals[i][1])).c_str());
		cat_cb.setLabel((to_string(vals[i][2])).c_str());
		cwww.setVal(vals[i][0]);
		ccw.setVal(vals[i][1]);
		cb.setVal(vals[i][2]);
		weight.setVal(rel_yield);
		data3D.add(parset,rel_yield);
		hist4fit->Fill(cwww.getVal(),ccw.getVal(),cb.getVal(),rel_yield);
	}

	//set initial parameter values
	double_t par[10];
	double_t parError[10];
	par[0]	= 1;
	par[1]	= 0.1;
	par[2]	= 0.1;
	par[3]	= 0.1;
	par[4]	= 0.1;
	par[5]	= 0.1;
	par[6]	= 0.1;
	par[7]	= 0.1;
	par[8]	= 0.1;
	par[9]	= 0.1;
/*
	double_t chi2_b4 = 0;
	chi2_minuit(npar,NULL,chi2_b4,par,0);
	std::cout<<chi2_b4<<std::endl;
*/

	TMinuit minuit(10);
	minuit.SetFCN(chi2_minuit);

	minuit.DefineParameter(0,"par0",1,0.1,0,10);
	minuit.DefineParameter(1,"par1",0.1,0.1,-1,1);
	minuit.DefineParameter(2,"par2",0.1,0.1,-1,1);
	minuit.DefineParameter(3,"par3",0.1,0.1,-1,1);
	minuit.DefineParameter(4,"par4",0.1,0.1,-1,1);
	minuit.DefineParameter(5,"par5",0.1,0.1,-1,1);
	minuit.DefineParameter(6,"par6",0.1,0.1,-1,1);
	minuit.DefineParameter(7,"par7",0.1,0.1,-1,1);
	minuit.DefineParameter(8,"par8",0.1,0.1,-1,1);
	minuit.DefineParameter(9,"par9",0.1,0.1,-1,1);

	minuit.mnderi();
	minuit.Migrad();

	for(unsigned int i=0; i<10;i++)
	{
		minuit.GetParameter(i,par[i],parError[i]);
		std::cout<<"par"<<i<<" : "<<par[i]<< " +- "<<parError[i]<<std::endl;
	}

	TF1 * func_cwww	= new TF1("func_cwww","[0]+[1]*x+[2]*x*x+[3]*[10]+[4]*[10]*[10]+[5]*[11]+[6]*[11]*[11]+[7]*x*[10]+[8]*x*[11]+[9]*[10]*[11]",-15,15);	
	TF1 * func_ccw	= new TF1("func_cwww","[0]+[1]*[10]+[2]*[10]*[10]+[3]*x+[4]*x*x+[5]*[11]+[6]*[11]*[11]+[7]*[10]*x+[8]*[10]*[11]+[9]*x*[11]",-25,25);
	TF1 * func_cb	= new TF1("func_cb","[0]+[1]*[10]+[2]*[10]*[10]+[3]*[11]+[4]*[11]*[11]+[5]*x+[6]*x*x+[7]*[10]*[11]+[8]*[10]*x+[9]*[11]*x",-75,75);
	for(int i=0; i<10; i++)
	{
		func_cwww->SetParameter(i,par[i]);
		func_ccw->SetParameter(i,par[i]);
		func_cb->SetParameter(i,par[i]);
	}

	
	TFile * fileOut	= TFile::Open("atgc_histos2.root","RECREATE");


	for(int i=0; i<7; i++)
		for(int j=0; j<7; j++)
		{
			func_cwww->SetParameter(10,ccws[i]);
			func_cwww->SetParameter(11,cbs[j]);
			func_cwww->SetName(("func_cwww"+to_string(i)+to_string(j)).c_str());
			func_cwww->Write();
		
			func_ccw->SetParameter(10,cwwws[i]);
			func_ccw->SetParameter(11,cbs[j]);
			func_ccw->SetName(("func_ccw"+to_string(i)+to_string(j)).c_str());
			func_ccw->Write();

			func_cb->SetParameter(10,cwwws[i]);
			func_cb->SetParameter(11,ccws[j]);
			func_cb->SetName(("func_cb"+to_string(i)+to_string(j)).c_str());
			func_cb->Write();
		}
	



	TF3 * fitfunc 	= new TF3("fitfunc","[0]+[1]*x+[2]*x*x+[3]*y+[4]*y*y+[5]*z+[6]*z*z+[7]*x*y+[8]*x*z+[9]*y*z");
	hist4fit->Fit(fitfunc);

	Float_t cwww_bins[8]	= {-15,-9,-3,-1,1,3,9,15};
	Float_t ccw_bins[8]	= {-20,-15,-5,-2,2,5,15,20};
	Float_t cb_bins[8]	= {-75,-45,-15,-5,5,15,45,75};

	RooBinning bins_cwww(-15,15);
	RooBinning bins_cwww2(-15,15);
	RooBinning bins_ccw(-25,25);
	RooBinning bins_ccw2(-25,25);
	RooBinning bins_cb(-75,75);
	RooBinning bins_cb2(-75,75);
	for(int i = 4;i<7;i++)
	{
		bins_cwww.addBoundaryPair(cwww_bins[i]);
		bins_ccw.addBoundaryPair(ccw_bins[i]);
		bins_cb.addBoundaryPair(cb_bins[i]);
		if(i>4)
		{
			bins_cwww2.addBoundaryPair(cwww_bins[i]);
			bins_ccw2.addBoundaryPair(ccw_bins[i]);
			bins_cb2.addBoundaryPair(cb_bins[i]);
		}
	}


	//make 2d histograms
	for(int i = 0; i<7; i++)
	{
		RooBinning binning_cwww(0,1);
		RooBinning binning_ccw(0,1);
		RooBinning binning_cb(0,1);
		if(i>=2 and i<=4)
		{
			binning_cwww	= bins_cwww;
			binning_ccw	= bins_ccw;
			binning_cb	= bins_cb;
		}else
		{
			binning_cwww	= bins_cwww2;
			binning_ccw	= bins_ccw2;
			binning_cb	= bins_cb2;
		}	
		TH1 * hist1	= data3D.reduce(RooArgSet(cwww,ccw),("cat_cb==cat_cb::"+to_string(i)).c_str())->createHistogram(("cb"+to_string(i)).c_str(),cwww,RooFit::Binning(binning_cwww),RooFit::YVar(ccw,RooFit::Binning(binning_ccw)));
		TH1 * hist2	= data3D.reduce(RooArgSet(cwww,cb),("cat_ccw==cat_ccw::"+to_string(i)).c_str())->createHistogram(("ccw"+to_string(i)).c_str(),cwww,RooFit::Binning(binning_cwww),RooFit::YVar(cb,RooFit::Binning(binning_cb)));
		TH1 * hist3	= data3D.reduce(RooArgSet(ccw,cb),("cat_cwww==cat_cwww::"+to_string(i)).c_str())->createHistogram(("cwww"+to_string(i)).c_str(),ccw,RooFit::Binning(binning_ccw),RooFit::YVar(cb,RooFit::Binning(binning_cb)));	
		hist1->GetZaxis()->SetRangeUser(1,4.5);
		hist2->GetZaxis()->SetRangeUser(1,4.5);
		hist3->GetZaxis()->SetRangeUser(1,4.5);
		if(i==2 or i==4)
		{
			hist1->GetZaxis()->SetRangeUser(1,1.15);
			hist1->GetXaxis()->SetRangeUser(-3,3);
			hist1->GetYaxis()->SetRangeUser(-5,5);
			
			hist2->GetZaxis()->SetRangeUser(1,1.15);
			hist2->GetXaxis()->SetRangeUser(-3,3);
			hist2->GetYaxis()->SetRangeUser(-15,15);
			
			hist3->GetZaxis()->SetRangeUser(1,1.15);
			hist3->GetXaxis()->SetRangeUser(-5,5);
			hist3->GetYaxis()->SetRangeUser(-15,15);
		}
		hist1->Write();
		hist2->Write();
		hist3->Write();
	}	


	//make 1d histograms
	for(int i=0; i<7; i++)
		for(int j=0; j<7; j++)
		{
			TH1 * hist1	= data3D.reduce(RooArgSet(cwww),("cat_ccw==cat_ccw::"+to_string(i)+"&&cat_cb==cat_cb::"+to_string(j)).c_str())->createHistogram(("cat"+to_string(i)+to_string(j)).c_str(),cwww,RooFit::Binning(bins_cwww));
			TH1 * hist2	= data3D.reduce(RooArgSet(ccw),("cat_cwww==cat_cwww::"+to_string(i)+"&&cat_cb==cat_cb::"+to_string(j)).c_str())->createHistogram(("cat"+to_string(i)+to_string(j)).c_str(),ccw,RooFit::Binning(bins_ccw));
			TH1 * hist3	= data3D.reduce(RooArgSet(cb),("cat_cwww==cat_cwww::"+to_string(i)+"&&cat_ccw==cat_ccw::"+to_string(j)).c_str())->createHistogram(("cat"+to_string(i)+to_string(j)).c_str(),cb,RooFit::Binning(bins_cb));
			hist1->Write();
			hist2->Write();
			hist3->Write();
		}

	w.import(data3D);
	std::cout<<tmp<<std::endl;
	w.Write();
	fileOut->Close();
}
