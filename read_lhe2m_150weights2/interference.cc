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
#include <TCanvas.h>
#include <TROOT.h>
#include <TPad.h>
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
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>


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



void interference(int channel)
{
	//channel: 1=el, 2=mu
	TString ch;	
	if(channel==1)
		ch = "el";
	if(channel==2)
		ch = "mu";

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

	TFile * fileIn 	= TFile::Open("root_files/atgc_tree.root");

	TTreeReader reader("tree",fileIn);
	TTreeReaderValue<std::vector<double>> weights(reader,"weight");
	TTreeReaderValue<double> MWW_tree(reader,"MWW");
	TTreeReaderValue<int> channel_tree(reader,"channel");

	RooRealVar cwww("cwww","cwww",0,-15,15);	cwww.setConstant(true);
	RooRealVar ccw("ccw","ccw",0,-25,25);		ccw.setConstant(true);
	RooRealVar cb("cb","cb",0,-70,70);		cb.setConstant(true);
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

	for(unsigned int i = 0; i<150; i++)
	{
		RooDataHist hist(("hist" + to_string(i)).c_str(),("hist" + to_string(i)).c_str(),RooArgSet(MWW));
		w.import(hist);
	}
	

	int tmp = 0;	
	while(reader.Next())
	{
		MWW.setVal(*MWW_tree);
		if(*channel_tree==channel)
			for (unsigned int i = 0; i<150; i++)
				w.data(("hist"+to_string(i)).c_str())->add(RooArgSet(MWW),(*weights)[i]);
		tmp++;
		if(tmp%5000==0)
			std::cout<<tmp<<std::endl;
	}

	w.factory("Exponential:SM_Pdf(MWW,a1[-0.1,0.1])");
	w.factory("Exponential:Pdf_cwww(MWW,a2[-0.1,0.1])");
	w.factory("Exponential:Pdf_cwww_lin(MWW,a22[-0.1,0.1])");
	w.factory("Exponential:Pdf_ccw(MWW,a3[-0.1,0.1])");
	w.factory("Exponential:Pdf_ccw_lin(MWW,a33[-0.1,0.1])");
	w.factory("Exponential:Pdf_cb(MWW,a4[-0.1,0.1])");
	w.factory("Exponential:Pdf_cb_lin(MWW,a44[-0.1,0.1])");
	//w.factory("Exponential:Int_cwww_ccw(MWW,a5[-0.1,0.1])");
	//w.factory("Exponential:Int_cwww_cb(MWW,a6[-0.1,0.1])");
	//w.factory("Exponential:Int_ccw_cb(MWW,a7[-0.1,0.1])");

	RooAbsPdf * SM_Pdf 	= w.pdf("SM_Pdf");
	RooAbsPdf * Pdf_cwww	= w.pdf("Pdf_cwww");
	RooAbsPdf * Pdf_cwww_lin= w.pdf("Pdf_cwww_lin");
	RooAbsPdf * Pdf_ccw	= w.pdf("Pdf_ccw");
	RooAbsPdf * Pdf_ccw_lin	= w.pdf("Pdf_ccw_lin");
	RooAbsPdf * Pdf_cb	= w.pdf("Pdf_cb");
	RooAbsPdf * Pdf_cb_lin	= w.pdf("Pdf_cb_lin");
	//RooAbsPdf * Int_cwww_ccw= w.pdf("Int_cwww_ccw");
	//RooAbsPdf * Int_cwww_cb	= w.pdf("Int_cwww_cb");
	//RooAbsPdf * Int_ccw_cb	= w.pdf("Int_ccw_cb");

	RooRealVar N_SM("N_SM","N_SM",w.data("hist0")->sumEntries());
	RooRealVar N__12("N__12","N__12",w.data("hist13")->sumEntries());
	RooRealVar N_12("N_12","N_12",w.data("hist112")->sumEntries());
	RooRealVar N__20("N__20","N__20",w.data("hist53")->sumEntries());
	RooRealVar N_20("N_20","N_20",w.data("hist72")->sumEntries());
	RooRealVar N__60("N__60","N__60",w.data("hist61")->sumEntries());
	RooRealVar N_60("N_60","N_60",w.data("hist64")->sumEntries());
	//RooRealVar N_12_20("N_12_20","N_12_20",w.data("hist122")->sumEntries());
	//RooRealVar N_12_60("N_12_60","N_12_60",w.data("hist114")->sumEntries());
	//RooRealVar N_20_60("N_20_60","N_20_60",w.data("hist74")->sumEntries());

	RooArgList N_norm_list(N_SM,N_12,N__12,N_20,N__20,N_60,N__60);
	N_norm_list.add(RooArgList(cwww,ccw,cb));
	//N_norm_list.add(RooArgList(N_12_20,N_12_60,N_20_60));

	TString cwww_f 		= "+((@1+@2)/2-@0)*(@7/12)**2+((@1-@2)/2)*(@7/12)";
	TString ccw_f 		= "+((@3+@4)/2-@0)*(@8/20)**2+((@3-@4)/2)*(@8/20)";
	TString cb_f		= "+((@5+@6)/2-@0)*(@9/60)**2+((@5-@6)/2)*(@9/60)";
	//TString cwww_ccw_f	= "+((@10+@0)-(@2+@4))*(@7/12)*(@8/20)";
	//TString cwww_cb_f	= "+((@11+@0)-(@2+@6))*(@7/12)*(@9/60)";
	//TString ccw_cb_f	= "+((@12+@0)-(@4+@6))*(@8/20)*(@9/60)";
	//RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+cb_f+cwww_ccw_f+cwww_cb_f+ccw_cb_f,N_norm_list);
	RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+cb_f,N_norm_list);

	RooFormulaVar N1("N1","@1/@0",RooArgList(Norm,N_SM));
	RooFormulaVar N2("N2","(((@1+@2)/2-@3)*(@4/12)**2)/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
	RooFormulaVar N3("N3","(((@1-@2)/2)*(@3/12))/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
	RooFormulaVar N4("N4","(((@1+@2)/2-@3)*(@4/20)**2)/@0",RooArgList(Norm,N_20,N__20,N_SM,ccw));
	RooFormulaVar N5("N5","(((@1-@2)/2)*(@3/20))/@0",RooArgList(Norm,N_20,N__20,ccw));
	RooFormulaVar N6("N6","(((@1+@2)/2-@3)*(@4/60)**2)/@0",RooArgList(Norm,N_60,N__60,N_SM,cb));
	RooFormulaVar N7("N7","(((@1-@2)/2)*(@3/60))/@0",RooArgList(Norm,N_60,N__60,cb));
	//RooFormulaVar N8("N8","(((@3+@4)-(@1+@2))*(@5/12)*(@6/20))/@0",RooArgList(Norm,N__12,N__20,N_12_20,N_SM,cwww,ccw));
	//RooFormulaVar N9("N9","(((@3+@4)-(@1+@2))*(@5/12)*(@6/60))/@0",RooArgList(Norm,N__12,N__60,N_12_60,N_SM,cwww,cb));
	//RooFormulaVar N10("N10","(((@3+@4)-(@1+@2))*(@5/20)*(@6/60))/@0",RooArgList(Norm,N__20,N__60,N_20_60,N_SM,ccw,cb));

	RooArgList N_list(N1,N2,N3,N4,N5,N6,N7);
	//N_list.add(RooArgList(N8,N9,N10));
	RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_cwww_lin,*Pdf_ccw,*Pdf_ccw_lin,*Pdf_cb,*Pdf_cb_lin);
	//Pdf_list.add(RooArgList(*Int_cwww_ccw,*Int_cwww_cb,*Int_ccw_cb));
	RooAddPdf model1("model","model",Pdf_list,N_list);
	model1.Print();


	RooRealVar * a1 	= w.var("a1");
	RooRealVar * a2		= w.var("a2");	a2->setConstant(true);
	RooRealVar * a22	= w.var("a22");	a22->setConstant(true);
	RooRealVar * a3		= w.var("a3");	a3->setConstant(true);
	RooRealVar * a33	= w.var("a33");	a33->setConstant(true);
	RooRealVar * a4		= w.var("a4");	a4->setConstant(true);
	RooRealVar * a44	= w.var("a44");	a44->setConstant(true);
	//RooRealVar * a5		= w.var("a5");	a5->setConstant(true);
	//RooRealVar * a6		= w.var("a6");	a6->setConstant(true);
	//RooRealVar * a7		= w.var("a7");	a7->setConstant(true);


	MWW.setRange(900,1500);
	//SM-fit
	cwww.setVal(0);	ccw.setVal(0); cb.setVal(0); 	
	model1.fitTo(*w.data("hist0"));
	a1->setConstant(true);
	//cwww-fit
	cwww.setVal(-12); ccw.setVal(0); cb.setVal(0); 
	a2->setConstant(false);
	a22->setConstant(false);
	model1.fitTo(*w.data("hist13"));
	a2->setConstant(true);
	a22->setConstant(true);
	//ccw-fit
	cwww.setVal(0);	ccw.setVal(-20); cb.setVal(0); 
	a3->setConstant(false);
	a33->setConstant(false);
	model1.fitTo(*w.data("hist53"));
	a3->setConstant(true);
	a33->setConstant(true);
	//cb-fit
	cwww.setVal(0); ccw.setVal(0); cb.setVal(-60);
	a4->setConstant(false);
	a44->setConstant(false);
	model1.fitTo(*w.data("hist61"));
	a4->setConstant(true);
	a44->setConstant(true);
/*	//int cwww-ccw-fit
	cwww.setVal(12); ccw.setVal(20); cb.setVal(0);
	w.var("a5")->setConstant(false);
	model1.fitTo(*w.data("hist122"));
	w.var("a5")->setConstant(true);
	//int cwww-cb-fit
	cwww.setVal(12); ccw.setVal(0); cb.setVal(60);
	w.var("a6")->setConstant(false);
	model1.fitTo(*w.data("hist114"));
	w.var("a6")->setConstant(true);
	//int cwww-ccw-fit
	cwww.setVal(0); ccw.setVal(20); cb.setVal(60);
	w.var("a7")->setConstant(false);
	model1.fitTo(*w.data("hist74"));
	w.var("a7")->setConstant(true);*/



	TCanvas c1("cwww0-ccw-10cb-30","cwww0-ccw-10cb-30",1);
	c1.cd();
	c1.SetLogy();
	RooPlot * plot1 = MWW.frame();
	cwww.setVal(0); ccw.setVal(-10); cb.setVal(-30);
	RooBinning bins(900,1500);
	bins.addUniform(20,900,1500);
	w.data("hist57")->plotOn(plot1,RooFit::Binning(bins));
	model1.plotOn(plot1);
	plot1->Draw();
	c1.Draw();
	c1.Update();

	int end;
	end = getchar();

	TCanvas c2("SM","SM",1);
	c2.cd();
	c2.SetLogy();
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	RooPlot * plot2 = MWW.frame();
	w.data("hist0")->plotOn(plot2,RooFit::Binning(bins));
	model1.plotOn(plot2);
	plot2->Draw();
	c2.Draw();
	c2.Update();
	end = getchar();/*

	TCanvas c3("ccw-20cb-60","ccw-20cb-60",1);
	c3.cd();
	c3.SetLogy();
	ccw.setVal(-20);
	cb.setVal(-60);
	RooPlot * plot3 = MWW.frame();
	w.data("hist51")->plotOn(plot3,RooFit::Binning(bins_ccw_cb));
	model1.plotOn(plot3);
	plot3->Draw();
	c3.Draw();
	c3.Update();
	end = getchar();

	TCanvas c4("c4","c4",1);
	c4.cd();
	c4.SetLogy();
	RooPlot * plot4 = MWW.frame();
	int i = 0;
	int j = 0;
	while(i<20 and j<60)
	{
		ccw.setVal(i);
		cb.setVal(j);
		model1.plotOn(plot4,RooFit::LineColor(kBlue+i));
		i++;
		j = j + 3;
		plot4->Draw();
		c4.Draw();
		std::cout<<i<<std::endl;	
	}
	c4.Update();
	end = getchar();



	RooBinning bins(900,3500);
	bins.addUniform(20,900,3500);

	RooRealVar * a1 = w.var("a1");
	RooRealVar * a2 = w.var("a2");
	RooRealVar * a22 = w.var("a22");

	ccw.setVal(0);
	a2->setConstant(true);
	a22->setConstant(true);
	model_lin.fitTo(*w.data("hist0"));
	a1->setConstant(true);

	a2->setConstant(false);
	a22->setConstant(false);
	ccw.setVal(20);
	model_lin.fitTo(*w.data("hist72"));
	a2->setConstant(true);
	a22->setConstant(true);

	TCanvas c1("c1","c1",1);
	c1.cd();
	c1.SetLogy();
	ccw.setVal(10);
	RooPlot * plot = MWW.frame();
	w.data("hist67")->plotOn(plot,RooFit::Binning(bins));
	model_lin.plotOn(plot);
	plot->Draw();
	c1.Draw();
	c1.Update();
	int end;
	end = getchar();

	TCanvas c2("c2","c2",1);
	c2.cd();
	c2.SetLogy();
	ccw.setVal(-3.5);
	RooPlot * plot2 = MWW.frame();
	w.data("hist134")->plotOn(plot2,RooFit::Binning(bins));
	model_lin.plotOn(plot2);
	plot2->Draw();
	c2.Draw();
	c2.Update();
	end = getchar();

	TCanvas c3("c3","c3",1);
	c3.cd();
	c3.SetLogy();
	ccw.setVal(-20);
	RooPlot * plot3 = MWW.frame();
	w.data("hist53")->plotOn(plot3,RooFit::Binning(bins));
	model_lin.plotOn(plot3);
	for(int i = 0; i<20; i++)
	{
		ccw.setVal(i);
		model_lin.plotOn(plot3,RooFit::LineColor(kBlue+i));
	}
	plot3->Draw();
	c3.Draw();
	c3.Update();
	end = getchar();*/

}


