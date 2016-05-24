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
double ccws[]		= {-20.,-10.,-3.5,0.,3.5,10.,20.};
double cbs[]		= {-60.,-30.,-10.,0.,10.,30.,60.};

double normSM;
RooWorkspace w("w","w");
RooWorkspace w2("w2","w2");


void interference(int channel)
{
	//channel: 1=el, 2=mu
	TString ch;	
	if(channel==1)
		ch = "el";
	if(channel==2)
		ch = "mu";

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
	
	double number_of_events = 0;
	double sum_of_weights[150] = {0};
	while(reader.Next())
		if(*channel_tree==channel)
		{
			for(int i=0; i<150; i++)
				sum_of_weights[i] += (*weights)[i];
			number_of_events += 1;
		}

	//for(int i=0; i<150; i++)
	//	std::cout<<sum_of_weights[i]<<std::endl;



	std::cout<<"number of events: "<<number_of_events<<std::endl;

	reader.SetEntry(0);
	int tmp = 0;	
	while(reader.Next())
	{
		MWW.setVal(*MWW_tree);
		if(*channel_tree==channel)
			for (int i = 0; i<150; i++)
			{
				//double weight_tmp = (*weights)[i] * (number_of_events/sum_of_weights[0]);
				double weight_tmp = (*weights)[i];
				w.data(("hist"+to_string(i)).c_str())->add(RooArgSet(MWW),weight_tmp);
			}
		tmp++;
		if(tmp%5000==0)
			std::cout<<tmp<<std::endl;
	}


	w.factory("Exponential:SM_Pdf(MWW,a1[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_cwww(MWW,a2[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_ccw(MWW,a3[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_cb(MWW,a4[-0.001,-0.01,0.])");
	w.factory("Exponential:Int_cwww_ccw(MWW,a5[-0.001,-0.01,0.])");
	w.factory("Exponential:Int_cwww_cb(MWW,a6[-0.001,-0.1,0.])");
	w.factory("Exponential:Int_ccw_cb(MWW,a7[-0.001,-0.01,0.])");

	RooAbsPdf * SM_Pdf 	= w.pdf("SM_Pdf");
	RooAbsPdf * Pdf_cwww	= w.pdf("Pdf_cwww");
	RooAbsPdf * Pdf_ccw	= w.pdf("Pdf_ccw");
	RooAbsPdf * Pdf_cb	= w.pdf("Pdf_cb");
	RooAbsPdf * Int_cwww_ccw= w.pdf("Int_cwww_ccw");
	RooAbsPdf * Int_cwww_cb	= w.pdf("Int_cwww_cb");
	RooAbsPdf * Int_ccw_cb	= w.pdf("Int_ccw_cb");

	RooRealVar N_SM("N_SM","N_SM",w.data("hist0")->sumEntries());			//hist0
	RooRealVar N__12("N_cwww__12","N_cwww__12",w.data("hist13")->sumEntries());	//hist13, hist128
	RooRealVar N_12("N_cwww_12","N_cwww_12",w.data("hist112")->sumEntries());	//hist112, hist145
	RooRealVar N__20("N_ccw__20","N_ccw__20",w.data("hist53")->sumEntries());	//hist53, hist134
	RooRealVar N_20("N_ccw_20","N_ccw_20",w.data("hist72")->sumEntries());		//hist72, hist139
	RooRealVar N__60("N_cb__60","N_cb__60",w.data("hist61")->sumEntries());		//hist61, hist136
	RooRealVar N_60("N_cb_60","N_cb_60",w.data("hist64")->sumEntries());		//hist64, hist137
	RooRealVar N__12_20("N_cwww_ccw__12_20","N_cwww_ccw__12_20",w.data("hist23")->sumEntries());
	RooRealVar N__12_60("N_cwww_cb__12_60","N_cwww_cb__12_60",w.data("hist15")->sumEntries());
	RooRealVar N__20_60("N_ccw_cb__20_60","N_ccw_cb__20_60",w.data("hist55")->sumEntries());

	RooRealVar N_4norm("N_4norm","N_4norm",w.data("hist1")->sumEntries());
	w2.import(N_4norm);

	w2.import(N_SM);
	w2.import(N__12);
	w2.import(N_12);
	w2.import(N__20);
	w2.import(N_20);
	w2.import(N__60);
	w2.import(N_60);
	w2.import(N__12_20);
	w2.import(N__12_60);
	w2.import(N__20_60);

	std::cout<<N_SM.GetName()<<": "<<N_SM.getVal()<<std::endl;
	std::cout<<N__12.GetName()<<": "<<N__12.getVal()<<std::endl;
	std::cout<<N_12.GetName()<<": "<<N_12.getVal()<<std::endl;
	std::cout<<N__20.GetName()<<": "<<N__20.getVal()<<std::endl;
	std::cout<<N_20.GetName()<<": "<<N_20.getVal()<<std::endl;
	std::cout<<N__60.GetName()<<": "<<N__60.getVal()<<std::endl;
	std::cout<<N_60.GetName()<<": "<<N_60.getVal()<<std::endl;
	std::cout<<N__12_20.GetName()<<": "<<N__12_20.getVal()<<std::endl;
	std::cout<<N__12_60.GetName()<<": "<<N__12_60.getVal()<<std::endl;
	std::cout<<N__20_60.GetName()<<": "<<N__20_60.getVal()<<std::endl;

	RooArgList N_norm_list(N_SM,N_12,N__12,N_20,N__20,N_60,N__60);
	N_norm_list.add(RooArgList(cwww,ccw,cb));
	N_norm_list.add(RooArgList(N__12_20,N__12_60,N__20_60));

	TString cwww_f 		= "+((@1+@2)/2-@0)*(@7/12)**2";
	TString ccw_f 		= "+((@3+@4)/2-@0)*(@8/20)**2";
	TString cb_f		= "+((@5+@6)/2-@0)*(@9/60)**2";
	TString cwww_ccw_f	= "+((@10+@0)-(@2+@4))*(@7/12)*(@8/20)";
	TString cwww_cb_f	= "+((@11+@0)-(@2+@6))*(@7/12)*(@9/60)";
	TString ccw_cb_f	= "+((@12+@0)-(@4+@6))*(@8/20)*(@9/60)";
	RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+cb_f+cwww_ccw_f+cwww_cb_f+ccw_cb_f,N_norm_list);

	RooFormulaVar N1("N1","@1/@0",RooArgList(Norm,N_SM));
	RooFormulaVar N2("N2","(((@1+@2)/2-@3)*(@4/12)**2)/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
	RooFormulaVar N3("N3","(((@1-@2)/2)*(@3/12))/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
	RooFormulaVar N4("N4","(((@1+@2)/2-@3)*(@4/20)**2)/@0",RooArgList(Norm,N_20,N__20,N_SM,ccw));
	RooFormulaVar N5("N5","(((@1-@2)/2)*(@3/20))/@0",RooArgList(Norm,N_20,N__20,ccw));
	RooFormulaVar N6("N6","(((@1+@2)/2-@3)*(@4/60)**2)/@0",RooArgList(Norm,N_60,N__60,N_SM,cb));
	RooFormulaVar N7("N7","(((@1-@2)/2)*(@3/60))/@0",RooArgList(Norm,N_60,N__60,cb));
	RooFormulaVar N8("N8","(((@3+@4)-(@1+@2))*(@5/12)*(@6/20))/@0",RooArgList(Norm,N__12,N__20,N__12_20,N_SM,cwww,ccw));
	RooFormulaVar N9("N9","(((@3+@4)-(@1+@2))*(@5/12)*(@6/60))/@0",RooArgList(Norm,N__12,N__60,N__12_60,N_SM,cwww,cb));
	RooFormulaVar N10("N10","(((@3+@4)-(@1+@2))*(@5/20)*(@6/60))/@0",RooArgList(Norm,N__20,N__60,N__20_60,N_SM,ccw,cb));



	RooArgList N_list(N1,N2,N4,N6);
	N_list.add(RooArgList(N8,N9,N10));
	RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_ccw,*Pdf_cb);
	Pdf_list.add(RooArgList(*Int_cwww_ccw,*Int_cwww_cb,*Int_ccw_cb));
	Pdf_list.Print();
	N_list.Print();
	RooAddPdf model1("model","model",Pdf_list,N_list);
	model1.Print();


	RooRealVar * a1 	= w.var("a1");
	RooRealVar * a2		= w.var("a2");	a2->setConstant(true);
	RooRealVar * a3		= w.var("a3");	a3->setConstant(true);
	RooRealVar * a4		= w.var("a4");	a4->setConstant(true);
	RooRealVar * a5		= w.var("a5");	a5->setConstant(true);
	RooRealVar * a6		= w.var("a6");	a6->setConstant(true);
	RooRealVar * a7		= w.var("a7");	a7->setConstant(true);


//SM-fit
	cwww.setVal(0);	ccw.setVal(0); cb.setVal(0); 	
	model1.fitTo(*w.data("hist0"));//hist0
	a1->setConstant(true);
//cwww-fit
	cwww.setVal(-12); ccw.setVal(0); cb.setVal(0); 
	a2->setConstant(false);
	model1.fitTo(*w.data("hist13"));//hist13, hist128
	a2->setConstant(true);
//ccw-fit
	cwww.setVal(0);	ccw.setVal(-20); cb.setVal(0); 
	a3->setConstant(false);
	model1.fitTo(*w.data("hist53"));//hist53, hist134
	a3->setConstant(true);/*
//cb-fit
	cwww.setVal(0); ccw.setVal(0); cb.setVal(-60);
	a4->setConstant(false);
	//a44->setConstant(false);
	//a44->setVal(0);
	model1.fitTo(*w.data("hist61"));//hist61, hist 136
	a4->setConstant(true);/*
//int cwww-ccw-fit
	cwww.setVal(-12); ccw.setVal(20); cb.setVal(0);
	w.var("a5")->setConstant(false);
	model1.fitTo(*w.data("hist23"));
	w.var("a5")->setConstant(true);
//int cwww-cb-fit
	cwww.setVal(-12); ccw.setVal(0); cb.setVal(60);
	w.var("a6")->setConstant(false);
	model1.fitTo(*w.data("hist15"));
	w.var("a6")->setConstant(true);
//int cwww-ccw-fit
	cwww.setVal(0); ccw.setVal(-20); cb.setVal(60);
	w.var("a7")->setConstant(false);
	model1.fitTo(*w.data("hist55"));
	w.var("a7")->setConstant(true);

	w2.import(*a1);
	w2.import(*a2);
	w2.import(*a3);
	w2.import(*a4);
	w2.import(*a5);
	w2.import(*a6);
	w2.import(*a7);

	TFile * fileOut = new TFile("genlevel_WZ_"+ch+".root","RECREATE");
	w2.Write();
	fileOut->Close();

	RooBinning bins(900,3500);
	bins.addUniform(20,900,3500);


	int end;

	TCanvas c7("1","1",1);
	c7.cd(); c7.SetLogy();
	RooPlot * plot7 = MWW.frame();
	cwww.setVal(-12); ccw.setVal(-20); cb.setVal(0);
	w.data("hist3")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(1),RooFit::LineWidth(2),RooFit::DrawOption("E"));
	model1.plotOn(plot7,RooFit::LineColor(1),RooFit::Normalization(w.data("hist3")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(-12); ccw.setVal(-10); cb.setVal(0);
	w.data("hist8")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(2),RooFit::LineWidth(2),RooFit::DrawOption("E"));
	model1.plotOn(plot7,RooFit::LineColor(2),RooFit::Normalization(w.data("hist8")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(-12); ccw.setVal(0); cb.setVal(0);
	w.data("hist13")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(3),RooFit::LineWidth(2),RooFit::DrawOption("E"));
	model1.plotOn(plot7,RooFit::LineColor(3),RooFit::Normalization(w.data("hist13")->sumEntries(),RooAbsReal::NumEvent));
	/*cwww.setVal(0); ccw.setVal(-10); cb.setVal(-60);
	w.data("hist56")->plotOn(plot7,RooFit::Binning(bins),RooFit::MarkerColor(4));
	model1.plotOn(plot7,RooFit::LineColor(4),RooFit::Normalization(w.data("hist56")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(-20); cb.setVal(-60);
	w.data("hist51")->plotOn(plot7,RooFit::Binning(bins),RooFit::MarkerColor(5));
	model1.plotOn(plot7,RooFit::LineColor(5),RooFit::Normalization(w.data("hist51")->sumEntries(),RooAbsReal::NumEvent));*/
	l7 = new TLegend(0.65,0.55,0.89,0.89);
	l7->SetFillColor(0);
	l7->SetBorderSize(0);
	l7->SetTextSize(0.03);
	l7->AddEntry(plot7->getObject(1),"#frac{c_{W}}{#Lambda^{2}} = -20 TeV^{-2}","le");
	l7->AddEntry(plot7->getObject(3),"#frac{c_{W}}{#Lambda^{2}} = -10 TeV^{-2}","le");
	l7->AddEntry(plot7->getObject(5),"#frac{c_{W}}{#Lambda^{2}} = 0 TeV^{-2}","le");
	l7->SetHeader("#frac{c_{WWW}}{#Lambda^{2}}=-12 TeV^{-2}, #frac{c_{B}}{#Lambda^{2}}=0 TeV^{-2}");
	plot7->SetTitle("");
	plot7->GetXaxis()->SetTitle("M_{WZ}");
	plot7->GetYaxis()->SetTitle("arb. units");
	plot7->SetMinimum(0.04);
	plot7->Draw();
	c7.Draw();
	l7->Draw();
	c7.Update();
	c7.SaveAs("aTGC_interference.png");
	end=getchar();


	for(int i = 1; i<8; i++)
	{
		std::cout<<("a"+to_string(i)).c_str()<<": "<<w.var(("a"+to_string(i)).c_str())->getVal()<<std::endl;
		//if(i>1 and i<5)
		//	std::cout<<("a"+to_string(i)+to_string(i)).c_str()<<": "<<w.var(("a"+to_string(i)+to_string(i)).c_str())->getVal()<<std::endl;
	}


	cwww.setVal(12); ccw.setVal(20); cb.setVal(60);
	std::cout<<"all params max:"<<std::endl;
	std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
	std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
	//std::cout<<"N3: "<<N3.getVal()<<", N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;

	cwww.setVal(1.2); ccw.setVal(2); cb.setVal(6);
	std::cout<<"all params 1/10:"<<std::endl;
	std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
	std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
	//std::cout<<"N3: "<<N3.getVal()<<", N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;
	

	exit(0);
}


