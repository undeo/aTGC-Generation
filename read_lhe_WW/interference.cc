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
#include <TSystem.h>


#include <HWWLVJRooPdfs.h>

double cwwws[]		= {-12.,-6.,-2.,0.,2.,6.,12.};
double ccws[]		= {-20.,-10.,-3.5,0.,3.5,10.,20.};
double cbs[]		= {-60.,-30.,-10.,0.,10.,30.,60.};

double normSM;
RooWorkspace w("w","w");
RooWorkspace w2("w2","w2");


void interference(int channel)
{
	gSystem->Load("HWWLVJRooPdfs.cxx");
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
	RooRealVar MWW("MWW","MWW",2500,600,3500);
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

	RooDataHist hist_diff_ccw("hist_diff_ccw","hist_diff_ccw",RooArgSet(MWW));
	RooDataHist hist_diff_cb("hist_diff_cb","hist_diff_cb",RooArgSet(MWW));
	RooDataHist hist_diff_cwww_ccw("hist_diff_cwww_ccw","hist_diff_cwww_ccw",RooArgSet(MWW));
	RooDataHist hist_diff_cwww_cb("hist_diff_cwww_cb","hist_diff_cwww_cb",RooArgSet(MWW));
	RooDataHist hist_diff_ccw_cb("hist_diff_ccw_cb","hist_diff_ccw_cb",RooArgSet(MWW));

	reader.SetEntry(0);
	while(reader.Next())
	{
		MWW.setVal(*MWW_tree);
		if(*channel_tree==channel)
		{
			double weight_cwww_tmp		= (*weights)[112]-(*weights)[13];
			double weight_ccw_tmp		= (*weights)[72]-(*weights)[53];
			hist_diff_ccw.add(RooArgSet(MWW),weight_ccw_tmp);
			double weight_cb_tmp		= (*weights)[64]-(*weights)[61];
			hist_diff_cb.add(RooArgSet(MWW),weight_cb_tmp);
			double weight_cwww_ccw_tmp	= ((*weights)[122]-(*weights)[23]) - weight_cwww_tmp;//([cwww12,ccw20]-[cwww-12,ccw20]) - ([cwww12]-[cwww-12])
			hist_diff_cwww_ccw.add(RooArgSet(MWW),weight_cwww_ccw_tmp);
			double weight_cwww_cb_tmp	= ((*weights)[114]-(*weights)[110]) - weight_cb_tmp;//([cwww12,cb60]-[cwww12,cb-60]) - ([cb60]-[cb-60])
			hist_diff_cwww_cb.add(RooArgSet(MWW),weight_cwww_cb_tmp);
			double weight_ccw_cb_tmp 	= ((*weights)[74]-(*weights)[70]) - weight_cb_tmp;//([ccw20,cb60]-[ccw20,cb-60]) - ([cb60]-[cb-60])
			hist_diff_ccw_cb.add(RooArgSet(MWW),weight_ccw_cb_tmp);
		}
	}
	
	TH1F* testhist = (TH1F*)hist_diff_cwww_ccw.createHistogram("hist_cwww_ccw",MWW);
	TCanvas dummy1("dummy1","dummy1",1);
	dummy1.cd(); dummy1.SetLogy();
	testhist->Draw();
	testhist->Fit("expo");
	dummy1.Draw();
	float slopeval = testhist->GetFunction("expo")->GetParameter(1);


	w.factory("Exponential:SM_Pdf(MWW,a1[-0.001,-0.01,0.])");
	//w.factory("Exponential:Pdf_cwww(MWW,a2[-0.001,-0.01,0.])");
	//w.factory("Exponential:Pdf_ccw(MWW,a3[-0.001,-0.01,0.])");
	//w.factory("Exponential:Pdf_cb(MWW,a4[-0.001,-0.01,0.])");
	//w.factory("Exponential:Pdf_cwww_lin(MWW,a22[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_ccw_lin(MWW,a33[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_cb_lin(MWW,a44[-0.001,-0.01,0.])");
	w.factory("Exponential:Int_cwww_ccw(MWW,a5[-0.0001,-0.01,0.01])");
	w.factory("Exponential:Int_cwww_cb(MWW,a6[-0.001,-0.01,0.01])");
	w.factory("Exponential:Int_ccw_cb(MWW,a7[-0.001,-0.01,0.01])");

	RooRealVar a2("a2","a2",-0.001,-0.01,0.);
	RooRealVar a3("a3","a3",-0.001,-0.01,0.);
	RooRealVar a4("a4","a4",-0.001,-0.01,0.);
	RooRealVar Erf_offset_cwww("Erf_offset_cwww","Erf_offset_cwww",700,200,2000);
	RooRealVar Erf_offset_ccw("Erf_offset_ccw","Erf_offset_ccw",700,200,2000);
	RooRealVar Erf_offset_cb("Erf_offset_cb","Erf_offset_cb",700,200,2000);
	RooRealVar Erf_width_cwww("Erf_width_cwww","Erf_width_cwww",200,10,10000);
	RooRealVar Erf_width_ccw("Erf_width_ccw","Erf_width_ccw",200,10,10000);
	RooRealVar Erf_width_cb("Erf_width_cb","Erf_width_cb",200,10,10000);
	Erf_offset_cwww.setConstant(true);
	Erf_offset_ccw.setConstant(true);
	Erf_offset_cb.setConstant(true);
	Erf_width_cwww.setConstant(true);
	Erf_width_ccw.setConstant(true);
	Erf_width_cb.setConstant(true);

	RooErfExpPdf Pdf_cwww("Pdf_cwww","Pdf_cwww",MWW,a2,Erf_offset_cwww,Erf_width_cwww);
	RooErfExpPdf Pdf_ccw("Pdf_ccw","Pdf_ccw",MWW,a3,Erf_offset_ccw,Erf_width_ccw);
	RooErfExpPdf Pdf_cb("Pdf_cb","Pdf_cb",MWW,a4,Erf_offset_cb,Erf_width_cb);

	RooAbsPdf * SM_Pdf 	= w.pdf("SM_Pdf");
	//RooAbsPdf * Pdf_cwww	= w.pdf("Pdf_cwww");
	//RooAbsPdf * Pdf_ccw	= w.pdf("Pdf_ccw");
	//RooAbsPdf * Pdf_cb	= w.pdf("Pdf_cb");
	//RooAbsPdf * Pdf_cwww_lin= w.pdf("Pdf_cwww_lin");
	RooAbsPdf * Pdf_ccw_lin	= w.pdf("Pdf_ccw_lin");
	RooAbsPdf * Pdf_cb_lin	= w.pdf("Pdf_cb_lin");
	RooAbsPdf * Int_cwww_ccw= w.pdf("Int_cwww_ccw");
	RooAbsPdf * Int_cwww_cb	= w.pdf("Int_cwww_cb");
	RooAbsPdf * Int_ccw_cb	= w.pdf("Int_ccw_cb");

	RooRealVar N_SM4fit("N_SM4fit","N_SM4fit",w.data("hist0")->sumEntries());			//hist0
	RooRealVar N_SM("N_SM","N_SM",w.data("hist0")->sumEntries("MWW>900"));
	RooRealVar N__124fit("N_cwww__124fit","N_cwww__124fit",w.data("hist13")->sumEntries());	//hist13, hist128
	RooRealVar N__12("N_cwww__12","N_cwww__12",w.data("hist13")->sumEntries("MWW>900"));
	RooRealVar N_124fit("N_cwww_124fit","N_cwww_124fit",w.data("hist112")->sumEntries());	//hist112, hist145
	RooRealVar N_12("N_cwww_12","N_cwww_12",w.data("hist112")->sumEntries("MWW>900"));
	RooRealVar N__204fit("N_ccw__204fit","N_ccw__204fit",w.data("hist53")->sumEntries());	//hist53, hist134
	RooRealVar N__20("N_ccw__20","N_ccw__20",w.data("hist53")->sumEntries("MWW>900"));
	RooRealVar N_204fit("N_ccw_204fit","N_ccw_204fit",w.data("hist72")->sumEntries());		//hist72, hist139
	RooRealVar N_20("N_ccw_20","N_ccw_20",w.data("hist72")->sumEntries("MWW>900"));
	RooRealVar N__604fit("N_cb__604fit","N_cb__604fit",w.data("hist61")->sumEntries());		//hist61, hist136
	RooRealVar N__60("N_cb__60","N_cb__60",w.data("hist61")->sumEntries("MWW>900"));
	RooRealVar N_604fit("N_cb_604fit","N_cb_604fit",w.data("hist64")->sumEntries());		//hist64, hist137
	RooRealVar N_60("N_cb_60","N_cb_60",w.data("hist64")->sumEntries("MWW>900"));
	//RooRealVar N__12_20("N_cwww_ccw__12_20","N_cwww_ccw__12_20",w.data("hist23")->sumEntries());
	//RooRealVar N__12_60("N_cwww_cb__12_60","N_cwww_cb__12_60",w.data("hist15")->sumEntries());
	//RooRealVar N__20_60("N_ccw_cb__20_60","N_ccw_cb__20_60",w.data("hist55")->sumEntries());
	RooRealVar N_12_204fit("N_cwww_ccw__12__204fit","N_cwww_ccw__12__204fit",w.data("hist122")->sumEntries());
	RooRealVar N_12_20("N_cwww_ccw__12__20","N_cwww_ccw__12__20",w.data("hist122")->sumEntries("MWW>900"));
	RooRealVar N_12_604fit("N_cwww_cb__12__604fit","N_cwww_cb__12__604fit",w.data("hist114")->sumEntries());
	RooRealVar N_12_60("N_cwww_cb__12__60","N_cwww_cb__12__60",w.data("hist114")->sumEntries("MWW>900"));
	RooRealVar N_20_604fit("N_ccw_cb__20__604fit","N_ccw_cb__20__604fit",w.data("hist74")->sumEntries());
	RooRealVar N_20_60("N_ccw_cb__20__60","N_ccw_cb__20__60",w.data("hist74")->sumEntries("MWW>900"));
	RooRealVar N_4norm4fit("N_4norm4fit","N_4norm4fit",w.data("hist1")->sumEntries());
	RooRealVar N_4norm("N_4norm","N_4norm",w.data("hist1")->sumEntries("MWW>900"));


	w2.import(N_4norm4fit);
	w2.import(N_SM4fit);
	w2.import(N__124fit);
	w2.import(N_124fit);
	w2.import(N__204fit);
	w2.import(N_204fit);
	w2.import(N__604fit);
	w2.import(N_604fit);
	w2.import(N_12_204fit);
	w2.import(N_12_604fit);
	w2.import(N_20_604fit);

	w2.import(N_4norm);
	w2.import(N_SM);
	w2.import(N__12);
	w2.import(N_12);
	w2.import(N__20);
	w2.import(N_20);
	w2.import(N__60);
	w2.import(N_60);
	w2.import(N_12_20);
	w2.import(N_12_60);
	w2.import(N_20_60);

	std::cout<<N_SM4fit.GetName()<<": "<<N_SM4fit.getVal()<<std::endl;
	std::cout<<N__124fit.GetName()<<": "<<N__124fit.getVal()<<std::endl;
	std::cout<<N_124fit.GetName()<<": "<<N_124fit.getVal()<<std::endl;
	std::cout<<N__204fit.GetName()<<": "<<N__204fit.getVal()<<std::endl;
	std::cout<<N_204fit.GetName()<<": "<<N_204fit.getVal()<<std::endl;
	std::cout<<N__604fit.GetName()<<": "<<N__604fit.getVal()<<std::endl;
	std::cout<<N_604fit.GetName()<<": "<<N_604fit.getVal()<<std::endl;
	std::cout<<N_12_204fit.GetName()<<": "<<N_12_204fit.getVal()<<std::endl;
	std::cout<<N_12_604fit.GetName()<<": "<<N_12_604fit.getVal()<<std::endl;
	std::cout<<N_20_604fit.GetName()<<": "<<N_20_604fit.getVal()<<std::endl;

	RooRealVar N2_tmp("N2_tmp","N2_tmp",(N_124fit.getVal()+N__124fit.getVal())/2-N_SM4fit.getVal());
	RooRealVar N4_tmp("N4_tmp","N4_tmp",(N_204fit.getVal()+N__204fit.getVal())/2-N_SM4fit.getVal());
	RooRealVar N5_tmp("N5_tmp","N5_tmp",(N_204fit.getVal()-N__204fit.getVal())/2);
	RooRealVar N6_tmp("N6_tmp","N6_tmp",(N_604fit.getVal()+N__604fit.getVal())/2-N_SM4fit.getVal());
	RooRealVar N7_tmp("N7_tmp","N7_tmp",(N_604fit.getVal()-N__604fit.getVal())/2);
	RooRealVar N8_tmp("N8_tmp","N8_tmp",(N_12_204fit.getVal()+N_SM4fit.getVal())-(N_124fit.getVal()+N_204fit.getVal()));
	RooRealVar N9_tmp("N9_tmp","N9_tmp",(N_12_604fit.getVal()+N_SM4fit.getVal())-(N_124fit.getVal()+N_604fit.getVal()));
	RooRealVar N10_tmp("N10_tmp","N10_tmp",(N_20_604fit.getVal()+N_SM4fit.getVal())-(N_204fit.getVal()+N_604fit.getVal()));

	RooArgList N_norm_list(N_SM4fit,N2_tmp,N4_tmp,N5_tmp,N6_tmp,N7_tmp);
	N_norm_list.add(RooArgList(cwww,ccw,cb));
	N_norm_list.add(RooArgList(N8_tmp,N9_tmp,N10_tmp));

	TString cwww_f 		= "+@1*(@6/12)**2";
	TString ccw_f 		= "+@2*(@7/20)**2";
	TString cb_f		= "+@4*(@8/60)**2";
	//TString cwww_lin_f	= "+((@1-@2)/2)*(@7/12)";
	TString ccw_lin_f 	= "+@3*(@7/20)";
	TString cb_lin_f	= "+@5*(@8/60)";
	TString cwww_ccw_f	= "+@9*(@6/12)*(@7/20)";
	TString cwww_cb_f	= "+@10*(@6/12)*(@8/60)";
	TString ccw_cb_f	= "+@11*(@7/20)*(@8/60)";
	RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+ccw_lin_f+cb_f+cb_lin_f+cwww_ccw_f+cwww_cb_f+ccw_cb_f,N_norm_list);

	RooFormulaVar N1("N1","@1/@0",RooArgList(Norm,N_SM4fit));
	RooFormulaVar N2("N2","(@2*(@1/12)**2)/@0",RooArgList(Norm,cwww,N2_tmp));
	//RooFormulaVar N3("N3","(((@1-@2)/2)*(@3/12))/@0",RooArgList(Norm,N_12,N__12,N_SM,cwww));
	RooFormulaVar N4("N4","(@2*(@1/20)**2)/@0",RooArgList(Norm,ccw,N4_tmp));
	RooFormulaVar N5("N5","(@2*(@1/20))/@0",RooArgList(Norm,ccw,N5_tmp));
	RooFormulaVar N6("N6","(@2*(@1/60)**2)/@0",RooArgList(Norm,cb,N6_tmp));
	RooFormulaVar N7("N7","(@2*(@1/60))/@0",RooArgList(Norm,cb,N7_tmp));
	RooFormulaVar N8("N8","(@3*(@1/12)*(@2/20))/@0",RooArgList(Norm,cwww,ccw,N8_tmp));
	RooFormulaVar N9("N9","(@3*(@1/12)*(@2/60))/@0",RooArgList(Norm,cwww,cb,N9_tmp));
	RooFormulaVar N10("N10","(@3*(@1/20)*(@2/60))/@0",RooArgList(Norm,ccw,cb,N10_tmp));



	//RooArgList N_list(N1,N2,N4,N6);
	RooArgList N_list(N1,N2,N4,N5,N6,N7);
	N_list.add(RooArgList(N8,N9,N10));
	//RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_ccw,*Pdf_cb);
	RooArgList Pdf_list(*SM_Pdf,Pdf_cwww,Pdf_ccw,*Pdf_ccw_lin,Pdf_cb,*Pdf_cb_lin);
	Pdf_list.add(RooArgList(*Int_cwww_ccw,*Int_cwww_cb,*Int_ccw_cb));
	Pdf_list.Print();
	N_list.Print();
	RooAddPdf model1("model","model",Pdf_list,N_list);


	RooRealVar * a1 	= w.var("a1");
	a2.setConstant(true);
	//RooRealVar * a22	= w.var("a22"); a22->setConstant(true);
	a3.setConstant(true);
	RooRealVar * a33	= w.var("a33");	a33->setConstant(true);
	a4.setConstant(true);
	RooRealVar * a44	= w.var("a44");	a44->setConstant(true);
	RooRealVar * a5		= w.var("a5");	a5->setConstant(true);
	RooRealVar * a6		= w.var("a6");	a6->setConstant(true);
	RooRealVar * a7		= w.var("a7");	a7->setConstant(true);

	RooBinning bins(600,3500);
	bins.addUniform(25,600,3500);
//SM-fit
	cwww.setVal(0);	ccw.setVal(0); cb.setVal(0); 	
	model1.fitTo(*w.data("hist0"));//hist0
	a1->setConstant(true);
//SM-interference-fits
	double N_SM_tmp_val = N_SM4fit.getVal();//SM
	double N2_tmp_val = N2_tmp.getVal();//cwww quad
	double N4_tmp_val = N4_tmp.getVal();//ccw quad
	double N6_tmp_val = N6_tmp.getVal();//cb quad
	double N5_tmp_val = N5_tmp.getVal();//ccw lin
	double N7_tmp_val = N7_tmp.getVal();//cb lin
	N_SM4fit.setVal(0);N2_tmp.setVal(0);N4_tmp.setVal(0);N6_tmp.setVal(0);	
	cwww.setVal(0);ccw.setVal(20); cb.setVal(0);
	a33->setConstant(false);
	model1.fitTo(hist_diff_ccw);
	a33->setConstant(true);
	cwww.setVal(0);ccw.setVal(0); cb.setVal(60);
	a44->setConstant(false);
	model1.fitTo(hist_diff_cb);
	a44->setConstant(true);

	TCanvas cc2("smint","smint",1);
	cc2.cd();cc2.SetLogy();
	RooPlot * plott2 = MWW.frame();
	cwww.setVal(0);ccw.setVal(20); cb.setVal(0);
	hist_diff_ccw.plotOn(plott2,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E"));
	model1.plotOn(plott2,RooFit::LineColor(kBlue),RooFit::Normalization(hist_diff_ccw.sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0);ccw.setVal(0); cb.setVal(60);
	hist_diff_cb.plotOn(plott2,RooFit::Binning(bins),RooFit::LineColor(kGreen),RooFit::DrawOption("E"));
	model1.plotOn(plott2,RooFit::LineColor(kGreen),RooFit::Normalization(hist_diff_cb.sumEntries(),RooAbsReal::NumEvent));
	plott2->GetYaxis()->SetRangeUser(0.001,5);
	plott2->Draw();
	cc2.Draw();cc2.Update();

	//int cwww-ccw-fit
	N_SM4fit.setVal(0);N2_tmp.setVal(0);N4_tmp.setVal(0);N5_tmp.setVal(0);N6_tmp.setVal(0);N7_tmp.setVal(0);
	cwww.setVal(12); ccw.setVal(20); cb.setVal(0);
	//a5->setConstant(false);
	//model1.fitTo(hist_diff_cwww_ccw);
	a5->setVal(slopeval);
	a5->setConstant(true);
/*	//int cwww-cb-fit
	cwww.setVal(12); ccw.setVal(0); cb.setVal(60);
	a6->setConstant(false);
	model1.fitTo(hist_diff_cwww_cb);
	a6->setConstant(true);*/
	//int ccw-cb-fit
	cwww.setVal(0); ccw.setVal(20); cb.setVal(60);//cb lin
	a7->setConstant(false);
	model1.fitTo(hist_diff_ccw_cb);
	a7->setConstant(true);
	int end2;

	TCanvas cc("atgcint","atgcint",1);
	cc.cd();cc.SetLogy();
	RooPlot * plott = MWW.frame();
	cwww.setVal(12);ccw.setVal(20);cb.setVal(0);
	hist_diff_cwww_ccw.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E"));
	model1.plotOn(plott,RooFit::LineColor(kBlue),RooFit::Normalization(hist_diff_cwww_ccw.sumEntries(),RooAbsReal::NumEvent));
	//cwww.setVal(12); ccw.setVal(0); cb.setVal(60);
	//hist_diff_cwww_cb.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kRed),RooFit::DrawOption("E"));
	//model1.plotOn(plott,RooFit::LineColor(kRed),RooFit::Normalization(hist_diff_cwww_cb.sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(20); cb.setVal(60);
	hist_diff_ccw_cb.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kGreen),RooFit::DrawOption("E"));
	model1.plotOn(plott,RooFit::LineColor(kGreen),RooFit::Normalization(hist_diff_ccw_cb.sumEntries(),RooAbsReal::NumEvent));
	plott->GetYaxis()->SetRangeUser(0.0001,1);
	plott->Draw();
	cc.Draw();
	cc.Update();

	N_SM4fit.setVal(N_SM_tmp_val);
	N2_tmp.setVal(N2_tmp_val);
	N4_tmp.setVal(N4_tmp_val);
	N5_tmp.setVal(N5_tmp_val);
	N6_tmp.setVal(N6_tmp_val);
	N7_tmp.setVal(N7_tmp_val);
//cwww-fit
	cwww.setVal(12); ccw.setVal(0); cb.setVal(0); 
	a2.setConstant(false);
	Erf_offset_cwww.setConstant(false);
	Erf_width_cwww.setConstant(false);
	model1.fitTo(*w.data("hist112"));//hist13, hist128
	a2.setConstant(true);
	Erf_offset_cwww.setConstant(true);
	Erf_width_cwww.setConstant(true);
//ccw-fit
	cwww.setVal(0);	ccw.setVal(20); cb.setVal(0); 
	a3.setConstant(false);
	Erf_offset_ccw.setConstant(false);
	Erf_width_ccw.setConstant(false);
	model1.fitTo(*w.data("hist72"));//hist53, hist134
	a3.setConstant(true);
	Erf_offset_ccw.setConstant(true);
	Erf_width_ccw.setConstant(true);
//cb-fit
	cwww.setVal(0); ccw.setVal(0); cb.setVal(60);
	a4.setConstant(false);
	Erf_offset_cb.setConstant(false);
	Erf_width_cb.setConstant(false);
	model1.fitTo(*w.data("hist64"));//hist61, hist 136
	a4.setConstant(true);
	Erf_offset_cb.setConstant(true);
	Erf_width_cb.setConstant(true);

	w2.import(*a1);
	w2.import(a2);
	w2.import(a3);
	w2.import(*a33);
	w2.import(a4);
	w2.import(*a44);
	w2.import(*a5);
	w2.import(*a6);
	w2.import(*a7);


	TFile * fileOut = new TFile("genlevel_WW_"+ch+".root","RECREATE");
	w2.Write();
	fileOut->Close();



	TCanvas c1("fits","fits0",1);
	c1.cd(); c1.SetLogy();
	RooPlot * plot1 = MWW.frame();
	cwww.setVal(12); ccw.setVal(20); cb.setVal(0);
	w.data("hist122")->plotOn(plot1,RooFit::Binning(bins),RooFit::MarkerColor(1));
	model1.plotOn(plot1,RooFit::LineColor(1),RooFit::Normalization(w.data("hist122")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(12); ccw.setVal(0); cb.setVal(60);
	w.data("hist114")->plotOn(plot1,RooFit::Binning(bins),RooFit::MarkerColor(2));
	model1.plotOn(plot1,RooFit::LineColor(2),RooFit::Normalization(w.data("hist114")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(20); cb.setVal(60);
	w.data("hist74")->plotOn(plot1,RooFit::Binning(bins),RooFit::MarkerColor(3));
	model1.plotOn(plot1,RooFit::LineColor(3),RooFit::Normalization(w.data("hist74")->sumEntries(),RooAbsReal::NumEvent));

	plot1->Draw();
	c1.Draw();
	c1.Update();

	int end;

/*
	TCanvas c2("SM","SM",1);
	c2.cd(); c2.SetLogy();
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	RooPlot * plot2 = MWW.frame();
	w.data("hist0")->plotOn(plot2,RooFit::Binning(bins));
	model1.plotOn(plot2);
	plot2->Draw();
	c2.Draw();
	c2.Update();
*/
/*
	TCanvas c6("allneg","allneg",1);
	c6.cd(); c6.SetLogy();
	RooPlot * plot6 = MWW.frame();
	cwww.setVal(-12); ccw.setVal(-20); cb.setVal(-60);
	w.data("hist1")->plotOn(plot6,RooFit::Binning(bins),RooFit::MarkerColor(kBlue));
	model1.plotOn(plot6);
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	w.data("hist0")->plotOn(plot6,RooFit::Binning(bins),RooFit::MarkerColor(kBlack));
	model1.plotOn(plot6);
	plot6->Draw();
	c6.Draw();
	c6.Update();

	TCanvas c5("all","all",1);
	c5.cd(); c5.SetLogy();
	RooPlot * plot5 = MWW.frame();
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	w.data("hist0")->plotOn(plot5,RooFit::Binning(bins),RooFit::MarkerColor(kBlue));
	model1.plotOn(plot5);
	cwww.setVal(2); ccw.setVal(3.5); cb.setVal(10);
	w.data("hist149")->plotOn(plot5,RooFit::Binning(bins),RooFit::MarkerColor(kRed));
	model1.plotOn(plot5,RooFit::LineColor(kRed),RooFit::Normalization(w.data("hist149")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(-6); ccw.setVal(0); cb.setVal(30);
	w.data("hist39")->plotOn(plot5,RooFit::Binning(bins),RooFit::MarkerColor(kMagenta));
	model1.plotOn(plot5,RooFit::LineColor(kMagenta),RooFit::Normalization(w.data("hist39")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(-6); ccw.setVal(20); cb.setVal(-60);
	w.data("hist46")->plotOn(plot5,RooFit::Binning(bins),RooFit::MarkerColor(kCyan));
	model1.plotOn(plot5,RooFit::LineColor(kCyan),RooFit::Normalization(w.data("hist46")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(12); ccw.setVal(-10); cb.setVal(-30);
	w.data("hist106")->plotOn(plot5,RooFit::Binning(bins),RooFit::MarkerColor(kGreen));
	model1.plotOn(plot5,RooFit::LineColor(kGreen),RooFit::Normalization(w.data("hist106")->sumEntries(),RooAbsReal::NumEvent));
	plot5->Draw();
	c5.Draw();
	c5.Update();
	c5.SaveAs("interference_plot.root","RECREATE");

	TCanvas c7("cwww","cwww",1);
	c7.cd(); c7.SetLogy();
	RooPlot * plot7 = MWW.frame();
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	w.data("hist0")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(kBlack),RooFit::DrawOption("E1"));
	model1.plotOn(plot7,RooFit::LineColor(kBlack),RooFit::LineWidth(1));
	cwww.setVal(0);	ccw.setVal(10); cb.setVal(-30);
	w.data("hist66")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E1"));
	model1.plotOn(plot7,RooFit::LineColor(kBlue),RooFit::LineWidth(1),RooFit::Normalization(w.data("hist66")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(1);	ccw.setVal(2); cb.setVal(-6);
	w.data("hist147")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(kGreen),RooFit::DrawOption("E1"));
	model1.plotOn(plot7,RooFit::LineColor(kGreen),RooFit::LineWidth(1),RooFit::Normalization(w.data("hist147")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(6);	ccw.setVal(10); cb.setVal(-30);
	w.data("hist91")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(kRed),RooFit::DrawOption("E1"));
	model1.plotOn(plot7,RooFit::LineColor(kRed),RooFit::LineWidth(1),RooFit::Normalization(w.data("hist91")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(12); ccw.setVal(10); cb.setVal(-30);
	w.data("hist116")->plotOn(plot7,RooFit::Binning(bins),RooFit::LineColor(kCyan),RooFit::DrawOption("E1"));
	model1.plotOn(plot7,RooFit::LineColor(kCyan),RooFit::LineWidth(1),RooFit::Normalization(w.data("hist116")->sumEntries(),RooAbsReal::NumEvent));
	plot7->Draw();
	c7.Draw();
	c7.Update();
*/
	//N8_tmp.setVal(0);
	//N9_tmp.setVal(0);
	//N10_tmp.setVal(0);
	TCanvas c4plot("c4plot","c4plot",1);
	c4plot.cd();
	c4plot.SetLogy();
	RooPlot * plot4plot = MWW.frame();
	cwww.setVal(0); ccw.setVal(0); cb.setVal(0);
	w.data("hist0")->plotOn(plot4plot,RooFit::Binning(bins),RooFit::LineColor(kBlack),RooFit::DrawOption("E1"));
	model1.plotOn(plot4plot,RooFit::LineColor(kBlack));
	cwww.setVal(0); ccw.setVal(-20); cb.setVal(0);
	w.data("hist53")->plotOn(plot4plot,RooFit::Binning(bins),RooFit::LineColor(kRed),RooFit::DrawOption("E1"));
	model1.plotOn(plot4plot,RooFit::LineColor(kRed),RooFit::Normalization(w.data("hist53")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(-20); cb.setVal(30);
	w.data("hist54")->plotOn(plot4plot,RooFit::Binning(bins),RooFit::LineColor(kGreen+1),RooFit::DrawOption("E1"));
	model1.plotOn(plot4plot,RooFit::LineColor(kGreen+1),RooFit::Normalization(w.data("hist54")->sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(-20); cb.setVal(60);
	w.data("hist55")->plotOn(plot4plot,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E1"));
	model1.plotOn(plot4plot,RooFit::LineColor(kBlue),RooFit::Normalization(w.data("hist55")->sumEntries(),RooAbsReal::NumEvent));
	plot4plot->GetXaxis()->SetRangeUser(900,3500);
	plot4plot->GetYaxis()->SetTitle("arbitrary units");
	plot4plot->GetXaxis()->SetTitle("M_{WW} (GeV)");
	plot4plot->SetTitle("");
	plot4plot->Draw();
	leg = new TLegend(0.6,0.6,0.89,0.89);
	leg->SetHeader("#frac{c_{WWW}}{#Lambda^{2}}=0 TeV^{-2}, #frac{c_{W}}{#Lambda^{2}}=-20 TeV^{-2}");
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->AddEntry(plot4plot->getObject(0),"SM","le");
	leg->AddEntry(plot4plot->getObject(2),"#frac{c_{B}}{#Lambda^{2}}=0 TeV^{-2}","le");
	leg->AddEntry(plot4plot->getObject(4),"#frac{c_{B}}{#Lambda^{2}}=30 TeV^{-2}","le");
	leg->AddEntry(plot4plot->getObject(6),"#frac{c_{B}}{#Lambda^{2}}=60 TeV^{-2}","le");
	leg->Draw("SAME");
	c4plot.Update();
	c4plot.SaveAs("atgcinterference_plot4talk.png");

	for(int i = 1; i<8; i++)
	{
		std::cout<<("a"+to_string(i)).c_str()<<": "<<w2.var(("a"+to_string(i)).c_str())->getVal()<<std::endl;
		if(i>2 and i<5)
			std::cout<<("a"+to_string(i)+to_string(i)).c_str()<<": "<<w2.var(("a"+to_string(i)+to_string(i)).c_str())->getVal()<<std::endl;
	}

	std::cout<<"offsets:"<<std::endl;
	std::cout<<"cwww: "<<Erf_offset_cwww.getVal()<<", ccw: "<<Erf_offset_ccw.getVal()<<", cb: "<<Erf_offset_cb.getVal()<<std::endl;
	std::cout<<"widths:"<<std::endl;
	std::cout<<"cwww: "<<Erf_width_cwww.getVal()<<", ccw: "<<Erf_width_ccw.getVal()<<", cb: "<<Erf_width_cb.getVal()<<std::endl;


	cwww.setVal(12); ccw.setVal(20); cb.setVal(60);
	std::cout<<"all params max:"<<std::endl;
	std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
	std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
	std::cout<<", N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;

	cwww.setVal(1.2); ccw.setVal(2); cb.setVal(6);
	std::cout<<"all params 1/10:"<<std::endl;
	std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
	std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
	std::cout<<", N5: "<<N5.getVal()<<", N7:"<<N7.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N9: "<<N9.getVal()<<", N10:"<<N10.getVal()<<std::endl;


	std::cout<<"cwww/ccw: "<<hist_diff_cwww_ccw.sumEntries()<<std::endl;
	std::cout<<"cwww/cb: "<<hist_diff_cwww_cb.sumEntries()<<std::endl;
	std::cout<<"ccw/cb: "<<hist_diff_ccw_cb.sumEntries()<<std::endl;	

	std::string raw_input;
	std::cin>>raw_input;
	exit(0);
}


