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

//usage: root -l, .L HWWLVJRooPdfs.cxx+, .L interference.cc, interference(1/2)


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


double normSM;
RooWorkspace w("w","w");
RooWorkspace w2("w2","w2");


void interference(int channel)
{
	//vals gives the atgc values, SM for [0][i]
	vals[0][0]	= 0.;
	vals[0][1]	= 0.;
	vals[0][2]	= 0.;
	//set other atgc values of vals
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
	vector<TString> histonames;
	for(unsigned int i = 0; i<150; i++)
	{
		TString cwww_tmp	= vals[i][0] >= 0 ? ("cwww"+to_string(int(vals[i][0]))).c_str() : ("cwww_"+to_string(int(abs(vals[i][0])))).c_str();
		TString ccw_tmp		= vals[i][1] >= 0 ? ("ccw"+to_string(int(vals[i][1]))).c_str() : ("ccw_"+to_string(int(abs(vals[i][1])))).c_str();
		TString cb_tmp		= vals[i][2] >= 0 ? ("cb"+to_string(int(vals[i][2]))).c_str() : ("cb_"+to_string(int(abs(vals[i][2])))).c_str();
		if(vals[i][1] - int(vals[i][1]) != 0)
			//if it looks stupid but it works, it ain't stupid
			ccw_tmp		= vals[i][1] >= 0 ? ("ccw"+to_string(int(vals[i][1]))+"_"+to_string(abs(int(10*(vals[i][1]-int(vals[i][1])))))).c_str() : ("ccw_"+to_string(int(abs(vals[i][1])))+"_"+to_string(abs(int(10*(vals[i][1]-int(vals[i][1])))))).c_str();
		TString name_tmp	= cwww_tmp + ccw_tmp + cb_tmp;
		histonames.push_back(name_tmp);
	}
	


	//gSystem->Load("HWWLVJRooPdfs.cxx");
	//channel: 1=el, 2=mu
	TString ch;	
	if(channel==1)
		ch = "el";
	if(channel==2)
		ch = "mu";
	//nbins	= (3500-900)/100;
	nbins	= 30;

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
		RooDataHist hist(histonames[i],histonames[i],RooArgSet(MWW));
		w.import(hist);
	}
	
	RooDataHist hist_diff_ccw("hist_diff_ccw","hist_diff_ccw",RooArgSet(MWW));
	RooDataHist hist_diff_cb("hist_diff_cb","hist_diff_cb",RooArgSet(MWW));
	RooDataHist hist_diff_cwww_ccw("hist_diff_cwww_ccw","hist_diff_cwww_ccw",RooArgSet(MWW));
	RooDataHist hist_diff_cwww_cb("hist_diff_cwww_cb","hist_diff_cwww_cb",RooArgSet(MWW));
	RooDataHist hist_diff_ccw_cb("hist_diff_ccw_cb","hist_diff_ccw_cb",RooArgSet(MWW));


	int tmp = 0;	
	while(reader.Next())
	{
		MWW.setVal(*MWW_tree);
		if(*channel_tree==channel)
			for (int i = 0; i<150; i++)
			{
				double weight_tmp = (*weights)[i];
				w.data(histonames[i])->add(RooArgSet(MWW),weight_tmp);
			}
			double weight_cwww_tmp		= (*weights)[112]-(*weights)[13];
			double weight_ccw_tmp		= (*weights)[72]-(*weights)[53];
			double weight_cb_tmp		= (*weights)[64]-(*weights)[61];
			double weight_cwww_ccw_tmp	= ((*weights)[122]-(*weights)[23]) - weight_cwww_tmp;//([cwww12,ccw20]-[cwww-12,ccw20]) - ([cwww12]-[cwww-12])
			double weight_cwww_cb_tmp	= (((*weights)[114]-(*weights)[110]) - weight_cb_tmp);//([cwww12,cb60]-[cwww12,cb-60]) - ([cb60]-[cb-60])
			double weight_ccw_cb_tmp 	= (((*weights)[74]-(*weights)[70]) - weight_cb_tmp);//([ccw20,cb60]-[ccw20,cb-60]) - ([cb60]-[cb-60])
			hist_diff_ccw.add(RooArgSet(MWW),weight_ccw_tmp);
			hist_diff_cb.add(RooArgSet(MWW),weight_cb_tmp);
			hist_diff_cwww_ccw.add(RooArgSet(MWW),weight_cwww_ccw_tmp);
			hist_diff_cwww_cb.add(RooArgSet(MWW),weight_cwww_cb_tmp);
			hist_diff_ccw_cb.add(RooArgSet(MWW),weight_ccw_cb_tmp);
		tmp++;
		if(tmp%5000==0)
			std::cout<<tmp<<std::endl;
	}


	int end;	
	//binning added after pre-approval-------------------------------------------------------------------------------------------------------------------------
	RooBinning bins2(900,3500);
	bins2.addUniform(26,900,3500);
	TH1F* testhist = (TH1F*)hist_diff_cwww_ccw.createHistogram("hist_cwww_ccw",MWW,RooFit::Binning(bins2));
	TH1F* testhist2 = (TH1F*)hist_diff_ccw_cb.createHistogram("hist_ccw_cb",MWW,RooFit::Binning(bins2));
	TCanvas dummy("dummy","dummy",1);
	TCanvas dummy2("dummy2","dummy2",1);
	dummy.cd();
	TF1 negexpo("negexpo","-exp([0]+[1]*x)",600,3500);
	testhist->Fit("expo");
	testhist2->Fit("negexpo");
	float slopeval = testhist->GetFunction("expo")->GetParameter(1);
	float slopeval2 = testhist2->GetFunction("negexpo")->GetParameter(1);
	testhist->Draw();
	dummy.Draw();
	dummy.Update();
	dummy2.cd();
	testhist2->Draw();
	dummy2.Draw();
	dummy2.Update();
	//end=getchar();
	//dummy.Close();


	w.factory("Exponential:SM_Pdf(MWW,a1[-0.001,-0.01,0.])");
	w.factory("Exponential:Pdf_ccw_lin(MWW,a33[-0.001,-0.01,0.])");
	//w.factory("Exponential:Pdf_cb_lin(MWW,a44[-0.001,-0.01,0.])");
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
	RooAbsPdf * Pdf_ccw_lin	= w.pdf("Pdf_ccw_lin");
	//RooAbsPdf * Pdf_cb_lin	= w.pdf("Pdf_cb_lin");
	RooAbsPdf * Int_cwww_ccw= w.pdf("Int_cwww_ccw");
	RooAbsPdf * Int_cwww_cb	= w.pdf("Int_cwww_cb");
	RooAbsPdf * Int_ccw_cb	= w.pdf("Int_ccw_cb");

	RooRealVar N_SM4fit("N_SM4fit","N_SM4fit",w.data("cwww0ccw0cb0")->sumEntries());			//hist0
	RooRealVar N__124fit("N_cwww__124fit","N_cwww__124fit",w.data("cwww_12ccw0cb0")->sumEntries());	//hist13, hist128
	RooRealVar N_124fit("N_cwww_124fit","N_cwww_124fit",w.data("cwww12ccw0cb0")->sumEntries());	//hist112, hist145
	RooRealVar N__204fit("N_ccw__204fit","N_ccw__204fit",w.data("cwww0ccw_20cb0")->sumEntries());	//hist53, hist134
	RooRealVar N_204fit("N_ccw_204fit","N_ccw_204fit",w.data("cwww0ccw20cb0")->sumEntries());		//hist72, hist139
	RooRealVar N__604fit("N_cb__604fit","N_cb__604fit",w.data("cwww0ccw0cb_60")->sumEntries());		//hist61, hist136
	RooRealVar N_604fit("N_cb_604fit","N_cb_604fit",w.data("cwww0ccw0cb60")->sumEntries());		//hist64, hist137
	RooRealVar N_12_204fit("N_cwww_ccw_12_204fit","N_cwww_ccw_12_204fit",w.data("cwww12ccw20cb0")->sumEntries());
	RooRealVar N_12_604fit("N_cwww_cb_12_604fit","N_cwww_cb_12_604fit",w.data("cwww12ccw0cb60")->sumEntries());
	RooRealVar N_20_604fit("N_ccw_cb_20_604fit","N_ccw_cb_20_604fit",w.data("cwww0ccw20cb60")->sumEntries());
	RooRealVar N_4norm4fit("N_4norm4fit","N_4norm4fit",w.data("cwww_12ccw_20cb_60")->sumEntries());

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

	//add for later, range 900-3500
	RooRealVar N_SM("N_SM","N_SM",w.data("cwww0ccw0cb0")->sumEntries("MWW>900"));			//hist0
	RooRealVar N__12("N_cwww__12","N_cwww__12",w.data("cwww_12ccw0cb0")->sumEntries("MWW>900"));	//hist13, hist128
	RooRealVar N_12("N_cwww_12","N_cwww_12",w.data("cwww12ccw0cb0")->sumEntries("MWW>900"));	//hist112, hist145
	RooRealVar N__20("N_ccw__20","N_ccw__20",w.data("cwww0ccw_20cb0")->sumEntries("MWW>900"));	//hist53, hist134
	RooRealVar N_20("N_ccw_20","N_ccw_20",w.data("cwww0ccw20cb0")->sumEntries("MWW>900"));		//hist72, hist139
	RooRealVar N__60("N_cb__60","N_cb__60",w.data("cwww0ccw0cb_60")->sumEntries("MWW>900"));		//hist61, hist136
	RooRealVar N_60("N_cb_60","N_cb_60",w.data("cwww0ccw0cb60")->sumEntries("MWW>900"));		//hist64, hist137
	RooRealVar N_12_20("N_cwww_ccw_12_20","N_cwww_ccw_12_20",w.data("cwww12ccw20cb0")->sumEntries("MWW>900"));
	RooRealVar N_12_60("N_cwww_cb_12_60","N_cwww_cb_12_60",w.data("cwww12ccw0cb60")->sumEntries("MWW>900"));
	RooRealVar N_20_60("N_ccw_cb_20_60","N_ccw_cb_20_60",w.data("cwww0ccw20cb60")->sumEntries("MWW>900"));
	RooRealVar N_4norm("N_4norm","N_4norm",w.data("cwww_12ccw_20cb_60")->sumEntries("MWW>900"));

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
	//RooRealVar N7_tmp("N7_tmp","N7_tmp",(N_60.getVal()-N__60.getVal())/2);
	RooRealVar N8_tmp("N8_tmp","N8_tmp",(N_12_204fit.getVal()+N_SM4fit.getVal())-(N_124fit.getVal()+N_204fit.getVal()));
	RooRealVar N9_tmp("N9_tmp","N9_tmp",(N_12_604fit.getVal()+N_SM4fit.getVal())-(N_124fit.getVal()+N_604fit.getVal()));
	RooRealVar N10_tmp("N10_tmp","N10_tmp",(N_20_604fit.getVal()+N_SM4fit.getVal())-(N_204fit.getVal()+N_604fit.getVal()));

	RooArgList N_norm_list(N_SM4fit,N2_tmp,N4_tmp,N5_tmp,N6_tmp);
	N_norm_list.add(RooArgList(cwww,ccw,cb));
	N_norm_list.add(RooArgList(N8_tmp,N9_tmp,N10_tmp));

	TString cwww_f 		= "+@1*(@5/12)**2";
	TString ccw_f 		= "+@2*(@6/20)**2";
	TString cb_f		= "+@4*(@7/60)**2";
	TString ccw_lin_f 	= "+@3*(@6/20)";
	//TString cb_lin_f	= "+@5*(@8/60)";
	TString cwww_ccw_f	= "+@8*(@5/12)*(@6/20)";
	TString cwww_cb_f	= "+@9*(@5/12)*(@7/60)";
	TString ccw_cb_f	= "+@10*(@6/20)*(@7/60)";
	RooFormulaVar Norm("Norm","@0"+cwww_f+ccw_f+ccw_lin_f+cb_f+cwww_ccw_f+cwww_cb_f+ccw_cb_f,N_norm_list);

	RooFormulaVar N1("N1","@1/@0",RooArgList(Norm,N_SM4fit));
	RooFormulaVar N2("N2","(@2*(@1/12)**2)/@0",RooArgList(Norm,cwww,N2_tmp));
	RooFormulaVar N4("N4","(@2*(@1/20)**2)/@0",RooArgList(Norm,ccw,N4_tmp));
	RooFormulaVar N5("N5","(@2*(@1/20))/@0",RooArgList(Norm,ccw,N5_tmp));
	RooFormulaVar N6("N6","(@2*(@1/60)**2)/@0",RooArgList(Norm,cb,N6_tmp));
	//RooFormulaVar N7("N7","(@2*(@1/60))/@0",RooArgList(Norm,cb,N7_tmp));
	RooFormulaVar N8("N8","(@3*(@1/12)*(@2/20))/@0",RooArgList(Norm,cwww,ccw,N8_tmp));
	//RooFormulaVar N9("N9","(@3*(@1/12)*(@2/60))/@0",RooArgList(Norm,cwww,cb,N9_tmp));
	RooFormulaVar N10("N10","(@3*(@1/20)*(@2/60))/@0",RooArgList(Norm,ccw,cb,N10_tmp));



	//RooArgList N_list(N1,N2,N4,N6);
	RooArgList N_list(N1,N2,N4,N5,N6);
	N_list.add(RooArgList(N8,N10));
	//RooArgList Pdf_list(*SM_Pdf,*Pdf_cwww,*Pdf_ccw,*Pdf_cb);
	RooArgList Pdf_list(*SM_Pdf,Pdf_cwww,Pdf_ccw,*Pdf_ccw_lin,Pdf_cb);
	Pdf_list.add(RooArgList(*Int_cwww_ccw,*Int_ccw_cb));
	Pdf_list.Print();
	N_list.Print();
	RooAddPdf model1("model","model",Pdf_list,N_list);
	
	RooRealVar * a1 	= w.var("a1");
	a2.setConstant(true);
	a3.setConstant(true);
	RooRealVar * a33	= w.var("a33");	a33->setConstant(true);
	a4.setConstant(true);
	//RooRealVar * a44	= w.var("a44");	a44->setConstant(true);
	RooRealVar * a5		= w.var("a5");	a5->setConstant(true);
	//RooRealVar * a6		= w.var("a6");	a6->setConstant(true);
	RooRealVar * a7		= w.var("a7");	a7->setConstant(true);

	RooBinning bins(600,3500);
	bins.addUniform(nbins-5,600,3500);
//SM-fit
	cwww.setVal(0);	ccw.setVal(0); cb.setVal(0); 	
	model1.fitTo(*w.data("cwww0ccw0cb0"));//hist0
	a1->setConstant(true);
//SM-interference-fits
	double N_SM_tmp_val = N_SM4fit.getVal();//SM
	double N2_tmp_val = N2_tmp.getVal();//cwww quad
	double N4_tmp_val = N4_tmp.getVal();//ccw quad
	double N6_tmp_val = N6_tmp.getVal();//cb quad
	double N5_tmp_val = N5_tmp.getVal();//ccw lin
	//double N7_tmp_val = N7_tmp.getVal();//cb lin
	N_SM4fit.setVal(0);N2_tmp.setVal(0);N4_tmp.setVal(0);N6_tmp.setVal(0);
	cwww.setVal(0);ccw.setVal(20); cb.setVal(0);
	a33->setConstant(false);
	model1.fitTo(hist_diff_ccw);
	a33->setConstant(true);
	//cwww.setVal(0);ccw.setVal(0); cb.setVal(60);
	//a44->setConstant(false);
	//model1.fitTo(hist_diff_cb);
	//a44->setConstant(true);
/*
	TCanvas cc2("smint","smint",1);
	cc2.cd();cc2.SetLogy();
	RooPlot * plott2 = MWW.frame();
	cwww.setVal(0);ccw.setVal(20); cb.setVal(0);
	hist_diff_ccw.plotOn(plott2,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E"));
	model1.plotOn(plott2,RooFit::LineColor(kBlue),RooFit::Normalization(hist_diff_ccw.sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0);ccw.setVal(0); cb.setVal(60);
	hist_diff_cb.plotOn(plott2,RooFit::Binning(bins),RooFit::LineColor(kGreen),RooFit::DrawOption("E"));
	//model1.plotOn(plott2,RooFit::LineColor(kGreen),RooFit::Normalization(hist_diff_cb.sumEntries(),RooAbsReal::NumEvent));
	plott2->GetYaxis()->SetRangeUser(0.001,5);
	plott2->Draw();
	cc2.Draw();cc2.Update();
*/
	//int cwww-ccw-fit

	a5->setVal(slopeval);
	a5->setConstant(true);

	a7->setVal(slopeval2);
	a7->setConstant(true);
/*
	TCanvas cc("atgcint","atgcint",1);
	cc.cd();cc.SetLogy();
	RooPlot * plott = MWW.frame();
	cwww.setVal(12);ccw.setVal(20);cb.setVal(0);
	//hist_diff_cwww_ccw.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kBlue),RooFit::DrawOption("E"));
	//model1.plotOn(plott,RooFit::LineColor(kBlue),RooFit::Normalization(hist_diff_cwww_ccw.sumEntries(),RooAbsReal::NumEvent));
	//cwww.setVal(12); ccw.setVal(0); cb.setVal(60);
	//hist_diff_cwww_cb.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kRed),RooFit::DrawOption("E"));
	//model1.plotOn(plott,RooFit::LineColor(kRed),RooFit::Normalization(hist_diff_cwww_cb.sumEntries(),RooAbsReal::NumEvent));
	cwww.setVal(0); ccw.setVal(20); cb.setVal(60);
	hist_diff_ccw_cb.plotOn(plott,RooFit::Binning(bins),RooFit::LineColor(kGreen),RooFit::DrawOption("E"));
	model1.plotOn(plott,RooFit::LineColor(kGreen),RooFit::Normalization(hist_diff_ccw_cb.sumEntries(),RooAbsReal::NumEvent));
	plott->GetYaxis()->SetRangeUser(-1,-0.0001);
	plott->Draw();
	cc.Draw();
	cc.Update();
*/

	N_SM4fit.setVal(N_SM_tmp_val);
	N2_tmp.setVal(N2_tmp_val);
	N4_tmp.setVal(N4_tmp_val);
	N5_tmp.setVal(N5_tmp_val);
	N6_tmp.setVal(N6_tmp_val);
	//N7_tmp.setVal(N7_tmp_val);
//cwww-fit
	cwww.setVal(12); ccw.setVal(0); cb.setVal(0); 
	a2.setConstant(false);
	Erf_offset_cwww.setConstant(false);
	Erf_width_cwww.setConstant(false);
	model1.fitTo(*w.data("cwww12ccw0cb0"));//hist13, hist128
	a2.setConstant(true);
	Erf_offset_cwww.setConstant(true);
	Erf_width_cwww.setConstant(true);
//ccw-fit
	cwww.setVal(0);	ccw.setVal(20); cb.setVal(0); 
	a3.setConstant(false);
	Erf_offset_ccw.setConstant(false);
	Erf_width_ccw.setConstant(false);
	model1.fitTo(*w.data("cwww0ccw20cb0"));//hist53, hist134
	a3.setConstant(true);
	Erf_offset_ccw.setConstant(true);
	Erf_width_ccw.setConstant(true);
//cb-fit
	cwww.setVal(0); ccw.setVal(0); cb.setVal(60);
	a4.setConstant(false);
	Erf_offset_cb.setConstant(false);
	Erf_width_cb.setConstant(false);
	model1.fitTo(*w.data("cwww0ccw0cb60"));//hist61, hist 136
	a4.setConstant(true);
	Erf_offset_cb.setConstant(true);
	Erf_width_cb.setConstant(true);

	w2.import(*a1);
	w2.import(a2);
	w2.import(a3);
	w2.import(*a33);
	w2.import(a4);
	//w2.import(*a44);
	w2.import(*a5);
	//w2.import(*a6);
	w2.import(*a7);


	TFile * fileOut = new TFile("genlevel_WZ_"+ch+".root","RECREATE");
	w2.Write();
	fileOut->Close();


	double sum_SM	= w.data("cwww0ccw0cb0")->sumEntries();
	//scaling for ccw
	TH1F hist4scale_ccw("hist4scale_ccw","hist4scale_ccw",3,-30,30);
	hist4scale_ccw.SetBinContent(1,w.data("cwww0ccw_20cb0")->sumEntries()/sum_SM);
	hist4scale_ccw.SetBinContent(2,w.data("cwww0ccw0cb0")->sumEntries()/sum_SM);
	hist4scale_ccw.SetBinContent(3,w.data("cwww0ccw20cb0")->sumEntries()/sum_SM);
	hist4scale_ccw.Fit("pol2");
	TF1 * fitfunc_ccw	= hist4scale_ccw.GetFunction("pol2");
	RooRealVar par0_ccw("par0_ccw","par0_ccw",fitfunc_ccw->GetParameter(0));
	RooRealVar par1_ccw("par1_ccw","par1_ccw",fitfunc_ccw->GetParameter(1));
	RooRealVar par2_ccw("par2_ccw","par2_ccw",fitfunc_ccw->GetParameter(2));
	RooFormulaVar scaleshape_ccw("scaleshape_ccw","scaleshape_ccw","(@0+@1*@3+@2*@3**2)",RooArgList(par0_ccw,par1_ccw,par2_ccw,ccw));
	//scaling for cb
	TH1F hist4scale_cb("hist4scale_cb","hist4scale_cb",3,-90,90);
	hist4scale_cb.SetBinContent(1,w.data("cwww0ccw0cb_60")->sumEntries()/sum_SM);
	hist4scale_cb.SetBinContent(2,w.data("cwww0ccw0cb0")->sumEntries()/sum_SM);
	hist4scale_cb.SetBinContent(3,w.data("cwww0ccw0cb60")->sumEntries()/sum_SM);
	hist4scale_cb.Fit("pol2");
	TCanvas ctmp("ctmp","ctmp",1);
	ctmp.cd();
	hist4scale_cb.Draw();
	ctmp.Draw();
	ctmp.Update();
	int abs;
	abs=getchar();
	TF1 * fitfunc_cb	= hist4scale_cb.GetFunction("pol2");
	RooRealVar par0_cb("par0_cb","par0_cb",fitfunc_cb->GetParameter(0));
	RooRealVar par1_cb("par1_cb","par1_cb",fitfunc_cb->GetParameter(1));
	RooRealVar par2_cb("par2_cb","par2_cb",fitfunc_cb->GetParameter(2));
	RooFormulaVar scaleshape_cb("scaleshape_cb","scaleshape_cb","(@0+@1*@3+@2*@3**2)",RooArgList(par0_cb,par1_cb,par2_cb,cb));


	TCanvas cc1("smintccwpos","smintccwpos",1);
	cc1.cd();cc1.SetLogy();
	RooPlot * pp = MWW.frame(900,3500);
	w.data("cwww0ccw20cb0")->plotOn(pp,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(20);cb.setVal(0);
	model1.plotOn(pp,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw3_5cb0")->plotOn(pp,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(3.5);cb.setVal(0);
	model1.plotOn(pp,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw10cb0")->plotOn(pp,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(10);cb.setVal(0);
	model1.plotOn(pp,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	for(unsigned int i=0; i<20; i++)
	{
		ccw.setVal(i);
		if(i==0)
			model1.plotOn(pp,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		else
			model1.plotOn(pp,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		cout<<i<<" : "<<scaleshape_ccw.getVal()*sum_SM<<endl;
	}
	w.data("cwww0ccw0cb0")->plotOn(pp,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));
	pp->GetYaxis()->SetRangeUser(0.00075,3);
	pp->Draw();
	cc1.Draw();cc1.Update();
	cc1.SaveAs("ARC_plots/ccw_SMint_pos.png");

	TCanvas ccc2("smintccwneg","smintccwneg",1);
	ccc2.cd();ccc2.SetLogy();
	RooPlot * pp2 = MWW.frame(600,3500);
	w.data("cwww0ccw_20cb0")->plotOn(pp2,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(-20);cb.setVal(0);
	model1.plotOn(pp2,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw_10cb0")->plotOn(pp2,RooFit::LineColor(kGreen),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(-10);cb.setVal(0);
	model1.plotOn(pp2,RooFit::LineColor(kGreen),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw_3_5cb0")->plotOn(pp2,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);ccw.setVal(-3.5);cb.setVal(0);
	model1.plotOn(pp2,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
	ccw.setVal(-20);
	for(int i=0; i<20; i++)
	{
		ccw.setVal(i*(-1));
		if(i==0)
			model1.plotOn(pp2,RooFit::LineWidth(3),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		else
			model1.plotOn(pp2,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		cout<<i<<" : "<<ccw.getVal()<<", "<<scaleshape_ccw.getVal()*sum_SM<<endl;
	}
	w.data("cwww0ccw0cb0")->plotOn(pp2,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));
	pp2->GetYaxis()->SetRangeUser(0.00075,3);
	pp2->Draw();
	ccc2.Draw();ccc2.Update();
	ccc2.SaveAs("ARC_plots/ccw_SMint_neg.png");

	TCanvas ccc3("smintcbneg","smintcbneg",1);
	ccc3.cd();ccc3.SetLogy();
	RooPlot * pp3 = MWW.frame(900,3500);
	w.data("cwww0ccw0cb_60")->plotOn(pp3,RooFit::LineColor(kBlue),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);cb.setVal(0);cb.setVal(-60);
	model1.plotOn(pp3,RooFit::LineColor(kBlue),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw0cb_30")->plotOn(pp3,RooFit::LineColor(kMagenta),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);cb.setVal(0);cb.setVal(-30);
	model1.plotOn(pp3,RooFit::LineColor(kMagenta),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
	w.data("cwww0ccw0cb_10")->plotOn(pp3,RooFit::LineColor(kCyan),RooFit::DrawOption("E"),RooFit::Binning(bins));
	cwww.setVal(0);cb.setVal(0);cb.setVal(-10);
	model1.plotOn(pp3,RooFit::LineColor(kCyan),RooFit::Normalization(scaleshape_cb.getVal()*sum_SM,RooAbsReal::NumEvent));
	cb.setVal(-20);
	for(int i=0; i<60; i++)
	{
		cb.setVal(i*(-3));
		if(i==0)
			model1.plotOn(pp3,RooFit::LineWidth(2),RooFit::LineColor(kBlack),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		else
			model1.plotOn(pp3,RooFit::LineWidth(1),RooFit::LineColor(kGray),RooFit::Normalization(scaleshape_ccw.getVal()*sum_SM,RooAbsReal::NumEvent));
		cout<<i<<" : "<<cb.getVal()<<", "<<scaleshape_cb.getVal()*sum_SM<<endl;
	}
	w.data("cwww0ccw0cb0")->plotOn(pp3,RooFit::LineColor(kBlack),RooFit::DrawOption("E"),RooFit::Binning(bins));
	pp3->GetYaxis()->SetRangeUser(0.00075,3);
	pp3->Draw();
	ccc3.Draw();ccc3.Update();














	end = getchar();

	for(int i = 1; i<8; i++)
	{
		if (i!=6)
			std::cout<<("a"+to_string(i)).c_str()<<": "<<w2.var(("a"+to_string(i)).c_str())->getVal()<<std::endl;
		if(i==3)
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
	std::cout<<", N5: "<<N5.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N10:"<<N10.getVal()<<std::endl;

	cwww.setVal(1.2); ccw.setVal(2); cb.setVal(6);
	std::cout<<"all params 1/10:"<<std::endl;
	std::cout<<"N_SM: "<<N1.getVal()<<std::endl;
	std::cout<<"N2: "<<N2.getVal()<<", N4: "<<N4.getVal()<<", N6:"<<N6.getVal()<<std::endl;
	std::cout<<", N5: "<<N5.getVal()<<std::endl;
	std::cout<<"N8: "<<N8.getVal()<<", N10:"<<N10.getVal()<<std::endl;


	std::cout<<"cwww/ccw: "<<hist_diff_cwww_ccw.sumEntries()<<std::endl;
	std::cout<<"cwww/cb: "<<hist_diff_cwww_cb.sumEntries()<<std::endl;
	std::cout<<"ccw/cb: "<<hist_diff_ccw_cb.sumEntries()<<std::endl;	

	exit(0);
}


