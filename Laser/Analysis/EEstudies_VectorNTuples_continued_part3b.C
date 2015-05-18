#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3F.h>
#include <TH3D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <THistPainter.h>
#include <TGraphPainter.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TF1.h>
#include <TString.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaletteAxis.h>
#include <TColor.h>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <vector>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <fstream>
//#include "LaserAnalysisReducedTrees.h"
#include "LaserAnalysisVectorTrees.C"
#include "RootMacros/Utilities.hh"


using namespace std;

//plot intercepts and slopes vs. integrated lumi
void EEstudies_VectorNTuples_continued_part3b(){

	bool startdiff         = false;//only true if normalized == false, always diff w.r.t. 2011 beginning
	bool normalized        = false;
	const int period       = 7;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11
				//5 periods for only2011, 2012 and 2011: 7 periods
	bool samelumi          = true;//if same lumi = true, use only last three periods with first 1.1fb-1 (2011), or 950 pb-1 (20112012)
	bool only2011          = false;

	bool plotscaled        = false;//plot also plots where VPT/PN has been scaled.

	bool plotslopeintercept= true; //you want this to be true
	bool VPTdiff_vs_mustd  = true; //you probably want this to be true -> includes fit
	bool plot_av_rms       = false; //plot does not help much
	bool print_av_rms      = true; //instead just print out the numbers
	bool plotdistancetofit = false; //distance to fit
	bool plotonlyAbsEta    = true; //if true don't print EE+ EE- separately


	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	string string_prod[2] = {"__russ", "__chin"};
	string string_aeta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_acap[1] = {"_"};
	string aname;

	TString outputdir;
	outputdir = "Plots/20121023/EEstudies/continued/diffwrt2011begin/InterceptSlopesVSintegrLumi";
	if(normalized) outputdir = "Plots/20121023/EEstudies/continued/normalized/InterceptSlopesVSintegrLumi";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
	else if(!startdiff) outputdir = "Plots/20121023/EEstudies/continued/notnormalized/InterceptSlopesVSintegrLumi";
	if(samelumi&& only2011) outputdir = outputdir + (TString)"_sameLumi1p1fb-1";
	if(samelumi&&!only2011) outputdir = outputdir + (TString)"_sameLumi950pb-1";
	Util::MakeOutputDir(outputdir);

	TString periodTS[period]; TString periods[period];
	TFile* oldfile[period];
	TH1D* a_intercept_vs_eta[period][2];
	TH1D* a_intercept_vs_eta_scaled[period][2];
	TH1D* a_slope_vs_eta[period][2];
	TH1D* a_slope_vs_eta_scaled[period][2];

	for(int p = 0; p<period; ++p){
		//p is periodnumber - 1
		string pe;
//		     if(p==0) {periodTS[p] = "/period150311_270311"; periods[p] = "_period150311_270311"; pe="_period1";}
//		else if(p==1) {periodTS[p] = "/period120411_050511"; periods[p] = "_period120411_050511"; pe="_period2";}
//		else if(p==2) {periodTS[p] = "/period130511_010711"; periods[p] = "_period130511_010711"; pe="_period3";}
//		else if(p==3) {periodTS[p] = "/period150711_250811"; periods[p] = "_period150711_250811"; pe="_period4";}
//		else if(p==4) {periodTS[p] = "/period040911_021111"; periods[p] = "_period040911_021111"; pe="_period5";}
//		else if(p==3&&samelumi){periodTS[p] = "/period150711_080811_1p1fb-1"; periods[p] = "_period150711_080811_1p1fb-1"; pe="_period4_1p1fb-1";}
//		else if(p==4&&samelumi){periodTS[p] = "/period040911_230911_1p1fb-1"; periods[p] = "_period040911_021111_1p1fb-1"; pe="_period5_1p1fb-1";}

		if     (p==0) {periodTS[p] = "/period110315_110327"; periods[p] = "_period110315_110327"; pe="_period1";}
		else if(p==1) {periodTS[p] = "/period110412_110505"; periods[p] = "_period110412_110505"; pe="_period2";}
		else if(p==2) {periodTS[p] = "/period110513_110701"; periods[p] = "_period110513_110701"; pe="_period3";}
		else if(p==3) {periodTS[p] = "/period110715_110825"; periods[p] = "_period110715_110825"; pe="_period4";}
		else if(p==4) {periodTS[p] = "/period110904_111102"; periods[p] = "_period110904_111102"; pe="_period5";}
		else if(p==5) {periodTS[p] = "/period120406_120420"; periods[p] = "_period120406_120420"; pe="_period6";}
		else if(p==6) {periodTS[p] = "/period120505_120618"; periods[p] = "_period120505_120618"; pe="_period7";}
		else if(p==3&&samelumi&& only2011){periodTS[p] = "/period110715_110808_1p1fb-1"; periods[p] = "_period110715_110808_1p1fb-1"; pe="_period4";}
		else if(p==4&&samelumi&& only2011){periodTS[p] = "/period110904_110923_1p1fb-1"; periods[p] = "_period110904_110923_1p1fb-1"; pe="_period5";}
		else if(p==2&&samelumi&&!only2011){periodTS[p] = "/period110513_110624_950pb-1"; periods[p] = "_period110513_110624_950pb-1"; pe="_period3";}
		else if(p==3&&samelumi&&!only2011){periodTS[p] = "/period110715_110807_950pb-1"; periods[p] = "_period110715_110807_950pb-1"; pe="_period4";}
		else if(p==4&&samelumi&&!only2011){periodTS[p] = "/period110904_110922_950pb-1"; periods[p] = "_period110904_110922_950pb-1"; pe="_period5";}
		else if(p==6&&samelumi&&!only2011){periodTS[p] = "/period120505_110517_950pb-1"; periods[p] = "_period120505_110517_950pb-1"; pe="_period7";}

		oldfile[p] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_diffwrt2011begin"+periods[p]+".root");
		if(normalized) oldfile[p] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted"+periods[p]+".root");//mustd <-> muSIC
		else if(!startdiff) oldfile[p] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_nonormalization"+periods[p]+".root");
		a_intercept_vs_eta[p][0]        = (TH1D*)oldfile[p]->Get("intercept_vs_AbsEta___russ");
		a_intercept_vs_eta[p][1]        = (TH1D*)oldfile[p]->Get("intercept_vs_AbsEta___chin");
		a_intercept_vs_eta_scaled[p][0] = (TH1D*)oldfile[p]->Get("intercept_vs_AbsEta__scaled___russ");
		a_intercept_vs_eta_scaled[p][1] = (TH1D*)oldfile[p]->Get("intercept_vs_AbsEta__scaled___chin");
		a_slope_vs_eta[p][0]            = (TH1D*)oldfile[p]->Get("slope_vs_AbsEta___russ");
		a_slope_vs_eta[p][1]            = (TH1D*)oldfile[p]->Get("slope_vs_AbsEta___chin");
		a_slope_vs_eta_scaled[p][0]     = (TH1D*)oldfile[p]->Get("slope_vs_AbsEta__scaled___russ");
		a_slope_vs_eta_scaled[p][1]     = (TH1D*)oldfile[p]->Get("slope_vs_AbsEta__scaled___chin");


		aname = "intercept_vs_AbsEta__russ" + pe;
		a_intercept_vs_eta[p][0]        ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "intercept_vs_AbsEta__chin" + pe;
		a_intercept_vs_eta[p][1]        ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "intercept_vs_AbsEta__scaled__russ" + pe;
		a_intercept_vs_eta_scaled[p][0] ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "intercept_vs_AbsEta__scaled__chin" + pe;
		a_intercept_vs_eta_scaled[p][1] ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "slope_vs_AbsEta__russ" + pe;
		a_slope_vs_eta[p][0]            ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "slope_vs_AbsEta__chin" + pe;
		a_slope_vs_eta[p][1]            ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "slope_vs_AbsEta__scaled__russ" + pe;
		a_slope_vs_eta_scaled[p][0]     ->SetNameTitle(aname.c_str(),aname.c_str());
		aname = "slope_vs_AbsEta__scaled__chin" + pe;
		a_slope_vs_eta_scaled[p][1]     ->SetNameTitle(aname.c_str(),aname.c_str());

	}

	//new definitions
	//const int Ntimebins = 11;
	//time between bins: beginning 2011, period1, period2, period3, period4, period5, end 2011
	//double timebins[Ntimebins+1] = {0., 30., 63., 64., 323., 324., 1430., 1431., 2910., 2911., 6123., 6124.};

	//lumi: 0: 2010, 1-5 beginning of period n, last bin: end of 2012
//	const int Ntimebins = 6;
//	double timebins[Ntimebins+1] = {0., 45., 78., 338., 1445., 2925., 6138.};
//	const int Ntimebins = 11;
//	double timebins[Ntimebins+1] = {0., 45., 77., 78., 337., 338., 1444., 1445., 2924., 2925., 6137., 6138.};//last ones are dummys for 2011
//	if(samelumi){ timebins[8] = 2538.; timebins[10] = 3986.; }
	const int Ntimebins = 15;
	double timebins[Ntimebins+1] = {0., 45., 77., 78., 337., 338., 1444., 1445., 2924., 2925., 6137., 6177.,7080.,7178.,12777.,15000.};
	if(only2011){//last ones are dummys for 2011
		timebins[11] = 6138.; timebins[12] = 6139.; timebins[13] = 6140.; timebins[14] = 6141.; timebins[15] = 7500.; }
	if( only2011 && samelumi){ timebins[8] = 2538.; timebins[10] = 3986.; }
	if(!only2011 && samelumi){ timebins[6] = 1300.; timebins[8] = 2400.; timebins[10] = 3875.; timebins[14] = 8130.;}


	TH1D *haxis = new TH1D("haxis", "", Ntimebins, timebins);

	TH1D * h_intercept_vs_time[8][2];
	TH1D * h_intercept_vs_time_scaled[8][2];
	TH1D * h_slope_vs_time[8][2];
	TH1D * h_slope_vs_time_scaled[8][2];
	for(int e = 0; e<8; ++e){
	   for(int nn = 0; nn<2; ++nn){
		aname = (string)"intercept_vs_integrLumi_" + string_aeta[e] + string_prod[nn];
	//	h_intercept_vs_time[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_intercept_vs_time[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_intercept_vs_time[e][nn]->Sumw2();
		aname = (string)"intercept_vs_integrLumi_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
	//	h_intercept_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_intercept_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_intercept_vs_time_scaled[e][nn]->Sumw2();
		aname = (string)"slope_vs_integrLumi_" + string_aeta[e] + string_prod[nn];
	//	h_slope_vs_time[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_slope_vs_time[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_slope_vs_time[e][nn]->Sumw2();
		aname = (string)"slope_vs_integrLumi_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
	//	h_slope_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_slope_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_slope_vs_time_scaled[e][nn]->Sumw2();

		h_intercept_vs_time[e][nn]->SetLineWidth(4); h_intercept_vs_time[e][nn]->SetMarkerStyle(20); h_intercept_vs_time[e][nn]->SetMarkerSize(2);
		h_intercept_vs_time_scaled[e][nn]->SetLineWidth(4); h_intercept_vs_time_scaled[e][nn]->SetMarkerStyle(20); h_intercept_vs_time_scaled[e][nn]->SetMarkerSize(2);
		h_slope_vs_time[e][nn]->SetLineWidth(4); h_slope_vs_time[e][nn]->SetMarkerStyle(20); h_slope_vs_time[e][nn]->SetMarkerSize(2);
		h_slope_vs_time_scaled[e][nn]->SetLineWidth(4); h_slope_vs_time_scaled[e][nn]->SetMarkerStyle(20); h_slope_vs_time_scaled[e][nn]->SetMarkerSize(2);
		if(nn==0){
		h_intercept_vs_time[e][nn]->SetLineColor(kRed); h_intercept_vs_time[e][nn]->SetMarkerColor(kRed);
		h_intercept_vs_time_scaled[e][nn]->SetLineColor(kRed); h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kRed);
		h_slope_vs_time[e][nn]->SetLineColor(kRed); h_slope_vs_time[e][nn]->SetMarkerColor(kRed);
		h_slope_vs_time_scaled[e][nn]->SetLineColor(kRed); h_slope_vs_time_scaled[e][nn]->SetMarkerColor(kRed);
		}
		else {
		h_intercept_vs_time[e][nn]->SetLineColor(kBlue); h_intercept_vs_time[e][nn]->SetMarkerColor(kBlue);
		h_intercept_vs_time_scaled[e][nn]->SetLineColor(kBlue); h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kBlue);
		h_slope_vs_time[e][nn]->SetLineColor(kBlue); h_slope_vs_time[e][nn]->SetMarkerColor(kBlue);
		h_slope_vs_time_scaled[e][nn]->SetLineColor(kBlue); h_slope_vs_time_scaled[e][nn]->SetMarkerColor(kBlue);
		}
	   }
	}

	//fill the plots
	for(int e = 0; e<8; ++e){//eta
	   for(int nn = 0; nn<2; ++nn){//produer
	      for(int p = 0; p<period; ++p){//period
		//p+2, since 1st bin starts at 1, not at 0, and first bin is empty, 2nd bin is period1, every second bin is a empty bin between the periods
		//e+1, since 1st bin starts at 1, not at 0
		//if(p==0) continue;
		if(samelumi && (p==0||p==1) && only2011) continue;
		if(samelumi && (p==0||p==1) && !only2011) continue;
		if(only2011 && p>=5) continue;
		h_intercept_vs_time[e][nn]       ->SetBinContent(2*p+2, a_intercept_vs_eta[p][nn]       ->GetBinContent(e+1));
		h_intercept_vs_time[e][nn]       ->SetBinError(  2*p+2, a_intercept_vs_eta[p][nn]       ->GetBinError(  e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinContent(2*p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinContent(e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinError(  2*p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinError(  e+1));
		h_slope_vs_time[e][nn]           ->SetBinContent(2*p+2, a_slope_vs_eta[p][nn]           ->GetBinContent(e+1));
		h_slope_vs_time[e][nn]           ->SetBinError(  2*p+2, a_slope_vs_eta[p][nn]           ->GetBinError(  e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinContent(2*p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinContent(e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinError(  2*p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinError(  e+1));
/*
		h_intercept_vs_time[e][nn]       ->SetBinContent(p+2, a_intercept_vs_eta[p][nn]       ->GetBinContent(e+1));
		h_intercept_vs_time[e][nn]       ->SetBinError(  p+2, a_intercept_vs_eta[p][nn]       ->GetBinError(  e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinContent(p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinContent(e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinError(  p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinError(  e+1));
		h_slope_vs_time[e][nn]           ->SetBinContent(p+2, a_slope_vs_eta[p][nn]           ->GetBinContent(e+1));
		h_slope_vs_time[e][nn]           ->SetBinError(  p+2, a_slope_vs_eta[p][nn]           ->GetBinError(  e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinContent(p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinContent(e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinError(  p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinError(  e+1));
*/
	      }
	   }
	}

	//plot plots
	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	TString outputfile;
	for(int e = 0; e<8; ++e){
	   for(int nn = 0; nn<2; ++nn){
	      if(plotslopeintercept){
		col->cd();
		col->SetName(h_intercept_vs_time[e][nn]->GetName() );
		col->SetTitle(h_intercept_vs_time[e][nn]->GetTitle() );
		h_intercept_vs_time[e][nn]->GetXaxis()->SetTitle("L_{int} (pb^{-1})");
		h_intercept_vs_time[e][nn]->GetYaxis()->SetTitle("intercept");
		if(only2011){
		if(startdiff)       {h_intercept_vs_time[e][nn]->SetMaximum(0.4);   h_intercept_vs_time[e][nn]->SetMinimum(-0.02); }
		if(normalized)      {h_intercept_vs_time[e][nn]->SetMaximum(0.25);  h_intercept_vs_time[e][nn]->SetMinimum(-0.02); }
		else if(!startdiff) {h_intercept_vs_time[e][nn]->SetMaximum(0.185); h_intercept_vs_time[e][nn]->SetMinimum(-0.01); }
		} else{
		if(startdiff)       {h_intercept_vs_time[e][nn]->SetMaximum(0.55);  h_intercept_vs_time[e][nn]->SetMinimum(-0.02); }
		if(normalized)      {h_intercept_vs_time[e][nn]->SetMaximum(0.27);  h_intercept_vs_time[e][nn]->SetMinimum(-0.025);}
		else if(!startdiff) {h_intercept_vs_time[e][nn]->SetMaximum(0.22);  h_intercept_vs_time[e][nn]->SetMinimum(-0.01); }
		}
		h_intercept_vs_time[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_intercept_vs_time[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_intercept_vs_time[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_intercept_vs_time[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_intercept_vs_time[e][nn]->GetXaxis()->SetLabelOffset(0.02); 
		h_intercept_vs_time[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_intercept_vs_time[e][nn]->Draw();
		outputfile = h_intercept_vs_time[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		if(plotscaled){
		col->cd();
		col->SetName(h_intercept_vs_time_scaled[e][nn]->GetName() );
		col->SetTitle(h_intercept_vs_time_scaled[e][nn]->GetTitle() );
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTitle("L_{int} (pb^{-1})");
		h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetTitle("intercept (scaled)");
		if(only2011){
		if(startdiff)       {h_intercept_vs_time_scaled[e][nn]->SetMaximum(5.);  h_intercept_vs_time_scaled[e][nn]->SetMinimum(-2.);}
		if(normalized)      {h_intercept_vs_time_scaled[e][nn]->SetMaximum(6.);  h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		else if(!startdiff) {h_intercept_vs_time_scaled[e][nn]->SetMaximum(7.);  h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		} else{
		if(startdiff)       {h_intercept_vs_time_scaled[e][nn]->SetMaximum(5.5); h_intercept_vs_time_scaled[e][nn]->SetMinimum(-2.);}
		if(normalized)      {h_intercept_vs_time_scaled[e][nn]->SetMaximum(6.);  h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		else if(!startdiff) {h_intercept_vs_time_scaled[e][nn]->SetMaximum(7.);  h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		}
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetLabelOffset(0.02); 
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_intercept_vs_time_scaled[e][nn]->Draw();
		outputfile = h_intercept_vs_time_scaled[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		}
		col->cd();
		col->SetName(h_slope_vs_time[e][nn]->GetName() );
		col->SetTitle(h_slope_vs_time[e][nn]->GetTitle() );
		h_slope_vs_time[e][nn]->GetXaxis()->SetTitle("L_{int} (pb^{-1})");
		h_slope_vs_time[e][nn]->GetYaxis()->SetTitle("slope");
		if(only2011){
		if(startdiff)       {h_slope_vs_time[e][nn]->SetMaximum(0.20); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time[e][nn]->SetMaximum(0.16); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time[e][nn]->SetMaximum(0.10); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		} else{
		if(startdiff)       {h_slope_vs_time[e][nn]->SetMaximum(0.24); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time[e][nn]->SetMaximum(0.17); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time[e][nn]->SetMaximum(0.10); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		}
		h_slope_vs_time[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_slope_vs_time[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_slope_vs_time[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_slope_vs_time[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_slope_vs_time[e][nn]->GetXaxis()->SetLabelOffset(0.02);  
		h_slope_vs_time[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_slope_vs_time[e][nn]->Draw();
		outputfile = h_slope_vs_time[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		if(plotscaled){
		col->cd();
		col->SetName(h_slope_vs_time_scaled[e][nn]->GetName() );
		col->SetTitle(h_slope_vs_time_scaled[e][nn]->GetTitle() );
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTitle("L_{int} (pb^{-1})");
		h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetTitle("slope (scaled)");
		if(only2011){
		if(startdiff)       {h_slope_vs_time_scaled[e][nn]->SetMaximum(15.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time_scaled[e][nn]->SetMaximum(10.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time_scaled[e][nn]->SetMaximum(10.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		} else{
		if(startdiff)       {h_slope_vs_time_scaled[e][nn]->SetMaximum(15.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time_scaled[e][nn]->SetMaximum(10.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time_scaled[e][nn]->SetMaximum(10.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		}
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetLabelOffset(0.02);
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_slope_vs_time_scaled[e][nn]->Draw();
		outputfile = h_slope_vs_time_scaled[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		}
	      }
	   }
	}

	//save plots
	TFile *newFile;
	if(only2011){
	if(samelumi){
	if(startdiff)  newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin_only2011_sameLumi1p1fb-1.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized_only2011_sameLumi1p1fb-1.root", "RECREATE");
	else if(!startdiff) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_sameLumi1p1fb-1.root", "RECREATE");
	} else{
	if(startdiff)  newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_diffwrt2011begin.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_normalized.root", "RECREATE");
	else if(!startdiff) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011.root", "RECREATE");
	}
	} else{
	if(samelumi){
	if(startdiff)  newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin_sameLumi950pb-1.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized_sameLumi950pb-1.root", "RECREATE");
	else if(!startdiff) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_sameLumi950pb-1.root", "RECREATE");
	} else{
	if(startdiff)  newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized.root", "RECREATE");//mustd <-> muSIC
	else if(!startdiff) newFile = new TFile("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples.root", "RECREATE");
	}
	}
	newFile->cd();
	for(int e = 0; e<8; ++e){
	   for(int nn = 0; nn<2; ++nn){
		h_intercept_vs_time[e][nn]->Write();
		h_intercept_vs_time_scaled[e][nn]->Write();
		h_slope_vs_time[e][nn]->Write();
		h_slope_vs_time_scaled[e][nn]->Write();
	   }
	}

}