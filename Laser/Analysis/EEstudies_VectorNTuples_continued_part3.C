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

//plot intercept and slopes vs. time
void EEstudies_VectorNTuples_continued_part3(){

	bool startdiff         = false;//only true if normalized == false, always diff w.r.t. beginning 2011
	bool normalized        = false;
	const int period       = 5;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11
	bool samelumi          = false;//if same lumi = true, use only last three periods with first 1.1fb-1

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
	outputdir = "Plots/20121008/EEstudies/continued/diffwrt2011begin/InterceptSlopesVStime";
	if(normalized) outputdir = "Plots/20121008/EEstudies/continued/normalized/InterceptSlopesVStime";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
	else if(!startdiff) outputdir = "Plots/20121008/EEstudies/continued/notnormalized/InterceptSlopesVStime";
	if(samelumi) outputdir = outputdir + (TString)"_sameLumi1p1fb-1";
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
		     if(p==0) {periodTS[p] = "/period150311_270311"; periods[p] = "_period150311_270311"; pe="_period1";}
		else if(p==1) {periodTS[p] = "/period120411_050511"; periods[p] = "_period120411_050511"; pe="_period2";}
		else if(p==2) {periodTS[p] = "/period130511_010711"; periods[p] = "_period130511_010711"; pe="_period3";}
		else if(p==3) {periodTS[p] = "/period150711_250811"; periods[p] = "_period150711_250811"; pe="_period4";}
		else if(p==4) {periodTS[p] = "/period040911_021111"; periods[p] = "_period040911_021111"; pe="_period5";}
		else if(p==3&&samelumi){periodTS[p] = "/period150711_080811_1p1fb-1"; periods[p] = "_period150711_080811_1p1fb-1"; pe="_period4_1p1fb-1";}
		else if(p==4&&samelumi){periodTS[p] = "/period040911_230911_1p1fb-1"; periods[p] = "_period040911_021111_1p1fb-1"; pe="_period5_1p1fb-1";}

		oldfile[p] = TFile::Open("FileOutputs/20121008_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_diffwrt2011begin"+periods[p]+".root");
		if(normalized) oldfile[p] = TFile::Open("FileOutputs/20121008_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted"+periods[p]+".root");//mustd <-> muSIC
		else if(!startdiff) oldfile[p] = TFile::Open("FileOutputs/20121008_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_nonormalization"+periods[p]+".root");

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
	const int Ntimebins = 11;
	//time between bins: beginning 2011, period1, period2, period3, period4, period5, end 2011
	double timebins[Ntimebins+1] = {1293840000., 1300200000., 1301210000., 1302590000., 1304620000., 1305260000., 1309540000., 1310710000., 1314270000., 1315120000., 1320250000., 1325376000.};
	if(samelumi){
		timebins[8] = 1312800923.; timebins[10] = 1316813137.;
	}

	TH1D *haxis = new TH1D("haxis", "", Ntimebins, timebins);

	TH1D * h_intercept_vs_time[8][2];
	TH1D * h_intercept_vs_time_scaled[8][2];
	TH1D * h_slope_vs_time[8][2];
	TH1D * h_slope_vs_time_scaled[8][2];
	for(int e = 0; e<8; ++e){
	   for(int nn = 0; nn<2; ++nn){
		aname = (string)"intercept_vs_time_" + string_aeta[e] + string_prod[nn];
	//	h_intercept_vs_time[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_intercept_vs_time[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_intercept_vs_time[e][nn]->Sumw2();
		aname = (string)"intercept_vs_time_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
	//	h_intercept_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_intercept_vs_time_scaled[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_intercept_vs_time_scaled[e][nn]->Sumw2();
		aname = (string)"slope_vs_time_" + string_aeta[e] + string_prod[nn];
	//	h_slope_vs_time[e][nn] = new TH1D(aname.c_str(), "", Ntimebins, timebins);
		h_slope_vs_time[e][nn] = new TH1D(aname.c_str(), aname.c_str(), Ntimebins, timebins);
		h_slope_vs_time[e][nn]->Sumw2();
		aname = (string)"slope_vs_time_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
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
	for(int e = 0; e<8; ++e){
	   for(int nn = 0; nn<2; ++nn){
	      for(int p = 0; p<period; ++p){
		//p+2, since 1st bin starts at 1, not at 0, and first bin is empty, 2nd bin is period1, every second bin is a empty bin between the periods
		//e+1, since 1st bin starts at 1, not at 0
		//if(p==0) continue;
		if(samelumi && (p==0||p==1)) continue;
		h_intercept_vs_time[e][nn]       ->SetBinContent(2*p+2, a_intercept_vs_eta[p][nn]       ->GetBinContent(e+1));
		h_intercept_vs_time[e][nn]       ->SetBinError(  2*p+2, a_intercept_vs_eta[p][nn]       ->GetBinError(  e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinContent(2*p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinContent(e+1));
		h_intercept_vs_time_scaled[e][nn]->SetBinError(  2*p+2, a_intercept_vs_eta_scaled[p][nn]->GetBinError(  e+1));
		h_slope_vs_time[e][nn]           ->SetBinContent(2*p+2, a_slope_vs_eta[p][nn]           ->GetBinContent(e+1));
		h_slope_vs_time[e][nn]           ->SetBinError(  2*p+2, a_slope_vs_eta[p][nn]           ->GetBinError(  e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinContent(2*p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinContent(e+1));
		h_slope_vs_time_scaled[e][nn]    ->SetBinError(  2*p+2, a_slope_vs_eta_scaled[p][nn]    ->GetBinError(  e+1));
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
		h_intercept_vs_time[e][nn]->GetXaxis()->SetTitle("time");
		h_intercept_vs_time[e][nn]->GetYaxis()->SetTitle("intercept");
		if(startdiff)       {h_intercept_vs_time[e][nn]->SetMaximum(0.38);  h_intercept_vs_time[e][nn]->SetMinimum(-0.02); }
		if(normalized)      {h_intercept_vs_time[e][nn]->SetMaximum(0.2);   h_intercept_vs_time[e][nn]->SetMinimum(-0.025);}
		else if(!startdiff) {h_intercept_vs_time[e][nn]->SetMaximum(0.175); h_intercept_vs_time[e][nn]->SetMinimum(-0.1);  }
		h_intercept_vs_time[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_intercept_vs_time[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_intercept_vs_time[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_intercept_vs_time[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_intercept_vs_time[e][nn]->GetXaxis()->SetTimeDisplay(1); h_intercept_vs_time[e][nn]->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
		h_intercept_vs_time[e][nn]->GetXaxis()->SetLabelOffset(0.02); h_intercept_vs_time[e][nn]->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); 
		h_intercept_vs_time[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_intercept_vs_time[e][nn]->Draw();
		outputfile = h_intercept_vs_time[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_intercept_vs_time_scaled[e][nn]->GetName() );
		col->SetTitle(h_intercept_vs_time_scaled[e][nn]->GetTitle() );
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTitle("time");
		h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetTitle("intercept (scaled)");
		if(startdiff)       {h_intercept_vs_time_scaled[e][nn]->SetMaximum(5.); h_intercept_vs_time_scaled[e][nn]->SetMinimum(-2.);}
		if(normalized)      {h_intercept_vs_time_scaled[e][nn]->SetMaximum(5.); h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		else if(!startdiff) {h_intercept_vs_time_scaled[e][nn]->SetMaximum(7.); h_intercept_vs_time_scaled[e][nn]->SetMinimum(0.); }
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_intercept_vs_time_scaled[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTimeDisplay(1); h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetLabelOffset(0.02); h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); 
		h_intercept_vs_time_scaled[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_intercept_vs_time_scaled[e][nn]->Draw();
		outputfile = h_intercept_vs_time_scaled[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_time[e][nn]->GetName() );
		col->SetTitle(h_slope_vs_time[e][nn]->GetTitle() );
		h_slope_vs_time[e][nn]->GetXaxis()->SetTitle("time");
		h_slope_vs_time[e][nn]->GetYaxis()->SetTitle("slope");
		if(startdiff)       {h_slope_vs_time[e][nn]->SetMaximum(0.20); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time[e][nn]->SetMaximum(0.15); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time[e][nn]->SetMaximum(0.10); h_slope_vs_time[e][nn]->SetMinimum(0.);}
		h_slope_vs_time[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_slope_vs_time[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_slope_vs_time[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_slope_vs_time[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_slope_vs_time[e][nn]->GetXaxis()->SetTimeDisplay(1); h_slope_vs_time[e][nn]->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
		h_slope_vs_time[e][nn]->GetXaxis()->SetLabelOffset(0.02); h_slope_vs_time[e][nn]->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); 
		h_slope_vs_time[e][nn]->GetXaxis()->SetNdivisions(510, true);
		h_slope_vs_time[e][nn]->Draw();
		outputfile = h_slope_vs_time[e][nn]->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_time_scaled[e][nn]->GetName() );
		col->SetTitle(h_slope_vs_time_scaled[e][nn]->GetTitle() );
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTitle("time");
		h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetTitle("slope (scaled)");
		if(startdiff)       {h_slope_vs_time_scaled[e][nn]->SetMaximum(12.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		if(normalized)      {h_slope_vs_time_scaled[e][nn]->SetMaximum(10.); h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		else if(!startdiff) {h_slope_vs_time_scaled[e][nn]->SetMaximum(9.);  h_slope_vs_time_scaled[e][nn]->SetMinimum(0.);}
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTitleOffset(1.25); h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetLabelSize(0.027);
		h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetTitleOffset(1.25); h_slope_vs_time_scaled[e][nn]->GetYaxis()->SetLabelSize(0.027);
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTimeDisplay(1); h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
		h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetLabelOffset(0.02); h_slope_vs_time_scaled[e][nn]->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); 
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

	//save plots
	TFile *newFile;
	if(samelumi){
	if(startdiff)  newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples_diffwrt2011begin_sameLumi1p1fb-1.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples_normalized_sameLumi1p1fb-1.root", "RECREATE");
	else if(!startdiff) newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples_sameLumi1p1fb-1.root", "RECREATE");
	}
	else{
	if(startdiff)  newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples_diffwrt2011begin.root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples_normalized.root", "RECREATE");//mustd <-> muSIC
	else if(!startdiff) newFile = new TFile("FileOutputs/20121008_intercept_slopes_VS_time_VectorNTuples.root", "RECREATE");
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