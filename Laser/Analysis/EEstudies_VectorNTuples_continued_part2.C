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
/*
Double_t fline(Double_t *x, Double_t *par)
{
    if (par[0] + par[1]*x[0] > 0.1) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}*/

//this function averages a vector<float>
float avg(vector<float> x){

  float sum = 0;
  for (unsigned int i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}

//this function fits the correlation of VPT/PN(max)-VPT/PN(min) vs. mu_std (in eta bins), and 
//              plots the fit parameters vs. eta (fit parameters are slope and intercept)
void EEstudies_VectorNTuples_continued_part2(){

	bool normalized        = true;
	bool startdiff         = false;//only true if normalized == false
	int period             = 0;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11

	bool fitlinepluserr    = true;	//also prints the error on the VPTdiff vs. XXX (only for burnin / mustd)
	bool plotslopeintercept= true; //you want this to be true
	bool VPTdiff_vs_mustd  = true; //you probably want this to be true -> includes fit
	bool plot_av_rms       = false; //plot does not help much
	bool print_av_rms      = false; //instead just print out the numbers
	bool plotdistancetofit = false; //distance to fit
	bool plotonlyAbsEta    = true; //if true don't print EE+ EE- separately


	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	TF1 *fitup;
	TF1 *fitdown;
	//define all Graphs, Histos, etc.
	TGraph* h_mustd_vs_VPTPNdiff[16][2];
	TGraph* h_mustd_vs_VPTPNdiff_scaled[16][2];
	TH1D* h_mustd_vs_VPTPNdiff_distancefit[16][2];
	TH1D* h_mustd_vs_VPTPNdiff_scaled_distancefit[16][2];
	TH1D* h_intercept_vs_eta[2][2];
	TH1D* h_intercept_vs_eta_scaled[2][2];
	TH1D* h_slope_vs_eta[2][2];
	TH1D* h_slope_vs_eta_scaled[2][2];
	TH1D* h_RMS_intercept[2][2];//Put RMS and Av in title here
	TH1D* h_RMS_intercept_scaled[2][2];
	TH1D* h_RMS_slope[2][2];
	TH1D* h_RMS_slope_scaled[2][2];
	TGraph* a_mustd_vs_VPTPNdiff[8][2];
	TGraph* a_mustd_vs_VPTPNdiff_scaled[8][2];
	TH1D* a_mustd_vs_VPTPNdiff_distancefit[8][2];
	TH1D* a_mustd_vs_VPTPNdiff_scaled_distancefit[8][2];
	TH1D* a_intercept_vs_eta[1][2];
	TH1D* a_intercept_vs_eta_scaled[1][2];
	TH1D* a_slope_vs_eta[1][2];
	TH1D* a_slope_vs_eta_scaled[1][2];
	TH1D* a_RMS_intercept[1][2];//Put RMS and Av in title here
	TH1D* a_RMS_intercept_scaled[1][2];
	TH1D* a_RMS_slope[1][2];
	TH1D* a_RMS_slope_scaled[1][2];
	vector<double> removedX[16][2];
	vector<double> removedY[16][2];
	vector<double> removedscX[16][2];
	vector<double> removedscY[16][2];
	vector<double> removedaX[8][2];
	vector<double> removedaY[8][2];
	vector<double> removedscaX[8][2];
	vector<double> removedscaY[8][2];
	for(int n=0;n<16;++n){
	for(int m=0;m<2; ++m){
		removedX[n][m].clear();
		removedY[n][m].clear();
		removedscX[n][m].clear();
		removedscY[n][m].clear();
		if(n>=8) continue;
		removedaX[n][m].clear();
		removedaY[n][m].clear();
		removedscaX[n][m].clear();
		removedscaY[n][m].clear();
	}}

	TString periodTS; TString periods;
//	if     (period==1) {periodTS = "/period150311_270311"; periods = "_period150311_270311";}
//	else if(period==2) {periodTS = "/period120411_050511"; periods = "_period120411_050511";}
//	else if(period==3) {periodTS = "/period130511_010711"; periods = "_period130511_010711";}
//	else if(period==4) {periodTS = "/period150711_250811"; periods = "_period150711_250811";}
//	else if(period==5) {periodTS = "/period040911_021111"; periods = "_period040911_021111";}
//	else if(period==41){periodTS = "/period150711_080811_1p1fb-1"; periods = "_period150711_080811_1p1fb-1";}
//	else if(period==51){periodTS = "/period040911_230911_1p1fb-1"; periods = "_period040911_021111_1p1fb-1";}
	if     (period==1) {periodTS = "/period110315_110327"; periods = "_period110315_110327";}
	else if(period==2) {periodTS = "/period110412_110505"; periods = "_period110412_110505";}
	else if(period==3) {periodTS = "/period110513_110701"; periods = "_period110513_110701";}
	else if(period==4) {periodTS = "/period110715_110825"; periods = "_period110715_110825";}
	else if(period==5) {periodTS = "/period110904_111102"; periods = "_period110904_111102";}
	else if(period==6) {periodTS = "/period120406_120420"; periods = "_period120406_120420";}
	else if(period==7) {periodTS = "/period120505_120618"; periods = "_period120505_120618";}
	else if(period==41){periodTS = "/period110715_110808_1p1fb-1"; periods = "_period110715_110808_1p1fb-1";}
	else if(period==51){periodTS = "/period110904_110923_1p1fb-1"; periods = "_period110904_110923_1p1fb-1";}
	else if(period==32){periodTS = "/period110513_110624_950pb-1"; periods = "_period110513_110624_950pb-1";}
	else if(period==42){periodTS = "/period110715_110807_950pb-1"; periods = "_period110715_110807_950pb-1";}
	else if(period==52){periodTS = "/period110904_110922_950pb-1"; periods = "_period110904_110922_950pb-1";}
	else if(period==72){periodTS = "/period120505_110517_950pb-1"; periods = "_period120505_110517_950pb-1";}
	else if(period==-1){periodTS = "/period110822_110907"; periods = "_period110822_110907";}//technical stop Sept.2011
	else               {periodTS = "/originalperiod";      periods = "_originalperiod";}//also defined as 0, which is chosen to be the default

	//define eta and abs(eta) (same for caps and producer)
	string string_eta[16] = {"m3p0", "m2p8", "m2p6", "m2p4", "m2p2", "m2p0", "m1p8", "m1p6", "1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	string string_cap[2] = {"__EEm", "__EEp"};
	string string_aeta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_acap[1] = {"_"};

	//load files and tgraphs
	//new scaling in 20121001
	//old scaling in 20120829
	TFile *oldFile;
	if(startdiff) oldFile =  TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_diffwrt2011begin_VectorNTuples"+periods+".root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples"+periods+"temp.root");//mustd <-> muSIC
//	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples"+periods+"temp.root");//mustd <-> muSIC//xxyyzz
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_nonormalization_VectorNTuples"+periods+".root");
	TString outputdir;
	if(startdiff) outputdir = "tempPlots/20130402/EEstudies/old/continued/diffwrt2011begin"+periodTS+"/fitted";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
	if(normalized) outputdir = "tempPlots/20130402/EEstudies/old/continued/normalized"+periodTS+"/fitted";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
//	if(normalized) outputdir = "Plots/20130402/EEstudies/oldtemp/continued/normalized"+periodTS+"/fitted";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized//xxyyzz
	else if(!startdiff) outputdir = "tempPlots/20130402/EEstudies/old/continued/notnormalized"+periodTS+"/fitted";
	Util::MakeOutputDir(outputdir);

	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		string hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn];
		string hname_s = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn] + (string)"__scaled";
		string aname;
		if(n<8) aname = (string)"mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_aeta[n] + string_prod[nn];
		string aname_s;
		if(n<8) aname_s = (string)"mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_aeta[n] + string_prod[nn] + (string)"__scaled";

		h_mustd_vs_VPTPNdiff_scaled[n][nn] = (TGraph*)oldFile->Get(hname_s.c_str() );
		h_mustd_vs_VPTPNdiff[n][nn] = (TGraph*)oldFile->Get(hname.c_str() );

		if(n<8) a_mustd_vs_VPTPNdiff_scaled[n][nn] = (TGraph*)oldFile->Get(aname_s.c_str() );
		if(n<8) a_mustd_vs_VPTPNdiff[n][nn] = (TGraph*)oldFile->Get(aname.c_str() );
		//remove points here -> excluding them in fit, add back after fit
		/*for(int m = h_mustd_vs_VPTPNdiff[n][nn]->GetN()-1; m>=0; --m){
			double x,y;
			h_mustd_vs_VPTPNdiff[n][nn]->GetPoint(m,x,y);
			if(y>0.4 || (n==7&&nn==0&&y>0.1&&x>0.2&&x<0.8)){
				removedX[n][nn].push_back(x);
				removedY[n][nn].push_back(y);
				h_mustd_vs_VPTPNdiff[n][nn]->RemovePoint(m);
			}
			if(y<=0.){
				h_mustd_vs_VPTPNdiff[n][nn]->RemovePoint(m);
			}
		}
		for(int m = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetN()-1; m>=0; --m){
			double x,y;
			h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetPoint(m,x,y);
			if(y>2.5){
				removedscX[n][nn].push_back(x);
				removedscY[n][nn].push_back(y);
				h_mustd_vs_VPTPNdiff_scaled[n][nn]->RemovePoint(m);
			}
			if(y<=0.){
				h_mustd_vs_VPTPNdiff_scaled[n][nn]->RemovePoint(m);
			}
		}
		if(n<8){
		for(int m = a_mustd_vs_VPTPNdiff[n][nn]->GetN()-1; m>=0; --m){
			double x,y;
			a_mustd_vs_VPTPNdiff[n][nn]->GetPoint(m,x,y);
			if(y>0.4 || (n==7&&nn==0&&y>0.1&&x>0.2&&x<0.8)){
				removedaX[n][nn].push_back(x);
				removedaY[n][nn].push_back(y);
				a_mustd_vs_VPTPNdiff[n][nn]->RemovePoint(m);
			}
			if(y<=0.){
				a_mustd_vs_VPTPNdiff[n][nn]->RemovePoint(m);
			}
		}
		for(int m = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetN()-1; m>=0; --m){
			double x,y;
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetPoint(m,x,y);
			if(y>2.5){
				removedscaX[n][nn].push_back(x);
				removedscaY[n][nn].push_back(y);
				a_mustd_vs_VPTPNdiff_scaled[n][nn]->RemovePoint(m);
			}
			if(y<=0.){
				a_mustd_vs_VPTPNdiff_scaled[n][nn]->RemovePoint(m);
			}
		}
		}*/
		//define distance to fit histos
		hname = (string)"distancefit__" + hname;
		hname_s = (string)"distancefit__" + hname_s;
		h_mustd_vs_VPTPNdiff_distancefit[n][nn] = new TH1D(hname.c_str(), hname.c_str(),20, 0., 0.2);
		h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(),50, 0., 0.5);
		if(n<8){
			aname = (string)"distancefit__" + aname;
			aname_s = (string)"distancefit__" + aname_s;
			a_mustd_vs_VPTPNdiff_distancefit[n][nn] = new TH1D(aname.c_str(), aname.c_str(),20, 0., 0.2);
			a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn] = new TH1D(aname_s.c_str(), aname_s.c_str(),50, 0., 0.5);
			/*if(hname==(string)"mustd__vs__VPTPNmax-VPTPNmin__eta_m1p6__chin") {
			for(int x=0; x<= h_mustd_vs_VPTPNdiff[n][nn]->GetNbinsX()+1 ;++x){
			for(int y=0; y<= h_mustd_vs_VPTPNdiff[n][nn]->GetNbinsY()+1 ;++y){
				double yy = h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->GetBinCenter(y);
				if(yy>0.1 && h_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(x,y)!=0) h_mustd_vs_VPTPNdiff[n][nn]->SetBinContent(x,y,0.);
			}}}*/
		}

	}
	}
	//define intercept, slopes and RMS histos
	for(int n = 0; n<2; ++n){
	for(int nn = 0; nn<2; ++ nn){
		string hname = (string)"intercept_vs_eta" + string_cap[n] + string_prod[nn];
		string hname_s = (string)"intercept_vs_eta__scaled" + string_cap[n] + string_prod[nn];
		double el, eu;
		if(n==0){el = -3.0; eu = -1.4;}
		else    {el =  1.4; eu =  3.0;}
		h_intercept_vs_eta[n][nn] = new TH1D(hname.c_str(), hname.c_str(), 8, el, eu);
		h_intercept_vs_eta_scaled[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(), 8, el, eu);
		hname = (string)"slope_vs_eta" + string_cap[n] + string_prod[nn];
		hname_s = (string)"slope_vs_eta__scaled" + string_cap[n] + string_prod[nn];
		h_slope_vs_eta[n][nn] = new TH1D(hname.c_str(), hname.c_str(), 8, el, eu);
		h_slope_vs_eta_scaled[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(), 8, el, eu);
		if(normalized){
			h_intercept_vs_eta[n][nn]->SetMinimum(0.); h_intercept_vs_eta[n][nn]->SetMaximum(0.25);
			h_intercept_vs_eta_scaled[n][nn]->SetMinimum(0.); h_intercept_vs_eta_scaled[n][nn]->SetMaximum(0.75);
			h_slope_vs_eta[n][nn]->SetMinimum(0.); h_slope_vs_eta[n][nn]->SetMaximum(0.15);
			h_slope_vs_eta_scaled[n][nn]->SetMinimum(0.); h_slope_vs_eta_scaled[n][nn]->SetMaximum(0.75);
		}
		else{
//			h_intercept_vs_eta[n][nn]->SetMinimum(-0.0125); h_intercept_vs_eta[n][nn]->SetMaximum(0.0875);
//			h_intercept_vs_eta_scaled[n][nn]->SetMinimum(-0.0125); h_intercept_vs_eta_scaled[n][nn]->SetMaximum(0.0875);
//			h_slope_vs_eta[n][nn]->SetMinimum(0.); h_slope_vs_eta[n][nn]->SetMaximum(0.115);
//			h_slope_vs_eta_scaled[n][nn]->SetMinimum(0.); h_slope_vs_eta_scaled[n][nn]->SetMaximum(0.115);
		}
		h_intercept_vs_eta[n][nn]->GetYaxis()->SetTitle("intercept"); h_intercept_vs_eta[n][nn]->GetXaxis()->SetTitle("#eta");
		h_intercept_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("intercept"); h_intercept_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("#eta");
		h_slope_vs_eta[n][nn]->GetYaxis()->SetTitle("slope"); h_slope_vs_eta[n][nn]->GetXaxis()->SetTitle("#eta");
		h_slope_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("slope"); h_slope_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("#eta");

		h_intercept_vs_eta[n][nn]->SetLineWidth(4); h_intercept_vs_eta[n][nn]->SetMarkerStyle(20); h_intercept_vs_eta[n][nn]->SetMarkerSize(2);
		h_intercept_vs_eta_scaled[n][nn]->SetLineWidth(4); h_intercept_vs_eta_scaled[n][nn]->SetMarkerStyle(20); h_intercept_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		h_slope_vs_eta[n][nn]->SetLineWidth(4); h_slope_vs_eta[n][nn]->SetMarkerStyle(20); h_slope_vs_eta[n][nn]->SetMarkerSize(2);
		h_slope_vs_eta_scaled[n][nn]->SetLineWidth(4); h_slope_vs_eta_scaled[n][nn]->SetMarkerStyle(20); h_slope_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		if(nn==0){
		h_intercept_vs_eta[n][nn]->SetLineColor(kRed); h_intercept_vs_eta[n][nn]->SetMarkerColor(kRed);
		h_intercept_vs_eta_scaled[n][nn]->SetLineColor(kRed); h_intercept_vs_eta_scaled[n][nn]->SetMarkerColor(kRed);
		h_slope_vs_eta[n][nn]->SetLineColor(kRed); h_slope_vs_eta[n][nn]->SetMarkerColor(kRed);
		h_slope_vs_eta_scaled[n][nn]->SetLineColor(kRed); h_slope_vs_eta_scaled[n][nn]->SetMarkerColor(kRed);
		}
		else {
		h_intercept_vs_eta[n][nn]->SetLineColor(kBlue); h_intercept_vs_eta[n][nn]->SetMarkerColor(kBlue);
		h_intercept_vs_eta_scaled[n][nn]->SetLineColor(kBlue); h_intercept_vs_eta_scaled[n][nn]->SetMarkerColor(kBlue);
		h_slope_vs_eta[n][nn]->SetLineColor(kBlue); h_slope_vs_eta[n][nn]->SetMarkerColor(kBlue);
		h_slope_vs_eta_scaled[n][nn]->SetLineColor(kBlue); h_slope_vs_eta_scaled[n][nn]->SetMarkerColor(kBlue);
		}
		//RMS
		hname = (string)"intercept_RMS" + string_cap[n] + string_prod[nn];
		hname_s = (string)"intercept_RMS__scaled" + string_cap[n] + string_prod[nn];
		h_RMS_intercept[n][nn] = new TH1D(hname.c_str(), hname.c_str(), 25, 0, 0.25);
		h_RMS_intercept_scaled[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(), 75, 0, 0.75);
		hname = (string)"slope_RMS" + string_cap[n] + string_prod[nn];
		hname_s = (string)"slope_RMS__scaled" + string_cap[n] + string_prod[nn];
		h_RMS_slope[n][nn] = new TH1D(hname.c_str(), hname.c_str(), 15, 0, 0.15);
		h_RMS_slope_scaled[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(), 75, 0, 0.75);
		h_RMS_intercept[n][nn]->SetLineWidth(4); h_RMS_intercept[n][nn]->SetMarkerStyle(20); h_RMS_intercept[n][nn]->SetMarkerSize(2);
		h_RMS_intercept_scaled[n][nn]->SetLineWidth(4); h_RMS_intercept_scaled[n][nn]->SetMarkerStyle(20); h_RMS_intercept_scaled[n][nn]->SetMarkerSize(2);
		h_slope_vs_eta[n][nn]->SetLineWidth(4); h_slope_vs_eta[n][nn]->SetMarkerStyle(20); h_slope_vs_eta[n][nn]->SetMarkerSize(2);
		h_slope_vs_eta_scaled[n][nn]->SetLineWidth(4); h_slope_vs_eta_scaled[n][nn]->SetMarkerStyle(20); h_slope_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		h_intercept_vs_eta[n][nn]->GetYaxis()->SetTitle("intercept"); h_intercept_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		h_intercept_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("intercept"); h_intercept_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");
		h_slope_vs_eta[n][nn]->GetYaxis()->SetTitle("slope"); h_slope_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		h_slope_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("slope"); h_slope_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");
		if(nn==0){
		h_RMS_intercept[n][nn]->SetLineColor(kRed); h_RMS_intercept[n][nn]->SetMarkerColor(kRed);
		h_RMS_intercept_scaled[n][nn]->SetLineColor(kRed); h_RMS_intercept_scaled[n][nn]->SetMarkerColor(kRed);
		h_RMS_slope[n][nn]->SetLineColor(kRed); h_RMS_slope[n][nn]->SetMarkerColor(kRed);
		h_RMS_slope_scaled[n][nn]->SetLineColor(kRed); h_RMS_slope_scaled[n][nn]->SetMarkerColor(kRed);
		}
		else {
		h_RMS_intercept[n][nn]->SetLineColor(kBlue); h_RMS_intercept[n][nn]->SetMarkerColor(kBlue);
		h_RMS_intercept_scaled[n][nn]->SetLineColor(kBlue); h_RMS_intercept_scaled[n][nn]->SetMarkerColor(kBlue);
		h_RMS_slope[n][nn]->SetLineColor(kBlue); h_RMS_slope[n][nn]->SetMarkerColor(kBlue);
		h_RMS_slope_scaled[n][nn]->SetLineColor(kBlue); h_RMS_slope_scaled[n][nn]->SetMarkerColor(kBlue);
		}
		
		if(n>=1) continue;
		string aname = (string)"intercept_vs_AbsEta" + string_acap[n] + string_prod[nn];
		string aname_s = (string)"intercept_vs_AbsEta__scaled" + string_acap[n] + string_prod[nn];
		el=1.4;eu=3.0;
		a_intercept_vs_eta[n][nn] = new TH1D(aname.c_str(), aname.c_str(), 8, el, eu);
		a_intercept_vs_eta_scaled[n][nn] = new TH1D(aname_s.c_str(), aname_s.c_str(), 8, el, eu);
		aname = (string)"slope_vs_AbsEta" + string_acap[n] + string_prod[nn];
		aname_s = (string)"slope_vs_AbsEta__scaled" + string_acap[n] + string_prod[nn];
		a_slope_vs_eta[n][nn] = new TH1D(aname.c_str(), aname.c_str(), 8, el, eu);
		a_slope_vs_eta_scaled[n][nn] = new TH1D(aname_s.c_str(), aname_s.c_str(), 8, el, eu);
		if(normalized){
//			a_intercept_vs_eta[n][nn]->SetMinimum(0.); a_intercept_vs_eta[n][nn]->SetMaximum(0.25);
//			a_intercept_vs_eta_scaled[n][nn]->SetMinimum(0.); a_intercept_vs_eta_scaled[n][nn]->SetMaximum(0.75);
//			a_slope_vs_eta[n][nn]->SetMinimum(0.); a_slope_vs_eta[n][nn]->SetMaximum(0.15);
//			a_slope_vs_eta_scaled[n][nn]->SetMinimum(0.); a_slope_vs_eta_scaled[n][nn]->SetMaximum(0.75);
		}
		else{
//			a_intercept_vs_eta[n][nn]->SetMinimum(-0.0125); a_intercept_vs_eta[n][nn]->SetMaximum(0.0875);
//			a_intercept_vs_eta_scaled[n][nn]->SetMinimum(-0.0125); a_intercept_vs_eta_scaled[n][nn]->SetMaximum(0.0875);
//			a_slope_vs_eta[n][nn]->SetMinimum(0.); a_slope_vs_eta[n][nn]->SetMaximum(0.115);
//			a_slope_vs_eta_scaled[n][nn]->SetMinimum(0.); a_slope_vs_eta_scaled[n][nn]->SetMaximum(0.115);
		}
		a_intercept_vs_eta[n][nn]->GetYaxis()->SetTitle("intercept"); a_intercept_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_intercept_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("intercept"); a_intercept_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_slope_vs_eta[n][nn]->GetYaxis()->SetTitle("slope"); a_slope_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_slope_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("slope"); a_slope_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");

		a_intercept_vs_eta[n][nn]->SetLineWidth(4); a_intercept_vs_eta[n][nn]->SetMarkerStyle(20); a_intercept_vs_eta[n][nn]->SetMarkerSize(2);
		a_intercept_vs_eta_scaled[n][nn]->SetLineWidth(4); a_intercept_vs_eta_scaled[n][nn]->SetMarkerStyle(20); a_intercept_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		a_slope_vs_eta[n][nn]->SetLineWidth(4); a_slope_vs_eta[n][nn]->SetMarkerStyle(20); a_slope_vs_eta[n][nn]->SetMarkerSize(2);
		a_slope_vs_eta_scaled[n][nn]->SetLineWidth(4); a_slope_vs_eta_scaled[n][nn]->SetMarkerStyle(20); a_slope_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		if(nn==0){
		a_intercept_vs_eta[n][nn]->SetLineColor(kRed); a_intercept_vs_eta[n][nn]->SetMarkerColor(kRed);
		a_intercept_vs_eta_scaled[n][nn]->SetLineColor(kRed); a_intercept_vs_eta_scaled[n][nn]->SetMarkerColor(kRed);
		a_slope_vs_eta[n][nn]->SetLineColor(kRed); a_slope_vs_eta[n][nn]->SetMarkerColor(kRed);
		a_slope_vs_eta_scaled[n][nn]->SetLineColor(kRed); a_slope_vs_eta_scaled[n][nn]->SetMarkerColor(kRed);
		}
		else {
		a_intercept_vs_eta[n][nn]->SetLineColor(kBlue); a_intercept_vs_eta[n][nn]->SetMarkerColor(kBlue);
		a_intercept_vs_eta_scaled[n][nn]->SetLineColor(kBlue); a_intercept_vs_eta_scaled[n][nn]->SetMarkerColor(kBlue);
		a_slope_vs_eta[n][nn]->SetLineColor(kBlue); a_slope_vs_eta[n][nn]->SetMarkerColor(kBlue);
		a_slope_vs_eta_scaled[n][nn]->SetLineColor(kBlue); a_slope_vs_eta_scaled[n][nn]->SetMarkerColor(kBlue);
		}
		aname = (string)"intercept_RMS_Abs" + string_acap[n] + string_prod[nn];
		aname_s = (string)"intercept_RMS_Abs__scaled" + string_acap[n] + string_prod[nn];
		a_RMS_intercept[n][nn] = new TH1D(aname.c_str(), aname.c_str(), 25, 0, 0.25);
		a_RMS_intercept_scaled[n][nn] = new TH1D(aname_s.c_str(), aname_s.c_str(), 75, 0, 0.75);
		aname = (string)"slope_RMS_Abs" + string_acap[n] + string_prod[nn];
		aname_s = (string)"slope_RMS_Abs__scaled" + string_acap[n] + string_prod[nn];
		a_RMS_slope[n][nn] = new TH1D(aname.c_str(), aname.c_str(), 15, 0, 0.15);
		a_RMS_slope_scaled[n][nn] = new TH1D(aname_s.c_str(), aname_s.c_str(), 75, 0, 0.75);
		a_RMS_intercept[n][nn]->SetLineWidth(4); a_RMS_intercept[n][nn]->SetMarkerStyle(20); a_RMS_intercept[n][nn]->SetMarkerSize(2);
		a_RMS_intercept_scaled[n][nn]->SetLineWidth(4); a_RMS_intercept_scaled[n][nn]->SetMarkerStyle(20); a_RMS_intercept_scaled[n][nn]->SetMarkerSize(2);
		a_slope_vs_eta[n][nn]->SetLineWidth(4); a_slope_vs_eta[n][nn]->SetMarkerStyle(20); a_slope_vs_eta[n][nn]->SetMarkerSize(2);
		a_slope_vs_eta_scaled[n][nn]->SetLineWidth(4); a_slope_vs_eta_scaled[n][nn]->SetMarkerStyle(20); a_slope_vs_eta_scaled[n][nn]->SetMarkerSize(2);
		a_intercept_vs_eta[n][nn]->GetYaxis()->SetTitle("intercept"); a_intercept_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_intercept_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("intercept"); a_intercept_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_slope_vs_eta[n][nn]->GetYaxis()->SetTitle("slope"); a_slope_vs_eta[n][nn]->GetXaxis()->SetTitle("|#eta|");
		a_slope_vs_eta_scaled[n][nn]->GetYaxis()->SetTitle("slope"); a_slope_vs_eta_scaled[n][nn]->GetXaxis()->SetTitle("|#eta|");
		if(nn==0){
		a_RMS_intercept[n][nn]->SetLineColor(kRed); a_RMS_intercept[n][nn]->SetMarkerColor(kRed);
		a_RMS_intercept_scaled[n][nn]->SetLineColor(kRed); a_RMS_intercept_scaled[n][nn]->SetMarkerColor(kRed);
		a_RMS_slope[n][nn]->SetLineColor(kRed); a_RMS_slope[n][nn]->SetMarkerColor(kRed);
		a_RMS_slope_scaled[n][nn]->SetLineColor(kRed); a_RMS_slope_scaled[n][nn]->SetMarkerColor(kRed);
		}
		else {
		a_RMS_intercept[n][nn]->SetLineColor(kBlue); a_RMS_intercept[n][nn]->SetMarkerColor(kBlue);
		a_RMS_intercept_scaled[n][nn]->SetLineColor(kBlue); a_RMS_intercept_scaled[n][nn]->SetMarkerColor(kBlue);
		a_RMS_slope[n][nn]->SetLineColor(kBlue); a_RMS_slope[n][nn]->SetMarkerColor(kBlue);
		a_RMS_slope_scaled[n][nn]->SetLineColor(kBlue); a_RMS_slope_scaled[n][nn]->SetMarkerColor(kBlue);
		}
	}}

	//define fit functions
	TF1 *fit_mustd_vs_VPTPNdiff_scaled[16][2];
	TF1 *fit_mustd_vs_VPTPNdiff[16][2];
	TF1 *afit_mustd_vs_VPTPNdiff_scaled[8][2];
	TF1 *afit_mustd_vs_VPTPNdiff[8][2];
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		string hname = (string)"fit__mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn];
		string hname_s = (string)"fit__mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn] + (string)"__scaled";
		//if(hname==(string)"mustd__vs__VPTPNmax-VPTPNmin__eta_m1p6__chin")
		//	fit_mustd_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), fline, 0.1, 1.9, 2);
		//else
		fit_mustd_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-0.0)", 0.1, 1.9);
		fit_mustd_vs_VPTPNdiff_scaled[n][nn] = new TF1(hname_s.c_str(), "[0]+[1]*(x-0.0)", 0.1, 1.9);
		fit_mustd_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);
		fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetLineColor(kRed);
		//set fit parameteres if necessary
		/*if(normalized){
			if(n<2||n>13)      fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.07);//eta >2.6: x0=0.0 ->0.07; x0=0.2 ->0.05
			else if(n<4||n>11) fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.04);//eta >2.2: x0=0.0 ->0.04; x0=0.2 ->0.03
			else               fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.01);//eta >1.4: x0=0.0 ->0.01; x0=0.2 ->0.01
			if(n<2||n>13)      fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.09);
			else if(n<4||n>11) fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.05);
			else               fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.01);
		//	if(n<5||n>10)      fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.01);//eta>2.0: x0=0.0 ->; for x0=0.2 -> 0.01
		//	else               fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.02);//eta>1.4: x0=0.0 ->; for x0=0.2 -> 0.01
			if(n<6||n>9)       fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.5);//eta>1.8: x0=0.0 ->; for x0=0.2 -> 
			else               fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.3);//eta>1.4: x0=0.0 ->; for x0=0.2 -> 
			if(n<5||n>10)      fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.07);
			else               fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.07);
		}
		else{
			fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.015);
			fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.015);
			fit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.03);
			fit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.03);
		}*/
		if(n>=8) continue;
		string aname = (string)"fit__mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_aeta[n] + string_prod[nn];
		string aname_s = (string)"fit__mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_aeta[n] + string_prod[nn] + (string)"__scaled";
		//if(aname==(string)"mustd__vs__VPTPNmax-VPTPNmin__eta_m1p6__chin")
		//	afit_mustd_vs_VPTPNdiff[n][nn] = new TF1(aname.c_str(), fline, 0.1, 1.9, 2);
		//else
		afit_mustd_vs_VPTPNdiff[n][nn] = new TF1(aname.c_str(), "[0]+[1]*(x-0.0)", 0.1, 1.9);
		afit_mustd_vs_VPTPNdiff_scaled[n][nn] = new TF1(aname_s.c_str(), "[0]+[1]*(x-0.0)", 0.1, 1.9);
		afit_mustd_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);
		afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetLineColor(kRed);

//		if(nn==0&&(n==0||n==1)&&normalized){
//			afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.02);
//			afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.01);
//			afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.2);
//			afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 1.0);
//		}
		/*if(normalized){
			if(n>5)            afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.07);//eta >2.6: x0=0.0 ->0.07; x0=0.2 ->0.05
			else if(n>3)       afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.04);//eta >2.2: x0=0.0 ->0.04; x0=0.2 ->0.03
			else               afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0, 0.01);//eta >1.4: x0=0.0 ->0.01; x0=0.2 ->0.01
			if(n>5)            afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.09);
			else if(n>3)       afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.05);
			else               afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.01);
		//	if(n<5||n>10)      afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.01);//eta>2.0: x0=0.0 ->; for x0=0.2 -> 0.01
		//	else               afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.02);//eta>1.4: x0=0.0 ->; for x0=0.2 -> 0.01
			if(n>2)            afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.5);//eta>1.8: x0=0.0 ->; for x0=0.2 -> 
			else               afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(0, 0.7);//eta>1.4: x0=0.0 ->; for x0=0.2 -> 
			if(n>2)            afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.07);
			else               afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.5);
		}
		else{
			afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.015);
			afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.015);
			afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(1, 0.03);
			afit_mustd_vs_VPTPNdiff_scaled[n][nn]->SetParameter(1, 0.03);
		}*/
	}
	}

	//perform the fit
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(h_mustd_vs_VPTPNdiff[n][nn])->Fit(fit_mustd_vs_VPTPNdiff[n][nn], "R");
		(h_mustd_vs_VPTPNdiff_scaled[n][nn])->Fit(fit_mustd_vs_VPTPNdiff_scaled[n][nn], "R");
		//add removed points back here (some points lay such far away that they might spoil the fit)
		for(unsigned int m = 0; m<(removedX[n][nn]).size();++m){
			h_mustd_vs_VPTPNdiff[n][nn]->SetPoint(h_mustd_vs_VPTPNdiff[n][nn]->GetN(), (removedX[n][nn])[m], (removedY[n][nn])[m]);
		}
		for(unsigned int m = 0; m<(removedscX[n][nn]).size();++m){
			h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetN(), (removedscX[n][nn])[m], (removedscY[n][nn])[m]);
		}
	}}
	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
//		if(nn==0&&(n==0||n==1)&&normalized){
	//		cout << endl << endl << endl << endl;
//			cout << "For n="<<n<<" nn="<<nn << ": Parameters before fit: [0]=" <<afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << " [1]=" << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << endl;
//			cout << "       Scaled Parameters before fit: [0]=" <<afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << " [1]=" << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << endl;
//		}
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(a_mustd_vs_VPTPNdiff[n][nn])->Fit(afit_mustd_vs_VPTPNdiff[n][nn], "R");
		(a_mustd_vs_VPTPNdiff_scaled[n][nn])->Fit(afit_mustd_vs_VPTPNdiff_scaled[n][nn], "R");
		
//		if(nn==0&&(n==0||n==1)&&normalized){
//			cout << "For n="<<n<<" nn="<<nn << ": Parameters after  fit: [0]=" <<afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << " [1]=" << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << endl;
//			cout << "       Scaled Parameters after  fit: [0]=" <<afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << " [1]=" << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << endl;
//		}
		//add removed points back here
		for(unsigned int m = 0; m<(removedaX[n][nn]).size();++m){
			a_mustd_vs_VPTPNdiff[n][nn]->SetPoint(a_mustd_vs_VPTPNdiff[n][nn]->GetN(), (removedaX[n][nn])[m], (removedaY[n][nn])[m]);
		}
		for(unsigned int m = 0; m<(removedscaX[n][nn]).size();++m){
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetN(), (removedscaX[n][nn])[m], (removedscaY[n][nn])[m]);
		}
	}}

	//print out full fit parameters (including goodness of fit)
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		cout << "Fit " << h_mustd_vs_VPTPNdiff[n][nn]->GetName() << endl;
		cout << "Fit constant = " << fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << fit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare() << "/" << fit_mustd_vs_VPTPNdiff[n][nn]->GetNDF() << " = " << (fit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare())/(fit_mustd_vs_VPTPNdiff[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << fit_mustd_vs_VPTPNdiff[n][nn]->GetProb() << endl;
		cout << "now scaled" << endl;
		cout << "Fit constant = " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << " +/- " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1) << " +/- " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetChisquare() << "/" << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetNDF() << " = " << (fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetChisquare())/(fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetProb() << endl << endl;
	}}
	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		cout << "Fit " << a_mustd_vs_VPTPNdiff[n][nn]->GetName() << endl;
		cout << "Fit constant = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare() << "/" << afit_mustd_vs_VPTPNdiff[n][nn]->GetNDF() << " = " << (afit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare())/(afit_mustd_vs_VPTPNdiff[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << afit_mustd_vs_VPTPNdiff[n][nn]->GetProb() << endl;
		cout << "now scaled" << endl;
		cout << "Fit constant = " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0) << " +/- " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1) << " +/- " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetChisquare() << "/" << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetNDF() << " = " << (afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetChisquare())/(afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetProb() << endl << endl;

	}}

	//correlate fit variables (slopes and intercepts) vs eta bins
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		int k; if(n<8) k=0; else k = 1;
		double slope     = fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double slopeerr  = fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1);
		double slopes    = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1);
		double slopeserr = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1);
		double incep     = fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double inceperr  = fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0);
		double inceps    = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0);
		double incepserr = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0);
		int bin;
		if(k==0) bin = n+1;
		else     bin = n+1 - 8;
		h_intercept_vs_eta[k][nn]->SetBinContent(bin, incep);
		h_intercept_vs_eta[k][nn]->SetBinError(  bin, inceperr);
		h_intercept_vs_eta_scaled[k][nn]->SetBinContent(bin, inceps);
		h_intercept_vs_eta_scaled[k][nn]->SetBinError(  bin, incepserr);
		h_slope_vs_eta[k][nn]->SetBinContent(bin, slope);
		h_slope_vs_eta[k][nn]->SetBinError(  bin, slopeerr);
		h_slope_vs_eta_scaled[k][nn]->SetBinContent(bin, slopes);
		h_slope_vs_eta_scaled[k][nn]->SetBinError(  bin, slopeserr);
		//RMS
		h_RMS_intercept[k][nn]->Fill(incep);
		h_RMS_intercept_scaled[k][nn]->Fill(inceps);
		h_RMS_slope[k][nn]->Fill(slope);
		h_RMS_slope_scaled[k][nn]->Fill(slopes);
	}}
	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		int k; if(n<8) k=0; else k = 1;//always true
		double slope     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double slopeerr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1);
		double slopes    = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1);
		double slopeserr = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1);
		double incep     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double inceperr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0);
		double inceps    = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0);
		double incepserr = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0);
		int bin;
		if(k==0) bin = n+1;
		else     bin = n+1 - 8;
		a_intercept_vs_eta[k][nn]->SetBinContent(bin, incep);
		a_intercept_vs_eta[k][nn]->SetBinError(  bin, inceperr);
		a_intercept_vs_eta_scaled[k][nn]->SetBinContent(bin, inceps);
		a_intercept_vs_eta_scaled[k][nn]->SetBinError(  bin, incepserr);
		a_slope_vs_eta[k][nn]->SetBinContent(bin, slope);
		a_slope_vs_eta[k][nn]->SetBinError(  bin, slopeerr);
		a_slope_vs_eta_scaled[k][nn]->SetBinContent(bin, slopes);
		a_slope_vs_eta_scaled[k][nn]->SetBinError(  bin, slopeserr);
		//RMS
		a_RMS_intercept[k][nn]->Fill(incep);
		a_RMS_intercept_scaled[k][nn]->Fill(inceps);
		a_RMS_slope[k][nn]->Fill(slope);
		a_RMS_slope_scaled[k][nn]->Fill(slopes);

	}}

	//calculate the distance to the fit line for all points and add this to histograms
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		double distance;
		double slope     = fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double slopes    = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1);
		double incep     = fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double inceps    = fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0);
		double x0, y0;
		for(int p = 0; p< h_mustd_vs_VPTPNdiff[n][nn]->GetN(); ++p){
			//not h_mustd_vs_VPTPNdiff_scaled and h_mustd_vs_VPTPNdiff have same number of entries

			h_mustd_vs_VPTPNdiff[n][nn]->GetPoint(p, x0, y0);
			distance = fabs(slope*x0 - y0 + incep)/sqrt(slope*slope + 1.);
			h_mustd_vs_VPTPNdiff_distancefit[n][nn]->Fill(distance);

			h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetPoint(p, x0, y0);
			distance = fabs(slopes*x0 - y0 + inceps)/sqrt(slopes*slopes + 1.);
			h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Fill(distance);
		}

		if(n>=8) continue;
		slope     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
		slopes    = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1);
		incep     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
		inceps    = afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0);
		for(int p = 0; p< a_mustd_vs_VPTPNdiff[n][nn]->GetN(); ++p){
			//not a_mustd_vs_VPTPNdiff_scaled and a_mustd_vs_VPTPNdiff have same number of entries
			a_mustd_vs_VPTPNdiff[n][nn]->GetPoint(p, x0, y0);
			distance = fabs(slope*x0 - y0 + incep)/sqrt(slope*slope + 1.);
			a_mustd_vs_VPTPNdiff_distancefit[n][nn]->Fill(distance);
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetPoint(p, x0, y0);
			distance = fabs(slopes*x0 - y0 + inceps)/sqrt(slopes*slopes + 1.);
			a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Fill(distance);
		}

	}}

	//plot VPT/PN(max) - VPT/PN(min) vs mu_std in eta bins including the fit
	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	for(int n = 0; n<16; ++n){
	for(int nn= 0; nn<2;++nn){
		TH1F *haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		TString outputfile;
	     if(!plotonlyAbsEta){
	      if(VPTdiff_vs_mustd){
		if(fitlinepluserr){
			fitup = new TF1("fitup", "[0]+[1]*(x-0.0)", 0.1, 1.9);
			fitdown = new TF1("fitdown", "[0]+[1]*(x-0.0)", 0.1, 1.9);
			fitup  ->SetLineColor(kRed); fitup  ->SetLineStyle(3); fitup  ->SetLineWidth(2);
			fitdown->SetLineColor(kRed); fitdown->SetLineStyle(3); fitdown->SetLineWidth(2);
			fitup  ->SetParameter(0, fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)+fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
			fitdown->SetParameter(0, fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)-fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
			fitup  ->SetParameter(1, fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)+fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
			fitdown->SetParameter(1, fit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)-fit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
		}
		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
		haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetName(), h_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
		else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		fit_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		if(fitlinepluserr){
			fitup  ->Draw("same");
			fitdown->Draw("same");
		}
		//h_mustd_vs_VPTPNdiff[n][nn]->Draw();
		//h_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		outputfile = h_mustd_vs_VPTPNdiff[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();


		if(fitlinepluserr){
			fitup  ->SetParameter(0, fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0)+fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0));
			fitdown->SetParameter(0, fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0)-fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0));
			fitup  ->SetParameter(1, fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1)+fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1));
			fitdown->SetParameter(1, fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1)-fit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1));
		}
		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
		haxis = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,2.5);
		else           haxis->GetYaxis()->SetRangeUser(0.,4.0);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		fit_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		if(fitlinepluserr){
			fitup  ->Draw("same");
			fitdown->Draw("same");
		}
		//h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw();
		//h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		outputfile = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
		if(fitlinepluserr){
			delete fitup;
			delete fitdown;
		}
	      }
	      if(plotdistancetofit){
		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_distancefit[n][nn]->Draw();
		outputfile = h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle();
		col->Update();
		if(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Draw();
		outputfile = h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle();
		col->Update();
		if(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	      }
	     }
		if(n>=8) continue;
	     if(VPTdiff_vs_mustd){
		if(fitlinepluserr){
			fitup = new TF1("fitup", "[0]+[1]*(x-0.0)", 0.1, 1.9);
			fitdown = new TF1("fitdown", "[0]+[1]*(x-0.0)", 0.1, 1.9);
			fitup  ->SetLineColor(kRed); fitup  ->SetLineStyle(3); fitup  ->SetLineWidth(2);
			fitdown->SetLineColor(kRed); fitdown->SetLineStyle(3); fitdown->SetLineWidth(2);
			fitup  ->SetParameter(0, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)+afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
			fitdown->SetParameter(0, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)-afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
			fitup  ->SetParameter(1, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)+afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
			fitdown->SetParameter(1, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)-afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
		}
		col->cd();
		col->SetName(a_mustd_vs_VPTPNdiff[n][nn]->GetName() );
		col->SetTitle(a_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
		a_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
		haxis = a_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		haxis->SetNameTitle(a_mustd_vs_VPTPNdiff[n][nn]->GetName(), a_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
		else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		a_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		afit_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		if(fitlinepluserr){
			fitup  ->Draw("same");
			fitdown->Draw("same");
		}
		//a_mustd_vs_VPTPNdiff[n][nn]->Draw();
		//a_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		outputfile = a_mustd_vs_VPTPNdiff[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		if(fitlinepluserr){
			fitup  ->SetParameter(0, afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0)+afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0));
			fitdown->SetParameter(0, afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(0)-afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(0));
			fitup  ->SetParameter(1, afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1)+afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1));
			fitdown->SetParameter(1, afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParameter(1)-afit_mustd_vs_VPTPNdiff_scaled[n][nn]->GetParError(1));
		}
		col->cd();
		col->SetName(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
		col->SetTitle(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
		haxis = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
		haxis->SetNameTitle(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,2.5);
		else           haxis->GetYaxis()->SetRangeUser(0.,4.0);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		afit_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		if(fitlinepluserr){
			fitup  ->Draw("same");
			fitdown->Draw("same");
		}
		//a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw();
		//a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		outputfile = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(a_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
		if(fitlinepluserr){
			delete fitup;
			delete fitdown;
		}
	     }
	     if(plotdistancetofit){
		col->cd();
		col->SetName(a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetName() );
		col->SetTitle(a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle() );
		a_mustd_vs_VPTPNdiff_distancefit[n][nn]->Draw();
		outputfile = a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle();
		col->Update();
		if(a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetName() );
		col->SetTitle(a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle() );
		a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Draw();
		outputfile = a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle();
		col->Update();
		if(a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	     }
	}
	}
	//plot intercept and slopes vs. eta histograms
	for(int n = 0; n< 2; ++n){
	for(int nn= 0; nn<2;++nn){
		TString outputfile;
		char buffer[100]; 
		string buf;
	     if(!plotonlyAbsEta){
	      if(plotslopeintercept){
		col->cd();
		col->SetName(h_intercept_vs_eta[n][nn]->GetName() );
		col->SetTitle(h_intercept_vs_eta[n][nn]->GetTitle() );
		h_intercept_vs_eta[n][nn]->Draw();
		//h_intercept_vs_eta[n][nn]->Draw("same");
		outputfile = h_intercept_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_intercept_vs_eta[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(h_intercept_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(h_intercept_vs_eta_scaled[n][nn]->GetTitle() );
		h_intercept_vs_eta_scaled[n][nn]->Draw();
		//h_intercept_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = h_intercept_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_eta[n][nn]->GetName() );
		col->SetTitle(h_slope_vs_eta[n][nn]->GetTitle() );
		h_slope_vs_eta[n][nn]->Draw();
		//h_slope_vs_eta[n][nn]->Draw("same");
		outputfile = h_slope_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_slope_vs_eta[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(h_slope_vs_eta_scaled[n][nn]->GetTitle() );
		h_slope_vs_eta_scaled[n][nn]->Draw();
		//h_slope_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = h_slope_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_slope_vs_eta_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	      }
	      if(plot_av_rms){
		//plot average and RMS histos
		string title = h_RMS_slope_scaled[n][nn]->GetTitle();
		double average = h_RMS_slope_scaled[n][nn]->GetMean();
		double RMS = h_RMS_slope_scaled[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		h_RMS_slope_scaled[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(h_RMS_slope_scaled[n][nn]->GetName() );
		col->SetTitle(h_RMS_slope_scaled[n][nn]->GetTitle() );
		h_RMS_slope_scaled[n][nn]->Draw();
		//h_RMS_slope_scaled[n][nn]->Draw("same");
		outputfile = h_RMS_slope_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_RMS_slope_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = h_RMS_slope[n][nn]->GetTitle();
		average = h_RMS_slope[n][nn]->GetMean();
		RMS = h_RMS_slope[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		h_RMS_slope[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(h_RMS_slope[n][nn]->GetName() );
		col->SetTitle(h_RMS_slope[n][nn]->GetTitle() );
		h_RMS_slope[n][nn]->Draw();
		//h_RMS_slope[n][nn]->Draw("same");
		outputfile = h_RMS_slope[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_RMS_slope[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_RMS_slope[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_RMS_slope[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = h_RMS_intercept[n][nn]->GetTitle();
		average = h_RMS_intercept[n][nn]->GetMean();
		RMS = h_RMS_intercept[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		h_RMS_intercept[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(h_RMS_intercept[n][nn]->GetName() );
		col->SetTitle(h_RMS_intercept[n][nn]->GetTitle() );
		h_RMS_intercept[n][nn]->Draw();
		//h_RMS_intercept[n][nn]->Draw("same");
		outputfile = h_RMS_intercept[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_RMS_intercept[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = h_RMS_intercept_scaled[n][nn]->GetTitle();
		average = h_RMS_intercept_scaled[n][nn]->GetMean();
		RMS = h_RMS_intercept_scaled[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		h_RMS_intercept_scaled[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(h_RMS_intercept_scaled[n][nn]->GetName() );
		col->SetTitle(h_RMS_intercept_scaled[n][nn]->GetTitle() );
		h_RMS_intercept_scaled[n][nn]->Draw();
		//h_RMS_intercept_scaled[n][nn]->Draw("same");
		outputfile = h_RMS_intercept_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_RMS_intercept_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	      }
	      if(print_av_rms){
		string title = h_RMS_slope[n][nn]->GetTitle();
		double average = h_RMS_slope[n][nn]->GetMean();
		double RMS = h_RMS_slope[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average <<  endl;
		title = h_RMS_slope_scaled[n][nn]->GetTitle();
		average = h_RMS_slope_scaled[n][nn]->GetMean();
		RMS = h_RMS_slope_scaled[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
		title = h_RMS_intercept[n][nn]->GetTitle();
		average = h_RMS_intercept[n][nn]->GetMean();
		RMS = h_RMS_intercept[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
		title = h_RMS_intercept_scaled[n][nn]->GetTitle();
		average = h_RMS_intercept_scaled[n][nn]->GetMean();
		RMS = h_RMS_intercept_scaled[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
	      }
	     }
		if(n>=1) continue;
	     if(plotslopeintercept){
		col->cd();
		col->SetName(a_intercept_vs_eta[n][nn]->GetName() );
		col->SetTitle(a_intercept_vs_eta[n][nn]->GetTitle() );
		a_intercept_vs_eta[n][nn]->Draw();
		//a_intercept_vs_eta[n][nn]->Draw("same");
		outputfile = a_intercept_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_intercept_vs_eta[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(a_intercept_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(a_intercept_vs_eta_scaled[n][nn]->GetTitle() );
		a_intercept_vs_eta_scaled[n][nn]->Draw();
		//a_intercept_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = a_intercept_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(a_slope_vs_eta[n][nn]->GetName() );
		col->SetTitle(a_slope_vs_eta[n][nn]->GetTitle() );
		a_slope_vs_eta[n][nn]->Draw();
		//a_slope_vs_eta[n][nn]->Draw("same");
		outputfile = a_slope_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_slope_vs_eta[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		col->cd();
		col->SetName(a_slope_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(a_slope_vs_eta_scaled[n][nn]->GetTitle() );
		a_slope_vs_eta_scaled[n][nn]->Draw();
		//a_slope_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = a_slope_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_slope_vs_eta_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	     }
	     if(plot_av_rms){
		string title = a_RMS_slope_scaled[n][nn]->GetTitle();
		double average = a_RMS_slope_scaled[n][nn]->GetMean();
		double RMS = a_RMS_slope_scaled[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		a_RMS_slope_scaled[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(a_RMS_slope_scaled[n][nn]->GetName() );
		col->SetTitle(a_RMS_slope_scaled[n][nn]->GetTitle() );
		a_RMS_slope_scaled[n][nn]->Draw();
		//a_RMS_slope_scaled[n][nn]->Draw("same");
		outputfile = a_RMS_slope_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_RMS_slope_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = a_RMS_slope[n][nn]->GetTitle();
		average = a_RMS_slope[n][nn]->GetMean();
		RMS = a_RMS_slope[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		a_RMS_slope[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(a_RMS_slope[n][nn]->GetName() );
		col->SetTitle(a_RMS_slope[n][nn]->GetTitle() );
		a_RMS_slope[n][nn]->Draw();
		//a_RMS_slope[n][nn]->Draw("same");
		outputfile = a_RMS_slope[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_RMS_slope[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_RMS_slope[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_RMS_slope[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = a_RMS_intercept[n][nn]->GetTitle();
		average = a_RMS_intercept[n][nn]->GetMean();
		RMS = a_RMS_intercept[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		a_RMS_intercept[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(a_RMS_intercept[n][nn]->GetName() );
		col->SetTitle(a_RMS_intercept[n][nn]->GetTitle() );
		a_RMS_intercept[n][nn]->Draw();
		//a_RMS_intercept[n][nn]->Draw("same");
		outputfile = a_RMS_intercept[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_RMS_intercept[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();

		title = a_RMS_intercept_scaled[n][nn]->GetTitle();
		average = a_RMS_intercept_scaled[n][nn]->GetMean();
		RMS = a_RMS_intercept_scaled[n][nn]->GetRMS();
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		buf = string(buffer);
		title = title + buf;
		a_RMS_intercept_scaled[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(a_RMS_intercept_scaled[n][nn]->GetName() );
		col->SetTitle(a_RMS_intercept_scaled[n][nn]->GetTitle() );
		a_RMS_intercept_scaled[n][nn]->Draw();
		//a_RMS_intercept_scaled[n][nn]->Draw("same");
		outputfile = a_RMS_intercept_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(a_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(a_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(a_RMS_intercept_scaled[n][nn]->GetEntries()>3) col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	     }
	     if(print_av_rms){
		string title = a_RMS_slope[n][nn]->GetTitle();
		double average = a_RMS_slope[n][nn]->GetMean();
		double RMS = a_RMS_slope[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average <<  endl;
		title = a_RMS_slope_scaled[n][nn]->GetTitle();
		average = a_RMS_slope_scaled[n][nn]->GetMean();
		RMS = a_RMS_slope_scaled[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
		title = a_RMS_intercept[n][nn]->GetTitle();
		average = a_RMS_intercept[n][nn]->GetMean();
		RMS = a_RMS_intercept[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
		title = a_RMS_intercept_scaled[n][nn]->GetTitle();
		average = a_RMS_intercept_scaled[n][nn]->GetMean();
		RMS = a_RMS_intercept_scaled[n][nn]->GetRMS();
		cout << "histogram " << title << " has average = " << average << " and rms = " << RMS << "--> rms/av =" << RMS/average << endl;
	     }
	}}
	//store all histos, tgraphs, fits, etc.
	TFile *newFile;
	newFile = new TFile("FileOutputs/20130402_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_diffwrt2011begin"+periods+".root", "RECREATE");
	if(normalized) newFile = new TFile("FileOutputs/20130402_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted"+periods+".root", "RECREATE");//mustd <-> muSIC
//	if(normalized) newFile = new TFile("FileOutputs/20130402_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted"+periods+"temp.root", "RECREATE");//mustd <-> muSIC//xxyyzz
	else if(!startdiff) newFile = new TFile("FileOutputs/20140402_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_nonormalization"+periods+".root", "RECREATE");
	newFile->cd();
	for(int n = 0; n<16; ++n){ for(int nn= 0; nn<2;++nn) {
	h_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	fit_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	fit_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	h_mustd_vs_VPTPNdiff_distancefit[n][nn]->Write();
	h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Write();
	if(n>=8) continue;
	a_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	a_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	afit_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	afit_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	a_mustd_vs_VPTPNdiff_distancefit[n][nn]->Write();
	a_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Write();
	}}

	for(int n = 0; n< 2; ++n){
	for(int nn= 0; nn<2;++nn){
	h_intercept_vs_eta[n][nn]->Write(); 
	h_intercept_vs_eta_scaled[n][nn]->Write(); 
	h_slope_vs_eta[n][nn]->Write(); 
	h_slope_vs_eta_scaled[n][nn]->Write(); 
	h_RMS_intercept[n][nn]->Write(); 
	h_RMS_intercept_scaled[n][nn]->Write(); 
	h_RMS_slope[n][nn]->Write(); 
	h_RMS_slope_scaled[n][nn]->Write(); 
	if(n>=1) continue;
	a_intercept_vs_eta[n][nn]->Write(); 
	a_intercept_vs_eta_scaled[n][nn]->Write(); 
	a_slope_vs_eta[n][nn]->Write(); 
	a_slope_vs_eta_scaled[n][nn]->Write(); 
	a_RMS_intercept[n][nn]->Write(); 
	a_RMS_intercept_scaled[n][nn]->Write(); 
	a_RMS_slope[n][nn]->Write(); 
	a_RMS_slope_scaled[n][nn]->Write(); 
	}}

	cout << "histograms for period " << period <<" called "<< periods << " (normalized/startdiff = " << normalized << "/" << startdiff << ") stored in " << newFile->GetName() << endl;


}

