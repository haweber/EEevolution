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
#include "LaserAnalysisReducedTrees.C"
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

float avg(vector<float> x){

  float sum = 0;
  for (Int_t i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}

void EEstudies_continued_v2_part2(){

	bool normalized = true;

	gROOT->SetStyle("Plain");

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
	vector<double> removedX[16][2];
	vector<double> removedY[16][2];
	vector<double> removedscX[16][2];
	vector<double> removedscY[16][2];
	for(int n=0;n<16;++n){
	for(int m=0;m<2; ++m){
		removedX[n][m].clear();
		removedY[n][m].clear();
		removedscX[n][m].clear();
		removedscY[n][m].clear();
	}}

	string string_eta[16] = {"m3p0", "m2p8", "m2p6", "m2p4", "m2p2", "m2p0", "m1p8", "m1p6", "1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	string string_cap[2] = {"__EEm", "__EEp"};
	TFile *oldFile;
	//if(normalized) oldFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120312_mustd_vs_VPTPNmaxMinusVPTPNmin.root");
	if(normalized) oldFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_muSIC_vs_VPTPNmaxMinusVPTPNmin_savecopy_TESTEWRIAC_withTDR.root");//mustd <-> muSIC
	else oldFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120312_mustd_vs_VPTPNmaxMinusVPTPNmin_nonorm.root");
	TString outputdir;
	if(normalized) outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2c/SIC/normalized/fitted";//normalized if normalize VPT/PN // v2b/normalized <-> v2b/SIC/normalized
	else outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2b/notnormalized/fitted";
	Util::MakeOutputDir(outputdir);
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		string hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn];
		string hname_s = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn] + (string)"__scaled";

		h_mustd_vs_VPTPNdiff_scaled[n][nn] = (TGraph*)oldFile->Get(hname_s.c_str() );
		h_mustd_vs_VPTPNdiff[n][nn] = (TGraph*)oldFile->Get(hname.c_str() );
		//remove points here -> excluding them in fit, add back after fit
		for(int m = h_mustd_vs_VPTPNdiff[n][nn]->GetN()-1; m>=0; --m){
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
		hname = (string)"distancefit__" + hname;
		hname_s = (string)"distancefit__" + hname_s;
		h_mustd_vs_VPTPNdiff_distancefit[n][nn] = new TH1D(hname.c_str(), hname.c_str(),20, 0., 0.2);
		h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn] = new TH1D(hname_s.c_str(), hname_s.c_str(),50, 0., 0.5);

		/*if(hname==(string)"mustd__vs__VPTPNmax-VPTPNmin__eta_m1p6__chin") {
		for(int x=0; x<= h_mustd_vs_VPTPNdiff[n][nn]->GetNbinsX()+1 ;++x){
		for(int y=0; y<= h_mustd_vs_VPTPNdiff[n][nn]->GetNbinsY()+1 ;++y){
			double yy = h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->GetBinCenter(y);
			if(yy>0.1 && h_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(x,y)!=0) h_mustd_vs_VPTPNdiff[n][nn]->SetBinContent(x,y,0.);
		}}}*/

	}
	}
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
		hname = (string)"slope_RMS" + string_cap[n] + string_prod[nn];
		hname_s = (string)"slope_RMS__scaled" + string_cap[n] + string_prod[nn];
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
	}}

	TF1 *fit_mustd_vs_VPTPNdiff_scaled[16][2];
	TF1 *fit_mustd_vs_VPTPNdiff[16][2];
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
	}
	}

	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()<=4) continue;
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(h_mustd_vs_VPTPNdiff[n][nn])->Fit(fit_mustd_vs_VPTPNdiff[n][nn], "R");
		(h_mustd_vs_VPTPNdiff_scaled[n][nn])->Fit(fit_mustd_vs_VPTPNdiff_scaled[n][nn], "R");
		//add removed points back here
		for(int m = 0; m<(removedX[n][nn]).size();++m){
			h_mustd_vs_VPTPNdiff[n][nn]->SetPoint(h_mustd_vs_VPTPNdiff[n][nn]->GetN(), (removedX[n][nn])[m], (removedY[n][nn])[m]);
		}
		for(int m = 0; m<(removedscX[n][nn]).size();++m){
			h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetN(), (removedscX[n][nn])[m], (removedscY[n][nn])[m]);
		}
		
	}}

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
	}}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	for(int n = 0; n<16; ++n){
	for(int nn= 0; nn<2;++nn){
		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
		TH1F *haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetName(), h_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
		else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		fit_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		//h_mustd_vs_VPTPNdiff[n][nn]->Draw();
		//h_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
		TString outputfile = h_mustd_vs_VPTPNdiff[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
		haxis = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,2.5);
		else if(n<=4 && n>=11) haxis->GetYaxis()->SetRangeUser(0.,0.05);
		else  haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->Draw();
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		fit_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		//h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw();
		//h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("same");
		outputfile = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		//if(h_mustd_vs_VPTPNdiff[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_distancefit[n][nn]->Draw();
		outputfile = h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetTitle();
		col->Update();
		if(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Draw();
		outputfile = h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetTitle();
		col->Update();
		if(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
	}
	}
	for(int n = 0; n< 2; ++n){
	for(int nn= 0; nn<2;++nn){
		col->cd();
		col->SetName(h_intercept_vs_eta[n][nn]->GetName() );
		col->SetTitle(h_intercept_vs_eta[n][nn]->GetTitle() );
		h_intercept_vs_eta[n][nn]->Draw();
		//h_intercept_vs_eta[n][nn]->Draw("same");
		TString outputfile = h_intercept_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_intercept_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_intercept_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(h_intercept_vs_eta_scaled[n][nn]->GetTitle() );
		h_intercept_vs_eta_scaled[n][nn]->Draw();
		//h_intercept_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = h_intercept_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_intercept_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_eta[n][nn]->GetName() );
		col->SetTitle(h_slope_vs_eta[n][nn]->GetTitle() );
		h_slope_vs_eta[n][nn]->Draw();
		//h_slope_vs_eta[n][nn]->Draw("same");
		outputfile = h_slope_vs_eta[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_slope_vs_eta[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_slope_vs_eta_scaled[n][nn]->GetName() );
		col->SetTitle(h_slope_vs_eta_scaled[n][nn]->GetTitle() );
		h_slope_vs_eta_scaled[n][nn]->Draw();
		//h_slope_vs_eta_scaled[n][nn]->Draw("same");
		outputfile = h_slope_vs_eta_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_slope_vs_eta_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		string title = h_RMS_slope_scaled[n][nn]->GetTitle();
		double average = h_RMS_slope_scaled[n][nn]->GetMean();
		double RMS = h_RMS_slope_scaled[n][nn]->GetRMS();
		char buffer[100]; 
		sprintf(buffer,"_Av_%f_RMS_%f",average,RMS);
		string buf = string(buffer);
		title = title + buf;
		h_RMS_slope_scaled[n][nn]->SetTitle(title.c_str());
		col->cd();
		col->SetName(h_RMS_slope_scaled[n][nn]->GetName() );
		col->SetTitle(h_RMS_slope_scaled[n][nn]->GetTitle() );
		h_RMS_slope_scaled[n][nn]->Draw();
		//h_RMS_slope_scaled[n][nn]->Draw("same");
		outputfile = h_RMS_slope_scaled[n][nn]->GetTitle();//TString("mustd__vs__VPTPNmax-VPTPNmin");
		col->Update();
		if(h_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_RMS_slope_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
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
		if(h_RMS_slope[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_RMS_slope[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
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
		if(h_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_RMS_intercept[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
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
		if(h_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_RMS_intercept_scaled[n][nn]->GetEntries()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

	}}
	TFile *newFile;
	if(normalized) newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_muSIC_vs_VPTPNmaxMinusVPTPNmin_savecopy_TESTEWRIAC_withTDR_withfits.root", "RECREATE");//mustd <-> muSIC
	else newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_mustd_vs_VPTPNmaxMinusVPTPNmin_withfits_nonorm.root", "RECREATE");
	newFile->cd();
	for(int n = 0; n<16; ++n){ for(int nn= 0; nn<2;++nn) {
	h_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	fit_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	fit_mustd_vs_VPTPNdiff_scaled[n][nn]->Write();
	h_mustd_vs_VPTPNdiff_distancefit[n][nn]->Write();
	h_mustd_vs_VPTPNdiff_scaled_distancefit[n][nn]->Write();
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
	}}

}

