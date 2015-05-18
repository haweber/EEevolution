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
#include <TGraphErrors.h>
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
#include <TGraphErrors.h>

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
//#include "LaserAnalysisVectorTrees.C"
#include "RootMacros/Utilities.hh"


using namespace std;

void EEstudies_VectorNTuples_burnin_continued_daily(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool plotoverview   = true;
	bool savefinalplots = true;


	TString outputdir = "Plots/20121119/EEstudies/continued/daily/burnin";//normalized";
	Util::MakeOutputDir(outputdir);


	TFile *stdmapsfile     = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *EEp_producer     = (TH2D*)stdmapsfile->Get("EEp_producer");//1 russian, 2 chinese
	TH2D *EEm_producer     = (TH2D*)stdmapsfile->Get("EEm_producer");

	TFile *burninfile      = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEburnin.root");
	TH2D *burnin_EEm       = (TH2D*)burninfile->Get("burnin_EEm");
	TH2D *burnin_EEp       = (TH2D*)burninfile->Get("burnin_EEp");

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	TGraphErrors* g_slopes[8][2];
	TGraphErrors* g_onevalues[8][2];
	bool usetgrapherrors = true;
	TGraphErrors* a_mustd_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	//TGraph* a_mustd_vs_VPTPNdiff[8][2];//usetgrapherrors = false;
	TF1 *afit_mustd_vs_VPTPNdiff[8][2];
	TLegend *Legend1 = new TLegend(.68,.5,.88,.88);
	Legend1 -> SetFillColor(0);
	Legend1 -> SetBorderSize(0);
   	TH1F* h[8]; //needed for legend

	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	string string_eta_leg[8] = {"1.4 #leq |#eta| < 1.6", "1.6 #leq |#eta| < 1.8", "1.8 #leq |#eta| < 2.0", "2.0 #leq |#eta| < 2.2", "2.2 #leq |#eta| < 2.4", "2.4 #leq |#eta| < 2.6", "2.6 #leq |#eta| < 2.8", "2.8 #leq |#eta|"};

	
	for(int n = 0; n<8; ++n){//eta bin
		h[n] = new TH1F("", (string_eta_leg[n]).c_str(), 1, 0, 1); 
	for(int nn = 0; nn<2; ++ nn){//producer
		string hname;
	
		hname = (string)"slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_slopes[n][nn] = new TGraphErrors();
		g_slopes[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
		g_slopes[n][nn]->GetYaxis()->SetTitle("burn-in slope"); g_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_slopes[n][nn]->SetMarkerSize(1); g_slopes[n][nn]->SetMarkerStyle(20);

		hname = (string)"onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_onevalues[n][nn] = new TGraphErrors();
		g_onevalues[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
		g_onevalues[n][nn]->GetYaxis()->SetTitle("burn-in=1 equivalent"); g_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_onevalues[n][nn]->SetMarkerSize(1); g_onevalues[n][nn]->SetMarkerStyle(20);

		//set colors for plotting
		Color_t color;
		if(n==0) color = kViolet;
		if(n==1) color = kBlue;
		if(n==2) color = kCyan+2;
		if(n==3) color = kGreen+2;
		if(n==4) color = kYellow+3;
		if(n==5) color = kOrange+5;
		if(n==6) color = kOrange;
		if(n==7) color = kRed;
		g_slopes[n][nn]     ->SetLineColor(color);   g_slopes[n][nn]     ->SetMarkerColor(color);
		g_onevalues[n][nn] ->SetLineColor(color);   g_onevalues[n][nn] ->SetMarkerColor(color);
		h[n]                ->SetFillColor(color);


	}}

	for(int n = 0; n<8; ++n){
		Legend1->AddEntry(h[n], h[n]->GetTitle(), "f");
	}

	//TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121023_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121116_EE_VPT_over_PN_all_VectorNTuples_2011and2012_pointperday.root");
	//TGraph* g_vptpn_las[101][101][2];
	TGraphErrors* g_vptpn_las[101][101][2];
	char gname[101];


	bool badsth[101][101][2];



	cout << "start computation" << endl;
	cout << "first cleaning and normalization for all crystals" << endl;
	int numberofdays = 0;
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		badsth[i][j][k] = false;
		bool risingcrystals = (k==1&&((i==95&&j==67)||(i==95&&j==66)||(i==94&&j==27)||(i==93&&j==34)||(i==92&&j==27)||(i==79&&j==92)||(i==79&&j==91)||(i==77&&j==91)||(i==69&&j==95)||(i==54&&j==99)||(i==54&&j==98)||(i==53&&j==100)||(i==52&&j==100)||(i==38&&j==97)||(i==36&&j==97)||(i==29&&j==95)||(i==27&&j==6)||(i==25&&j==92)||(i==24&&j==91)||(i==19&&j==83)||(i==17&&j==84)||(i==17&&j==83)||(i==16&&j==82)||(i==10&&j==30)||(i==8&&j==27)||(i==5&&j==62)||(i==4&&j==64)||(i==4&&j==63)||(i==4&&j==50)))||(i==87&&j==25) || (k==0&&((i==99&&j==57)||(i==97&&j==59)||(i==97&&j==38)||(i==96&&j==42)||(i==96&&j==40)||(i==96&&j==37)||(i==94&&j==73)||(i==92&&j==27)||(i==84&&j==15)||(i==80&&j==15)||(i==80&&j==13)||(i==80&&j==12)||(i==80&&j==9)||(i==78&&j==91)||(i==78&&j==10)||(i==77&&j==9)||(i==76&&j==11)||(i==76&&j==9)||(i==63&&j==96)||(i==61&&j==5)||(i==51&&j==3)||(i==29&&j==8)||(i==24&&j==9)||(i==15&&j==25)||(i==15&&j==19)||(i==15&&j==16)||(i==14&&j==85)||(i==14&&j==82)||(i==14&&j==19)||(i==10&&j==77)||(i==10&&j==27)||(i==10&&j==24)||(i==8&&j==28)||(i==8&&j==27)||(i==7&&j==28)));
		//these crystals are crystals whose VPT/PN rises even during some time of radiation, but not always
		bool partiallyrising= (k==1&&((i==94&&j==70)||(j==91&&j==71)||(i==91&&j==66)||(i==91&&j==27)||(i==80&&j==91)||(i==76&&j==92)||(i==76&&j==91)||(i==70&&j==91)||(i==65&&j==97)||(i==55&&j==97)||(i==52&&j==97)||(i==36&&j==96)||(i==30&&j==8)||(i==28&&j==10)||(i==28&&j==9)||(i==25&&j==91)||(i==22&&j==91)||(i==19&&j==39))) || (k==0&&((i==96&&j==60)||(i==92&&j==70)||(i==85&&j==14)||(i==84&&j==14)||(i==76&&j==12)||(i==72&&j==93)||(i==70&&j==95)||(i==64&&j==5)||(i==61&&j==4)||(i==55&&j==15)||(i==55&&j==3)||(i==31&&j==10)||(i==26&&j==7)||(i==25&&j==18)||(i==15&&j==20)||(i==12&&j==21)||(i==8&&j==34)));
		//these crystals are crystals whose VPT/PN rises even during beginning 2011/2012, but then loss dominates
		bool partiallyrising2 = (k==1&&((i==41&&j==65)||(i==21&&j==91))) || (k==0&&((i==93&&j==69)||(i==32&&j==73)||(i==14&&j==83)));
		bool bad2012 = (k==0&&((i==14&&j==50)||(i==18&&j==21)||(i==24&&j==51)||(i==25&&j==64)||(i==26&&j==50)||(i==31&&j==38)||(i==31&&j==72)||(i==32&&j==32)||(i==32&&j==56)||(i==32&&j==61)||(i==32&&j==62)||(i==32&&j==73)||(i==33&&j==57)||(i==33&&j==75)||(i==38&&j==15)||(i==38&&j==40)||(i==39&&j==82)||(i==44&&j==13)||(i==49&&j==82)||(i==52&&j==2)||(i==70&&j==50)||(i==91&&j==76))) || (k==1&&((i==54&&j==36)||(i==56&&j==2)||(i==58&&j==38)||(i==59&&j==38)||(i==60&&j==39)||(i==62&&j==51)||(i==65&&j==27)||(i==70&&j==32)||(i==70&&j==33)||(i==81&&j==35)));
		bool badsthsingle = (k==0&&((i==9&&j==33)||(i==47&&j==78)))||(k==1&&( (i==39&&j==46)||(i==52&&j==4)||(i==78&&j==76)||(i==83&&j==21)||(i==84&&j==24)||(i==86&&j==24)||(i==94&&j==64)||(i==19&&j==84)||(i==24&&j==77) ) )||(k==1&&(i==91||i==92)&&(j>=21&&j<=25) )||bad2012||risingcrystals || partiallyrising || partiallyrising2;

		badsth[i][j][k] = badsthsingle;

		//if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEm", i, j);
		//else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEp", i, j);
		//g_vptpn_las[i][j][k] = (TGraph*)vpt_values->Get(gname);
		if (k==0) sprintf(gname,"tg_vptpn_las_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"tg_vptpn_las_ix%d_iy%d_EEp", i, j);
		g_vptpn_las[i][j][k] = (TGraphErrors*)vpt_values->Get(gname);

		if(g_vptpn_las[i][j][k]->GetN()>0){ // start normalization
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y, xprev,yprev, xpost,ypost, xprev2,yprev2;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(n==g_vptpn_las[i][j][k]->GetN()-1) g_vptpn_las[i][j][k]->GetPoint(n,xpost,ypost);
				else g_vptpn_las[i][j][k]->GetPoint(n+1,xpost,ypost);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev,yprev);
				else g_vptpn_las[i][j][k]->GetPoint(n-1,xprev,yprev);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev2,yprev2);
				else if(n==1) g_vptpn_las[i][j][k]->GetPoint(n-1,xprev2,yprev2);
				else g_vptpn_las[i][j][k]->GetPoint(n-2,xprev2,yprev2);
				double ctu = 1.15; double ctd = 0.85;//c-leaning t-hreshold u-p/d-own
				if((yprev>y*ctu && ypost>y*ctu) || (yprev<y*ctd && ypost<y*ctd)) g_vptpn_las[i][j][k]->SetPoint(n,x,yprev);
				//upper 2cleaning
				if(ypost<y*1.1 && yprev2<y*ctu && ypost<yprev*ctu && yprev2<yprev*ctu){
					g_vptpn_las[i][j][k]->SetPoint(n,x,yprev2);//clean only n
					if(n>0) g_vptpn_las[i][j][k]->SetPoint(n-1,xprev,yprev2);//clean only n
				}
			}
			double norm = -1;
			for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(y==0) continue;
				if(x<1295000000) continue;//safety cut
				norm = y;
				if(norm>0) break;
			}
			if(norm<0) {cout << "error - no norm " << norm << " of " << g_vptpn_las[i][j][k]->GetN() << " point at ijk "<<i<<" "<<j<<" " <<k << endl; norm = 1.; }
			if(norm!=1. && norm!=0.){
			for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(y==0) continue;
				if(x<1295000000) continue;//safety cut
				g_vptpn_las[i][j][k]->SetPoint(n,x,y/norm);
				double ex = g_vptpn_las[i][j][k]->GetErrorX(n);
				double ey = g_vptpn_las[i][j][k]->GetErrorY(n);
				g_vptpn_las[i][j][k]->SetPointError(n,ex/norm,ey/norm);

			}
			}
			double yprev = -1.; int sameyprev = 0;
			for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(y==0) continue;
				if(x<1295000000) continue;//safety cut
				if(fabs(y-yprev)<=10e-10) ++sameyprev;
				else {sameyprev=0; yprev = y;}
				if(sameyprev>29 && yprev!=1.){
					cout << "veto crystal " << i<<":"<<j<<":"<<k<<" because it has " << sameyprev << " times the same value " << yprev << " in a row, last at t= " << int(x) << endl;
					badsth[i][j][k] = true;
					break;
				}
			}

		}

		if(badsth[i][j][k]) continue;
		//cout << "x y z " << i << " " << j << " " << k << " has " << g_vptpn_las[i][j][k]->GetN() << " entries " << endl;
		if(g_vptpn_las[i][j][k]->GetN()>numberofdays) numberofdays = g_vptpn_las[i][j][k]->GetN();
	}}}//ijk
	//now start the fitting
	cout << "now start the fitting of burnin vs VPTdiff --> slopes onevalues " << endl;
	int oneday = 86400;
	for(int d = /*1298688352-100*/1293836400; d<1356994800; d+=86400){//from 1st Jan. of 2011 up to 31st of Dec. of 2012
	if((d-1293836400+1293840000)%(25*86400)==0) cout << "day " << int((d-1293836400)/(86400)) << endl;

	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		gStyle->SetOptStat(0);
		string hname;
		//create here
		hname = (string)"burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_mustd_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
		a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.7,1.2);
		a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in ratio");

		//create here
		hname = (string)"fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_mustd_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-1.0)", 0.75, 1.15);
		afit_mustd_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()!=0) cout << "Error graph " << n << " " << nn << "  has " << a_mustd_vs_VPTPNdiff[n][nn]->GetN() << " entries, but should have 0" << endl;

	}}
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		if(badsth[i][j][k]) continue;//kill bad crystals
		if(g_vptpn_las[i][j][k]->GetN()==0 ) continue;
		double x0,y0; bool runfurther = false; double x(-1), y(-1); double ex,ey;
		for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
			g_vptpn_las[i][j][k]->GetPoint(n,x0,y0);
			if(x0<d) continue;
			if(x0>d && x0<(d+oneday)){ 
				ex = g_vptpn_las[i][j][k]->GetErrorX(n);
				ey = g_vptpn_las[i][j][k]->GetErrorY(n);
				x = x0; y = y0; runfurther = true; }
			if(runfurther) break;
		}
		if(!runfurther) continue;
		//now point passed the day selection

		//get mustd
		double mustd = 0; int prod = -1;
		int bin; float eta; string seta; stringstream russchin; string rs; stringstream etass; string es;
		//get mu_std and eta for each crystals
		if(k==0){//choose xx_EEm_yy histos
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = EEm_producer->FindBin(i,j);
			prod = EEm_producer->GetBinContent(bin);
			if(prod==0) continue;
			if(prod==1) mustd = burnin_EEm->GetBinContent(bin);
			if(prod==2) mustd = burnin_EEm->GetBinContent(bin);
			russchin << mustd;
			rs = russchin.str();
			if(prod==1) seta = (string)"_eta_" + es + (string)"_russian_burnin_" + rs;
			if(prod==2) seta = (string)"_eta_" + es + (string)"_chinese_burnin_" + rs;
		}
		else{
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = EEp_producer->FindBin(i,j);
			prod = EEp_producer->GetBinContent(bin);
			if(prod==0) continue;
			if(prod==1) mustd = burnin_EEp->GetBinContent(bin);
			if(prod==2) mustd = burnin_EEp->GetBinContent(bin);
			russchin << mustd;
			rs = russchin.str();
			if(prod==1) seta = (string)"_eta_" + es + (string)"_russian_burnin_" + rs;
			if(prod==2) seta = (string)"_eta_" + es + (string)"_chinese_burnin_" + rs;
		}

		
		double vptdiff = 1.-y;
		int etaindex = -1;
		if(fabs(eta)<1.6&&fabs(eta)>=1.4) etaindex = 0;
		else if(fabs(eta)<( 1.8))         etaindex = 1;
		else if(fabs(eta)<( 2.0))         etaindex = 2;
		else if(fabs(eta)<( 2.2))         etaindex = 3;
		else if(fabs(eta)<( 2.4))         etaindex = 4;
		else if(fabs(eta)<( 2.6))         etaindex = 5;
		else if(fabs(eta)<( 2.8))         etaindex = 6;
		else if(fabs(eta)<=(3.0))         etaindex = 7;
		if(etaindex>=0 && etaindex <=7 && mustd!=0){
			//if(i>95&&k==1) cout<<"ijk "<<i<<" "<<j<<" "<<k<<" has eta "<<etaindex<<" prod-1 "<<prod-1<<" mustd "<<mustd<<" vptdiff "<<vptdiff<<endl;
			int graphnbins = a_mustd_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
			a_mustd_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, mustd, vptdiff);
			a_mustd_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
		}

	}}}//ijk
	//now run over all endcaps for one day
	//now make fit
	
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		//fit now
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(a_mustd_vs_VPTPNdiff[n][nn])->Fit(afit_mustd_vs_VPTPNdiff[n][nn], "RQ");//quiet
		/*
		cout << "Fit " << a_mustd_vs_VPTPNdiff[n][nn]->GetName() << endl;
		cout << "Fit constant = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << afit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare() << "/" << afit_mustd_vs_VPTPNdiff[n][nn]->GetNDF() << " = " << (afit_mustd_vs_VPTPNdiff[n][nn]->GetChisquare())/(afit_mustd_vs_VPTPNdiff[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << afit_mustd_vs_VPTPNdiff[n][nn]->GetProb() << endl;
		*/
		double slope     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double slopeerr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1);
		double incep     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double inceperr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0);
		int slopebin = g_slopes[n][nn]    ->GetN();
		int incepbin = g_onevalues[n][nn]->GetN();
		//check
		if(slopebin!=incepbin) cout << __LINE__ << "  slopebin " << slopebin << "  incepbin " << incepbin << endl;
		g_slopes[n][nn]    ->SetPoint(     slopebin, d, slope   );
		g_slopes[n][nn]    ->SetPointError(slopebin, 0, slopeerr);
		g_onevalues[n][nn]->SetPoint(     incepbin, d, incep   );
		g_onevalues[n][nn]->SetPointError(incepbin, 0, inceperr);
		//if((d-1293836400+1293840000)%(25*86400)==0) cout << "d " << d << " has slope " << slope << "+/-" << slopeerr<< " onevalue " << incep<<"+/-"<<inceperr<< " at etaind/prodind " << n << "/" << nn << " check " << g_slopes[n][nn]->GetN() << " " << g_onevalues[n][nn]->GetN() << " " << endl;

	}}

	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		//delete now
		delete a_mustd_vs_VPTPNdiff[n][nn];
		delete afit_mustd_vs_VPTPNdiff[n][nn];

	}}
	
	}// int d

	if(plotoverview){
	//now do plotting and saving
	TString outname;
	double slopemax(-99), incepmax(-99), slopemin(999), incepmin(999);
	double timemin(10e16), timemax(-1);
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		for(int nnn=0; nnn<g_slopes[n][nn]->GetN();++nnn){
			double dx,dy;
			g_slopes[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>slopemax) slopemax = dy;
			if(dy<slopemin) slopemin = dy;
			if(dx>timemax ) timemax  = dx;
			if(dx<timemin ) timemin  = dx;
			//cout << "slope " << n << " " << nn << " point " << nnn << " has time " << (int)dx <<  " and slope " << dy << endl;
		}
		for(int nnn=0; nnn<g_onevalues[n][nn]->GetN();++nnn){
			double dx,dy;
			g_onevalues[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>incepmax) incepmax = dy;
			if(dy<incepmin) incepmin = dy;
			if(dx>timemax ) timemax  = dx;
			if(dx<timemin ) timemin  = dx;
			//cout << "incep " << n << " " << nn << " point " << nnn << " has time " << (int)dx << " and onevalue " << dy << endl;
		}
	}}
	/*
	timemin = timemin - oneday*7; //minus 1 week
	timemax = timemax + oneday*75;//plus  75 days --> make Legend not overlap with the rest
	slopemax = slopemax * 1.2;
	incepmax = incepmax * 1.2;
	if(slopemin<0) slopemin = slopemin * 1.2;
	else slopemin = 0;
	if(incepmin<0) incepmin = incepmin * 1.2;
	else incepmin = 0;
	cout << "timemin,timemax " << (int)timemin << " " << (int)timemax << endl;
	cout << "slopemin,slopemax " << slopemin << " " << slopemax << endl;
	cout << "incepmin,incepmax " << incepmin << " " << incepmax << endl;
	*/
	timemin  = 1297983600;
	timemax  = 1369044400;
	slopemin = -2.;
	slopemax = 0.5;
	incepmin = -0.1;
	incepmax = 1.1;

	//legend

	TH1F *haxis = new TH1F("haxis", "haxis", 25, timemin, timemax);
	haxis->GetXaxis()->SetTitle("time"); haxis->GetXaxis()->SetTitleOffset(1.25); haxis->GetXaxis()->SetLabelSize(0.027);
	haxis->GetXaxis()->SetLabelOffset(0.02);haxis->GetYaxis()->SetTitleOffset(1.25); haxis->GetYaxis()->SetLabelSize(0.027);
	haxis->GetXaxis()->SetTimeDisplay(1); haxis->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
	haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); haxis->GetXaxis()->SetNdivisions(510, true);

	haxis->SetNameTitle("slopes_vs_time", "slopes vs. time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("burn-in slope"); 
	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	col->cd();
	col->Clear();

	col->cd();
	Legend1->SetHeader("russian crystals");
	col->SetName( "slopes_vs_time_daily_russian_burnin" );
	col->SetTitle("slopes_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	for(int n = 0; n<8; ++n) {if(g_slopes[n][0]->GetN()>0) g_slopes[n][0]->Draw("P");}
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	Util::PrintNoEPS(col, outname, outputdir);
	Util::PrintEPS(col, outname, outputdir);

	col->Clear();
	col->cd();
	Legend1->SetHeader("chinese crystals");
	col->SetName( "slopes_vs_time_daily_chinese_burnin" );
	col->SetTitle("slopes_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	for(int n = 0; n<8; ++n) {if(g_slopes[n][1]->GetN()>0) g_slopes[n][1]->Draw("P");}
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	Util::PrintNoEPS(col, outname, outputdir);
	Util::PrintEPS(col, outname, outputdir);

	haxis->SetNameTitle("onevalues_vs_time", "onevalues vs. time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("burn-in=1 equivalent"); 

	col->Clear();
	col->cd();
	Legend1->SetHeader("russian crystals");
	col->SetName( "onevalues_vs_time_daily_russian_burnin" );
	col->SetTitle("onevalues_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	for(int n = 0; n<8; ++n) {if(g_onevalues[n][0]->GetN()>0) g_onevalues[n][0]->Draw("P");}
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	Util::PrintNoEPS(col, outname, outputdir);
	Util::PrintEPS(col, outname, outputdir);

	col->Clear();
	col->cd();
	Legend1->SetHeader("chinese crystals");
	col->SetName( "onevalues_vs_time_daily_chinese_burnin" );
	col->SetTitle("onevalues_vs_time_daily_chinese_burnin" );
	haxis->Draw("hist");
	for(int n = 0; n<8; ++n) {if(g_onevalues[n][1]->GetN()>0) g_onevalues[n][1]->Draw("P");}
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	Util::PrintNoEPS(col, outname, outputdir);
	Util::PrintEPS(col, outname, outputdir);	

	}//plotoverview

	if(savefinalplots){
	TFile *newfile = new TFile ("FileOutputs/20121120_burnin_vs_VPTPNdiff_VectorNTuples_fit_slopes_onevalues_daily.root", "RECREATE");
	newfile->cd();
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		g_slopes[n][nn]->Write();
		g_onevalues[n][nn]->Write();
	}}
	newfile->Close();
	cout << "file saved as " << newfile->GetName() << endl;
	}
	
}