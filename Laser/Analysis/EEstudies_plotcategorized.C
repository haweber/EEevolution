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

void EEstudies_plotcategorized(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool plotoverview   = true;
	bool savefinalplots = true;


	TString outputdir               = "Plots/20121207/EEstudies/allcrystalsVPTPN/";

	Util::MakeOutputDir(outputdir);


	TFile *stdmapsfile     = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *mu_ECAL_EEp_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_russ");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEp_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_chin");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEm_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_russ");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEm_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_chin");//in Get: SIC <-> ECAL
	TH2D *EEp_producer     = (TH2D*)stdmapsfile->Get("EEp_producer");//1 russian, 2 chinese
	TH2D *EEm_producer     = (TH2D*)stdmapsfile->Get("EEm_producer");

	TFile *burninfile      = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEburnin.root");
	TH2D *burnin_EEm       = (TH2D*)burninfile->Get("burnin_EEm");
	TH2D *burnin_EEp       = (TH2D*)burninfile->Get("burnin_EEp");

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	//TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121023_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121116_EE_VPT_over_PN_all_VectorNTuples_2011and2012_pointperday.root");
	//TGraph* g_vptpn_las[101][101][2];
	TGraphErrors* g_vptpn_las[101][101][2];
	char gname[101];


	bool badsth[101][101][2];


	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
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
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				g_vptpn_las[i][j][k]->SetPoint(n,x,y/norm);
				double ex = g_vptpn_las[i][j][k]->GetErrorX(n);
				double ey = g_vptpn_las[i][j][k]->GetErrorY(n);
				g_vptpn_las[i][j][k]->SetPointError(n,ex/norm,ey/norm);

			}
			}
		}

		if(badsth[i][j][k]) continue;
		//cout << "x y z " << i << " " << j << " " << k << " has " << g_vptpn_las[i][j][k]->GetN() << " entries " << endl;
		if(g_vptpn_las[i][j][k]->GetN()>numberofdays) numberofdays = g_vptpn_las[i][j][k]->GetN();


		double mustd = 0; double burnin = 0; int prod = -1;
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
			if(prod==1) burnin = burnin_EEm->GetBinContent(bin);
			if(prod==2) burnin = burnin_EEm->GetBinContent(bin);
			if(prod==1) mustd = mu_ECAL_EEm_russ->GetBinContent(bin);
			if(prod==2) mustd = mu_ECAL_EEm_chin->GetBinContent(bin);
			if(mustd==0 || (mustd>=1.-10e-6&& mustd<=1.+10e-6) || mustd>=2.) continue;
			if(prod==1) russchin << "_eta_" <<eta << "_russian_" << "mustd_" << mustd << "_burnin_" << burnin;
			if(prod==2) russchin << "_eta_" <<eta << "_chinese_" << "mustd_" << mustd << "_burnin_" << burnin;
			rs = russchin.str();
			seta = rs;
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
			if(prod==1) burnin = burnin_EEp->GetBinContent(bin);
			if(prod==2) burnin = burnin_EEp->GetBinContent(bin);
			if(prod==1) mustd = mu_ECAL_EEp_russ->GetBinContent(bin);
			if(prod==2) mustd = mu_ECAL_EEp_chin->GetBinContent(bin);
			if(mustd==0 || (mustd>=1.-10e-6&& mustd<=1.+10e-6) || mustd>=2.) continue;
			if(prod==1) russchin << "_eta_" <<eta << "_russian_" << "mustd_" << mustd << "_burnin_" << burnin;
			if(prod==2) russchin << "_eta_" <<eta << "_chinese_" << "mustd_" << mustd << "_burnin_" << burnin;			rs = russchin.str();
			seta = rs;
		}
		bool highmustd  = false;
		bool  lowmustd  = false;
		bool highburnin = false;
		bool  lowburnin = false;

		if(burnin>=0.95 && burnin<=1.0)        highburnin = true;
		if(burnin<=0.9  && burnin>=0.8)         lowburnin = true;
		if(prod==2 && mustd>=1.6  && mustd<2.)  highmustd = true;
		if(prod==2 && mustd<=0.45 && mustd>0.)   lowmustd = true;
		if(prod==1 && mustd>=1.4  && mustd<2.)  highmustd = true;
		if(prod==1 && mustd<=0.45 && mustd>0.)   lowmustd = true;


		if(prod==2&&fabs(eta)>=1.8&&fabs(eta)<2.0) {
		stringstream charstring;
		charstring << gname;
		string cs = charstring.str();
		cs = cs + (string)"  " + seta;
		TString cs2 = (TString)cs;
		g_vptpn_las[i][j][k]->SetNameTitle(cs2, cs2);//oled instead of las
		col->SetName(cs2);
		col->SetTitle(cs2);
		col->cd();
		g_vptpn_las[i][j][k]->Draw("AP");
		TH1F *haxis = g_vptpn_las[i][j][k]->GetHistogram();//oled instead of las
		haxis->SetNameTitle(cs2, cs2);
	        haxis->SetMinimum(0.25);
		haxis->SetMaximum(1.25);
				haxis->GetXaxis()->SetTitle("time"); haxis->GetXaxis()->SetTitleOffset(1.25); haxis->GetXaxis()->SetLabelSize(0.027);
				haxis->GetYaxis()->SetTitle("VPT/PN"); haxis->GetYaxis()->SetTitleOffset(1.25); haxis->GetYaxis()->SetLabelSize(0.027);
				haxis->GetXaxis()->SetTimeDisplay(1); haxis->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
				haxis->GetXaxis()->SetLabelOffset(0.02); haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); haxis->GetXaxis()->SetNdivisions(510, true);
		g_vptpn_las[i][j][k]->SetHistogram(haxis);//oled instead of las
		col->Clear();
		haxis->Draw("hist");
		g_vptpn_las[i][j][k]->Draw("P");//oled instead of las
		TString outputfile = "";
		TString temp = "";
	//	if(highburnin) outputfile = "BurninHigh";
	//	if( lowburnin) outputfile = "BurninLow";
	//	if(highmustd ) outputfile = outputfile + "_MustdHigh_";
	//	if( lowmustd ) outputfile = outputfile + "_mustdLow_";
 		temp = TString::Format("VPTPN_EE%i_IX%i_IY%i",k, i, j);//with orange
		outputfile = outputfile + temp + seta;
//		temp = TString::Format("_Burnin_%g",burnin);
//		outputfile = outputfile + temp;
//		temp = TString::Format("_Mustd_%g",mustd);
//		outputfile = outputfile + temp;
		col->Update();
		if(g_vptpn_las[i][j][k]->GetN()>10) Util::PrintNoEPS(col, outputfile, outputdir);//oled instead of las
		if(g_vptpn_las[i][j][k]->GetN()>10) Util::PrintEPS(col, outputfile, outputdir);//oled instead of las
		}
	}}}//ijk

	
}