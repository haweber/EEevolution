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

void EEstudies_continued_orange_v2(){

	bool plottingall_VPTPN = false;
	bool VPTdiff_vs_mustd  = false;//plot VPTdiff vs mustd --> you probably want this to be false as later you plot this including the fit
	bool plotonlyAbsEta    = true; //if true don't print EE+ EE- separately
	bool normalized        = false;
	bool startdiff         = false;//only true if normalized == false
	int period             = -1;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	//NOTE: Also change possibly newFileName at the very bottom!!!!!
	//new scaling in 20121001, old scaling in 20120829
	char gname[101];
	TString outputdir = "Plots/20121023/EEstudies/continued/orange/notnormalized";//normalized";
	if(startdiff) outputdir = "Plots/20121023/EEstudies/continued/orange/diffwrt2011begin";
	if(normalized) 
	        outputdir = "Plots/20121023/EEstudies/continued/orange/normalized";
	TString outputdirac = "Plots/20121023/EEstudies/continued/orange/allcrystalswithmustd/notnormalized";
	if(normalized) 
	        outputdirac = "Plots/20121023/EEstudies/continued/orange/allcrystalswithmustd/normalized";
	TString periodTS; TString periods;

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
	outputdir   = outputdir   + periodTS;
	//outputdirac = outputdirac + periodTS;
	Util::MakeOutputDir(outputdir);
	if(plottingall_VPTPN)	Util::MakeOutputDir(outputdirac);

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	TFile *stdmapsfile     = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *mu_ECAL_EEp_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_russ");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEp_chin = (TH2D*)stdmapsfile->Get("mu_SIC_EEp_chin");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEm_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_russ");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEm_chin = (TH2D*)stdmapsfile->Get("mu_SIC_EEm_chin");//in Get: SIC <-> ECAL
	//as all crystals chinese crystals should have SIC information, all crystals with mu_SIC==0 are russians
	TGraph* h_mustd_vs_VPTPNdiff[16][2];
	TGraph* a_mustd_vs_VPTPNdiff[8][2];
	TH1D *h_VPTPNdiff_vs_eta_EEm_russ = new TH1D("h_VPTPNdiff_vs_eta_EEm_russ", "h_VPTPNdiff_vs_eta_EEm_russ", 8, -3.0, -1.4);
	TH1D *h_VPTPNdiff_vs_eta_EEp_russ = new TH1D("h_VPTPNdiff_vs_eta_EEp_russ", "h_VPTPNdiff_vs_eta_EEp_russ", 8,  1.4,  3.0);
	TH1D *h_VPTPNdiff_vs_eta_russ     = new TH1D("h_VPTPNdiff_vs_eta_russ"    , "h_VPTPNdiff_vs_eta_russ"    , 8,  1.4,  3.0);
	TH1D *h_VPTPNdiff_vs_eta_EEm_chin = new TH1D("h_VPTPNdiff_vs_eta_EEm_chin", "h_VPTPNdiff_vs_eta_EEm_chin", 8, -3.0, -1.4);
	TH1D *h_VPTPNdiff_vs_eta_EEp_chin = new TH1D("h_VPTPNdiff_vs_eta_EEp_chin", "h_VPTPNdiff_vs_eta_EEp_chin", 8,  1.4,  3.0);
	TH1D *h_VPTPNdiff_vs_eta_chin     = new TH1D("h_VPTPNdiff_vs_eta_chin"    , "h_VPTPNdiff_vs_eta_chin"    , 8,  1.4,  3.0);
	string string_eta[16] = {"m3p0", "m2p8", "m2p6", "m2p4", "m2p2", "m2p0", "m1p8", "m1p6", "1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	for(int n = 0; n<16; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer

	string hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn];
	h_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
	h_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
//	h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
	//if(normalized) h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.3);
	//else           h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.22);
	h_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //h_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
	h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPTPNmax-VPTPNmin"); h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("mu_std (m^-1)");

	if(n<8){
		hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_feta[n] + string_prod[nn];
		a_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
		a_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
	//	a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
	}//if(n<8)

	}//for(int nn = 0; nn<2; ++ nn) // producer
	}//for(int n = 0; n<16; ++n)    // eta bin


	h_VPTPNdiff_vs_eta_EEm_russ->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_EEm_russ->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_EEm_russ->SetLineWidth(4);       h_VPTPNdiff_vs_eta_EEm_russ->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_EEm_russ->SetMarkerColor(kRed);  h_VPTPNdiff_vs_eta_EEm_russ->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_EEm_russ->SetLineColor(kRed);    h_VPTPNdiff_vs_eta_EEm_russ->GetXaxis()->SetTitle("#eta");
	h_VPTPNdiff_vs_eta_EEm_russ->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
	h_VPTPNdiff_vs_eta_EEp_russ->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_EEp_russ->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_EEp_russ->SetLineWidth(4);       h_VPTPNdiff_vs_eta_EEp_russ->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_EEp_russ->SetMarkerColor(kRed);  h_VPTPNdiff_vs_eta_EEp_russ->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_EEp_russ->SetLineColor(kRed);    h_VPTPNdiff_vs_eta_EEp_russ->GetXaxis()->SetTitle("#eta");
	h_VPTPNdiff_vs_eta_EEp_russ->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
	h_VPTPNdiff_vs_eta_russ    ->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_russ    ->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_russ    ->SetLineWidth(4);       h_VPTPNdiff_vs_eta_russ    ->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_russ    ->SetMarkerColor(kRed);  h_VPTPNdiff_vs_eta_russ    ->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_russ    ->SetLineColor(kRed);    h_VPTPNdiff_vs_eta_russ    ->GetXaxis()->SetTitle("|#eta|");
	h_VPTPNdiff_vs_eta_russ    ->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
	h_VPTPNdiff_vs_eta_EEm_chin->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_EEm_chin->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_EEm_chin->SetLineWidth(4);       h_VPTPNdiff_vs_eta_EEm_chin->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_EEm_chin->SetMarkerColor(kBlue); h_VPTPNdiff_vs_eta_EEm_chin->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_EEm_chin->SetLineColor(kBlue);   h_VPTPNdiff_vs_eta_EEm_chin->GetXaxis()->SetTitle("#eta");
	h_VPTPNdiff_vs_eta_EEm_chin->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
	h_VPTPNdiff_vs_eta_EEp_chin->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_EEp_chin->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_EEp_chin->SetLineWidth(4);       h_VPTPNdiff_vs_eta_EEp_chin->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_EEp_chin->SetMarkerColor(kBlue); h_VPTPNdiff_vs_eta_EEp_chin->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_EEp_chin->SetLineColor(kBlue);   h_VPTPNdiff_vs_eta_EEp_chin->GetXaxis()->SetTitle("#eta");
	h_VPTPNdiff_vs_eta_EEp_chin->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
	h_VPTPNdiff_vs_eta_chin    ->SetMaximum(0.4);       h_VPTPNdiff_vs_eta_chin    ->SetMinimum(-0.05);
	h_VPTPNdiff_vs_eta_chin    ->SetLineWidth(4);       h_VPTPNdiff_vs_eta_chin    ->SetMarkerStyle(20); 
	h_VPTPNdiff_vs_eta_chin    ->SetMarkerColor(kBlue); h_VPTPNdiff_vs_eta_chin    ->SetMarkerSize(2);
	h_VPTPNdiff_vs_eta_chin    ->SetLineColor(kBlue);   h_VPTPNdiff_vs_eta_chin    ->GetXaxis()->SetTitle("|#eta|");
	h_VPTPNdiff_vs_eta_chin    ->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");

	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121102_EE_VPT_over_PN_WG1_oled.root");
	TGraph* g_vptpn_las[101][101][2];
	double mustd_russ = 0; double mustd_chin = 0;

	//define and reset all values
	double vptpnnorm[101][101][2];
	double vptpnmax[101][101][2];
	double vptpnmin[101][101][2];
	double timemax[101][101][2];
	double timemin[101][101][2];
	int    vptpnmaxindex[101][101][2];
	int    vptpnminindex[101][101][2];
	double vptpnmmustd[101][101][2];
	bool   vptpnruss[101][101][2];
	int    lastsameentry[101][101][2];
	double lastentry[101][101][2];
	bool   hasbadperiod[101][101][2];
	int    lastsameentryI[101][101][2];//check only period of study
	double lastentryI[101][101][2];//check only period of study
	bool   hasbadperiodI[101][101][2];//check only period of study
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		vptpnnorm[i][j][k] = 0; 
		vptpnmax[i][j][k] = 0;
		vptpnmin[i][j][k] = 99; 
		timemax[i][j][k] = 0;
		timemin[i][j][k] = 0; 
		vptpnmaxindex[i][j][k] = -1;
		vptpnminindex[i][j][k] = -1; 
		vptpnmmustd[i][j][k] = 0; 
		vptpnruss[i][j][k] = false;
		lastentry[i][j][k] = -1;
		lastsameentry[i][j][k] = 0;
		hasbadperiod[i][j][k] = false;
		lastentryI[i][j][k] = -1;
		lastsameentryI[i][j][k] = 0;
		hasbadperiodI[i][j][k] = false;
	}}}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "start computation" << endl;
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		if (k==0) sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEp", i, j);

		g_vptpn_las[i][j][k] = (TGraph*)vpt_values->Get(gname);
		if(g_vptpn_las[i][j][k]->GetN()==0) continue;

		if(g_vptpn_las[i][j][k]->GetN()>0){ // start normalization
			double norm = 0; int countnorm = 0;
			//first clean up
			vector<int> removeindex; removeindex.clear();//has to go from end to start
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y, xprev,yprev, xpost,ypost, xprev2,yprev2;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(hasbadperiod[i][j][k]==false){
					if(x<1320120000) continue;//veto 2012 (>) or 2011(<)
					if(y==lastentry[i][j][k] && y > 0) lastsameentry[i][j][k] = lastsameentry[i][j][k] + 1;
					else if(y!=lastentry[i][j][k] && y>0) lastsameentry[i][j][k] = 0;
					lastentry[i][j][k] = y;
					if(lastsameentry[i][j][k]>=1000) hasbadperiod[i][j][k] = true;
					if(period==6 && x>1333530000 && x<1334240000 && lastsameentry[i][j][k]>=100) hasbadperiod[i][j][k] = true;
				}
				//1st period
				if(period==1 && !plottingall_VPTPN && x<1300200000) continue;//15/03/2011 	
				if(period==1 && !plottingall_VPTPN && x>1301210000) continue;//27/03/2011	
				//2nd period
				if(period==2 && !plottingall_VPTPN && x<1302590000) continue;//12/04/2011	
				if(period==2 && !plottingall_VPTPN && x>1304620000) continue;//05/05/2011	
				//3rd period
				if(period==3 && !plottingall_VPTPN && x<1305280000) continue;//13/05/2011	
				if(period==3 && !plottingall_VPTPN && x>1309540000) continue;//01/07/2011	
				//4th period
				if(period==4 && !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
				if(period==4 && !plottingall_VPTPN && x>1314270000) continue;//25/08/2011	
				//5th period
				if(period==5 && !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
				if(period==5 && !plottingall_VPTPN && x>1320250000) continue;//02/11/2011	
				//6th period
				if(period==6 && !plottingall_VPTPN && x<1333740000) continue;//06/04/2012
				if(period==6 && !plottingall_VPTPN && x>1334910000) continue;//20/04/2012
				//7th period
				if(period==7 && !plottingall_VPTPN && x<1336230000) continue;//05/05/2012
				if(period==7 && !plottingall_VPTPN && x>1340000000) continue;//18/06/2012
				//original period
				if(period==0 && !plottingall_VPTPN && x<1314920000) continue;//02/09/2011
				if(period==0 && !plottingall_VPTPN && x>1317690000) continue;//04/10/2011
				//4th period - first 1.1 fb-1
				if(period==41&& !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
				if(period==41&& !plottingall_VPTPN && x>1312800923) continue;//08/08/2011	
				//5th period - first 1.1 fb-1
				if(period==51&& !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
				if(period==51&& !plottingall_VPTPN && x>1316813137) continue;//23/09/2011	
				//3rd period - first 950 pb-1
				if(period==32&& !plottingall_VPTPN && x<1305280000) continue;//13/05/2011	
				if(period==32&& !plottingall_VPTPN && x>1308900000) continue;//24/06/2011
				//4th period - first 950 pb-1
				if(period==42&& !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
				if(period==42&& !plottingall_VPTPN && x>1312720000) continue;//07/08/2011
				//5th period - first 950 pb-1
				if(period==52&& !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
				if(period==52&& !plottingall_VPTPN && x>1316695000) continue;//22/09/2011
				//7th period - first 950 pb-1
				if(period==72&& !plottingall_VPTPN && x<1336230000) continue;//05/05/2012
				if(period==72&& !plottingall_VPTPN && x>1337260000) continue;//17/05/2012
				//-1th period - technical stop Sept.2011
				if(period==-1&& !plottingall_VPTPN && x<1313972000) continue;//22/08/2012
				if(period==-1&& !plottingall_VPTPN && x>1315405000) continue;//07/09/2012
				//cout << "t=" << x << " VPT/PN=" << y << endl;
				if(hasbadperiodI[i][j][k]==false){
					if(y==lastentryI[i][j][k] && y > 0) lastsameentryI[i][j][k] = lastsameentryI[i][j][k] + 1;
					else if(y!=lastentryI[i][j][k] && y>0) lastsameentryI[i][j][k] = 0;
					lastentryI[i][j][k] = y;
					if(lastsameentryI[i][j][k]==200) hasbadperiodI[i][j][k] = true;
				}
				if(n==g_vptpn_las[i][j][k]->GetN()-1) g_vptpn_las[i][j][k]->GetPoint(n,xpost,ypost);
				else g_vptpn_las[i][j][k]->GetPoint(n+1,xpost,ypost);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev,yprev);
				else g_vptpn_las[i][j][k]->GetPoint(n-1,xprev,yprev);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev2,yprev2);
				else if(n==1) g_vptpn_las[i][j][k]->GetPoint(n-1,xprev2,yprev2);
				else g_vptpn_las[i][j][k]->GetPoint(n-2,xprev2,yprev2);

				if(y==0) continue;
				if((yprev>y*1.1 && ypost>y*1.1) || (yprev<y*0.9 && ypost<y*0.9)) g_vptpn_las[i][j][k]->SetPoint(n,x,yprev);
				if(ypost<y*1.1 && yprev2<y*1.1 && ypost<yprev*1.1 && yprev2<yprev*1.1){
					g_vptpn_las[i][j][k]->SetPoint(n,x,yprev2);//clean only n
					if(n>0) g_vptpn_las[i][j][k]->SetPoint(n-1,xprev,yprev2);
				}
				//add extra cleaning
				double meany = g_vptpn_las[i][j][k]->GetMean(2);
				if(y>1.67*meany)  removeindex.push_back(n);
				else if(y<0.67*meany) removeindex.push_back(n);
			}//cleaning
			if(removeindex.size()>0){
			for(unsigned int ri = 0; ri<removeindex.size(); ++ri){
				g_vptpn_las[i][j][k]->RemovePoint(removeindex[ri]);
			}
			}

		//start normalized
			if(normalized){//i.e. normalize to end of technical stop 2011 in september
				for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
					double x,y;
					g_vptpn_las[i][j][k]->GetPoint(n,x,y);
					//1st period
					if(period==1 && !plottingall_VPTPN && x<1300200000) continue;//15/03/2011 	
					if(period==1 && !plottingall_VPTPN && x>1301210000) continue;//27/03/2011	
					//2nd period
					if(period==2 && !plottingall_VPTPN && x<1302590000) continue;//12/04/2011	
					if(period==2 && !plottingall_VPTPN && x>1304620000) continue;//05/05/2011	
					//3rd period
					if(period==3 && !plottingall_VPTPN && x<1305280000) continue;//13/05/2011	
					if(period==3 && !plottingall_VPTPN && x>1309540000) continue;//01/07/2011	
					//4th period
					if(period==4 && !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
					if(period==4 && !plottingall_VPTPN && x>1314270000) continue;//25/08/2011	
					//5th period
					if(period==5 && !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
					if(period==5 && !plottingall_VPTPN && x>1320250000) continue;//02/11/2011	
					//6th period
					if(period==6 && !plottingall_VPTPN && x<1333740000) continue;//06/04/2012
					if(period==6 && !plottingall_VPTPN && x>1334910000) continue;//20/04/2012
					//7th period
					if(period==7 && !plottingall_VPTPN && x<1336230000) continue;//05/05/2012
					if(period==7 && !plottingall_VPTPN && x>1340000000) continue;//18/06/2012
					//original period
					if(period==0 && !plottingall_VPTPN && x<1314920000) continue;//02/09/2011
					if(period==0 && !plottingall_VPTPN && x>1317690000) continue;//04/10/2011
					//4th period - first 1.1 fb-1
					if(period==41&& !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
					if(period==41&& !plottingall_VPTPN && x>1312800923) continue;//08/08/2011	
					//5th period - first 1.1 fb-1
					if(period==51&& !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
					if(period==51&& !plottingall_VPTPN && x>1316813137) continue;//23/09/2011	
					//3rd period - first 950 pb-1
					if(period==32&& !plottingall_VPTPN && x<1305280000) continue;//13/05/2011	
					if(period==32&& !plottingall_VPTPN && x>1308900000) continue;//24/06/2011
					//4th period - first 950 pb-1
					if(period==42&& !plottingall_VPTPN && x<1310710000) continue;//15/07/2011	
					if(period==42&& !plottingall_VPTPN && x>1312720000) continue;//07/08/2011
					//5th period - first 950 pb-1
					if(period==52&& !plottingall_VPTPN && x<1315120000) continue;//04/09/2011	
					if(period==52&& !plottingall_VPTPN && x>1316695000) continue;//22/09/2011
					//7th period - first 950 pb-1
					if(period==72&& !plottingall_VPTPN && x<1336230000) continue;//05/05/2012
					if(period==72&& !plottingall_VPTPN && x>1337260000) continue;//17/05/2012
					//-1th period - technical stop Sept.2011
					if(period==-1&& !plottingall_VPTPN && x<1313972000) continue;//22/08/2012
					if(period==-1&& !plottingall_VPTPN && x>1315405000) continue;//07/09/2012
					//cout << "t=" << x << " VPT/PN=" << y << endl;
					if(y==0) continue;
					norm += y;
					++countnorm;
					if(countnorm>5) break;
				}
				//if(countnorm==0) cout << "wfk i.j.k = " << i << ":"<< j << ":" << k << endl;
				if(countnorm>0) vptpnnorm[i][j][k] = norm/double(countnorm);
				//do normalization
				if(countnorm!=0){//normalize only if flag is set
				for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
					double x,y;
					g_vptpn_las[i][j][k]->GetPoint(n,x,y);
					g_vptpn_las[i][j][k]->SetPoint(n,x,y/(vptpnnorm[i][j][k]));//normalization
					if(y/(vptpnnorm[i][j][k])>1.5){//a crazy jump??
						g_vptpn_las[i][j][k]->GetPoint(n,x,y);
						double xx,yy;
						g_vptpn_las[i][j][k]->GetPoint(n+1,xx,yy);//already normalized
						g_vptpn_las[i][j][k]->SetPoint(n,x,yy);
					}
				}
				}
				else {
					cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" could not normalized " << countnorm << " number of entries " << g_vptpn_las[i][j][k]->GetN() << endl;}
			}
			else { //new 2012/08/29: notnormalize means normalizing at beginning of 2011 data taking
				for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
					double x,y;
					g_vptpn_las[i][j][k]->GetPoint(n,x,y);
					//cout << "t=" << x << " VPT/PN=" << y << endl;
					if(y==0) continue;
					if(x<1296518400) continue;//01/02/2012 --> i.e. take first ten measurements from year 2011
					norm += y;
					++countnorm;
					if(countnorm>10) break;
				}
				//if(countnorm==0) cout << "wfk i.j.k = " << i << ":"<< j << ":" << k << endl;
				if(countnorm>0) vptpnnorm[i][j][k] = norm/double(countnorm);
				//do normalization
				if(countnorm!=0){//normalize only if flag is set
				for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
					double x,y;
					g_vptpn_las[i][j][k]->GetPoint(n,x,y);
					g_vptpn_las[i][j][k]->SetPoint(n,x,y/(vptpnnorm[i][j][k]));//normalization
					if(y/(vptpnnorm[i][j][k])>1.5){// a crazy jump
						g_vptpn_las[i][j][k]->GetPoint(n,x,y);
						double xx,yy;
						g_vptpn_las[i][j][k]->GetPoint(n+1,xx,yy);//already normalized
						g_vptpn_las[i][j][k]->SetPoint(n,x,yy);
					}
				}
				}
				else {cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" could not set startvalue " << countnorm << " number of entries " << g_vptpn_las[i][j][k]->GetN() << endl;}
			}
		}//if(g_vptpn_las[i][j][k]->GetN()>0)  // end of normalizing
		int bin; float eta; string seta; stringstream russchin; string rs; stringstream etass; string es;
		if(k==0){//choose xx_EEm_yy histos
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEm_russ->FindBin(i,j);
			mustd_russ = mu_ECAL_EEm_russ->GetBinContent(bin);
			mustd_chin = mu_ECAL_EEm_chin->GetBinContent(bin);
		} else{ 
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEp_russ->FindBin(i,j);
			mustd_russ = mu_ECAL_EEp_russ->GetBinContent(bin);
			mustd_chin = mu_ECAL_EEp_chin->GetBinContent(bin);
		}
	//		if(mustd_russ==0 && mustd_chin==0) continue;
			/*else*/ if(mustd_chin!=0){
			//	if(mustd_chin>=1.-10e-6&& mustd_chin<=1.+10e-6) continue;
			//	if(mustd_chin>=2.) continue;
				russchin << mustd_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_chin;
				vptpnruss[i][j][k] = false;//is not russian thus it is chinese
			}
			else /*if(mustd_russ!=0)*/{
			//	if(mustd_russ>=1.-10e-6&& mustd_russ<=1.+10e-6) continue;
			//	if(mustd_russ>=2.) continue;
				russchin << mustd_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_russ;
				vptpnruss[i][j][k] = true;//is russian
			}
		//clean out noisy crystals and select time window of study

//		everything below commented at col.1 (i.e. without 'tab') is due to new ntuples
		//rising crystals are crystals whose VPT/PN rises even during radiation
		bool risingcrystals = (k==1&&((i==95&&j==67)||(i==95&&j==66)||(i==94&&j==27)||(i==93&&j==34)||(i==92&&j==27)||(i==79&&j==92)||(i==79&&j==91)||(i==77&&j==91)||(i==69&&j==95)||(i==54&&j==99)||(i==54&&j==98)||(i==53&&j==100)||(i==52&&j==100)||(i==38&&j==97)||(i==36&&j==97)||(i==29&&j==95)||(i==27&&j==6)||(i==25&&j==92)||(i==24&&j==91)||(i==19&&j==83)||(i==17&&j==84)||(i==17&&j==83)||(i==16&&j==82)||(i==10&&j==30)||(i==8&&j==27)||(i==5&&j==62)||(i==4&&j==64)||(i==4&&j==63)||(i==4&&j==50)))||(i==87&&j==25) || (k==0&&((i==99&&j==57)||(i==97&&j==59)||(i==97&&j==38)||(i==96&&j==42)||(i==96&&j==40)||(i==96&&j==37)||(i==94&&j==73)||(i==92&&j==27)||(i==84&&j==15)||(i==80&&j==15)||(i==80&&j==13)||(i==80&&j==12)||(i==80&&j==9)||(i==78&&j==91)||(i==78&&j==10)||(i==77&&j==9)||(i==76&&j==11)||(i==76&&j==9)||(i==63&&j==96)||(i==61&&j==5)||(i==51&&j==3)||(i==29&&j==8)||(i==24&&j==9)||(i==15&&j==25)||(i==15&&j==19)||(i==15&&j==16)||(i==14&&j==85)||(i==14&&j==82)||(i==14&&j==19)||(i==10&&j==77)||(i==10&&j==27)||(i==10&&j==24)||(i==8&&j==28)||(i==8&&j==27)||(i==7&&j==28)));
		//these crystals are crystals whose VPT/PN rises even during some time of radiation, but not always
		bool partiallyrising= (k==1&&((i==94&&j==70)||(j==91&&j==71)||(i==91&&j==66)||(i==91&&j==27)||(i==80&&j==91)||(i==76&&j==92)||(i==76&&j==91)||(i==70&&j==91)||(i==65&&j==97)||(i==55&&j==97)||(i==52&&j==97)||(i==36&&j==96)||(i==30&&j==8)||(i==28&&j==10)||(i==28&&j==9)||(i==25&&j==91)||(i==22&&j==91)||(i==19&&j==39))) || (k==0&&((i==96&&j==60)||(i==92&&j==70)||(i==85&&j==14)||(i==84&&j==14)||(i==76&&j==12)||(i==72&&j==93)||(i==70&&j==95)||(i==64&&j==5)||(i==61&&j==4)||(i==55&&j==15)||(i==55&&j==3)||(i==31&&j==10)||(i==26&&j==7)||(i==25&&j==18)||(i==15&&j==20)||(i==12&&j==21)||(i==8&&j==34)));
		//these crystals are crystals whose VPT/PN rises even during beginning 2011/2012, but then loss dominates
		bool partiallyrising2 = (k==1&&((i==41&&j==65)||(i==21&&j==91))) || (k==0&&((i==93&&j==69)||(i==32&&j==73)||(i==14&&j==83)));
		bool bad2012 = (k==0&&((i==14&&j==50)||(i==18&&j==21)||(i==24&&j==51)||(i==25&&j==64)||(i==26&&j==50)||(i==31&&j==38)||(i==31&&j==72)||(i==32&&j==32)||(i==32&&j==56)||(i==32&&j==61)||(i==32&&j==62)||(i==32&&j==73)||(i==33&&j==57)||(i==33&&j==75)||(i==38&&j==15)||(i==38&&j==40)||(i==39&&j==82)||(i==44&&j==13)||(i==49&&j==82)||(i==52&&j==2)||(i==70&&j==50)||(i==91&&j==76))) || (k==1&&((i==54&&j==36)||(i==56&&j==2)||(i==58&&j==38)||(i==59&&j==38)||(i==60&&j==39)||(i==62&&j==51)||(i==65&&j==27)||(i==70&&j==32)||(i==70&&j==33)||(i==81&&j==35)));

		bool Ospecialcrystals = ((k==0&&( (i==24&&j==51)||(i==26&&j==50)||(i==32&&j==56)||(i==33&&j==57)||(i==34&&j==48)||(i==36&&j==40)||(i==36&&j==48)||(i==38&&j==50)||(i==38&&j==63)||(i==39&&j==49)||(i==39&&j==50)||(i==42&&j==39)||(i==43&&j==41)||(i==47&&j==36)||(i==49&&j==36)||(i==49&&j==39)||(i==52&&j==11)||(i==53&&j==63)||(i>=58&&i<=59&&j>=60&&j<=61)||(i==61&&j==57)||(i<=70&&i>=62&&j>=46&&j<=58) ) ) || (k==1&&((i==24&&j==50)||(i>=36&&i<=37&&j>=45&&j<=60)||(i==38&&j==64)||(i==39&&j>=55&&j<=59)||(i==39&&j==64)||(i==40&&j==63)||(i==45&&j==39)||(i==45&&j==39)||(i==55&&j==63)||(i==56&&j>=38&&j<=39)||(i>=56&&i<=57&&j>=28&&j<=62)||(i>=59&&i<=60&&j>=39&&j<=40)||(i==61&&j==65)||(i==62&&j==51)||(i==62&&j>=54&&j<=55)||(i==65&&j==51) ) ) );
		bool Onoisecrystals2 = (k==0&&((i==26&&j==50)||(i==36&&j==48)))||(k==1&&((i==24&&j==50)))||((k==0&&((i==24&&j==51))))||((k==0&&((i==32&&j==56)||(i==38&&j==63)||(i==52&&j==11)||(i==53&&j>=62&&j<=68)||(i==55&&j==63)||(i>=58&&i<=59&&j==60)||(i==58&&j==61)||(i==60&&j==63)||(i==61&&j==57)))||(k==1&&((i==39&&j>=55&&j<=59)||(i==39&&j==64)||(i==39&&j==64)||(i>=59&&i<=60&&j>=39&&j<=40))))||(k==1&&((i==45&&j==39)||(i==56&&j==38)))||((k==0&&((i==26&&j==50))))||(k==1&&((i==45&&j==39)||(i==56&&j==38)))||((k==0&&((i==24&&j==51)||(i==34&&j==48)||(i==38&&j==50)||(i==47&&j==36)||(i==58&&j==61)||(i==60&&j==63)||(i==61&&j==57)||(i>=62&&i<=70&&j>=49&&j<=58)))||(k==1&&((i==24&&j==50)||(i==38&&j==64)||(i==39&&j>=55&&j<=59)||(i==62&&j==51)||(i==56&&j==39)||(i==39&&j==64))))||((k==0&&((i==36&&j==48)||(i==49&&(j==36||j==39)))));
		bool Onoisecrystals = ((k==0&&((i==40&&j==40)||(i==50&&j==39)||(i>=54&&i<=55&&j>=31&&j<=35)||(i==51&&j==63)||(i==54&&(j==65&&j==86))||(i>=61&&i<=74&&j>=91&&j<=97)||(i>=76&&i<=85&&j>=56&&j<=71)||(i>=92&&i<=97&&j>=27&&j<=42)))||(k==1&&((i>=1&&i<=4&&j>=41&&j<=50)||(i>=6&&i<=10&&j>=27&&j<=30)||(i==4&&j==63)||(i>=7&&i<=13&&j>=74&&j<=80)||(i==11&&j==65)||(i==15&&j==57)||(i==16&&j==82)||(i==17&&j==84)||(i==18&&j==82)||(i==19&&j==84)||(i>=27&&i<=29&&j>=6&&j<=10)||(i>=27&&i<=28&&j>=46&&j<=49)||(i==29&&j==95)||(i==30&&(j==9||j==95))||(i>=31&&i<=35&&j>=41&&j<=55)||(i==31&&j==94)||(i==34&&j>=17&&j<=29)||(i==36&&j>=48)||(i>=36&&i<=42&&j>=92&&j<=97)||(i==39&&j==49)||(i==41&&j==18)||(i>=41&&i<=50&&j>=62&&j<=65)||(i>=48&&i<=49&&j>=18&&j<=35)||(i==52&&j==4)||(i>=78&&i<=85&&j>=12&&j<=18)||(i==78&&j==76))));
		int ijk = i + j*1000;
		if(k==0) ijk = -ijk;
bool Onoisecrystalsp5 = (ijk==94030)||(k==1&&((i>=26&&i<=28)&&(j>=91&&j<=95)))||(k==1&&((i>=29&&i<=30)&&(j>=91&&j<=92)))||(ijk==93029);
bool Onoisecrystals3 = ijk==-30010||ijk==-23012||ijk==-81015||ijk==-83015||ijk==-52016||ijk==-86016||ijk==-86017||ijk==-58018||ijk==-87018||ijk==-53020||ijk==-66020||ijk==-70020||ijk==-86020||ijk==-54025||ijk==-64026||ijk==-88031||ijk==-46032||ijk==-72033||ijk==-57034||ijk==-37036||ijk==-40036||ijk==-43036||ijk==-47036||(k==0&&i==36&&j>=36&&j<=53)||ijk==-62036||ijk==-80036||ijk==-41037||ijk==-54037||ijk==-61037||ijk==-62037||ijk==-39038||ijk==-47038||ijk==-48038||ijk==-51038||ijk==-53038||ijk==-55038||ijk==-43039||(k==0&&i==39&&j>=48&&j<=54)||ijk==-65039||ijk==-76039||ijk==-28040||ijk==-44040||ijk==-42041||ijk==-58041||ijk==-60041||ijk==-59042||ijk==-60042||ijk==-39044||ijk==-41043||ijk==-37045||ijk==-36046||ijk==-65046||ijk==-39047||ijk==-63047||ijk==-39048||ijk==-62048||ijk==-65048||ijk==-37049||ijk==-38048||ijk==-63049||ijk==-1051||(k==0&&i==51&&j>=31&&j<=39)||ijk==-62051||ijk==-64051||ijk==-66051||(k==0&&i==52&&j>=1&&j<=62)||(k==0&&i==53&&j>=3&&j<=39)||(k==0&&i==54&&j>=1&&j<=63)||ijk==-65054||ijk==-86054||(k==0&&i==55&&j>=36&&j<=39)||ijk==-64055||ijk==-66055||ijk==-67055||ijk==-69055||ijk==-6056||ijk==-62056||(k==0&&i==56&&j>=36&&j<=40)||ijk==-64056||ijk==-65056||(k==0&&i==57&&j>=36&&j<=40)||ijk==-63057||ijk==-64057||(k==0&&i==58&&j>=36&&j<=41)||ijk==-63058||(k==0&&i==59&&j>=36&&j<=42)||(k==0&&i==60&&j>=36&&j<=62)||ijk==-21061||(k==0&&i==61&&j>=36&&j<=56)||ijk==-6062||ijk==-21062||(k==0&&i==62&&j>=36&&j<=48)||ijk==-64062||ijk==-13063||ijk==-24063||(k==0&&i==63&&j>=36&&j<=59)||(k==0&&i==64&&j>=36&&j<=47)||(k==0&&i==65&&j>=21&&j<=48)||ijk==-61065||ijk==-63065||ijk==-48066||ijk==-75066||(k==0&&i>=67&&i<=70&&j>=46&&j<=48)||ijk==-8075||ijk==-55075||ijk==-53081||(k==0&&i==75&&j>=91&&j<=95)||ijk==-91077||ijk==-92077||ijk==-75085||ijk==-81085||ijk==-36086||ijk==-76089||(k==0&&i==91&&j>=26&&j<=30)||(k==0&&i>=92&&i<=95&&j==26)||(k==0&&i>=96&&i<=97&&j>=43&&j<=45)||(k==0&&i>=98&&i<=100&&j>=41&&j<=45) || (k==1&&i>=4&&i<=5&&j>=61&&j<=65)||(k==1&&i==6&&j>=71&&j<=74)||ijk==7207||ijk==72008||(k==1&&i==9&&j>=21&&j<=25)||ijk==72009||ijk==23010||ijk==72010||ijk==31012||(k==1&&i==13&&j>=21&&j<=23)||(k==1&&i==14&&j>=16&&j<=20)||(k==0&&i==14&&j>=83&&j<=85)||ijk==85015||(k==0&&i==16&&j>=81&&j<=83)||ijk==82017||ijk==25018||ijk==26018||ijk==39018||ijk==83018||(k==1&&i==19&&j>=92&&j<=85)||ijk==87019||ijk==37020||ijk==40020||ijk==75020||ijk==82020||ijk==84020||ijk==86020||ijk==35023||ijk==35024||(k==1&&i==29&&j>=82&&j<=84)||ijk==94029||ijk==93030||(k==1&&i==30&&j>=7&&j<=10)||ijk==59031||ijk==77034||ijk==78034||ijk==89035||(k==1&&i==36&&j>=44&&j<=47)||ijk==36037||ijk==37037||ijk==44037||(k==1&&j==37&&j>=46&&j<=62)||(k==1&&i>=38&&i<=39&&j>=46&&j<=60)||ijk==65039||ijk==77038||ijk==77039||ijk==39040||ijk==45040||(k==1&&i==40&&j>=56&&j<=60)||(k==0&&i==40&&j>=63&&j<=65)||ijk==36041||ijk==41041||ijk==42041||(k==1&&i==41&&j>=58&&j<=61)||ijk==38042||ijk==39042||ijk==40042||ijk==41042||ijk==59042||ijk==60042||ijk==61042||ijk==38043||ijk==60043||ijk==61043||ijk==76043||ijk==40044||ijk==61044||ijk==40045||ijk==61045||(k==1&&i>=46&&i<=54&&j>=31&&j<=39)||(k==1&&i==47&&j>=71&&j<=79)||(k==1&&i>=46&&i<=50&&j>=66&&j<=70)||ijk==80046||ijk==37052||ijk==67052||(k==1&&i==54&&j>=26&&j<=30)||ijk==31055||ijk==67055||ijk==69055||ijk==61056||ijk==65056||ijk==64057||ijk==60058||ijk==5059||ijk==63058||ijk==64058||ijk==59059||ijk==60059||ijk==65059||ijk==36060||ijk==58060||ijk==59060||ijk==60060||ijk==63060||ijk==44061||ijk==59061||ijk==65061||ijk==42062||(k==1&&(i>=62&&i<=70)&&(j>=46&&j<=55))||ijk==57062||ijk==59062||ijk==60062||ijk==57063||(k==1&&i==64&&j>=57&&j<=60)||ijk==63064||ijk==63065||ijk==79065||ijk==36072||ijk==21073||(k==1&&i==73&&j>=40&&j<=50)||ijk==6074||ijk==9074||ijk==84074||ijk==3576||ijk==27077||ijk==28077||ijk==29077||ijk==28078||ijk==29079||ijk==46080||ijk==83080||ijk==31081||ijk==82081||ijk==31082||ijk==82082||ijk==83082||ijk==85082||ijk==81083||ijk==86083||ijk==34084||ijk==58084||ijk==81084||ijk==82084||ijk==84084||ijk==85084||ijk==81085||ijk==84085||ijk==52087||ijk==85088||ijk==50089||ijk==44090||(k==1&&i==91&&j>=26&&j<=29)||(k==1&&i==92&&j>=28&&j<=29)||(k==1&&i==94&&j>=29&&j<=40)||(k==1&&i==95&&j>=26&&j<=30)||ijk==33091||ijk==32092||ijk==74092||ijk==28093||ijk==30093||ijk==26094||ijk==45094||ijk==49094||ijk==35095||ijk==56095||ijk==43096||ijk==50100;

		//e.g. k==0, j==23, i==94 --> ijk = -23094;
		//e.g. k==1, j==46, i==76 --> ijk = +46076;

//		bool noisecrystals = ((k==1&&((i==78&&j==76)||(i==52&&j==4)||(i==34&&(j==29||j==17))||(i<=30&&is>=26&&(j<=10&&j>=6))||(i==25&&(j==92))||(i<=24&&i>=21&&j<=55&&j>=49)||(i==19&&(j==84||j==83))||(i==4&&j==63)||(i==4&&j==63)||(i==16&&j==63)||(i==24&&j==65)||(i==33&&j==69)||(i==35&&j==72)||(i==38&&j==19)||(i==38&&j==96)||(i==42&&j==75)||(i==49&&j==6)||(i==52&&j==27)||(i==66&&j==52)||(i==80&&j==22)||(i==96&&j==57)||(i==100&&j==59)))||(k==0&&((i<=97&&i>=92&&j<=42&&j>=29)||(i==86&&j==36)||(i<=55&&i>=51&&j<=35&&j>=31)||(i==52&&j==11)||(i==41&&j==60)||(i==46&&j==67)||(i==54&&j==94)||(i==55&&j==36)||(i==64&&j==41)||(i==71&&j==67)||(i==76&&j==61)||(i==82&&j==36)||(i==92&&j==27))));
			//kill crazy crystals for Vector NTuples
			if(k==0&&((i==9&&j==33)||(i==47&&j==78))) continue;
			if(k==1&&( (i==39&&j==46)||(i==52&&j==4)||(i==78&&j==76)||(i==83&&j==21)||(i==84&&j==24)||(i==86&&j==24)||(i==94&&j==64)||(i==19&&j==84)||(i==24&&j==77) ) ) continue;
			//no values
			if(k==1&&(i==91||i==92)&&(j>=21&&j<=25) ) continue;
			if(bad2012) continue;
			if(Onoisecrystals || Onoisecrystals2) continue;
			if(Onoisecrystals3) continue;//first orange kill
			if((period==52||period==51||period==5)&& Onoisecrystalsp5) continue;
			//kill crystals where we have constant period
		//	if(hasbadperiod[i][j][k] && !(risingcrystals || partiallyrising || partiallyrising2)) {
		//	cout << "crystal ix:iy:cap " << i << ":" << j << ":" << k << " has period of constant VPT/PN --> VETO" << endl;
		//	continue; }
			if(hasbadperiod[i][j][k]) continue;
		for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){

			double x,y;
			g_vptpn_las[i][j][k]->GetPoint(n,x,y);
			if(y==0) continue;
			if(y<0.55) continue;
			if(y>1.2)  continue;
			//1st period
			if(period==1  && x<1300200000) continue;//15/03/2011 	
			if(period==1  && x>1301210000) continue;//27/03/2011	
			//2nd period
			if(period==2  && x<1302590000) continue;//12/04/2011	
			if(period==2  && x>1304620000) continue;//05/05/2011	
			//3rd period
			if(period==3  && x<1305280000) continue;//13/05/2011	
			if(period==3  && x>1309540000) continue;//01/07/2011	
			//4th period
			if(period==4  && x<1310710000) continue;//15/07/2011	
			if(period==4  && x>1314270000) continue;//25/08/2011	
			//5th period
			if(period==5  && x<1315120000) continue;//04/09/2011	
			if(period==5  && x>1320250000) continue;//02/11/2011	
			//6th period
			if(period==6  && x<1333740000) continue;//06/04/2012
			if(period==6  && x>1334910000) continue;//20/04/2012
			//7th period
			if(period==7  && x<1336230000) continue;//05/05/2012
			if(period==7  && x>1340000000) continue;//18/06/2012
			//original period
			if(period==0  && x<1314920000) continue;//04/10/2011
			if(period==0  && x>1317690000) continue;//02/09/2011
			//4th period - first 1.1 fb-1
			if(period==41 && x<1310710000) continue;//15/07/2011	
			if(period==41 && x>1312800923) continue;//08/08/2011	
			//5th period - first 1.1 fb-1
			if(period==51 && x<1315120000) continue;//04/09/2011	
			if(period==51 && x>1316813137) continue;//23/09/2011	
			//3rd period - first 950 pb-1
			if(period==32 && x<1305280000) continue;//13/05/2011	
			if(period==32 && x>1308900000) continue;//24/06/2011
			//4th period - first 950 pb-1
			if(period==42 && x<1310710000) continue;//15/07/2011	
			if(period==42 && x>1312680000) continue;//07/08/2011
			//5th period - first 950 pb-1
			if(period==52 && x<1315120000) continue;//04/09/2011	
			if(period==52 && x>1316695000) continue;//22/09/2011
			//7th period - first 950 pb-1
			if(period==72 && x<1336230000) continue;//05/05/2012
			if(period==72 && x>1337260000) continue;//17/05/2012
			//-1th period - technical stop Sept.2011
			if(period==-1 && x<1313972000) continue;//22/08/2012
			if(period==-1 && x>1315405000) continue;//07/09/2012
			//if(risingcrystals || partiallyrising) continue;
			//these should be dummies now
			if(k==0&&((i==9&&j==33)||(i==47&&j==78))) continue;
			if(k==1&&( (i==39&&j==46)||(i==52&&j==4)||(i==78&&j==76)||(i==83&&j==21)||(i==84&&j==24)||(i==86&&j==24)||(i==94&&j==64)||(i==19&&j==84)||(i==24&&j==77) ) ) continue;
			//no values
			if(k==1&&(i==91||i==92)&&(j>=21&&j<=25) ) continue;
			if(bad2012) continue;
			if(Onoisecrystals || Onoisecrystals2) continue;
			if(Onoisecrystals3) continue;//first orange kill
			if((period==52||period==51||period==5)&& Onoisecrystalsp5) continue;
//			if(y>1.1 && normalized) continue;//veto outliers	 --> check this
//			if(y<0.565 && normalized) continue;			 --> check this

			if(y>vptpnmax[i][j][k]) { vptpnmax[i][j][k] = y; vptpnmaxindex[i][j][k] = n; timemax[i][j][k] = x; }
			if(y<vptpnmin[i][j][k]) { vptpnmin[i][j][k] = y; vptpnminindex[i][j][k] = n; timemin[i][j][k] = x; }

		}

		//take maximal value and average 5 entries around it (same for minimum)
		//new 2012/08/29 average 5 entries to get minimum/maximum value to be save against little fluctuations
		if(vptpnmaxindex[i][j][k]>=0){
			double y1, x1; int num1 = 0; int n1;
			double sum = 0;
			if(vptpnmaxindex[i][j][k]<2){//maximal value is at 0th entries
				n1 = 0;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==g_vptpn_las[i][j][k]->GetN()-1) break;
					++n1;
				}
				if(num1!=0) vptpnmax[i][j][k] = sum/double(num1);
			}
			else if(vptpnmaxindex[i][j][k]>g_vptpn_las[i][j][k]->GetN()-3){// maximal value is at last entries
				n1 = g_vptpn_las[i][j][k]->GetN()-1;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==0) break;
					--n1;
				}
				if(num1!=0) vptpnmax[i][j][k] = sum/double(num1);
			}
			else{ //start 2 entries before maximum value --> maximum value = average of 5 entries around maximum value
				n1 = vptpnmaxindex[i][j][k]-2;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==g_vptpn_las[i][j][k]->GetN()-1) break;
					++n1;
				}
				if(num1!=0) vptpnmax[i][j][k] = sum/double(num1);
			}
		}
		if(vptpnminindex[i][j][k]>=0){
			double y1, x1; int num1 = 0; int n1;
			double sum = 0;
			if(vptpnminindex[i][j][k]<2){// minimum value is at 0th entries
				n1 = 0;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==g_vptpn_las[i][j][k]->GetN()-1) break;
					++n1;
				}
				if(num1!=0) vptpnmin[i][j][k] = sum/double(num1);
			}
			else if(vptpnminindex[i][j][k]>g_vptpn_las[i][j][k]->GetN()-3){// minimum value is at last entries
				n1 = g_vptpn_las[i][j][k]->GetN()-1;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==0) break;
					--n1;
				}
				if(num1!=0) vptpnmin[i][j][k] = sum/double(num1);
			}
			else{ //start 2 entries before minimum value --> minimum value = average of 5 entries minimum maximum value
				n1 = vptpnminindex[i][j][k]-2;
				while(num1<15){
					g_vptpn_las[i][j][k]->GetPoint(n1,x1,y1);
					if(y1==0) continue;
					sum += y1;
					++num1;
					if(n1==g_vptpn_las[i][j][k]->GetN()-1) break;
					++n1;
				}
				if(num1!=0) vptpnmin[i][j][k] = sum/double(num1);
			}
		}

		//plot each crystal VPT/PN vs. time (each meaning each crystal having amu_std value)
		if(plottingall_VPTPN && g_vptpn_las[i][j][k]->GetN()>10/* && vptpnmmustd[i][j][k]>0*/ ){//with mustd oled instead of las
			if(k==0&&((i==9&&j==33)||(i==47&&j==78))) continue;
			if(k==1&&( (i==39&&j==46)||(i==52&&j==4)||(i==78&&j==76)||(i==83&&j==21)||(i==84&&j==24)||(i==86&&j==24)||(i==94&&j==64)||(i==19&&j==84)||(i==24&&j==77) ) ) continue;
			//no values
			if(k==1&&(i==91||i==92)&&(j>=21&&j<=25) ) continue;
			if(bad2012) continue;
			if(Onoisecrystals || Onoisecrystals2 ||Onoisecrystals3) continue;//first orange kill
			//kill crystals where we have constant period
			if((risingcrystals || partiallyrising || partiallyrising2)) continue; 
			if((period==52||period==51||period==5)&& Onoisecrystalsp5) continue;
			if(hasbadperiod[i][j][k]) continue;
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
	//	haxis->SetMinimum(0.25);//not normalized
	//	haxis->SetMaximum(0.75);//not normalized
//		if(normalized) haxis->SetMinimum(0.565);//normalized, las 0.565, orange 0.55
//		else           haxis->SetMinimum(0.25);
//		if(normalized) haxis->SetMaximum(1.165);//normalized, las 1.165, orange 1.30
//		else           haxis->SetMaximum(1.5);
	//	haxis->SetMinimum(0.0);//test
	//	haxis->SetMaximum(3.0);//test
				haxis->GetXaxis()->SetTitle("time"); haxis->GetXaxis()->SetTitleOffset(1.25); haxis->GetXaxis()->SetLabelSize(0.027);
				haxis->GetYaxis()->SetTitle("VPT/PN"); haxis->GetYaxis()->SetTitleOffset(1.25); haxis->GetYaxis()->SetLabelSize(0.027);
				haxis->GetXaxis()->SetTimeDisplay(1); haxis->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
				haxis->GetXaxis()->SetLabelOffset(0.02); haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); haxis->GetXaxis()->SetNdivisions(510, true);
		g_vptpn_las[i][j][k]->SetHistogram(haxis);//oled instead of las
		col->Clear();
		haxis->Draw("hist");
	//	TLine *l1 = new TLine(haxis->GetXaxis()->GetBinLowEdge(1), 0.65, haxis->GetXaxis()->GetBinUpEdge(haxis->GetXaxis()->GetNbins() ), 0.65 );
	//	 l1->SetLineStyle(3); l1->SetLineColor(kRed); l1->SetLineWidth(2);
	//	TLine *l2 = new TLine(haxis->GetXaxis()->GetBinLowEdge(1), 1.25, haxis->GetXaxis()->GetBinUpEdge(haxis->GetXaxis()->GetNbins() ), 1.25 );
	//	 l2->SetLineStyle(3); l2->SetLineColor(kRed); l2->SetLineWidth(2);
		g_vptpn_las[i][j][k]->Draw("P");//oled instead of las
	//	l1->Draw();
	//	l2->Draw();
		TString outputfile = TString::Format("VPTPN_EE%i_IX%i_IY%i",k, i, j);//with orange
		col->Update();
		if(vptpnmmustd[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintNoEPS(col, outputfile, outputdirac);//oled instead of las
//		if(vptpnmmustd[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintEPS(col, outputfile, outputdirac);//oled instead of las
	//	delete l1; delete l2;
		}

		//make correlation plots VPT/PN(max)-VPT/PN(min) vs mu_std in eta bins
		col->Clear();
		if(vptpnmax[i][j][k]-vptpnmin[i][j][k] <= 0 && period!=-1) continue;//not for techincal stop Sept.2011
		if(vptpnmax[i][j][k]==0) continue;
		if(vptpnmin[i][j][k]==99) continue;
		if(timemax[i][j][k]>timemin[i][j][k] && period!=-1) continue;//not for techincal stop Sept.2011
		if(risingcrystals || partiallyrising || partiallyrising2){
		//	cout << "crystal ix:iy:iz " << i << ":" << j << ":" << k << " has max value " << vptpnmax[i][j][k] << " at " << int(timemax[i][j][k]) << " and min value " <<  vptpnmin[i][j][k] << " at " << int(timemin[i][j][k]) << endl;
			continue;
		}
		int nn = -1;
		if(vptpnruss[i][j][k]) nn = 0;
		else nn=1;
		double mustd_ = vptpnmmustd[i][j][k];
		if(startdiff) vptpnmax[i][j][k] = 1.;//compare to start of 2011 which is by definition 1
		double vptdiff = vptpnmax[i][j][k]-vptpnmin[i][j][k];
		//if(!(normalized) && vptdiff > 0.3) cout << "vptdiff " << vptdiff << " mustd " << mustd_ << " eta " << eta << " nn " << nn << " i " << i << " j " << j << " k " << k << endl;
		if(mustd_==0) mustd_= 1;
		//check for crazy crystals
	//	if(nn==0&&mustd_>0.9&&mustd_<1.1&&fabs(eta)>=1.6&&fabs(eta)<=1.8&&vptdiff>0.2)cout<<"i:j:k="<<i<<":"<<j<<":"<<k<< " eta " << eta << " mustd " << mustd_ << " russian=0? " << nn << " vptdiff " << vptdiff << " (" << vptpnmax[i][j][k] << "-" << vptpnmin[i][j][k] << ")" << endl;
	//	cout << "eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << "  " << nn << endl;

		if(eta<(-2.8) && mustd_!=0)           h_mustd_vs_VPTPNdiff[0][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[0][nn]->GetN(), mustd_,vptdiff);//m3p0
		else if(eta<(-2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[1][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[1][nn]->GetN(), mustd_,vptdiff);
		else if(eta<(-2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff[2][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[2][nn]->GetN(), mustd_,vptdiff);//m2p6
		else if(eta<(-2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff[3][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[3][nn]->GetN(), mustd_,vptdiff);
		else if(eta<(-2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[4][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[4][nn]->GetN(), mustd_,vptdiff);//m2p2
		else if(eta<(-1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[5][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[5][nn]->GetN(), mustd_,vptdiff);
		else if(eta<(-1.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[6][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[6][nn]->GetN(), mustd_,vptdiff);//m1p8
		else if(eta<=(-1.4)&& mustd_!=0)      h_mustd_vs_VPTPNdiff[7][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[7][nn]->GetN(), mustd_,vptdiff);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_mustd_vs_VPTPNdiff[8][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[8][nn]->GetN(), mustd_,vptdiff);//1p4
		else if(eta<( 1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[9][nn] ->SetPoint(h_mustd_vs_VPTPNdiff[9][nn]->GetN(), mustd_,vptdiff);
		else if(eta<( 2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[10][nn]->SetPoint(h_mustd_vs_VPTPNdiff[10][nn]->GetN(), mustd_,vptdiff);
		else if(eta<( 2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff[11][nn]->SetPoint(h_mustd_vs_VPTPNdiff[11][nn]->GetN(), mustd_,vptdiff);
		else if(eta<( 2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff[12][nn]->SetPoint(h_mustd_vs_VPTPNdiff[12][nn]->GetN(), mustd_,vptdiff);
		else if(eta<( 2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[13][nn]->SetPoint(h_mustd_vs_VPTPNdiff[13][nn]->GetN(), mustd_,vptdiff);
		else if(eta<( 2.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[14][nn]->SetPoint(h_mustd_vs_VPTPNdiff[14][nn]->GetN(), mustd_,vptdiff);
		else if(eta<=(3.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[15][nn]->SetPoint(h_mustd_vs_VPTPNdiff[15][nn]->GetN(), mustd_,vptdiff);

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&mustd_!=0)a_mustd_vs_VPTPNdiff[0][nn]->SetPoint(a_mustd_vs_VPTPNdiff[0][nn]->GetN(), mustd_,vptdiff);//1p4
		else if(fabs(eta)<( 1.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff[1][nn]->SetPoint(a_mustd_vs_VPTPNdiff[1][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<( 2.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff[2][nn]->SetPoint(a_mustd_vs_VPTPNdiff[2][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<( 2.2) && mustd_!=0)      a_mustd_vs_VPTPNdiff[3][nn]->SetPoint(a_mustd_vs_VPTPNdiff[3][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<( 2.4) && mustd_!=0)      a_mustd_vs_VPTPNdiff[4][nn]->SetPoint(a_mustd_vs_VPTPNdiff[4][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<( 2.6) && mustd_!=0)      a_mustd_vs_VPTPNdiff[5][nn]->SetPoint(a_mustd_vs_VPTPNdiff[5][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<( 2.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff[6][nn]->SetPoint(a_mustd_vs_VPTPNdiff[6][nn]->GetN(), mustd_,vptdiff);
		else if(fabs(eta)<=(3.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff[7][nn]->SetPoint(a_mustd_vs_VPTPNdiff[7][nn]->GetN(), mustd_,vptdiff);

	}}}

	TString outputfile;
	//draw correlation plots VPT/PN vs mu_std in eta bins
	for(int n = 0; n<16; ++n){//eta bins
	for(int nn= 0; nn<2;++nn){//producer: 0: russian, 1: chinese
		TH1F *haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
	     if(!plotonlyAbsEta){
	      if(VPTdiff_vs_mustd){
		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
		haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff[n][nn]->GetName(), h_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
		else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
		haxis->Draw("hist");
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		outputfile = h_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

	      }
	     }
		if(n<8){//abs eta
		     if(VPTdiff_vs_mustd){
			col->cd();
			col->SetName(a_mustd_vs_VPTPNdiff[n][nn]->GetName() );
			col->SetTitle(a_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
			a_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
			haxis = a_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_mustd_vs_VPTPNdiff[n][nn]->GetName(), a_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.,2.);
			haxis->GetXaxis()->SetLimits(0.,2.);
			haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
			haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("#mu_{std} (m^{-1})");
			haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
			else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
			haxis->Draw("hist");
			a_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
			outputfile = a_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
		     }
		}//if(n<8)

	     double meanyEEpm = h_mustd_vs_VPTPNdiff[n][nn]->GetMean(2);
	     double   rmsEEpm = h_mustd_vs_VPTPNdiff[n][nn]->GetRMS(2);
	     if(n<8) { 
		if(nn==0) { h_VPTPNdiff_vs_eta_EEm_russ->SetBinContent(n+1  , meanyEEpm); h_VPTPNdiff_vs_eta_EEm_russ->SetBinError(n+1  , rmsEEpm); }
		else      { h_VPTPNdiff_vs_eta_EEm_chin->SetBinContent(n+1  , meanyEEpm); h_VPTPNdiff_vs_eta_EEm_chin->SetBinError(n+1  , rmsEEpm); }
	     } else { 
		if(nn==0) { h_VPTPNdiff_vs_eta_EEp_russ->SetBinContent(n+1-8, meanyEEpm); h_VPTPNdiff_vs_eta_EEp_russ->SetBinError(n+1-8, rmsEEpm); }
		else      { h_VPTPNdiff_vs_eta_EEp_chin->SetBinContent(n+1  , meanyEEpm); h_VPTPNdiff_vs_eta_EEp_chin->SetBinError(n+1  , rmsEEpm); }
	     }
	     if(n<8){
		double meanyEE = a_mustd_vs_VPTPNdiff[n][nn]->GetMean(2);
		double   rmsEE = a_mustd_vs_VPTPNdiff[n][nn]->GetRMS(2);
		if(nn==0) { h_VPTPNdiff_vs_eta_russ->SetBinContent(n+1  , meanyEE); h_VPTPNdiff_vs_eta_russ->SetBinError(n+1  , rmsEE); }
		else      { h_VPTPNdiff_vs_eta_chin->SetBinContent(n+1  , meanyEE); h_VPTPNdiff_vs_eta_chin->SetBinError(n+1  , rmsEE); }
	     }
	}
	}
	if(!plotonlyAbsEta){
		col->cd();
		col->SetName(h_VPTPNdiff_vs_eta_EEm_russ->GetName() );
		col->SetTitle(h_VPTPNdiff_vs_eta_EEm_russ->GetTitle() );
		h_VPTPNdiff_vs_eta_EEm_russ->Draw();
		outputfile = h_VPTPNdiff_vs_eta_EEm_russ->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		col->cd();
		col->SetName(h_VPTPNdiff_vs_eta_EEm_chin->GetName() );
		col->SetTitle(h_VPTPNdiff_vs_eta_EEm_chin->GetTitle() );
		h_VPTPNdiff_vs_eta_EEm_chin->Draw();
		outputfile = h_VPTPNdiff_vs_eta_EEm_chin->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_VPTPNdiff_vs_eta_EEp_russ->GetName() );
		col->SetTitle(h_VPTPNdiff_vs_eta_EEp_russ->GetTitle() );
		h_VPTPNdiff_vs_eta_EEp_russ->Draw();
		outputfile = h_VPTPNdiff_vs_eta_EEp_russ->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
		col->cd();
		col->SetName(h_VPTPNdiff_vs_eta_EEp_chin->GetName() );
		col->SetTitle(h_VPTPNdiff_vs_eta_EEp_chin->GetTitle() );
		h_VPTPNdiff_vs_eta_EEp_chin->Draw();
		outputfile = h_VPTPNdiff_vs_eta_EEp_chin->GetName();
		col->Update();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();
	}
	col->cd();
	col->SetName(h_VPTPNdiff_vs_eta_russ->GetName() );
	col->SetTitle(h_VPTPNdiff_vs_eta_russ->GetTitle() );
	h_VPTPNdiff_vs_eta_russ->Draw();
	outputfile = h_VPTPNdiff_vs_eta_russ->GetName();
	col->Update();
	Util::PrintNoEPS(col, outputfile, outputdir);
	Util::PrintEPS(col, outputfile, outputdir);
	col->Clear();
	col->cd();
	col->SetName(h_VPTPNdiff_vs_eta_chin->GetName() );
	col->SetTitle(h_VPTPNdiff_vs_eta_chin->GetTitle() );
	h_VPTPNdiff_vs_eta_chin->Draw();
	outputfile = h_VPTPNdiff_vs_eta_chin->GetName();
	col->Update();
	Util::PrintNoEPS(col, outputfile, outputdir);
	Util::PrintEPS(col, outputfile, outputdir);
	col->Clear();

	//store histos into file
	TFile *newFile;
	if(startdiff) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_diffwrt2011begin_"+periods+".root", "RECREATE");
	if(normalized) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_"+periods+".root", "RECREATE");//0312->0410; //muSIC <-> mustd
	else if(!startdiff)  newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_nonormalization_"+periods+".root", "RECREATE");//0312->0410
	newFile->cd();
	for(int n = 0; n<16; ++n){ for(int nn= 0; nn<2;++nn) {
	h_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	if(n<8){
	a_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	}
	h_VPTPNdiff_vs_eta_EEm_russ->Write();
	h_VPTPNdiff_vs_eta_EEp_russ->Write();
	h_VPTPNdiff_vs_eta_russ    ->Write();
	h_VPTPNdiff_vs_eta_EEm_chin->Write();
	h_VPTPNdiff_vs_eta_EEp_chin->Write();
	h_VPTPNdiff_vs_eta_chin    ->Write();
	}}
	cout << "histograms for period " << period << " called " << periods << " (normalized/startdiff = " << normalized << "/" << startdiff << ") stored in " << newFile->GetName() << endl;

}

