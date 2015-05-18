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
float avg(vector<float> x){

  float sum = 0;
  for (Int_t i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}
*/
void EEstudies_continued_SIC_v2(){

	bool plottingall_VPTPN = false;
	bool normalized = true;

	gROOT->SetStyle("Plain");

	//NOTE: Also change possibly newFileName at the very bottom!!!!!
	char gname[101];
	TString outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2c/SIC/notnormalized";//normalized";//normalized if normalize VPT/PN
	if(normalized) outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2c/SIC/normalized/withSashaUnScaled";
//	cout << __LINE__<<endl;
//	Util::MakeOutputDir(outputdir);
//	cout << __LINE__<<endl;
	TString outputdirac = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2c/SIC/normalized/allcrystalswithmustd";//normalized if normalize VPT/PN
	if(normalized) outputdirac = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/v2c/SIC/notnormalized/allcrystalswithmustd";
//	if(plottingall_VPTPN)	Util::MakeOutputDir(outputdirac);

//	TFile *etaFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_eta_in_ix_iy.root");
	TFile *etaFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp                = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm                = (TH2F*)etaFile->Get("eta_EEm");

	//TFile *stdmapsfile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20111219_EEplots.root");
	TFile *stdmapsfile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *mu_ECAL_EEp_chin = (TH2D*)stdmapsfile->Get("mu_SIC_EEp_chin");//in Get: SIC <-> ECAL
	TH2D *mu_ECAL_EEm_chin = (TH2D*)stdmapsfile->Get("mu_SIC_EEm_chin");//in Get: SIC <-> ECAL

	TGraph* h_mustd_vs_VPTPNdiff[16][2];
	TGraph* h_mustd_vs_VPTPNdiff_scaled[16][2];
	TGraph* h_orange_mustd_vs_VPTPNdiff[16][2];
	TGraph* h_orange_mustd_vs_VPTPNdiff_scaled[16][2];
	TGraph* a_mustd_vs_VPTPNdiff[8][2];
	TGraph* a_mustd_vs_VPTPNdiff_scaled[8][2];
	TGraph* a_orange_mustd_vs_VPTPNdiff[8][2];
	TGraph* a_orange_mustd_vs_VPTPNdiff_scaled[8][2];
	string string_eta[16] = {"m3p0", "m2p8", "m2p6", "m2p4", "m2p2", "m2p0", "m1p8", "m1p6", "1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	for(int n = 0; n<16; ++n){
	for(int nn = 0; nn<2; ++ nn){
	string hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn];
	string hname_s = (string)"mustd__vs__VPTPNmax-VPTPNmin__eta_" + string_eta[n] + string_prod[nn] + (string)"__scaled";

	h_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
	h_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
	h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
	//if(normalized) h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.3);
	//else           h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.22);
	h_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //h_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
	h_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPTPNmax-VPTPNmin"); h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("mu_std (m^-1)");

	h_mustd_vs_VPTPNdiff_scaled[n][nn] = new TGraph();
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetNameTitle(hname_s.c_str(), hname_s.c_str());
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
	//if(normalized) h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetRangeUser(0.,0.15);
	//else           h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetRangeUser(0.,0.22);
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerStyle(4);// h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerSize(4);
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetTitle("VPTPNmax-VPTPNmin (scaled)"); h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetTitle("mu_std (m^-1)");


	h_orange_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
	h_orange_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(((string)"Orange__" + hname).c_str(), ((string)"Orange__" + hname).c_str());
	h_orange_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
//	if(normalized) h_orange_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.3);
//	else           h_orange_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetRangeUser(0.,0.22);
	h_orange_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //h_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
	h_orange_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPTPNmax-VPTPNmin"); h_orange_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("mu_std (m^-1)");

	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn] = new TGraph();
	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetNameTitle(((string)"Orange__" + hname_s).c_str(), ((string)"Orange__" + hname_s).c_str());
	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
//	if(normalized) h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetRangeUser(0.,0.15);
//	else           h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetRangeUser(0.,0.22);
	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerStyle(4);// h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerSize(4);
	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetTitle("VPTPNmax-VPTPNmin (scaled)"); h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetTitle("mu_std (m^-1)");

	if(n<8){
		hname = (string)"mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_feta[n] + string_prod[nn];
		hname_s = (string)"mustd__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_feta[n] + string_prod[nn] + (string)"__scaled";
		a_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
		a_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
		a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
	
		a_mustd_vs_VPTPNdiff_scaled[n][nn] = new TGraph();
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetNameTitle(hname_s.c_str(), hname_s.c_str());
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerStyle(4);// a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerSize(4);
		a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min} (scaled)"); a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
	
	
		a_orange_mustd_vs_VPTPNdiff[n][nn] = new TGraph();
		a_orange_mustd_vs_VPTPNdiff[n][nn]->SetNameTitle(((string)"Orange__" + hname).c_str(), ((string)"Orange__" + hname).c_str());
		a_orange_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_orange_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_orange_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_orange_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
	
		a_orange_mustd_vs_VPTPNdiff_scaled[n][nn] = new TGraph();
		a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetNameTitle(((string)"Orange__" + hname_s).c_str(), ((string)"Orange__" + hname_s).c_str());
		a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerStyle(4);// a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMarkerSize(4);
		a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min} (scaled)"); a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
	}

	}
	}

//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120202_EE_VPT_over_PN.root");
	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN.root");
//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN_new2sorted.root");
	TGraph* g_vptpn_las[101][101][2];
	TGraph* g_vptpn_oled[101][101][2];
	double mustd_russ = 0; double mustd_chin = 0;
	//start normalized
	double vptpnnorm[101][101][2];
	double vptpnnormorange[101][101][2];

	double vptpnmax[101][101][2];
	double vptpnmin[101][101][2];
	double vptpnmaxorange[101][101][2];
	double vptpnminorange[101][101][2];
	double vptpnmmustd[101][101][2];
	bool   vptpnruss[101][101][2];
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		vptpnnorm[i][j][k] = 0; 
		vptpnnormorange[i][j][k] = 0;
		vptpnmax[i][j][k] = 0;
		vptpnmin[i][j][k] = 99; 
		vptpnmaxorange[i][j][k] = 0;
		vptpnminorange[i][j][k] = 99; 
		vptpnmmustd[i][j][k] = 0; 
		vptpnruss[i][j][k] = false;
	}}}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "start computation" << endl;
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		if (k==0) sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEp", i, j);
		g_vptpn_oled[i][j][k] = (TGraph*)vpt_values->Get(gname);
		if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEp", i, j);
		g_vptpn_las[i][j][k] = (TGraph*)vpt_values->Get(gname);

		if(g_vptpn_las[i][j][k]->GetN()>0){
			double norm = 0; int countnorm = 0;
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y, xprev,yprev, xpost,ypost, xprev2,yprev2;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
			//	if(y>1000e6) { g_vptpn_las[i][j][k]->RemovePoint(n); continue;}
				if(!plottingall_VPTPN && x>1317690000) continue;//04/10 daystart	//timesaving
				if(!plottingall_VPTPN && x<1314920000) continue;//02/09 daystart	//timesaving
				if(n==g_vptpn_las[i][j][k]->GetN()-1) g_vptpn_las[i][j][k]->GetPoint(n,xpost,ypost);
				else g_vptpn_las[i][j][k]->GetPoint(n+1,xpost,ypost);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev,yprev);
				else g_vptpn_las[i][j][k]->GetPoint(n-1,xprev,yprev);
				if(n==0) g_vptpn_las[i][j][k]->GetPoint(n,xprev2,yprev2);
				else if(n==1) g_vptpn_las[i][j][k]->GetPoint(n-1,xprev2,yprev2);
				else g_vptpn_las[i][j][k]->GetPoint(n-2,xprev2,yprev2);

				if(y==0) continue;
				//if(x>1315398150) continue;
				if((yprev>y*1.1 && ypost>y*1.1) || (yprev<y*0.9 && ypost<y*0.9)) g_vptpn_las[i][j][k]->SetPoint(n,x,yprev);
				//upper 2cleaning
				if(ypost<y*1.1 && yprev2<y*1.1 && ypost<yprev*1.1 && yprev2<yprev*1.1){
					g_vptpn_las[i][j][k]->SetPoint(n,x,yprev2);//clean only n
					if(n>0) g_vptpn_las[i][j][k]->SetPoint(n-1,xprev,yprev2);//clean only n
				}
			}//cleaning
		//start normalized
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(!plottingall_VPTPN && x>1317690000) continue;//04/10 daystart	//timesaving
				if(!plottingall_VPTPN && x<1314920000) continue;//02/09 daystart	//timesaving
				//cout << "t=" << x << " VPT/PN=" << y << endl;
				if(y==0) continue;
				if(x>1315398150) continue;
				norm += y;
				++countnorm;
				if(countnorm>5) break;
			}
			//if(countnorm==0) cout << "wfk i.j.k = " << i << ":"<< j << ":" << k << endl;
			vptpnnorm[i][j][k] = norm/countnorm;
			//do normalization
			if(countnorm!=0){
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				//cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" and t=" << x << " had VPT/PN = " << y << " which is now normalized to ";
				g_vptpn_las[i][j][k]->SetPoint(n,x,y/(vptpnnorm[i][j][k]));
				//cout << y << endl;
				if(y/(vptpnnorm[i][j][k])>1.5){
					g_vptpn_las[i][j][k]->GetPoint(n,x,y);
					double xx,yy;
					g_vptpn_las[i][j][k]->GetPoint(n+1,xx,yy);//already normalized
					g_vptpn_las[i][j][k]->SetPoint(n,x,yy);
				}
			}
			}
			else {cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" could not normalized " << countnorm << " number of entries " << g_vptpn_las[i][j][k]->GetN() << endl;}
		}
		if(g_vptpn_oled[i][j][k]->GetN()>0){
			double norm = 0; int countnorm = 0;
			for(int n = g_vptpn_oled[i][j][k]->GetN()-1; n>=0; --n){
				double x,y, xprev,yprev, xpost,ypost, xprev2,yprev2;
				g_vptpn_oled[i][j][k]->GetPoint(n,x,y);
				if(!plottingall_VPTPN && x>1317690000) continue;//04/10 daystart	//timesaving
				if(!plottingall_VPTPN && x<1314920000) continue;//02/09 daystart	//timesaving
				if(n==g_vptpn_oled[i][j][k]->GetN()-1) g_vptpn_oled[i][j][k]->GetPoint(n,xpost,ypost);
				else g_vptpn_oled[i][j][k]->GetPoint(n+1,xpost,ypost);
				if(n==0) g_vptpn_oled[i][j][k]->GetPoint(n,xprev,yprev);
				else g_vptpn_oled[i][j][k]->GetPoint(n-1,xprev,yprev);
				if(n==0) g_vptpn_oled[i][j][k]->GetPoint(n,xprev2,yprev2);
				else if(n==1) g_vptpn_oled[i][j][k]->GetPoint(n-1,xprev2,yprev2);
				else g_vptpn_oled[i][j][k]->GetPoint(n-2,xprev2,yprev2);

				if(y==0) continue;
				//if(x>1315398150) continue;
				if((yprev>y*1.1 && ypost>y*1.1) || (yprev<y*0.9 && ypost<y*0.9)) g_vptpn_oled[i][j][k]->SetPoint(n,x,yprev);
				//upper 2cleaning
				if(ypost<y*1.1 && yprev2<y*1.1 && ypost<yprev*1.1 && yprev2<yprev*1.1){
					g_vptpn_oled[i][j][k]->SetPoint(n,x,yprev2);//clean only n
					if(n>0) g_vptpn_oled[i][j][k]->SetPoint(n-1,xprev,yprev2);//clean only n
				}
			}//cleaning
		//start normalized
			for(int n = g_vptpn_oled[i][j][k]->GetN()-1; n>=0; --n){
				double x,y;
				g_vptpn_oled[i][j][k]->GetPoint(n,x,y);
				if(!plottingall_VPTPN && x>1317690000) continue;//04/10 daystart	//timesaving
				if(!plottingall_VPTPN && x<1314920000) continue;//02/09 daystart	//timesaving
				if(y==0) continue;
				if(x>1315008150) continue;//00 <->39
				norm += y;
				++countnorm;
				if(countnorm>5) break;
			}
			//if(countnorm==0) cout << "wfk i.j.k = " << i << ":"<< j << ":" << k << endl;
			vptpnnormorange[i][j][k] = norm/countnorm;
			//cout << "i.j.k = " << i << ":"<< j << ":" << k << " norm " << vptpnnormorange[i][j][k] << endl;
			//do normalization
			if(countnorm!=0){
			for(int n = g_vptpn_oled[i][j][k]->GetN()-1; n>=0; --n){
				double x,y;
				g_vptpn_oled[i][j][k]->GetPoint(n,x,y);
				g_vptpn_oled[i][j][k]->SetPoint(n,x,y/(vptpnnormorange[i][j][k]));
				//cout << " time " << x << " yprev:ypost " << y << " : " << y/(vptpnnormorange[i][j][k]);
				/*if(y/(vptpnnormorange[i][j][k])>1.5 && n !=g_vptpn_oled[i][j][k]->GetN()-1){
					double xx,yy;
					g_vptpn_oled[i][j][k]->GetPoint(n,x,y);
					g_vptpn_oled[i][j][k]->GetPoint(n+1,xx,yy);//already normalized
					if(yy<1.5) g_vptpn_oled[i][j][k]->SetPoint(n,x,yy);
					cout << " - ypostpost " << yy;
					if(yy>10e6) g_vptpn_oled[i][j][k]->RemovePoint(n);
				}*/
				//cout << endl;
			}
			}
		}



		int bin; float eta; string seta; stringstream russchin; string rs; stringstream etass; string es;
		if(k==0){
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEm_chin->FindBin(i,j);
			mustd_russ = 0.;
			mustd_chin = mu_ECAL_EEm_chin->GetBinContent(bin);

			if(mustd_russ==0 && mustd_chin==0) continue;
			else if(mustd_russ!=0){
				russchin << mustd_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_russ;
				vptpnruss[i][j][k] = true;
			}
			else{
				russchin << mustd_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_chin;
				vptpnruss[i][j][k] = false;
			}
		}
		else{
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEp_chin->FindBin(i,j);
			mustd_russ = 0.;
			mustd_chin = mu_ECAL_EEp_chin->GetBinContent(bin);
			if(mustd_russ==0 && mustd_chin==0) continue;
			else if(mustd_russ!=0){
				russchin << mustd_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_russ;
				vptpnruss[i][j][k] = true;
			}
			else{
				russchin << mustd_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_chin;
				vptpnruss[i][j][k] = false;
			}
		}
		bool noisecrystals = ( (k==1&&((i==78&&j==76)||(i==52&&j==4)||(i==34&&(j==29||j==17))||(i<=30&&i>=26&&(j<=10&&j>=6))||(i==25&&(j==92))||(i<=24&&i>=21&&j<=55&&j>=49)||(i==19&&(j==84||j==83))||(i==4&&j==63)||(i==4&&j==63)||(i==16&&j==63)||(i==24&&j==65)||(i==33&&j==69)||(i==35&&j==72)||(i==38&&j==19)||(i==38&&j==96)||(i==42&&j==75)||(i==49&&j==6)||(i==52&&j==27)||(i==66&&j==52)||(i==80&&j==22)||(i==96&&j==57)||(i==100&&j==59)||(i==4&&j==40)||(i>=4&&i<=5&&j>=61&&j<=65)||(i==7&&j==28)||(i==8&&j==27)||(i==9&&j==23)||(i==9&&j==29)||(i==16&&j==84)||(i==16&&j>=85&&j<=86)||(i==17&&j==81)||(i==18&&j>=81&&j<=83)||(i==18&&j==85)||(i==19&&j>=81&&j<=85)||(i==20&&(j==81||j==83))||(i==21&&j==92)||(i==25&&j==92)||(i==26&&j==92)||(i==27&&j==91)||(i==29&&j==94)||(i==30&&j>=93&&j<=94)||(i==32&&j==53)||(i==37&&j>=39&&j<=40)||(i==39&&j==46)||(i==39&&j==96)||(i==52&&j==96)||(i==53&&j==35)||(i==55&&j==100)||(i==58&&j==38)||(i==62&&j==4)||(i==66&&j==95)||(i==68&&j==93)||(i==68&&j==95)||(i==72&&j==6)||(i==72&&j==95)||(i==73&&j==6)||(i==73&&j==8)||(i==73&&j==10)||(i==73&&j<=91&&j<=95)||(i==75&&j==91)||(i==75&&j>=94&&j<=95)||(i==76&&j==91)||(i==78&&j==10)||(i==79&&j==9)||(i==80&&j==10)||(i==81&&j==17)||(i==81&&j==84)||(i==81&&j==86)||(i==82&&j==18)||(i==82&&j>=82&&j<=86)||(i>=83&&i<=84&&j==19)||(i==83&&j==82)||(i==84&&j>=83&&j<=84)||(i>=84&&i<=85&&j==86)||(i>=86&&i<=87&&j==81)||(i==87&&j>=82&&j<=84)||(i==91&&j==31)||(i==91&&j==66)||(i==91&&j==69)||(i==91&&j>=71&&j<=73)||(i==91&&j==75)||(i==92&&j==28)||(i==92&&j==71)||(i==92&&j==75)||(i==93&&j==69)||(i==93&&j==71)||(i==94&&j>=27&&j<=28)||(i==94&&j>=70&&j<=71)||(i==95&&j>=28&&j<=29)||(i==95&&j>=73&&j<=74) ) ) || (k==0&&((i<=100&&i>=91&&j<=45&&j>=26)||(i==86&&j==36)||(i<=59&&i>=51&&j<=40&&j>=31)||(i==52&&j==11)||(i==41&&j==60)||(i==46&&j==67)||(i==54&&j==94)||(i==55&&j==36)||(i==64&&j==41)||(i==71&&j==67)||(i==76&&j==61)||(i==82&&j==36)||(i==92&&j==27)||(i==10&&j==30)||(i==12&&j>=22&&j<=23)||(i==16&&j==52)||(i==20&&j==53)||(i==29&&j==10)||(i==33&&j==50)||(i==34&&j>=46&&j<=48)||(i==41&&j==58)||(i==42&&j==60)||(i==61&&j==96)||(i==69&&j==93)||(i==72&&j==6)||(i==72&&j==92)||(i==73&&j==93)||(i==74&&j>=9&&j<=10)||(i==75&&j>=6&&j<=7)||(i==76&&j==91)||(i==77&&j==91)||(i==78&&j>=15&&j<=91)||(i==79&&j==14)||(i==81&&j==15)||(i==81&&j==87)||(i==82&&j==86)||(i==86&&j==19)||(i==90&&j==22)||(i==92&&j==78)||(i==93&&j==74)||(i>=94&&i<=95&&j>=73&&j<=74)||(i==97&&j==62)||(i==99&&j>=58&&j<=60)||(i==100&&j==56) ) ) );
		bool specialcrystals = false;
		for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){
			double x,y;
			g_vptpn_las[i][j][k]->GetPoint(n,x,y);
			if(y==0) continue;
			if(x<1314920000) continue;//02/09 daystart
			if(x>1317690000) continue;//04/10 daystart
			if(y>1.1 && normalized && !(specialcrystals) ) continue;//veto outliers
			if(y<0.565 && normalized && !(specialcrystals) ) continue;
			//special crystals
			if(specialcrystals){

			}
			//kill crazy crystals
		//wide seems to have two levels parallel
			if(i==17&&j==16&&k==1&&y<0.9) continue;
			if(noisecrystals) continue;
			if(i==4&&j==64&&k==0&&y>1.08&&x<1315000000.) continue;
			if(i==6&&j==38&&k==0&&y>1.05&&x<1315000000.) continue;
		//	if(i==5&&j==65&&k==0) continue;
		//	if(i==13&&j==25&&k==0) continue;
		//	if(i==16&&j==14&&k==0) continue;
		//	if(i==17&&j==65&&k==0) continue;
		//	if(i==18&&j==87&&k==0) continue;
		//	if(i==26&&j==7&&k==0) continue;
		//	if(i==35&&j==50&&k==0) continue;
		//	if(i==37&&j>=46&&j<=65&&k==0) continue;
		//	if(i==50&&j>=38&&j<=39&&k==0) continue;
		//	if(i>=52&&i<=55&&j>=62&&j<=66&&k==0) continue;
		//	if(i>=60&&i<=61&&j<=36&&j>=41&&k==0) continue;
		//	if(i>=61&&i<=62&&j<=44&&j>=63&&k==0) continue;
		//	if(i>=64&&i<=64&&j>=46&&j<=60&&k==0) continue;
		//	if(i==65&&j<=63&&j>=64&&k==0) continue;
			//if(i==72&&j==91&&k==0) continue;
		//	if(i==32&&j>=49&&j>=55&&k==1) continue;
		//	if(i==53&&j>=36&&j<=70&&k==1) continue;
		//	if(i==67&&j==92&&k==1) continue;
		//	if(i==70&&j==92&&k==1) continue;
		//	if(i==70&&j==93&&k==1) continue;




/*			if(i==78&&j==76&&k==1) continue;//noise
			//if(i==73&&(j<=50&&j>=40)&&k==1) continue;
			if(i==52&&j==4&&k==1) continue;//noise
			if(i==39&&j==46&&k==1) continue;
			if(i==34&&(j==29||j==17)&&k==1) continue;//very strange peak
			if(i<=30&&i>=26&&(j<=10&&j>=6)&&k==1) continue;//jumps
			if(i==25&&(j==92)&&k==1) continue;//rising up
		//	if(i==24&&j==77&&k==1) continue;//very strange wideness
			if(i<=24&&i>=21&&j<=55&&j>=49&&k==1) continue;//strange peak
			if(i==19&&(j==84||j==83)&&k==1) continue;
		//	if(i==17&&j<=39&&j>=37&&k==1) continue;
		//	if(i==9&&j==75&&k==1) continue;
		//	if(i==8&&j==27&&k==1) continue;
		//	if(i<=8&&i>=6&&j==48&&k==1) continue;
			if(i==4&&j==63&&k==1) continue;
			if(i<=97&&i>=92&&j<=42&&j>=29&&k==0) continue;//noise
		//	if(i==92&&j==78&&k==1) continue;
			if(i==86&&j==36&&k==0) continue;
		//	if(i==84&&j==15&&k==0) continue;
		//	if(i==83&&j==30&&k==0) continue;
		//	if(i==81&&j==87&&k==0) continue;
			//if(i=79&&j==80&&k==0) continue;
			//if(i=78&&j==91&&k==0) continue;
			//if(i=75&&j==77&&k==0) continue;
			//if(i=72&&j==93&&k==0) continue;
			//if(i=77&&j==8&&k==0) continue;
			if(i<=55&&i>=51&&j<=35&&j>=31&&k==0) continue;//noise
			if(i==52&&j==11&&k==0) continue;
		//	if(i==31&&j==10&&k==0) continue;
		//	if(i==26&&j==7&&k==0) continue;
			//if(i==25&&j==18&&k==0) continue;
		//	if(i==21&&j==58&&k==0) continue;
		//	if(i==18&&j==21&&k==0) continue;
		//	if(i==13&&j==25&&k==0) continue;
		//	if(i==3&&j==44&&k==0) continue;
			//no normalization possible
			if(i==41&&j==60&&k==0) continue;
			if(i==46&&j==67&&k==0) continue;
			if(i==54&&j==94&&k==0) continue;
			if(i==55&&j==36&&k==0) continue;
			if(i==64&&j==41&&k==0) continue;
			if(i==71&&j==67&&k==0) continue;
			if(i==76&&j==61&&k==0) continue;
			if(i==82&&j==36&&k==0) continue;
			if(i==92&&j==27&&k==0) continue;
			if(i==4 &&j==63&&k==1) continue;
			if(i==16&&j==63&&k==1) continue;
			if(i==24&&j==65&&k==1) continue;
			if(i==33&&j==69&&k==1) continue;
			if(i==35&&j==72&&k==1) continue;
			if(i==38&&j==19&&k==1) continue;
			if(i==38&&j==96&&k==1) continue;
			if(i==42&&j==75&&k==1) continue;
			if(i==49&&j==6 &&k==1) continue;
			if(i==52&&j==27&&k==1) continue;
			if(i==66&&j==52&&k==1) continue;
			if(i==80&&j==22&&k==1) continue;
			if(i==96&&j==57&&k==1) continue;
			if(i==100&&j==59&&k==1) continue;
*/
			if(y>vptpnmax[i][j][k]) vptpnmax[i][j][k] = y;
			if(y<vptpnmin[i][j][k]) vptpnmin[i][j][k] = y;

		}
		bool Ospecialcrystals = ((k==0&&( (i==24&&j==51)||(i==26&&j==50)||(i==32&&j==56)||(i==33&&j==57)||(i==34&&j==48)||(i==36&&j==40)||(i==36&&j==48)||(i==38&&j==50)||(i==38&&j==63)||(i==39&&j==49)||(i==39&&j==50)||(i==42&&j==39)||(i==43&&j==41)||(i==47&&j==36)||(i==49&&j==36)||(i==49&&j==39)||(i==52&&j==11)||(i==53&&j==63)||(i>=58&&i<=59&&j>=60&&j<=61)||(i==61&&j==57)||(i<=70&&i>=62&&j>=46&&j<=58) ) ) || (k==1&&((i==24&&j==50)||(i>=36&&i<=37&&j>=45&&j<=60)||(i==38&&j==64)||(i==39&&j>=55&&j<=59)||(i==39&&j==64)||(i==40&&j==63)||(i==45&&j==39)||(i==45&&j==39)||(i==55&&j==63)||(i==56&&j>=38&&j<=39)||(i>=56&&i<=57&&j>=28&&j<=62)||(i>=59&&i<=60&&j>=39&&j<=40)||(i==61&&j==65)||(i==62&&j==51)||(i==62&&j>=54&&j<=55)||(i==65&&j==51) ) ) );
		bool O15 = (k==0&&((i==26&&j==50)||(i==36&&j==48)));
		bool O14 = (k==1&&((i==24&&j==50)));
		bool O135= ((k==0&&((i==24&&j==51))));
		bool O13 = ((k==0&&((i==32&&j==56)||(i==38&&j==63)||(i==52&&j==11)||(i==53&&j>=62&&j<=68)||(i==55&&j==63)||(i>=58&&i<=59&&j==60)||(i==58&&j==61)||(i==60&&j==63)||(i==61&&j==57)))||(k==1&&((i==39&&j>=55&&j<=59)||(i==39&&j==64)||(i==39&&j==64)||(i>=59&&i<=60&&j>=39&&j<=40))));
		bool O125= (k==1&&((i==45&&j==39)||(i==56&&j==38)));
		bool O12 = Ospecialcrystals && !(O15)&& !(O14)&& !(O135) && !(O13)&& !(O125);
		bool O75 = ((k==0&&((i==26&&j==50))));
		bool O65 = (k==1&&((i==45&&j==39)||(i==56&&j==38)));
		bool O6  = ((k==0&&((i==24&&j==51)||(i==34&&j==48)||(i==38&&j==50)||(i==47&&j==36)||(i==58&&j==61)||(i==60&&j==63)||(i==61&&j==57)||(i>=62&&i<=70&&j>=49&&j<=58)))||(k==1&&((i==24&&j==50)||(i==38&&j==64)||(i==39&&j>=55&&j<=59)||(i==62&&j==51)||(i==56&&j==39)||(i==39&&j==64))));
		bool O5  = ((k==0&&((i==36&&j==48)||(i==49&&(j==36||j==39)))));
		bool O7  = Ospecialcrystals && !(O75)&& !(O65)&& !(O6) && !(O5);
		bool Onoisecrystals = ((k==0&&((i==40&&j==40)||(i==50&&j==39)||(i>=54&&i<=55&&j>=31&&j<=35)||(i==51&&j==63)||(i==54&&(j==65&&j==86))||(i>=61&&i<=74&&j>=91&&j<=97)||(i>=76&&i<=85&&j>=56&&j<=71)||(i>=92&&i<=97&&j>=27&&j<=42)))||(k==1&&((i>=1&&i<=4&&j>=41&&j<=50)||(i>=6&&i<=10&&j>=27&&j<=30)||(i==4&&j==63)||(i>=7&&i<=13&&j>=74&&j<=80)||(i==11&&j==65)||(i==15&&j==57)||(i==16&&j==82)||(i==17&&j==84)||(i==18&&j==82)||(i==19&&j==84)||(i>=27&&i<=29&&j>=6&&j<=10)||(i>=27&&i<=28&&j>=46&&j<=49)||(i==29&&j==95)||(i==30&&(j==9||j==95))||(i>=31&&i<=35&&j>=41&&j<=55)||(i==31&&j==94)||(i==34&&j>=17&&j<=29)||(i==36&&j>=48)||(i>=36&&i<=42&&j>=92&&j<=97)||(i==39&&j==49)||(i==41&&j==18)||(i>=41&&i<=50&&j>=62&&j<=65)||(i>=48&&i<=49&&j>=18&&j<=35)||(i==52&&j==4)||(i>=78&&i<=85&&j>=12&&j<=18)||(i==78&&j==76))));
		for(int n = 0; n<g_vptpn_oled[i][j][k]->GetN(); ++n){
			double x,y;
			g_vptpn_oled[i][j][k]->GetPoint(n,x,y);
			if(y==0) continue;
			if(x<1314920000) continue;//02/09 daystart
			if(x>1317690000) continue;//04/10 daystart
				//special treatment
			if(!(Ospecialcrystals) ){
				if(y>1.2 && normalized)  continue;//veto outliers
				if(y<0.8 && normalized) continue;
			}
			//now special treatment crystals
			//upperborder
			if(y>1.5 && normalized && (O15)) continue;
			if(y>1.4 && normalized && (O14)) continue;
			if(y>1.35&& normalized && (O135)) continue;
			if(y>1.3 && normalized && (O13)) continue;
			if(y>1.25&& normalized && (O125)) continue;
			if(y>1.2 && normalized && (O12)) continue;
			//lower border
			if(y<0.75&& normalized && (O75)) continue;
			if(y<0.7 && normalized && (O7)) continue;
			if(y<0.65&& normalized && (O65)) continue;
			if(y<0.6 && normalized && (O6)) continue;
			if(y<0.5 && normalized && (O5)) continue;
			//single crystals
			if(Onoisecrystals) continue;
			if(y<0.975&&normalized&&k==0&&i==14&&j==16) continue;
			if(y>1.125&&normalized&&k==0&&i==14&&j==18) continue;
			if(y<0.999&&x>1316478000&&x<1319070000&&normalized&&k==0&&i==15&&j==19) continue;
			if(y>1.025&&x>1314400000&&x<1315100000&&normalized&&k==0&&i==16&&j==67) continue;
			if(y>1.05 &&normalized&&k==0&&i==16&&j==69) continue;
			if(y>1.1  &&normalized&&k==0&&i==24&&j==22) continue;
			if(y>1.05 &&x<1315000000&&normalized&&k==0&&i==24&&j==56) continue;
			if(y<0.8  &&x<1316000000&&normalized&&k==0&&i==34&&j==48) continue;
			if((y<0.8||y>1.2)&&x<1316000000&&x>1315000000&&normalized&&k==0&&i==49&&(j==36||j==39)) continue;
			if(y>1.1  &&normalized&&k==0&&i==53&&j==3) continue;
			if(y<0.8  &&x<1315900000&&normalized&&k==0&&i==60&&j==63) continue;
			if(y<0.76 &&x<1316000000&&x>1315000000&&normalized&&k==0&&i==66&&j==49) continue;
			if(y<0.9  &&k==0&&((i==86&&j==36)||(i==91&&j==68))) continue;
			if(y>1.1  &&normalized&&k==1&&i==23&&j==35) continue;
			if(y>0.83 &&normalized&&k==1&&i==24&&j==35) continue;
			if(y<0.8  &&x<1316000000&&x>1315000000&&normalized&&k==1&&i==39&&j==64) continue;
			if(y<0.7  &&x<1316000000&&x>1315000000&&normalized&&k==1&&i==56&&j==39) continue;
			if(y>1.1  &&normalized&&k==1&&i==60&&j==99) continue;
			if(y<0.9  &&normalized&&k==1&&i==96&&j==43) continue;
			if(y<0.9  &&normalized&&k==1&&i==100&&j==50) continue;

			if(y>vptpnmaxorange[i][j][k]) vptpnmaxorange[i][j][k] = y;
			if(y<vptpnminorange[i][j][k]) vptpnminorange[i][j][k] = y;
		}
		if(plottingall_VPTPN && g_vptpn_oled[i][j][k]->GetN()>10 && vptpnmmustd[i][j][k]>0 ){//with mustd oled instead of las
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
		if(normalized) haxis->SetMinimum(0.565);//normalized, las 0.565, orange 0.55
		if(normalized) haxis->SetMaximum(1.165);//normalized, las 1.165, orange 1.30
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
		if(vptpnmmustd[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintNoEPS(col, outputfile, outputdirac, false);//oled instead of las
		if(vptpnmmustd[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintEPS(col, outputfile, outputdirac);//oled instead of las
	//	delete l1; delete l2;
		}

		col->Clear();
		if(vptpnmax[i][j][k]-vptpnmin[i][j][k] <= 0) continue;
		if(vptpnmax[i][j][k]==0) continue;
		if(vptpnmin[i][j][k]==99) continue;
		int nn = -1;
		if(vptpnruss[i][j][k]) nn = 0;
		else nn=1;
		double mustd_ = vptpnmmustd[i][j][k];
		double vptdiff = vptpnmax[i][j][k]-vptpnmin[i][j][k];
		double vptdifforange = vptpnmaxorange[i][j][k]-vptpnminorange[i][j][k];
		//if(!(normalized) && vptdiff > 0.3) cout << "vptdiff " << vptdiff << " mustd " << mustd_ << " eta " << eta << " nn " << nn << " i " << i << " j " << j << " k " << k << endl;
		/* //histograms
		if(eta<(-2.8) && mustd_!=0)           h_mustd_vs_VPTPNdiff[0][nn] ->Fill(mustd_, vptdiff);//m3p0
		else if(eta<(-2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[1][nn] ->Fill(mustd_, vptdiff);
		else if(eta<(-2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff[2][nn] ->Fill(mustd_, vptdiff);//m2p6
		else if(eta<(-2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff[3][nn] ->Fill(mustd_, vptdiff);
		else if(eta<(-2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[4][nn] ->Fill(mustd_, vptdiff);//m2p2
		else if(eta<(-1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[5][nn] ->Fill(mustd_, vptdiff);
		else if(eta<(-1.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[6][nn] ->Fill(mustd_, vptdiff);//m1p8
		else if(eta<=(-1.4)&& mustd_!=0)      h_mustd_vs_VPTPNdiff[7][nn] ->Fill(mustd_, vptdiff);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_mustd_vs_VPTPNdiff[8][nn] ->Fill(mustd_, vptdiff);//1p4
		else if(eta<( 1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[9][nn] ->Fill(mustd_, vptdiff);
		else if(eta<( 2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[10][nn]->Fill(mustd_, vptdiff);
		else if(eta<( 2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff[11][nn]->Fill(mustd_, vptdiff);
		else if(eta<( 2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff[12][nn]->Fill(mustd_, vptdiff);
		else if(eta<( 2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff[13][nn]->Fill(mustd_, vptdiff);
		else if(eta<( 2.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff[14][nn]->Fill(mustd_, vptdiff);
		else if(eta<=(3.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff[15][nn]->Fill(mustd_, vptdiff);
		//scaled
		if(eta<(-2.8) && mustd_!=0)           h_mustd_vs_VPTPNdiff_scaled[0][nn] ->Fill(mustd_, vptdiff/15.);
		else if(eta<(-2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[1][nn] ->Fill(mustd_, vptdiff/10.);
		else if(eta<(-2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[2][nn] ->Fill(mustd_, vptdiff/7.);
		else if(eta<(-2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[3][nn] ->Fill(mustd_, vptdiff/4.5);
		else if(eta<(-2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[4][nn] ->Fill(mustd_, vptdiff/2.2);
		else if(eta<(-1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[5][nn] ->Fill(mustd_, vptdiff/0.85);
		else if(eta<(-1.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[6][nn] ->Fill(mustd_, vptdiff/0.5);
		else if(eta<=(-1.4)&& mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[7][nn] ->Fill(mustd_, vptdiff/0.3);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_mustd_vs_VPTPNdiff_scaled[8][nn] ->Fill(mustd_, vptdiff/0.3);
		else if(eta<( 1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[9][nn] ->Fill(mustd_, vptdiff/0.5);
		else if(eta<( 2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[10][nn]->Fill(mustd_, vptdiff/0.85);
		else if(eta<( 2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[11][nn]->Fill(mustd_, vptdiff/2.2);
		else if(eta<( 2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[12][nn]->Fill(mustd_, vptdiff/4.5);
		else if(eta<( 2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[13][nn]->Fill(mustd_, vptdiff/7.);
		else if(eta<( 2.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[14][nn]->Fill(mustd_, vptdiff/10.);
		else if(eta<=(3.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[15][nn]->Fill(mustd_, vptdiff/15.);		//if(fabs(eta)<2.8) continue;
		//if(vptpnmmustd[i][j][k]!=0){ h_mustd_vs_VPTPNdiff->Fill(vptpnmmustd[i][j][k], vptpnmax[i][j][k]-vptpnmin[i][j][k]);
		//cout << "mustd " << vptpnmmustd[i][j][k] << " diff " << vptpnmax[i][j][k]-vptpnmin[i][j][k] << endl;
		//}
		*/

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

		if(vptdifforange>0 && vptpnmaxorange[i][j][k]!=0 && vptpnminorange[i][j][k]!=99){
		if(eta<(-2.8) && mustd_!=0)           h_orange_mustd_vs_VPTPNdiff[0][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[0][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-2.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[1][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[1][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-2.4) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[2][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[2][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-2.2) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[3][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[3][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-2.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[4][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[4][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-1.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[5][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[5][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<(-1.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[6][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[6][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<=(-1.4)&& mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[7][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[7][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_orange_mustd_vs_VPTPNdiff[8][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[8][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 1.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[9][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff[9][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 2.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[10][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[10][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 2.2) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[11][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[11][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 2.4) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[12][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[12][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 2.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[13][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[13][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<( 2.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[14][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[14][nn]->GetN(), mustd_,vptdifforange);
		else if(eta<=(3.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff[15][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff[15][nn]->GetN(), mustd_,vptdifforange);

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&mustd_!=0)a_orange_mustd_vs_VPTPNdiff[0][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[0][nn]->GetN(), mustd_,vptdifforange);//1p4
		else if(fabs(eta)<( 1.8) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[1][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[1][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<( 2.0) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[2][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[2][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<( 2.2) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[3][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[3][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<( 2.4) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[4][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[4][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<( 2.6) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[5][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[5][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<( 2.8) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[6][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[6][nn]->GetN(), mustd_,vptdifforange);
		else if(fabs(eta)<=(3.0) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff[7][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff[7][nn]->GetN(), mustd_,vptdifforange);
		}
		//scaled
		//scaled
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./20.)*log(1.*pow((100.)*exp(-13.5112 + 7.91386*x - 0.998649*x*x),2))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./20.)*log(1.*pow((100.)*exp(-4.11065 + 0.258478*x - 0.*x*x),2))", 0., 1.6);//tried forfactor 1/20
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./20.)*log(100.*exp(-13.5112 + 7.91386*x - 0.998649*x*x))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./20.)*log(100.*exp(-4.11065 + 0.258478*x - 0.*x*x))", 0., 1.6);//tried forfactor 1/20
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./1.)*sqrt(1.*exp(-13.5112 + 7.91386*x - 0.998649*x*x))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./1.)*sqrt(1.*exp(-4.11065 + 0.258478*x - 0.*x*x))", 0., 1.6);//tried forfactor 1/20
		//TF1 *dosefunction = new TF1("dosefunction", "0.072+0.00813887*pow(x,6)",0.,3.3);
		//TF1 *dosefunction = new TF1("dosefunction", "0.072+0.0.00833127*pow(x,5.97824)",0.,3.3);
		TF1 *scalefunc = new TF1("scalefunc", "exp(-13.5112 + 7.91386*x - 0.998649*x*x)", 1.3, 3.3);//original function --> NORMALIZED TO eta==2.5
		TF1 *scalefuncEB = new TF1("scalefuncEB", "exp(-4.11065 + 0.258478*x - 0.*x*x)", 0., 1.6);//original function --> NORMALIZED TO eta==2.5
		//TF1 *scalefunc = new TF1("scalefunc", "(6.0168/3.464398)* exp(-13.5112 + 7.91386*x - 0.998649*x*x)", 1.3, 3.3);//original function --> NORMALIZED TO dosefunction->Eval(3.)
		//TF1 *scalefuncEB = new TF1("scalefuncEB", "(6.0168/3.464398) *exp(-4.11065 + 0.258478*x - 0.*x*x)", 0., 1.6);//original function --> NORMALIZED TO dosefunction->Eval(3.)
		TF1 *dosefunction = new TF1("dosefunction", "0.072+0.00826968*exp(2.19256*x)",0.,3.3);//in Gy
		TF1 *ewriacfunctionRuss = new TF1("ewriacfunctionRuss", "0.208617*log(1+0.0281834*(x*100))",0.,3.3);// /100 as in rad and 1 rad = 0.01 Grey
		TF1 *ewriacfunctionChin = new TF1("ewriacfunctionChin", "0.231791*log(1+0.0512092*(x*100))",0.,3.3);
		double dose = dosefunction->Eval(fabs(eta));
		double doseAlt = scalefunc->Eval(fabs(eta));
		if(fabs(eta)<1.5) doseAlt = scalefuncEB->Eval(fabs(eta));
		double ewriac;
		//DO HERE if russ / chin
		if(nn==0) ewriac = ewriacfunctionRuss->Eval(doseAlt);
		else      ewriac = ewriacfunctionChin->Eval(doseAlt);
		double scale = ewriac;
		if(scale<=0) cout << "ERROR!!!!!!!! scale " << scale << endl;
		if(eta<(-2.8) && mustd_!=0)           h_mustd_vs_VPTPNdiff_scaled[0][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdiff/scale/*/15.*/);
		else if(eta<(-2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[1][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdiff/scale/*/10.*/);
		else if(eta<(-2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[2][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdiff/scale/*/7.*/);
		else if(eta<(-2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[3][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdiff/scale/*/4.5*/);
		else if(eta<(-2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[4][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdiff/scale/*/2.2*/);
		else if(eta<(-1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[5][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdiff/scale/*/0.85*/);
		else if(eta<(-1.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[6][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdiff/scale/*/0.5*/);
		else if(eta<=(-1.4)&& mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[7][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdiff/scale/*/0.3*/);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_mustd_vs_VPTPNdiff_scaled[8][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[8][nn]->GetN(), mustd_,vptdiff/scale/*/0.3*/);
		else if(eta<( 1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[9][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[9][nn]->GetN(), mustd_,vptdiff/scale/*/0.5*/);
		else if(eta<( 2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[10][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[10][nn]->GetN(), mustd_,vptdiff/scale/*/0.85*/);
		else if(eta<( 2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[11][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[11][nn]->GetN(), mustd_,vptdiff/scale/*/2.2*/);
		else if(eta<( 2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[12][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[12][nn]->GetN(), mustd_,vptdiff/scale/*/4.5*/);
		else if(eta<( 2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[13][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[13][nn]->GetN(), mustd_,vptdiff/scale/*/7.*/);
		else if(eta<( 2.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[14][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[14][nn]->GetN(), mustd_,vptdiff/scale/*/10.*/);
		else if(eta<=(3.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[15][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[15][nn]->GetN(), mustd_,vptdiff/scale/*/15.*/);

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&mustd_!=0)a_mustd_vs_VPTPNdiff_scaled[0][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdiff/scale);//1p4
		else if(fabs(eta)<( 1.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[1][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[2][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.2) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[3][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.4) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[4][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.6) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[5][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[6][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<=(3.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[7][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdiff/scale);

		if(vptdifforange>0 && vptpnmaxorange[i][j][k]!=0 && vptpnminorange[i][j][k]!=99){
		if(eta<(-2.8) && mustd_!=0)           h_orange_mustd_vs_VPTPNdiff_scaled[0][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdifforange/scale/*/15.*/);
		else if(eta<(-2.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[1][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdifforange/scale/*/10.*/);
		else if(eta<(-2.4) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[2][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdifforange/scale/*/7.*/);
		else if(eta<(-2.2) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[3][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdifforange/scale/*/4.5*/);
		else if(eta<(-2.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[4][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdifforange/scale/*/2.2*/);
		else if(eta<(-1.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[5][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdifforange/scale/*/0.85*/);
		else if(eta<(-1.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[6][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdifforange/scale/*/0.5*/);
		else if(eta<=(-1.4)&& mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[7][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdifforange/scale/*/0.3*/);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_orange_mustd_vs_VPTPNdiff_scaled[8][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[8][nn]->GetN(), mustd_,vptdifforange/scale/*/0.3*/);
		else if(eta<( 1.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[9][nn] ->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[9][nn]->GetN(), mustd_,vptdifforange/scale/*/0.5*/);
		else if(eta<( 2.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[10][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[10][nn]->GetN(), mustd_,vptdifforange/scale/*/0.85*/);
		else if(eta<( 2.2) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[11][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[11][nn]->GetN(), mustd_,vptdifforange/scale/*/2.2*/);
		else if(eta<( 2.4) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[12][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[12][nn]->GetN(), mustd_,vptdifforange/scale/*/4.5*/);
		else if(eta<( 2.6) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[13][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[13][nn]->GetN(), mustd_,vptdifforange/scale/*/7.*/);
		else if(eta<( 2.8) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[14][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[14][nn]->GetN(), mustd_,vptdifforange/scale/*/10.*/);
		else if(eta<=(3.0) && mustd_!=0)      h_orange_mustd_vs_VPTPNdiff_scaled[15][nn]->SetPoint(h_orange_mustd_vs_VPTPNdiff_scaled[15][nn]->GetN(), mustd_,vptdifforange/scale/*/15.*/);

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&mustd_!=0)a_orange_mustd_vs_VPTPNdiff_scaled[0][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdifforange/scale);//1p4
		else if(fabs(eta)<( 1.8) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[1][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<( 2.0) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[2][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<( 2.2) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[3][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<( 2.4) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[4][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<( 2.6) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[5][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<( 2.8) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[6][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdifforange/scale);
		else if(fabs(eta)<=(3.0) && mustd_!=0)      a_orange_mustd_vs_VPTPNdiff_scaled[7][nn]->SetPoint(a_orange_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdifforange/scale);
		}

		if(vptdiff>0.12&&vptdiff<999.&&mustd_>0.45&&mustd_<0.5&&nn==0&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << endl;
		if(vptdiff/scale>40.&&vptdiff/scale<999.&&mustd_>0.75&&mustd_<0.8&&nn==1&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;
		if(vptdiff>35.&&vptdiff<999.&&mustd_>0.95&&mustd_<1.05&&nn==0&&eta<=1.6&&eta>1.4) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << endl;
		if(vptdiff/scale>4.&&vptdiff/scale<999.&&mustd_>0.45&&mustd_<0.5&&nn==0&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;
		if(vptdiff/scale>0.2&&vptdiff/scale<999.&&mustd_>1.&&mustd_<1.3&&nn==1&&eta<=-2.6&&eta>-2.8) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;

	}}}

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
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		//if(n==1&&(nn==4||nn==13)){
		//haxis->GetXaxis()->SetRangeUser(0.,2.);
		//haxis->GetXaxis()->SetLimits(0.,2.);
		//haxis->SetAxisRange(0.,2.);
		//h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.);
		//h_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetLimits(0.,2.);
		//cout << col->GetTitle() << " has xaxis boundaries of " << haxis->GetXaxis()->GetXmin() << " to " << haxis->GetXaxis()->GetXmax() << endl;
		//}
		//else{
		//cout << col->GetTitle() << " has xaxis boundaries of " << haxis->GetXaxis()->GetXmin() << " to " << haxis->GetXaxis()->GetXmax() << endl;
		//}
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
//		else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->Draw("hist");
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		TString outputfile = h_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
		//col->Update();
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
		col->SetTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
	//	if(n<=4 && n>=11) h_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMaximum(0.05);//maybe uncomment this
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
		haxis = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		//cout << "n " << n << " nn " << nn << " xtitle " << h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->GetTitle() << endl;
	//	if(n==1&&(nn==4||nn==13)){
	//	haxis->GetXaxis()->SetRangeUser(0.,2.);
	//	haxis->GetXaxis()->SetLimits(0.,2.);
	//	haxis->SetAxisRange(0.,2.);
	//	h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.);
	//	h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetLimits(0.,2.);
	//	}
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,2.5);
//		else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->Draw("hist");
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		outputfile = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
		//col->Update();
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();


		col->cd();
		col->SetName(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetName() );
		col->SetTitle(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
		h_orange_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
		haxis = h_orange_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetName(), h_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.55);
		//else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->Draw("hist");
		h_orange_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		outputfile = h_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
		//col->Update();
		if(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		col->cd();
		col->SetName(h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
		col->SetTitle(h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
		//if(n<=4 && n>=11) h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMaximum(0.05);//maybe uncomment this
		h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
		haxis = h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
		haxis->SetNameTitle(h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
		haxis->GetXaxis()->SetRangeUser(0.,2.);
		haxis->GetXaxis()->SetLimits(0.,2.);
		haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		haxis->SetXTitle("#mu_{std} (m^{-1})");
		haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,5.5);
		//else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
		haxis->Draw("hist");
		h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		outputfile = h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
		//col->Update();
		if(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
		if(h_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
		col->Clear();

		if(n<8){
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
			//if(n==1&&(nn==4||nn==13)){
			//haxis->GetXaxis()->SetRangeUser(0.,2.);
			//haxis->GetXaxis()->SetLimits(0.,2.);
			//haxis->SetAxisRange(0.,2.);
			//a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.);
			//a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetLimits(0.,2.);
			//cout << col->GetTitle() << " has xaxis boundaries of " << haxis->GetXaxis()->GetXmin() << " to " << haxis->GetXaxis()->GetXmax() << endl;
			//}
			//else{
			//cout << col->GetTitle() << " has xaxis boundaries of " << haxis->GetXaxis()->GetXmin() << " to " << haxis->GetXaxis()->GetXmax() << endl;
			//}
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
	//		else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
			haxis->Draw("hist");
			a_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
			outputfile = a_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
			//col->Update();
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
	
			col->cd();
			col->SetName(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
			col->SetTitle(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
		//	if(n<=4 && n>=11) a_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMaximum(0.05);//maybe uncomment this
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
			haxis = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.,2.);
			haxis->GetXaxis()->SetLimits(0.,2.);
			haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
			haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("#mu_{std} (m^{-1})");
			haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
		//	if(n==1&&(nn==4||nn==13)){
		//	haxis->GetXaxis()->SetRangeUser(0.,2.);
		//	haxis->GetXaxis()->SetLimits(0.,2.);
		//	haxis->SetAxisRange(0.,2.);
		//	a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetRangeUser(0.,2.);
		//	a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetXaxis()->SetLimits(0.,2.);
		//	}
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,2.5);
	//		else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
			haxis->Draw("hist");
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
			outputfile = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
			//col->Update();
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
	
	
			col->cd();
			col->SetName(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetName() );
			col->SetTitle(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle() );
			a_orange_mustd_vs_VPTPNdiff[n][nn]->Draw("AP");
			haxis = a_orange_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetName(), a_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.,2.);
			haxis->GetXaxis()->SetLimits(0.,2.);
			haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
			haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("#mu_{std} (m^{-1})");
			haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.55);
			//else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
			haxis->Draw("hist");
			a_orange_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
			outputfile = a_orange_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
			//col->Update();
			if(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
			if(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
	
			col->cd();
			col->SetName(a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName() );
			col->SetTitle(a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle() );
			//if(n<=4 && n>=11) a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->SetMaximum(0.05);//maybe uncomment this
			a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("AP");
			haxis = a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetName(), a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.,2.);
			haxis->GetXaxis()->SetLimits(0.,2.);
			haxis->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
			haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("#mu_{std} (m^{-1})");
			haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,5.5);
			//else           haxis->GetYaxis()->SetRangeUser(0.,0.22);
			haxis->Draw("hist");
			a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
			outputfile = a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
			//col->Update();
			if(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir, false);
			if(a_orange_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
		}//if(n<8)

	}
	}
	cout << "saving ..." << endl;
	TFile *newFile;
//	if(normalized) newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_muSIC_vs_VPTPNmaxMinusVPTPNmin_savecopy.root", "RECREATE");//0312->0410; //muSIC <-> mustd
	if(normalized) newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_muSIC_vs_VPTPNmaxMinusVPTPNmin_savecopy_TESTEWRIAC_withSasha_unscaled.root", "RECREATE");//0312->0410; //muSIC <-> mustd

	else newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120410_mustd_vs_VPTPNmaxMinusVPTPNmin_nonorm.root", "RECREATE");//0312->0410
	newFile->cd();
	for(int n = 0; n<16; ++n){ for(int nn= 0; nn<2;++nn) {
	h_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	h_orange_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	h_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	if(n<8){
	a_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	a_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	a_orange_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	a_orange_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	}
	}}
	cout << "saved" << endl;

}
