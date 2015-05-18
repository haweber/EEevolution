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
float avg(vector<float> x){

  float sum = 0;
  for (Int_t i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}
*/

// THIS FUNCTION HAS TO BE CALLED IN THE DIRECTORY   /shome/haweber/ECAL/Laser/Analysis   SO THAT IT WORKS!!!!!!!!!!!!!!
//this function a) plots VPT/PN vs. time for all crystals that have mu_std information
//              b) computes VPT/PN(max)-VPT/PN(min) of a given time period and plots this vs. mu_std --> done in eta bins
//              c) this can be done with different normalization (false means normalizing to beginning of 2011 data running, 
//                 otherwise often choose beginning of time period of study
void EEstudies_VectorNTuples_continued(){

	bool plottingall_VPTPN = false;
	bool VPTdiff_vs_mustd  = false;//plot VPTdiff vs mustd --> you probably want this to be false as later you plot this including the fit
	bool plotonlyAbsEta    = true; //if true don't print EE+ EE- separately
	bool normalized        = true;
	bool startdiff         = false;//only true if normalized == false
	int period             = 0;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11
	bool burnin            = false;

	gROOT->SetStyle("Plain");

	//NOTE: Also change possibly newFileName at the very bottom!!!!!
	//new scaling in 20121001, old scaling in 20120829
	char gname[101];
	TString outputdir = "Plots/20121023/EEstudies/continued/notnormalized";//normalized";
	if(startdiff) outputdir = "Plots/20121023/EEstudies/continued/diffwrt2011begin";
	if(normalized) 
	        outputdir = "Plots/20121023/EEstudies/continued/normalized";
	TString outputdirac = "Plots/20121023/EEstudies/continued/allcrystalswithmustd/notnormalized";
	if(normalized) 
	        outputdirac = "Plots/20121023/EEstudies/continued/allcrystalswithmustd/normalized";
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
	if(burnin) outputdir = outputdir + "/burnin" + periodTS;
	else outputdir   = outputdir   + periodTS;
	outputdirac = outputdirac + periodTS;
	Util::MakeOutputDir(outputdir);
	if(plottingall_VPTPN)	Util::MakeOutputDir(outputdirac);

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

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

	TGraph* h_mustd_vs_VPTPNdiff[16][2];
	TGraph* h_mustd_vs_VPTPNdiff_scaled[16][2];
	TGraph* a_mustd_vs_VPTPNdiff[8][2];
	TGraph* a_mustd_vs_VPTPNdiff_scaled[8][2];
	string string_eta[16] = {"m3p0", "m2p8", "m2p6", "m2p4", "m2p2", "m2p0", "m1p8", "m1p6", "1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	for(int n = 0; n<16; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer

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
	}//if(n<8)

	}//for(int nn = 0; nn<2; ++ nn) // producer
	}//for(int n = 0; n<16; ++n)    // eta bin

//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120829_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121023_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121116_EE_VPT_over_PN_all_VectorNTuples_2011and2012_pointperday.root");//xxyyzz
	TGraph* g_vptpn_las[101][101][2];
//	TGraphErrors* g_vptpn_las[101][101][2];//xxyyzz
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
//	bool badsth[101][101][2];//xxyyzz
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
	//	badsth[i][j][k] = false;//xxyyzz
	}}}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	cout << "start computation" << endl;
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEm", i, j);//remove t
		else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEp", i, j);//remove t
		g_vptpn_las[i][j][k] = (TGraph*)vpt_values->Get(gname);
	//	g_vptpn_las[i][j][k] = (TGraphErrors*)vpt_values->Get(gname);//xxyyzz

		if(g_vptpn_las[i][j][k]->GetN()>0){ // start normalization
			double norm = 0; int countnorm = 0;
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){
				double x,y, xprev,yprev, xpost,ypost, xprev2,yprev2;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(hasbadperiod[i][j][k]==false){
				//	bool testcheck = false;
					//if(x<1202500000) continue;//before 2011 running
					//if(x>1320120000 && x<1333700000) continue;//between 2011 and 2012
					//if(x<1302590000) continue;//veto before period 2
					if(x<1320120000) continue;//veto 2012 (>) or 2011(<)
					//if(x>1340000000) continue;
					if(y==lastentry[i][j][k] && y > 0) lastsameentry[i][j][k] = lastsameentry[i][j][k] + 1;
					else if(y!=lastentry[i][j][k] && y>0) lastsameentry[i][j][k] = 0;
					lastentry[i][j][k] = y;
					if(lastsameentry[i][j][k]>=1000) {hasbadperiod[i][j][k] = true;
					//cout << "i:j:k " << i<<":"<<j<<":"<<k <<" y " << y << " time " << int(x) << endl;
					}
					if(period==6 && x>1333530000 && x<1334240000 && lastsameentry[i][j][k]>=100) hasbadperiod[i][j][k] = true;
				}
			//	if(y>1000e6) { g_vptpn_las[i][j][k]->RemovePoint(n); continue;}
				//1st period
				if(period==1 && !plottingall_VPTPN && x<1300200000) continue;//15/03/2011 	
				if(period==1 && !plottingall_VPTPN && x>1301210000) continue;//27/03/2011	
				//2nd period
				if(period==2 && !plottingall_VPTPN && x<1302590000) continue;//12/04/2011	
				if(period==2 && !plottingall_VPTPN && x>1304620000) continue;//05/05/2011	
				//3rd period
				if(period==3 && !plottingall_VPTPN && x<1305260000) continue;//13/05/2011	
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
				if(period==32&& !plottingall_VPTPN && x<1305260000) continue;//13/05/2011	
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
				//if(x>1315398150) continue;
				if((yprev>y*1.1 && ypost>y*1.1) || (yprev<y*0.9 && ypost<y*0.9)) g_vptpn_las[i][j][k]->SetPoint(n,x,yprev);
				//upper 2cleaning
				if(ypost<y*1.1 && yprev2<y*1.1 && ypost<yprev*1.1 && yprev2<yprev*1.1){
					g_vptpn_las[i][j][k]->SetPoint(n,x,yprev2);//clean only n
					if(n>0) g_vptpn_las[i][j][k]->SetPoint(n-1,xprev,yprev2);//clean only n
				}
			}//cleaning
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
					if(period==3 && !plottingall_VPTPN && x<1305260000) continue;//13/05/2011	
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
					if(period==32&& !plottingall_VPTPN && x<1305260000) continue;//13/05/2011	
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
					//1st period
					//if(period==1 && x>1300450000) continue;
					//2nd period
					//if(period==2 && x>1302840000) continue;
					//3rd period
					//if(period==3 && x>1305510000) continue;
					//4th period
					//if(period==4 && x>1310840000) continue;
					//5th period
					//if(period==5 && x>1315370000) continue;
					//original
					//if(period==0 && x>1315398150) continue;// 07/09/2011
					//as we start from 0 entry take the first entries of the period is just fine

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
					//cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" and t=" << x << " had VPT/PN = " << y << " which is now normalized to ";
					g_vptpn_las[i][j][k]->SetPoint(n,x,y/(vptpnnorm[i][j][k]));//normalization
					//cout << y << endl;
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
					//cout << "at i:j:k=" << i << ":"<<j<<":"<<k<<" and t=" << x << " had VPT/PN = " << y << " which is now normalized to ";
					g_vptpn_las[i][j][k]->SetPoint(n,x,y/(vptpnnorm[i][j][k]));//normalization
					//cout << y << endl;
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
		//get mu_std and eta for each crystals
		if(k==0){//choose xx_EEm_yy histos
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEm_russ->FindBin(i,j);
			mustd_russ = mu_ECAL_EEm_russ->GetBinContent(bin);
			mustd_chin = mu_ECAL_EEm_chin->GetBinContent(bin);

			if(mustd_russ==0 && mustd_chin==0) continue;
			else if(mustd_russ!=0){
				russchin << mustd_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_russ;
				vptpnruss[i][j][k] = true;//is russian
			}
			else{
				russchin << mustd_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_chin;
				vptpnruss[i][j][k] = false;//is not russian thus it is chinese
			}
		}
		else{//choose xx_EEp_yy histos
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = mu_ECAL_EEp_russ->FindBin(i,j);
			mustd_russ = mu_ECAL_EEp_russ->GetBinContent(bin);
			mustd_chin = mu_ECAL_EEp_chin->GetBinContent(bin);
			if(mustd_russ==0 && mustd_chin==0) continue;
			
			else if(mustd_russ!=0){
				if(mustd_russ>=1.-10e-6&& mustd_russ<=1.+10e-6) continue;
				if(mustd_russ>=2.) continue;
				russchin << mustd_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_russ;
				vptpnruss[i][j][k] = true;
			}
			else{
				if(mustd_chin>=1.-10e-6&& mustd_chin<=1.+10e-6) continue;
				if(mustd_chin>=2.) continue;
				russchin << mustd_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__mustd_" + rs;
				vptpnmmustd[i][j][k] = mustd_chin;
				vptpnruss[i][j][k] = false;
			}
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


//		bool noisecrystals = ((k==1&&((i==78&&j==76)||(i==52&&j==4)||(i==34&&(j==29||j==17))||(i<=30&&i>=26&&(j<=10&&j>=6))||(i==25&&(j==92))||(i<=24&&i>=21&&j<=55&&j>=49)||(i==19&&(j==84||j==83))||(i==4&&j==63)||(i==4&&j==63)||(i==16&&j==63)||(i==24&&j==65)||(i==33&&j==69)||(i==35&&j==72)||(i==38&&j==19)||(i==38&&j==96)||(i==42&&j==75)||(i==49&&j==6)||(i==52&&j==27)||(i==66&&j==52)||(i==80&&j==22)||(i==96&&j==57)||(i==100&&j==59)))||(k==0&&((i<=97&&i>=92&&j<=42&&j>=29)||(i==86&&j==36)||(i<=55&&i>=51&&j<=35&&j>=31)||(i==52&&j==11)||(i==41&&j==60)||(i==46&&j==67)||(i==54&&j==94)||(i==55&&j==36)||(i==64&&j==41)||(i==71&&j==67)||(i==76&&j==61)||(i==82&&j==36)||(i==92&&j==27))));
			//kill crazy crystals for Vector NTuples
			if(k==0&&((i==9&&j==33)||(i==47&&j==78))) continue;
			if(k==1&&( (i==39&&j==46)||(i==52&&j==4)||(i==78&&j==76)||(i==83&&j==21)||(i==84&&j==24)||(i==86&&j==24)||(i==94&&j==64)||(i==19&&j==84)||(i==24&&j==77) ) ) continue;
			//no values
			if(k==1&&(i==91||i==92)&&(j>=21&&j<=25) ) continue;
			if(bad2012) continue;
			//kill crystals where we have constant period
		//	if(hasbadperiod[i][j][k] && !(risingcrystals || partiallyrising || partiallyrising2)) {
		//	cout << "crystal ix:iy:cap " << i << ":" << j << ":" << k << " has period of constant VPT/PN --> VETO" << endl;
		//	continue; }
			if(hasbadperiod[i][j][k]) continue;
/*
			double yprev = -1.; int sameyprev = 0;//xxyyzz start
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
		if(badsth[i][j][k]) continue; //xxyyzz end
*/
		for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){

			double x,y;
			g_vptpn_las[i][j][k]->GetPoint(n,x,y);
			if(y==0) continue;

			//1st period
			if(period==1  && x<1300200000) continue;//15/03/2011 	
			if(period==1  && x>1301210000) continue;//27/03/2011	
			//2nd period
			if(period==2  && x<1302590000) continue;//12/04/2011	
			if(period==2  && x>1304620000) continue;//05/05/2011	
			//3rd period
			if(period==3  && x<1305260000) continue;//13/05/2011	
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
			if(period==32 && x<1305260000) continue;//13/05/2011	
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

//			if(y>1.1 && normalized) continue;//veto outliers	 --> check this
//			if(y<0.565 && normalized) continue;			 --> check this
			//kill crazy crystals
		//wide seems two have to levels parallel
//			if(i==17&&j==16&&k==1&&y<0.9) continue;
//			if(noisecrystals) continue;
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
				while(num1<5){
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
				while(num1<5){
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
				while(num1<5){
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
				while(num1<5){
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
				while(num1<5){
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
				while(num1<5){
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
		if(plottingall_VPTPN && g_vptpn_las[i][j][k]->GetN()>10 && vptpnmmustd[i][j][k]>0 ){//with mustd oled instead of las
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
		else           haxis->SetMinimum(0.25);
		if(normalized) haxis->SetMaximum(1.165);//normalized, las 1.165, orange 1.30
		else           haxis->SetMaximum(1.5);
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
		if(vptpnmmustd[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintEPS(col, outputfile, outputdirac);//oled instead of las
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

		//check for crazy crystals
		if(nn==0&&mustd_>0.9&&mustd_<1.1&&fabs(eta)>=1.6&&fabs(eta)<=1.8&&vptdiff>0.2)cout<<"i:j:k="<<i<<":"<<j<<":"<<k<< " eta " << eta << " mustd " << mustd_ << " russian=0? " << nn << " vptdiff " << vptdiff << " (" << vptpnmax[i][j][k] << "-" << vptpnmin[i][j][k] << ")" << endl;

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

		//scaled
		//all scale functions
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./20.)*log(1.*pow((100.)*exp(-13.5112 + 7.91386*x - 0.998649*x*x),2))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./20.)*log(1.*pow((100.)*exp(-4.11065 + 0.258478*x - 0.*x*x),2))", 0., 1.6);//tried forfactor 1/20
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./20.)*log(100.*exp(-13.5112 + 7.91386*x - 0.998649*x*x))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./20.)*log(100.*exp(-4.11065 + 0.258478*x - 0.*x*x))", 0., 1.6);//tried forfactor 1/20
	//	TF1 *scalefunc = new TF1("scalefunc", "(1./1.)*sqrt(1.*exp(-13.5112 + 7.91386*x - 0.998649*x*x))", 1.3, 3.3);//
	//	TF1 *scalefuncEB = new TF1("scalefuncEB", "(1./1.)*sqrt(1.*exp(-4.11065 + 0.258478*x - 0.*x*x))", 0., 1.6);//tried forfactor 1/20
		//TF1 *dosefunction = new TF1("dosefunction", "0.072+0.148382*exp(0.410485*x)",0.,3.3);// ECAL TDR Fig. 2.7
		//TF1 *dosefunction = new TF1("dosefunction", "0.072+0.00813887*pow(x,6)",0.,3.3);
		//TF1 *dosefunction = new TF1("dosefunction", "0.072+0.0.00833127*pow(x,5.97824)",0.,3.3);
		//TF1 *scalefunc = new TF1("scalefunc", "(6.0168/3.464398)* exp(-13.5112 + 7.91386*x - 0.998649*x*x)", 1.3, 3.3);//original function --> NORMALIZED TO dosefunction->Eval(3.) --> this is Sasha's function normalized to TDR dose function fit of ECAL plot
		//TF1 *scalefuncEB = new TF1("scalefuncEB", "(6.0168/3.464398) *exp(-4.11065 + 0.258478*x - 0.*x*x)", 0., 1.6);//original function --> NORMALIZED TO dosefunction->Eval(3.)
		TF1 *scalefunc = new TF1("scalefunc", "exp(-13.5112 + 7.91386*x - 0.998649*x*x)", 1.3, 3.3);//original function --> NORMALIZED TO eta==2.5 --> this is Sasha's function from the indico slides
		TF1 *scalefuncEB = new TF1("scalefuncEB", "exp(-4.11065 + 0.258478*x - 0.*x*x)", 0., 1.6);//original function --> NORMALIZED TO eta==2.5
	//x	TF1 *scalefunc = new TF1("scalefunc", "(6.0168/3.464398)* exp(-13.5112 + 7.91386*x - 0.998649*x*x)", 1.3, 3.3);//original function --> NORMALIZED TO dosefunction->Eval(3.)
	//x	TF1 *scalefuncEB = new TF1("scalefuncEB", "(6.0168/3.464398) *exp(-4.11065 + 0.258478*x - 0.*x*x)", 0., 1.6);//original function --> NORMALIZED TO dosefunction->Eval(3.)
		//factor 4.329756 = ECAL TDR Fig.2.7 / sasha's function at eta = 3.0 (now it is Gy/h
		//factor 0.8 = ECAL dose rate at 7 TeV / ECAL dose rate at 14 TeV (as in ECAL TDR Fig.2.7) due to xsection difference.
		//factor luminosity: 1 Hz nb = 10^33 cm-2 s-1 = 1/10 L_design = 1/10 10^34 L_design (--> factor 10 as basis)
		// divide by 1/e to get a better estimate of <L_inst> instead of <L_inst^peak>
		//period 1: <L_inst^peak> ~ 0.015 Hz nb (increasing) --> <L_inst> ~ 0.0055 Hz nb
		//period 2: <L_inst^peak> ~ 0.575 Hz nb (strong increasing) --> <L_inst> ~ 0.2115 Hz nb
		//period 3: <L_inst^peak> ~ 1.2 Hz nb (stable) --> <L_inst> ~ 0.4415 Hz nb
		//period 4: <L_inst^peak> ~ 2.025 Hz nb (slight increasing) --> <L_inst> ~ 0.7496 Hz nb
		//period 5: <L_inst^peak> ~ 3.35 Hz nb (stable) --> <L_inst> ~ 1.2324 Hz nb
		//period 0: <L_inst^peak> ~ 3.3 Hz nb --> <L_inst> ~ 1.2140 Hz nb
	//	TF1 *dosefunction = new TF1("dosefunction", "0.072+0.00826968*exp(2.19256*x)",0.,3.3);//in Gy --> ECAL TDR Fig. A.11
		TF1 *ewriacfunctionRuss = new TF1("ewriacfunctionRuss", "0.208617*log(1+0.0281834*(x*100))",0.,3.3);// /100 as in rad and 1 rad = 0.01 Grey
		TF1 *ewriacfunctionChin = new TF1("ewriacfunctionChin", "0.231791*log(1+0.0512092*(x*100))",0.,3.3);// this function is in rad/h --> *100
//		double dose = dosefunction->Eval(fabs(eta));
		double doseAlt = scalefunc->Eval(fabs(eta));
		if(fabs(eta)<1.5) doseAlt = scalefuncEB->Eval(fabs(eta));
		//factor = 4.329756 * 0.8 *(1/10) * <L_inst>
		double factor = 4.329756 * 0.8 * 0.1;
		if(     period==1) factor = factor * 0.0055;
		else if(period==2) factor = factor * 0.2115;
		else if(period==3) factor = factor * 0.4415;
		else if(period==4) factor = factor * 0.7496;
		else if(period==5) factor = factor * 1.2324;
		else               factor = factor * 1.2140;//period 0
		doseAlt = factor*doseAlt;
		//for period 0 also use //x from above 
		double ewriac;
		//DO HERE if russ / chin
		if(nn==0) ewriac = ewriacfunctionRuss->Eval(doseAlt);
		else      ewriac = ewriacfunctionChin->Eval(doseAlt);
		double scale = ewriac;
		if(scale<=0) cout << "ERROR!!!!!!!! scale " << scale << endl;
		if(eta<(-2.8) && mustd_!=0)           h_mustd_vs_VPTPNdiff_scaled[0][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[1][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[2][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[3][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[4][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[5][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<(-1.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[6][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<=(-1.4)&& mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[7][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<1.6&&eta>=1.4&&mustd_!=0) h_mustd_vs_VPTPNdiff_scaled[8][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[8][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 1.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[9][nn] ->SetPoint(h_mustd_vs_VPTPNdiff_scaled[9][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 2.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[10][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[10][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 2.2) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[11][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[11][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 2.4) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[12][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[12][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 2.6) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[13][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[13][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<( 2.8) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[14][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[14][nn]->GetN(), mustd_,vptdiff/scale);
		else if(eta<=(3.0) && mustd_!=0)      h_mustd_vs_VPTPNdiff_scaled[15][nn]->SetPoint(h_mustd_vs_VPTPNdiff_scaled[15][nn]->GetN(), mustd_,vptdiff/scale);

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&mustd_!=0)a_mustd_vs_VPTPNdiff_scaled[0][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[0][nn]->GetN(), mustd_,vptdiff/scale);//1p4
		else if(fabs(eta)<( 1.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[1][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[1][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[2][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[2][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.2) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[3][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[3][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.4) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[4][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[4][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.6) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[5][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[5][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<( 2.8) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[6][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[6][nn]->GetN(), mustd_,vptdiff/scale);
		else if(fabs(eta)<=(3.0) && mustd_!=0)      a_mustd_vs_VPTPNdiff_scaled[7][nn]->SetPoint(a_mustd_vs_VPTPNdiff_scaled[7][nn]->GetN(), mustd_,vptdiff/scale);
		//check some crazy crystals
//		if(vptdiff>0.12&&vptdiff<999.&&mustd_>0.45&&mustd_<0.5&&nn==0&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << endl;
//		if(vptdiff/scale>40.&&vptdiff/scale<999.&&mustd_>0.75&&mustd_<0.8&&nn==1&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;
//		if(vptdiff>35.&&vptdiff<999.&&mustd_>0.95&&mustd_<1.05&&nn==0&&eta<=1.6&&eta>1.4) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << endl;
//		if(vptdiff/scale>4.&&vptdiff/scale<999.&&mustd_>0.45&&mustd_<0.5&&nn==0&&eta<=-1.4&&eta>-1.6) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;
//		if(vptdiff/scale>0.2&&vptdiff/scale<999.&&mustd_>1.&&mustd_<1.3&&nn==1&&eta<=-2.6&&eta>-2.8) cout << "i:j:k="<<i<<":"<<j<<":"<<k<<" producer " << nn << " eta " << eta << " mustd " << mustd_ << " vptdiff " << vptdiff << " scale " << scale << endl;

	}}}

	//draw correlation plots VPT/PN vs mu_std in eta bins
	for(int n = 0; n<16; ++n){//eta bins
	for(int nn= 0; nn<2;++nn){//producer: 0: russian, 1: chinese
		TH1F *haxis = h_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
		TString outputfile;
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
		else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
		haxis->Draw("hist");
		h_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
		outputfile = h_mustd_vs_VPTPNdiff[n][nn]->GetTitle();
		//col->Update();
		if(h_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
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
		if(normalized) haxis->GetYaxis()->SetRangeUser(0.,3.5);
		else           haxis->GetYaxis()->SetRangeUser(0.,4.);
		haxis->Draw("hist");
		h_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
		outputfile = h_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
		//col->Update();
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
			//col->Update();
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
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
			if(normalized) haxis->GetYaxis()->SetRangeUser(0.,3.5);
			else           haxis->GetYaxis()->SetRangeUser(0.,4.);
			haxis->Draw("hist");
			a_mustd_vs_VPTPNdiff_scaled[n][nn]->Draw("P");
			outputfile = a_mustd_vs_VPTPNdiff_scaled[n][nn]->GetTitle();
			//col->Update();
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
		     }
		}//if(n<8)

	}
	}
	//store histos into file
	TFile *newFile;
	if(startdiff) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_diffwrt2011begin_VectorNTuples"+periods+".root", "RECREATE");
	if(normalized) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples"+periods+".root", "RECREATE");//0312->0410; //muSIC <-> mustd
	//if(normalized) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples"+periods+"temp.root", "RECREATE");//0312->0410; //muSIC <-> mustd //xxyyzz
	else if(!startdiff)  newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_nonormalization_VectorNTuples"+periods+".root", "RECREATE");//0312->0410
	newFile->cd();
	for(int n = 0; n<16; ++n){ for(int nn= 0; nn<2;++nn) {
	h_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	h_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	if(n<8){
	a_mustd_vs_VPTPNdiff[n][nn]->Write(); 
	a_mustd_vs_VPTPNdiff_scaled[n][nn]->Write(); 
	}
	}}
	cout << "histograms for period " << period << " called " << periods << " (normalized/startdiff = " << normalized << "/" << startdiff << ") stored in " << newFile->GetName() << endl;

}
