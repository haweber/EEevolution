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
void EEstudies_VectorNTuples_burnin(){

	bool plottingall_VPTPN = false;
	bool VPTdiff_vs_burnin = false;//plot VPTdiff vs burnin --> you probably want this to be false as later you plot this including the fit
	bool normalized        = true;
	bool startdiff         = false;//only true if normalized == false
	int period             = -1;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11

	gROOT->SetStyle("Plain");
	gStyle->SetOptStat(0);

	//NOTE: Also change possibly newFileName at the very bottom!!!!!
	//new scaling in 20121001, old scaling in 20120829
	char gname[101];
	TString outputdir = "Plots/20121120/EEstudies/continued/notnormalized";//normalized";
	if(startdiff) outputdir = "Plots/20121120/EEstudies/continued/diffwrt2011begin";
	if(normalized) 
	        outputdir = "Plots/20121120/EEstudies/continued/normalized";
	TString outputdirac = "Plots/20121120/EEstudies/continued/allcrystalswithburnin/notnormalized";
	if(normalized) 
	        outputdirac = "Plots/20121120/EEstudies/continued/allcrystalswithburnin/normalized";
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
	outputdir = outputdir + "/burnin" + periodTS;
	outputdirac = outputdirac + periodTS;
	outputdir = "Plots/20130402/EEstudies/tempburnin";//xxyyzz
	if(startdiff) outputdir = outputdir + "/diffwrt2011begin";
	if(normalized) outputdir = outputdir + "/normalized";
	Util::MakeOutputDir(outputdir);
	if(plottingall_VPTPN)	Util::MakeOutputDir(outputdirac);

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	TFile *stdmapsfile     = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *EEp_producer     = (TH2D*)stdmapsfile->Get("EEp_producer");//1 russian, 2 chinese
	TH2D *EEm_producer     = (TH2D*)stdmapsfile->Get("EEm_producer");

	TFile *burninfile      = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEburnin.root");
	TH2D *burnin_EEm       = (TH2D*)burninfile->Get("burnin_EEm");
	TH2D *burnin_EEp       = (TH2D*)burninfile->Get("burnin_EEp");

	TGraph* a_burnin_vs_VPTPNdiff[8][2];
	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		string hname;
		hname = (string)"burnin__vs__VPTPNmax-VPTPNmin__AbsEta_" + string_feta[n] + string_prod[nn];
		a_burnin_vs_VPTPNdiff[n][nn] = new TGraph();
		a_burnin_vs_VPTPNdiff[n][nn]->SetNameTitle(hname.c_str(), hname.c_str());
		//a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_burnin_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in ratio");
	
	}//for(int nn = 0; nn<2; ++ nn) // producer
	}//for(int n = 0; n<8; ++n)    // eta bin

//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120829_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
//	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121023_EE_VPT_over_PN_all_VectorNTuples_2011and2012.root");
	TFile *vpt_values = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121116_EE_VPT_over_PN_all_VectorNTuples_2011and2012_pointperday.root");//xxyyzz
//	TGraph* g_vptpn_las[101][101][2];
	TGraphErrors* g_vptpn_las[101][101][2];//xxyyzz
	double burnin_russ = 0; double burnin_chin = 0;

	//define and reset all values
	double vptpnnorm[101][101][2];
	double vptpnmax[101][101][2];
	double vptpnmin[101][101][2];
	double timemax[101][101][2];
	double timemin[101][101][2];
	int    vptpnmaxindex[101][101][2];
	int    vptpnminindex[101][101][2];
	double vptpnmburnin[101][101][2];
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
		vptpnmburnin[i][j][k] = 0; 
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
		if (k==0) sprintf(gname,"tg_vptpn_las_ix%d_iy%d_EEm", i, j);//remove t
		else sprintf(gname,"tg_vptpn_las_ix%d_iy%d_EEp", i, j);//remove t
	//	g_vptpn_las[i][j][k] = (TGraph*)vpt_values->Get(gname);
		g_vptpn_las[i][j][k] = (TGraphErrors*)vpt_values->Get(gname);//xxyyzz

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
			bin = EEm_producer->FindBin(i,j);
			int prod = EEm_producer->GetBinContent(bin);
			if(prod==1){ burnin_russ = burnin_EEm->GetBinContent(bin);
			             burnin_chin = 0; }
			if(prod==2){ burnin_chin = burnin_EEm->GetBinContent(bin);
			             burnin_russ = 0; }

			if(burnin_russ==0 && burnin_chin==0) continue;
			else if(burnin_russ!=0){
				russchin << burnin_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__burnin_" + rs;
				vptpnmburnin[i][j][k] = burnin_russ;
				vptpnruss[i][j][k] = true;//is russian
			}
			else{
				russchin << burnin_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__burnin_" + rs;
				vptpnmburnin[i][j][k] = burnin_chin;
				vptpnruss[i][j][k] = false;//is not russian thus it is chinese
			}
		}
		else{//choose xx_EEp_yy histos
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
	//		if(fabs(eta)<2.8) continue;
			etass << eta;
			es = etass.str();
			bin = EEp_producer->FindBin(i,j);
			int prod = EEp_producer->GetBinContent(bin);
			if(prod==1){ burnin_russ = burnin_EEp->GetBinContent(bin);
			             burnin_chin = 0; }
			if(prod==2){ burnin_chin = burnin_EEp->GetBinContent(bin);
			             burnin_russ = 0; }
			if(burnin_russ==0 && burnin_chin==0) continue;
			else if(burnin_russ!=0){

				russchin << burnin_russ;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__russian__burnin_" + rs;
				vptpnmburnin[i][j][k] = burnin_russ;
				vptpnruss[i][j][k] = true;
			}
			else{
				russchin << burnin_chin;
				rs = russchin.str();
				seta = (string)"__eta_" + es + (string)"__chinese__burnin_" + rs;
				vptpnmburnin[i][j][k] = burnin_chin;
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
		if(plottingall_VPTPN && g_vptpn_las[i][j][k]->GetN()>10 && vptpnmburnin[i][j][k]>0 ){//with burnin oled instead of las
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
		if(normalized) haxis->SetMinimum(0.565);//normalized, las 0.565, orange 0.55
		else           haxis->SetMinimum(0.25);
		if(normalized) haxis->SetMaximum(1.165);//normalized, las 1.165, orange 1.30
		else           haxis->SetMaximum(1.5);
				haxis->GetXaxis()->SetTitle("time"); haxis->GetXaxis()->SetTitleOffset(1.25); haxis->GetXaxis()->SetLabelSize(0.027);
				haxis->GetYaxis()->SetTitle("VPT/PN"); haxis->GetYaxis()->SetTitleOffset(1.25); haxis->GetYaxis()->SetLabelSize(0.027);
				haxis->GetXaxis()->SetTimeDisplay(1); haxis->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
				haxis->GetXaxis()->SetLabelOffset(0.02); haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); haxis->GetXaxis()->SetNdivisions(510, true);
		g_vptpn_las[i][j][k]->SetHistogram(haxis);//oled instead of las
		col->Clear();
		haxis->Draw("hist");
		g_vptpn_las[i][j][k]->Draw("P");//oled instead of las
		TString outputfile = TString::Format("VPTPN_EE%i_IX%i_IY%i",k, i, j);//with orange
		col->Update();
		if(vptpnmburnin[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintNoEPS(col, outputfile, outputdirac);//oled instead of las
		if(vptpnmburnin[i][j][k]>0 && g_vptpn_las[i][j][k]->GetN()>10) Util::PrintEPS(col, outputfile, outputdirac);//oled instead of las
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
		double burnin_ = vptpnmburnin[i][j][k];
		if(startdiff) vptpnmax[i][j][k] = 1.;//compare to start of 2011 which is by definition 1
		double vptdiff = vptpnmax[i][j][k]-vptpnmin[i][j][k];

		if(fabs(eta)<1.6&&fabs(eta)>=1.4&&burnin_!=0)a_burnin_vs_VPTPNdiff[0][nn]->SetPoint(a_burnin_vs_VPTPNdiff[0][nn]->GetN(), burnin_,vptdiff);//1p4
		else if(fabs(eta)<( 1.8) && burnin_!=0)      a_burnin_vs_VPTPNdiff[1][nn]->SetPoint(a_burnin_vs_VPTPNdiff[1][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<( 2.0) && burnin_!=0)      a_burnin_vs_VPTPNdiff[2][nn]->SetPoint(a_burnin_vs_VPTPNdiff[2][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<( 2.2) && burnin_!=0)      a_burnin_vs_VPTPNdiff[3][nn]->SetPoint(a_burnin_vs_VPTPNdiff[3][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<( 2.4) && burnin_!=0)      a_burnin_vs_VPTPNdiff[4][nn]->SetPoint(a_burnin_vs_VPTPNdiff[4][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<( 2.6) && burnin_!=0)      a_burnin_vs_VPTPNdiff[5][nn]->SetPoint(a_burnin_vs_VPTPNdiff[5][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<( 2.8) && burnin_!=0)      a_burnin_vs_VPTPNdiff[6][nn]->SetPoint(a_burnin_vs_VPTPNdiff[6][nn]->GetN(), burnin_,vptdiff);
		else if(fabs(eta)<=(3.0) && burnin_!=0)      a_burnin_vs_VPTPNdiff[7][nn]->SetPoint(a_burnin_vs_VPTPNdiff[7][nn]->GetN(), burnin_,vptdiff);
	}}}

	TF1 *afit_burnin_vs_VPTPNdiff[8][2];
	TH1D *a_slope_vs_eta_burnin[2];
	TH1D *a_onevalue_vs_eta_burnin[2];
	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
		if(n>=8) continue;
		string aname = (string)"fit_burnin_vs_VPTdiff_AbsEta_" + string_feta[n] + string_prod[nn];
		afit_burnin_vs_VPTPNdiff[n][nn] = new TF1(aname.c_str(), "[0]+[1]*(x-1.0)", 0.75, 1.25);
		afit_burnin_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);
	}
	}
	for(int nn = 0; nn<2; ++ nn){
		string aname = (string)"slope_burnin_AbsEta_" + string_prod[nn];
		a_slope_vs_eta_burnin[nn] = new TH1D(aname.c_str(), aname.c_str(), 8, 1.4, 3.);
		a_slope_vs_eta_burnin[nn]->SetMinimum(-0.9); a_slope_vs_eta_burnin[nn]->SetMaximum( 0.2);
		a_slope_vs_eta_burnin[nn]->Sumw2();
		a_slope_vs_eta_burnin[nn]->GetYaxis()->SetTitle("burn-in slope"); a_slope_vs_eta_burnin[nn]->GetXaxis()->SetTitle("|#eta|");
		a_slope_vs_eta_burnin[nn]->SetLineWidth(4); a_slope_vs_eta_burnin[nn]->SetMarkerStyle(20); a_slope_vs_eta_burnin[nn]->SetMarkerSize(2);
		if(nn==0){ a_slope_vs_eta_burnin[nn]->SetLineColor(kRed ); a_slope_vs_eta_burnin[nn]->SetMarkerColor(kRed );}
		else     { a_slope_vs_eta_burnin[nn]->SetLineColor(kBlue); a_slope_vs_eta_burnin[nn]->SetMarkerColor(kBlue);}
		aname = (string)"onevalue_burnin_AbsEta_" + string_prod[nn];
		a_onevalue_vs_eta_burnin[nn] = new TH1D(aname.c_str(), aname.c_str(), 8, 1.4, 3.);
		a_onevalue_vs_eta_burnin[nn]->SetMinimum(0.); a_onevalue_vs_eta_burnin[nn]->SetMaximum( 0.45);
		a_onevalue_vs_eta_burnin[nn]->Sumw2();
		a_onevalue_vs_eta_burnin[nn]->GetYaxis()->SetTitle("burn-in=1 equivalent"); a_onevalue_vs_eta_burnin[nn]->GetXaxis()->SetTitle("|#eta|");
		a_onevalue_vs_eta_burnin[nn]->SetLineWidth(4); a_onevalue_vs_eta_burnin[nn]->SetMarkerStyle(20); a_onevalue_vs_eta_burnin[nn]->SetMarkerSize(2);
		if(nn==0){ a_onevalue_vs_eta_burnin[nn]->SetLineColor(kRed ); a_onevalue_vs_eta_burnin[nn]->SetMarkerColor(kRed );}
		else     { a_onevalue_vs_eta_burnin[nn]->SetLineColor(kBlue); a_onevalue_vs_eta_burnin[nn]->SetMarkerColor(kBlue);}
	}

	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
		if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(a_burnin_vs_VPTPNdiff[n][nn])->Fit(afit_burnin_vs_VPTPNdiff[n][nn], "R");
	}}

	//print out full fit parameters (including goodness of fit)
	for(int n = 0; n<8; ++n){
	for(int nn = 0; nn<2; ++ nn){
		if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		cout << "Fit " << a_burnin_vs_VPTPNdiff[n][nn]->GetName() << endl;
		cout << "Fit value at 1 = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0) << endl;
		cout << "Fit slope      = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetChisquare() << "/" << afit_burnin_vs_VPTPNdiff[n][nn]->GetNDF() << " = " << (afit_burnin_vs_VPTPNdiff[n][nn]->GetChisquare())/(afit_burnin_vs_VPTPNdiff[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << afit_burnin_vs_VPTPNdiff[n][nn]->GetProb() << endl;
		a_slope_vs_eta_burnin[nn]->SetBinContent(n+1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1) );
		a_slope_vs_eta_burnin[nn]->SetBinError  (n+1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParError( 1) );
		a_onevalue_vs_eta_burnin[nn]->SetBinContent(n+1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0) );
		a_onevalue_vs_eta_burnin[nn]->SetBinError  (n+1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParError( 0) );
	}}


	//draw correlation plots VPT/PN vs mu_std in eta bins
	for(int n = 0; n<8; ++n){//eta bins
	for(int nn= 0; nn<2;++nn){//producer: 0: russian, 1: chinese
		TH1F *haxis = a_burnin_vs_VPTPNdiff[n][nn]->GetHistogram();
		TString outputfile;
		     if(VPTdiff_vs_burnin){
			col->cd();
			col->SetName(a_burnin_vs_VPTPNdiff[n][nn]->GetName() );
			col->SetTitle(a_burnin_vs_VPTPNdiff[n][nn]->GetTitle() );
			a_burnin_vs_VPTPNdiff[n][nn]->Draw("AP");
			haxis = a_burnin_vs_VPTPNdiff[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_burnin_vs_VPTPNdiff[n][nn]->GetName(), a_burnin_vs_VPTPNdiff[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.8,1.2);
			haxis->GetXaxis()->SetLimits(0.8,1.2);
			haxis->GetXaxis()->SetTitle("burn-in ratio");
			haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("burn-in ratio");
			haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			//if(normalized) haxis->GetYaxis()->SetRangeUser(0.,0.45);
			//else           haxis->GetYaxis()->SetRangeUser(0.,0.25);
			haxis->SetMinimum(-0.1);
			haxis->SetMaximum( 0.7);
			haxis->Draw("hist");
			a_burnin_vs_VPTPNdiff[n][nn]->Draw("P");
			outputfile = a_burnin_vs_VPTPNdiff[n][nn]->GetTitle();
			//col->Update();
			if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintNoEPS(col, outputfile, outputdir);
			if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3) Util::PrintEPS(col, outputfile, outputdir);
			col->Clear();
		     }

	}
	}
	for(int nn= 0; nn<2;++nn){//producer: 0: russian, 1: chinese
		TString outputfile;
		col->cd();
		col->SetName(a_slope_vs_eta_burnin[nn]->GetName() );
		col->SetTitle(a_slope_vs_eta_burnin[nn]->GetTitle() );
		a_slope_vs_eta_burnin[nn]->Draw();
		outputfile = a_slope_vs_eta_burnin[nn]->GetName();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
		col->cd();
		col->SetName(a_onevalue_vs_eta_burnin[nn]->GetName() );
		col->SetTitle(a_onevalue_vs_eta_burnin[nn]->GetTitle() );
		a_onevalue_vs_eta_burnin[nn]->Draw();
		outputfile = a_onevalue_vs_eta_burnin[nn]->GetName();
		Util::PrintNoEPS(col, outputfile, outputdir);
		Util::PrintEPS(col, outputfile, outputdir);
		col->SaveAs(outputdir + "/" + outputfile + ".C");
		col->Clear();
	}
//xxyyzz
/*	//store histos into file
	TFile *newFile;
	if(startdiff) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121120_burnin_vs_VPTPNmaxMinusVPTPNmin_diffwrt2011begin_VectorNTuples"+periods+".root", "RECREATE");
	if(normalized) newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121120_burnin_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples"+periods+".root", "RECREATE");//0312->0410; //muSIC <-> burnin
	else if(!startdiff)  newFile = new TFile("~/ECAL/Laser/Analysis/FileOutputs/20121120_burnin_vs_VPTPNmaxMinusVPTPNmin_nonormalization_VectorNTuples"+periods+".root", "RECREATE");//0312->0410
	newFile->cd();
	for(int n = 0; n<8; ++n){ for(int nn= 0; nn<2;++nn) {
	a_burnin_vs_VPTPNdiff[n][nn]->Write(); 
	afit_burnin_vs_VPTPNdiff[n][nn]->Write();
	}}
	for(int nn= 0; nn<2;++nn) {
	a_slope_vs_eta_burnin[nn]->Write();
	a_onevalue_vs_eta_burnin[nn]->Write();
	}

	cout << "histograms for period " << period << " called " << periods << " (normalized/startdiff = " << normalized << "/" << startdiff << ") stored in " << newFile->GetName() << endl;
*/
}
