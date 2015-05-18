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
#include <TProfile.h>
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

void EEstudies_VectorNTuples_everything_continued_daily(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool p_png          = true;
	bool p_eps          = true;
	bool p_c            = true;

	//default trues: fitwithouterror, ploterrshaded. plotoverview, fprinttime, savefinalplots, makecorrcoeffplots
	bool fitprofile     = false;//fit profile plot of mustd/burnin vs. VPT diff correlation instead of scatter plot
	bool fitprofile2    = false;//fit profile plot going through 2d histo first, only true if fitprofile==true, only for mustd and burnin

	bool fitwithouterror= true; //set the error in the scatter plot to 0. / not for profile right now (although easy to implement)

	bool ploterrshaded  = true; // plot error as shaded band
	bool dontploterrer  = false;//do not plot error
	bool cleanviaSIC    = false;//clean mustd for chinese by requiring muSIC>0
	bool cleanviaSICo   = false;//clean also other variables with upper requirement // only true if cleanviaSIC==true

	bool plotoverview   = true;
	bool fprinttime     = true;	//prints out VPTdiff vs. XXX for a given time (which is after technical stop 09/2011)
	bool fitlinepluserr = true;	//also prints the error on the VPTdiff vs. XXX (only for burnin / mustd)
	bool savefinalplots = true;
	bool makecorrcoeffplots      = true;//only for mustd/burnin, correlation plot

	bool useMuSICinsteafofMuStd  = false;

	bool makemustdplots          = true;
	bool makeburninplots         = true;//make true xxyyzz
	bool makeburnincorrplots     = false;
	bool makeburninvsmustdplots  = true;//make true xxyyzz
	bool makesplittedburninplots = false;
	bool makegainqeplots         = false;
	bool makegainplots           = false;
	bool makeqeplots             = false;

	bool rotatesc          = false;
	vector<int>  whichscvec; 	whichscvec.clear();
	vector<bool> rotatesc_m90vec; 	rotatesc_m90vec.clear();
	vector<bool> rotatesc_p90vec; 	rotatesc_p90vec.clear();
	vector<bool> rotatesc_180vec; 	rotatesc_180vec.clear();

	whichscvec.push_back(3);	//	whichscvec.push_back(17);		whichscvec.push_back(15);		whichscvec.push_back(9);
	rotatesc_m90vec.push_back(false );//	rotatesc_m90vec.push_back(true);	rotatesc_m90vec.push_back(false);	rotatesc_m90vec.push_back(false);
	rotatesc_p90vec.push_back(false);//	rotatesc_p90vec.push_back(false);	rotatesc_p90vec.push_back(true);	rotatesc_p90vec.push_back(false);
	rotatesc_180vec.push_back(false);//	rotatesc_180vec.push_back(false);	rotatesc_180vec.push_back(false);	rotatesc_180vec.push_back(true);
	whichscvec.push_back(17);
	rotatesc_m90vec.push_back(false);
	rotatesc_p90vec.push_back(false);
	rotatesc_180vec.push_back(false);
	whichscvec.push_back(15);
	rotatesc_m90vec.push_back(false);
	rotatesc_p90vec.push_back(false);
	rotatesc_180vec.push_back(false);
	whichscvec.push_back(9);
	rotatesc_m90vec.push_back(false);
	rotatesc_p90vec.push_back(false);
	rotatesc_180vec.push_back(false);
//	int  whichsc           = 3;
//	bool rotatesc_m90      = false;
//	bool rotatesc_p90      = false;
//	bool rotatesc_180      = true;


	string daily = "/daily";
	if(fitwithouterror) daily = "/tempdaily";
	if(cleanviaSIC) daily = daily + "/SICcleaned";
	if(cleanviaSIC&&cleanviaSICo) daily = daily + "All";
	if(dontploterrer) daily = daily + "/errorsurpressed";
	string strf = "continued";
	if(fitprofile) strf = "profilefit/continued";
	if(fitprofile&&fitprofile2) strf = "profilefit2/continued";
	TString outputdir               = "tempPlots/20130402/EEstudies/"+strf+daily;
	TString outdirburnin            = "tempPlots/20130402/EEstudies/"+strf+daily+"/burnin";
	TString outdirburninmustdcorr   = "tempPlots/20130402/EEstudies/"+strf+daily+"/burnin/mustdslopecorrected";
	TString outdirburninsplitgainqe = "tempPlots/20130402/EEstudies/"+strf+daily+"/burnin/splittedgainqe";
	TString outdirmustd             = "tempPlots/20130402/EEstudies/"+strf+daily+"/mustd";
	TString outdirmustdcorrburnin   = "tempPlots/20130402/EEstudies/"+strf+daily+"/mustd/correlatewithburnin";
	if(useMuSICinsteafofMuStd) outdirmustd             = "tempPlots/20130402/EEstudies/"+strf+daily+"/muSIC";
	if(useMuSICinsteafofMuStd) outdirmustdcorrburnin   = "tempPlots/20130402/EEstudies/"+strf+daily+"/muSIC/correlatewithburnin";
	TString outdirgain              = "tempPlots/20130402/EEstudies/"+strf+daily+"/gain";
	TString outdirqe                = "tempPlots/20130402/EEstudies/"+strf+daily+"/qe";
	TString outdirgainqe            = "tempPlots/20130402/EEstudies/"+strf+daily+"/gainqe";
	if(rotatesc){
		TString wsc = "_SC/SC";
		for(unsigned int w = 0; w<whichscvec.size(); ++w){
		if(rotatesc_m90vec[w] || rotatesc_p90vec[w] || rotatesc_180vec[w]) wsc = wsc + TString::Format("%d", whichscvec[w]);
		if(rotatesc_m90vec[w]) wsc = wsc + TString::Format("m");
		if(rotatesc_p90vec[w]) wsc = wsc + TString::Format("p");
		if(rotatesc_180vec[w]) wsc = wsc + TString::Format("pp");
		}
		string rotatestring = (string)"rotate"+wsc.Data()+(string)"/";
		if(fitprofile) rotatestring = (string)"rotate"+wsc.Data()+(string)"/profilefit/";
		if(fitprofile&&fitprofile2) rotatestring = (string)"rotate"+wsc.Data()+(string)"/profilefit2/";
		outputdir               = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily;
		outdirburnin            = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/burnin";
		outdirburninmustdcorr   = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/burnin/mustdslopecorrected";
		outdirburninsplitgainqe = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/burnin/splittedgainqe";
		outdirmustd             = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/mustd";
		outdirmustdcorrburnin   = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/mustd/correlatewithburnin";
		if(useMuSICinsteafofMuStd) outdirmustd             = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/muSIC";
		if(useMuSICinsteafofMuStd) outdirmustdcorrburnin   = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/muSIC/correlatewithburnin";
		outdirgain              = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/gain";
		outdirqe                = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/qe";
		outdirgainqe            = "Plots/20130402/EEstudies/"+rotatestring+"continued"+daily+"/gainqe";
	}
	Util::MakeOutputDir(outputdir);
	if(makeburninplots)         Util::MakeOutputDir(outdirburnin);
	if(makeburnincorrplots)     Util::MakeOutputDir(outdirburninmustdcorr);
	if(makesplittedburninplots) Util::MakeOutputDir(outdirburninsplitgainqe);
	if(makemustdplots)          Util::MakeOutputDir(outdirmustd);
	if(makeburninvsmustdplots)  Util::MakeOutputDir(outdirmustdcorrburnin);
	if(makegainplots)           Util::MakeOutputDir(outdirgain);
	if(makeqeplots)             Util::MakeOutputDir(outdirqe);
	if(makegainqeplots)         Util::MakeOutputDir(outdirgainqe);


	TFile *stdmapsfile     = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEmustd.root");
	TH2D *EEp_producer     = (TH2D*)stdmapsfile->Get("EEp_producer");//1 russian, 2 chinese
	TH2D *EEm_producer     = (TH2D*)stdmapsfile->Get("EEm_producer");
	TH2D *mu_std_EEp       = (TH2D*)stdmapsfile->Get("mu_std_EEp");
	TH2D *mu_std_EEm       = (TH2D*)stdmapsfile->Get("mu_std_EEm");
	TH2D *mu_SIC_EEp       = (TH2D*)stdmapsfile->Get("mu_SIC_EEp");
	TH2D *mu_SIC_EEm       = (TH2D*)stdmapsfile->Get("mu_SIC_EEm");
	TH2D *VPT_gainqe_EEp   = (TH2D*)stdmapsfile->Get("VPT_gainqe_EEp");
	TH2D *VPT_gainqe_EEm   = (TH2D*)stdmapsfile->Get("VPT_gainqe_EEm");
	TH2D *VPT_gain_EEp     = (TH2D*)stdmapsfile->Get("VPT_gain_EEp");
	TH2D *VPT_gain_EEm     = (TH2D*)stdmapsfile->Get("VPT_gain_EEm");
	TH2D *VPT_qe_EEp       = (TH2D*)stdmapsfile->Get("VPT_qe_EEp");
	TH2D *VPT_qe_EEm       = (TH2D*)stdmapsfile->Get("VPT_qe_EEm");

	TFile *burninfile      = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEburnin.root");
	TH2D *burnin_EEm       = (TH2D*)burninfile->Get("burnin_EEm");
	TH2D *burnin_EEp       = (TH2D*)burninfile->Get("burnin_EEp");

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	TFile *scfile          = TFile::Open("/shome/haweber/ECAL/DataLaser/20130118_EESC_center_alt.root");
	TH2D *SCcenter_EEm    = (TH2D*)scfile->Get("scbarcode_center_EEm");
	TH2D *SCcenter_EEp    = (TH2D*)scfile->Get("scbarcode_center_EEp");

	//mustd
	TGraphErrors* g_mustd_corrcoeff[8][2];
	TGraphErrors* g_mustd_slopes[8][2];
	TGraphErrors* g_mustd_intercepts[8][2];
	TGraphErrors* a_mustd_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TProfile* p_mustd_vs_VPTPNdiff[8][2];
//	TH2D* h_mustd_vs_VPTPNdiff[8][2];
	TH1D* pc_mustd_vs_VPTPNdiff[8][2];
	TF1 *afit_mustd_vs_VPTPNdiff[8][2];
	//burnin
	TGraphErrors* g_burnin_corrcoeff[8][2];
	TGraphErrors* g_burnin_slopes[8][2];
	TGraphErrors* g_burnin_onevalues[8][2];
	TGraphErrors* a_burnin_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TProfile* p_burnin_vs_VPTPNdiff[8][2];
//	TH2D* h_burnin_vs_VPTPNdiff[8][2];
	TH1D* pc_burnin_vs_VPTPNdiff[8][2];
	TF1 *afit_burnin_vs_VPTPNdiff[8][2];
	//burnin - mustd_corr
	TGraphErrors* g_burnin_corrected_slopes[8][2];
	TGraphErrors* g_burnin_corrected_onevalues[8][2];
	TGraphErrors* a_burnin_corrected_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TF1 *afit_burnin_corrected_vs_VPTPNdiff[8][2];
	TProfile* p_burnin_corrected_vs_VPTPNdiff[8][2];
	//mustd burnin fit parameter correlations
	TGraphErrors* g_burnin_slope_vs_mustd_intercept[8][2];
	TGraphErrors* g_burnin_onevalue_vs_mustd_slope[8][2];
	//burnin in 2 bins of gain, gain*qe or qe
	TGraphErrors* g_splittedburnin_slopes[8][2][6];
	TGraphErrors* g_splittedburnin_onevalues[8][2][6];
	TGraphErrors* a_splittedburnin_vs_VPTPNdiff[8][2][6];//usetgrapherrors = true;
	TF1 *afit_splittedburnin_vs_VPTPNdiff[8][2][6];
	//qe*gain
	TGraphErrors* g_gainqe_slopes[8][2];
	TGraphErrors* g_gainqe_intercepts[8][2];
	TGraphErrors* a_gainqe_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TF1 *afit_gainqe_vs_VPTPNdiff[8][2];
	//gain
	TGraphErrors* g_gain_slopes[8][2];
	TGraphErrors* g_gain_intercepts[8][2];
	TGraphErrors* a_gain_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TF1 *afit_gain_vs_VPTPNdiff[8][2];
	//quantum efficiency
	TGraphErrors* g_qe_slopes[8][2];
	TGraphErrors* g_qe_intercepts[8][2];
	TGraphErrors* a_qe_vs_VPTPNdiff[8][2];//usetgrapherrors = true;
	TF1 *afit_qe_vs_VPTPNdiff[8][2];

	TF1 *fitup;
	TF1 *fitdown;

	TH1D *h_burnin_eq1[2];
	TH1D *h_burnin_slope[2];
	TH1D *h_mustd_intercept[2];
	TH1D *h_mustd_slope[2];

	TF1 *fmustdtemp  = new TF1("fmustdtemp", "[0]+[1]*(x-0.0)", 0.0, 2.0);
	fmustdtemp ->SetLineColor(kRed);
	TF1 *fburnintemp = new TF1("fburnintemp", "[0]+[1]*(x-1.0)", 0.8, 1.1);
	fburnintemp->SetLineColor(kRed);

	TLegend *Legend1  = new TLegend(0.7455357,0.5600649,0.9598214,0.9577922,NULL,"brNDC");
	Legend1->SetBorderSize(0);
	Legend1->SetTextSize(0.04575163);
	Legend1->SetLineColor(1);
	Legend1->SetLineStyle(1);
	Legend1->SetLineWidth(1);
	Legend1->SetFillColor(0);
	Legend1->SetFillStyle(1001);
   	TH1F* h[8]; //needed for legend

	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	string string_eta_leg[8] = {"1.4 #leq |#eta| < 1.6", "1.6 #leq |#eta| < 1.8", "1.8 #leq |#eta| < 2.0", "2.0 #leq |#eta| < 2.2", "2.2 #leq |#eta| < 2.4", "2.4 #leq |#eta| < 2.6", "2.6 #leq |#eta| < 2.8", "2.8 #leq |#eta|"};
	
	for(int n = 0; n<8; ++n){//eta bin
		h[n] = new TH1F("", (string_eta_leg[n]).c_str(), 1, 0, 1); 
	for(int nn = 0; nn<2; ++ nn){//producer
		string hname;

		hname = (string)"mustd_corrcoeff_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_mustd_corrcoeff[n][nn] = new TGraphErrors();
		g_mustd_corrcoeff[n][nn]->SetName(hname.c_str());
		g_mustd_corrcoeff[n][nn]->GetYaxis()->SetTitle("#mu-#rho"); g_mustd_corrcoeff[n][nn]->GetXaxis()->SetTitle("time");
		g_mustd_corrcoeff[n][nn]->SetMarkerSize(1); g_mustd_corrcoeff[n][nn]->SetMarkerStyle(20); 
		if(ploterrshaded) g_mustd_corrcoeff[n][nn]->SetFillStyle(3002);
		hname = (string)"mustd_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_mustd_slopes[n][nn] = new TGraphErrors();
		g_mustd_slopes[n][nn]->SetName(hname.c_str());
		g_mustd_slopes[n][nn]->GetYaxis()->SetTitle("#mu-slope (m)"); g_mustd_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_mustd_slopes[n][nn]->SetMarkerSize(1); g_mustd_slopes[n][nn]->SetMarkerStyle(20); 
		if(ploterrshaded) g_mustd_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"mustd_intercept_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_mustd_intercepts[n][nn] = new TGraphErrors();
		g_mustd_intercepts[n][nn]->SetName(hname.c_str());
		g_mustd_intercepts[n][nn]->GetYaxis()->SetTitle("#mu-intercept"); g_mustd_intercepts[n][nn]->GetXaxis()->SetTitle("time");
		g_mustd_intercepts[n][nn]->SetMarkerSize(1); g_mustd_intercepts[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_mustd_intercepts[n][nn]->SetFillStyle(3002);

		hname = (string)"burnin_corrcoeff_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_corrcoeff[n][nn] = new TGraphErrors();
		g_burnin_corrcoeff[n][nn]->SetName(hname.c_str());
		g_burnin_corrcoeff[n][nn]->GetYaxis()->SetTitle("burn-in #rho"); g_burnin_corrcoeff[n][nn]->GetXaxis()->SetTitle("time");
		g_burnin_corrcoeff[n][nn]->SetMarkerSize(1); g_burnin_corrcoeff[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_corrcoeff[n][nn]->SetFillStyle(3002);
		hname = (string)"burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_slopes[n][nn] = new TGraphErrors();
		g_burnin_slopes[n][nn]->SetName(hname.c_str());
		g_burnin_slopes[n][nn]->GetYaxis()->SetTitle("burn-in slope"); g_burnin_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_burnin_slopes[n][nn]->SetMarkerSize(1); g_burnin_slopes[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_onevalues[n][nn] = new TGraphErrors();
		g_burnin_onevalues[n][nn]->SetName(hname.c_str());
		g_burnin_onevalues[n][nn]->GetYaxis()->SetTitle("burn-in = 1"); g_burnin_onevalues[n][nn]->GetXaxis()->SetTitle("time");
		g_burnin_onevalues[n][nn]->SetMarkerSize(1); g_burnin_onevalues[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_onevalues[n][nn]->SetFillStyle(3002);

		hname = (string)"burnin_corrected_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_corrected_slopes[n][nn] = new TGraphErrors();
		g_burnin_corrected_slopes[n][nn]->SetName(hname.c_str());
		g_burnin_corrected_slopes[n][nn]->GetYaxis()->SetTitle("burn-in slope"); g_burnin_corrected_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_burnin_corrected_slopes[n][nn]->SetMarkerSize(1); g_burnin_corrected_slopes[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_corrected_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"burnin_corrected_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_corrected_onevalues[n][nn] = new TGraphErrors();
		g_burnin_corrected_onevalues[n][nn]->SetName(hname.c_str());
		g_burnin_corrected_onevalues[n][nn]->GetYaxis()->SetTitle("burn-in = 1"); g_burnin_corrected_onevalues[n][nn]->GetXaxis()->SetTitle("time");
		g_burnin_corrected_onevalues[n][nn]->SetMarkerSize(1); g_burnin_corrected_onevalues[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_corrected_onevalues[n][nn]->SetFillStyle(3002);

		hname = (string)"burnin_slope_vs_mustd_intercept_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_slope_vs_mustd_intercept[n][nn] = new TGraphErrors();
		g_burnin_slope_vs_mustd_intercept[n][nn]->SetName(hname.c_str());
		g_burnin_slope_vs_mustd_intercept[n][nn]->GetYaxis()->SetTitle("burn-in slope"); g_burnin_slope_vs_mustd_intercept[n][nn]->GetXaxis()->SetTitle("#mu-intercept");
		g_burnin_slope_vs_mustd_intercept[n][nn]->SetMarkerSize(1); g_burnin_slope_vs_mustd_intercept[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_slope_vs_mustd_intercept[n][nn]->SetFillStyle(3002);
		hname = (string)"burnin_onevalue_vs_mustd_slope_Eta_" + string_feta[n] + string_prod[nn];
		g_burnin_onevalue_vs_mustd_slope[n][nn] = new TGraphErrors();
		g_burnin_onevalue_vs_mustd_slope[n][nn]->SetName(hname.c_str());
		g_burnin_onevalue_vs_mustd_slope[n][nn]->GetYaxis()->SetTitle("burn-in = 1"); g_burnin_onevalue_vs_mustd_slope[n][nn]->GetXaxis()->SetTitle("#mu-slope");
		g_burnin_onevalue_vs_mustd_slope[n][nn]->SetMarkerSize(1); g_burnin_onevalue_vs_mustd_slope[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_burnin_onevalue_vs_mustd_slope[n][nn]->SetFillStyle(3002);

		hname = (string)"gainqe_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_gainqe_slopes[n][nn] = new TGraphErrors();
		g_gainqe_slopes[n][nn]->SetName(hname.c_str());
		g_gainqe_slopes[n][nn]->GetYaxis()->SetTitle("q.e. #cdot gain-slope"); g_gainqe_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_gainqe_slopes[n][nn]->SetMarkerSize(1); g_gainqe_slopes[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_gainqe_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"gainqe_intercept_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_gainqe_intercepts[n][nn] = new TGraphErrors();
		g_gainqe_intercepts[n][nn]->SetName(hname.c_str());
		g_gainqe_intercepts[n][nn]->GetYaxis()->SetTitle("q.e. #cdot gain-intercept"); g_gainqe_intercepts[n][nn]->GetXaxis()->SetTitle("time");
		g_gainqe_intercepts[n][nn]->SetMarkerSize(1); g_gainqe_intercepts[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_gainqe_intercepts[n][nn]->SetFillStyle(3002);
	
		hname = (string)"gain_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_gain_slopes[n][nn] = new TGraphErrors();
		g_gain_slopes[n][nn]->SetName(hname.c_str());
		g_gain_slopes[n][nn]->GetYaxis()->SetTitle("gain-slope"); g_gain_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_gain_slopes[n][nn]->SetMarkerSize(1); g_gain_slopes[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_gain_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"gain_intercept_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_gain_intercepts[n][nn] = new TGraphErrors();
		g_gain_intercepts[n][nn]->SetName(hname.c_str());
		g_gain_intercepts[n][nn]->GetYaxis()->SetTitle("gain-intercept"); g_gain_intercepts[n][nn]->GetXaxis()->SetTitle("time");
		g_gain_intercepts[n][nn]->SetMarkerSize(1); g_gain_intercepts[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_gain_intercepts[n][nn]->SetFillStyle(3002);

		hname = (string)"qe_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_qe_slopes[n][nn] = new TGraphErrors();
		g_qe_slopes[n][nn]->SetName(hname.c_str());
		g_qe_slopes[n][nn]->GetYaxis()->SetTitle("q.e.-slope"); g_qe_slopes[n][nn]->GetXaxis()->SetTitle("time");
		g_qe_slopes[n][nn]->SetMarkerSize(1); g_qe_slopes[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_qe_slopes[n][nn]->SetFillStyle(3002);
		hname = (string)"qe_intercept_vs_time_Eta_" + string_feta[n] + string_prod[nn];
		g_qe_intercepts[n][nn] = new TGraphErrors();
		g_qe_intercepts[n][nn]->SetName(hname.c_str());
		g_qe_intercepts[n][nn]->GetYaxis()->SetTitle("q.e.-intercept"); g_qe_intercepts[n][nn]->GetXaxis()->SetTitle("time");
		g_qe_intercepts[n][nn]->SetMarkerSize(1); g_qe_intercepts[n][nn]->SetMarkerStyle(20);
		if(ploterrshaded) g_qe_intercepts[n][nn]->SetFillStyle(3002);

		for(int nnn = 0; nnn<6; ++nnn){
			if(nnn==0) hname = (string) "lowqe_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==1) hname = (string)"highqe_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==2) hname = (string) "lowgain_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==3) hname = (string)"highgain_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==4) hname = (string) "lowgainqe_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==5) hname = (string)"highgainqe_burnin_slopes_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			g_splittedburnin_slopes[n][nn][nnn] = new TGraphErrors();
			g_splittedburnin_slopes[n][nn][nnn]->SetName(hname.c_str());
			g_splittedburnin_slopes[n][nn][nnn]->GetYaxis()->SetTitle("burn-in slope"); 
			g_splittedburnin_slopes[n][nn][nnn]->GetXaxis()->SetTitle("time");
			g_splittedburnin_slopes[n][nn][nnn]->SetMarkerSize(1); g_splittedburnin_slopes[n][nn][nnn]->SetMarkerStyle(20);
			if(ploterrshaded) g_splittedburnin_slopes[n][nn][nnn]->SetFillStyle(3002);
			if(nnn==0) hname = (string) "lowqe_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==1) hname = (string)"highqe_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==2) hname = (string) "lowgain_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==3) hname = (string)"highgain_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==4) hname = (string) "lowgainqe_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==5) hname = (string)"highgainqe_burnin_onevalue_vs_time_Eta_" + string_feta[n] + string_prod[nn];
			g_splittedburnin_onevalues[n][nn][nnn] = new TGraphErrors();
			g_splittedburnin_onevalues[n][nn][nnn]->SetName(hname.c_str());
			g_splittedburnin_onevalues[n][nn][nnn]->GetYaxis()->SetTitle("burn-in = 1"); 
			g_splittedburnin_onevalues[n][nn][nnn]->GetXaxis()->SetTitle("time");
			g_splittedburnin_onevalues[n][nn][nnn]->SetMarkerSize(1); g_splittedburnin_onevalues[n][nn][nnn]->SetMarkerStyle(20);
			if(ploterrshaded) g_splittedburnin_onevalues[n][nn][nnn]->SetFillStyle(3002);
		}

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
		g_mustd_corrcoeff[n][nn]                ->SetLineColor(color);   g_mustd_corrcoeff[n][nn]                ->SetMarkerColor(color);
		g_mustd_slopes[n][nn]                   ->SetLineColor(color);   g_mustd_slopes[n][nn]                   ->SetMarkerColor(color);
		g_mustd_intercepts[n][nn]               ->SetLineColor(color);   g_mustd_intercepts[n][nn]               ->SetMarkerColor(color);
		g_burnin_corrcoeff[n][nn]               ->SetLineColor(color);   g_burnin_corrcoeff[n][nn]               ->SetMarkerColor(color);
		g_burnin_slopes[n][nn]                  ->SetLineColor(color);   g_burnin_slopes[n][nn]                  ->SetMarkerColor(color);
		g_burnin_onevalues[n][nn]               ->SetLineColor(color);   g_burnin_onevalues[n][nn]               ->SetMarkerColor(color);
		g_burnin_corrected_slopes[n][nn]        ->SetLineColor(color);   g_burnin_corrected_slopes[n][nn]        ->SetMarkerColor(color);
		g_burnin_corrected_onevalues[n][nn]     ->SetLineColor(color);   g_burnin_corrected_onevalues[n][nn]     ->SetMarkerColor(color);
		g_burnin_slope_vs_mustd_intercept[n][nn]->SetLineColor(color);   g_burnin_slope_vs_mustd_intercept[n][nn]->SetMarkerColor(color);
		g_burnin_onevalue_vs_mustd_slope[n][nn] ->SetLineColor(color);   g_burnin_onevalue_vs_mustd_slope[n][nn] ->SetMarkerColor(color);
		g_gainqe_slopes[n][nn]                  ->SetLineColor(color);   g_gainqe_slopes[n][nn]                  ->SetMarkerColor(color);
		g_gainqe_intercepts[n][nn]              ->SetLineColor(color);   g_gainqe_intercepts[n][nn]              ->SetMarkerColor(color);
		g_gain_slopes[n][nn]                    ->SetLineColor(color);   g_gain_slopes[n][nn]                    ->SetMarkerColor(color);
		g_gain_intercepts[n][nn]                ->SetLineColor(color);   g_gain_intercepts[n][nn]                ->SetMarkerColor(color);
		g_qe_slopes[n][nn]                      ->SetLineColor(color);   g_qe_slopes[n][nn]                      ->SetMarkerColor(color);
		g_qe_intercepts[n][nn]                  ->SetLineColor(color);   g_qe_intercepts[n][nn]                  ->SetMarkerColor(color);
		if(ploterrshaded){
			g_mustd_corrcoeff[n][nn]                ->SetFillColor(color);
			g_mustd_slopes[n][nn]                   ->SetFillColor(color);
			g_mustd_intercepts[n][nn]               ->SetFillColor(color);
			g_burnin_corrcoeff[n][nn]               ->SetFillColor(color);
			g_burnin_slopes[n][nn]                  ->SetFillColor(color);
			g_burnin_onevalues[n][nn]               ->SetFillColor(color);
			g_burnin_corrected_slopes[n][nn]        ->SetFillColor(color);
			g_burnin_corrected_onevalues[n][nn]     ->SetFillColor(color);
			g_burnin_slope_vs_mustd_intercept[n][nn]->SetFillColor(color);
			g_burnin_onevalue_vs_mustd_slope[n][nn] ->SetFillColor(color);
			g_gainqe_slopes[n][nn]                  ->SetFillColor(color);
			g_gainqe_intercepts[n][nn]              ->SetFillColor(color);
			g_gain_slopes[n][nn]                    ->SetFillColor(color);
			g_gain_intercepts[n][nn]                ->SetFillColor(color);
			g_qe_slopes[n][nn]                      ->SetFillColor(color);
			g_qe_intercepts[n][nn]                  ->SetFillColor(color);
		}
		for(int nnn = 0; nnn<6; ++nnn){
			g_splittedburnin_slopes[n][nn][nnn]        ->SetLineColor(color); g_splittedburnin_slopes[n][nn][nnn]    ->SetMarkerColor(color);
			g_splittedburnin_onevalues[n][nn][nnn]     ->SetLineColor(color); g_splittedburnin_onevalues[n][nn][nnn] ->SetMarkerColor(color);
			if(ploterrshaded){
				g_splittedburnin_slopes[n][nn][nnn]->SetFillColor(color); g_splittedburnin_onevalues[n][nn][nnn] ->SetFillColor(color);
			}
		}
		h[n]                ->SetFillColor(color);

		if(n==0){
			hname = (string)"mustd_intercept_vs_eta" + string_prod[nn];
			h_mustd_intercept[nn] = new TH1D(hname.c_str(), "", 8, 1.4, 3.0);
			hname = (string)"mustd_slope_vs_eta" + string_prod[nn];
			h_mustd_slope[nn]    = new TH1D(hname.c_str(), "", 8, 1.4, 3.0);
			h_mustd_intercept[nn]->SetMinimum(0.); h_mustd_intercept[nn]->SetMaximum(0.175);
			h_mustd_slope[nn]    ->SetMinimum(0.); h_mustd_slope[nn]    ->SetMaximum(0.13);
			hname = (string)"burnin_eq1_vs_eta" + string_prod[nn];
			h_burnin_eq1[nn] = new TH1D(hname.c_str(), "", 8, 1.4, 3.0);
			hname = (string)"burnin_slope_vs_eta" + string_prod[nn];
			h_burnin_slope[nn] = new TH1D(hname.c_str(), "", 8, 1.4, 3.0);
			h_burnin_eq1[nn]  ->SetMinimum(0.);   h_burnin_eq1[nn]  ->SetMaximum(0.45);
			h_burnin_slope[nn]->SetMinimum(-0.9); h_burnin_slope[nn]->SetMaximum(0.25);
			h_mustd_intercept[nn]->GetYaxis()->SetTitle("#mu-intercept"); h_mustd_intercept[nn]->GetXaxis()->SetTitle("|#eta|");
			h_mustd_slope[nn]    ->GetYaxis()->SetTitle("#mu-slope (m)"); h_mustd_slope[nn]    ->GetXaxis()->SetTitle("|#eta|");
			h_burnin_eq1[nn]     ->GetYaxis()->SetTitle("burn-in = 1");   h_burnin_eq1[nn]     ->GetXaxis()->SetTitle("|#eta|");
			h_burnin_slope[nn]   ->GetYaxis()->SetTitle("burn-in slope"); h_burnin_slope[nn]   ->GetXaxis()->SetTitle("|#eta|");
			h_mustd_intercept[nn]->SetLineWidth(4); h_mustd_intercept[nn]->SetMarkerStyle(20); h_mustd_intercept[nn]->SetMarkerSize(2);
			h_mustd_slope[nn]    ->SetLineWidth(4); h_mustd_slope[nn]    ->SetMarkerStyle(20); h_mustd_slope[nn]    ->SetMarkerSize(2);
			h_burnin_eq1[nn]     ->SetLineWidth(4); h_burnin_eq1[nn]     ->SetMarkerStyle(20); h_burnin_eq1[nn]     ->SetMarkerSize(2);
			h_burnin_slope[nn]   ->SetLineWidth(4); h_burnin_slope[nn]   ->SetMarkerStyle(20); h_burnin_slope[nn]   ->SetMarkerSize(2);
			if(nn==0){
			h_mustd_intercept[nn]->SetLineColor(kRed);  h_mustd_intercept[nn]->SetMarkerColor(kRed);
			h_mustd_slope[nn]    ->SetLineColor(kRed);  h_mustd_slope[nn]    ->SetMarkerColor(kRed);
			h_burnin_eq1[nn]     ->SetLineColor(kRed);  h_burnin_eq1[nn]     ->SetMarkerColor(kRed);
			h_burnin_slope[nn]   ->SetLineColor(kRed);  h_burnin_slope[nn]   ->SetMarkerColor(kRed);
			}
			else {
			h_mustd_intercept[nn]->SetLineColor(kBlue); h_mustd_intercept[nn]->SetMarkerColor(kBlue);
			h_mustd_slope[nn]    ->SetLineColor(kBlue); h_mustd_slope[nn]    ->SetMarkerColor(kBlue);
			h_burnin_eq1[nn]     ->SetLineColor(kBlue); h_burnin_eq1[nn]     ->SetMarkerColor(kBlue);
			h_burnin_slope[nn]   ->SetLineColor(kBlue); h_burnin_slope[nn]   ->SetMarkerColor(kBlue);
			}
		}

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

//	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	TCanvas *col = new TCanvas("col", "",404,289,900,646);
	gStyle->SetOptStat(0);
	col->Range(1.287727e+09,-0.3061225,1.370522e+09,0.6306122);
	col->SetFillColor(0);
	col->SetBorderMode(0);
	col->SetBorderSize(2);
	col->SetLeftMargin(0.1238839);
	col->SetRightMargin(0.01785714);
	col->SetTopMargin(0.03267974);
	col->SetBottomMargin(0.1666667);
	col->SetFrameBorderMode(0);
	col->SetFrameBorderMode(0);

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
		else      sprintf(gname,"tg_vptpn_las_ix%d_iy%d_EEp", i, j);
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
	if((d-1293836400+1293840000)%(25*86400)==0) cout << "day " << int((d-1293836400)/(86400)) << ", UTC time " << int(d) << endl;
	bool printtime = false;
	TString daydir = "/dummy";
	if(d<1316813137 && d>(1316813137-100000) && fprinttime) {printtime = true; daydir = "/day20110923";}//23/09/2011
	if(d<1349215200 && d>(1349215200-100000) && fprinttime) {printtime = true; daydir = "/day20121002";}//02/10/2012

	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		gStyle->SetOptStat(0);
		string hname;
		//create here
		//"VPT/PN_{max}-VPT/PN_{min}" --> #Delta #frac{VPT}{PN}
		hname = (string)"mustd_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_mustd_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_mustd_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.,2.0);
		a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_mustd_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		hname = (string)"fit_mustd_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_mustd_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-0.0)", 0.0, 2.0);
		afit_mustd_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

		hname = (string)"profile_mustd_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		if(fitprofile2) p_mustd_vs_VPTPNdiff[n][nn] = new TProfile(hname.c_str(),"", 20, 0., 2.0,"s");//xxyyzz 20-->40
		else            p_mustd_vs_VPTPNdiff[n][nn] = new TProfile(hname.c_str(),"", 20, 0., 2.0);
		p_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		p_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); p_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");
		hname = (string)"profile_clean_mustd_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		pc_mustd_vs_VPTPNdiff[n][nn] = new TH1D(hname.c_str(),"", 20, 0., 2.0);
		pc_mustd_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		pc_mustd_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); pc_mustd_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("#mu_{std} (m^{-1})");

		hname = (string)"burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_burnin_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_burnin_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.85,1.05);
		a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_burnin_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in");
		hname = (string)"fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_burnin_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-1.0)", 0.8, 1.1);
		afit_burnin_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);
		if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()!=0) cout << "Error graph " << n << " " << nn << "  has " << a_burnin_vs_VPTPNdiff[n][nn]->GetN() << " entries, but should have 0" << endl;

		hname = (string)"profile_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		if(fitprofile2) p_burnin_vs_VPTPNdiff[n][nn] = new TProfile(hname.c_str(),"", 20, 0.85, 1.05,"s");//xxyyzz 30-->60
		else            p_burnin_vs_VPTPNdiff[n][nn] = new TProfile(hname.c_str(),"", 20, 0.85, 1.05);
		p_burnin_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		p_burnin_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); p_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in");
		hname = (string)"profile_clean_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		pc_burnin_vs_VPTPNdiff[n][nn] = new TH1D(hname.c_str(),"", 20, 0.85, 1.05);
		pc_burnin_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		pc_burnin_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); pc_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in");

		for(int nnn=0;nnn<6;++nnn){
			if(nnn==0) hname = (string) "lowqe_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==1) hname = (string)"highqe_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==2) hname = (string) "lowgain_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==3) hname = (string)"highgain_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==4) hname = (string) "lowgainqe_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==5) hname = (string)"highgainqe_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			a_splittedburnin_vs_VPTPNdiff[n][nn][nnn] = new TGraphErrors();
			a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->SetName(hname.c_str());
			a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetXaxis()->SetRangeUser(0.85,1.05);
			a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->SetMarkerStyle(4);
			a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetXaxis()->SetTitle("burn-in");
			if(nnn==0) hname = (string) "lowqe_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==1) hname = (string)"highqe_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==2) hname = (string) "lowgain_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==3) hname = (string)"highgain_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==4) hname = (string) "lowgainqe_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			if(nnn==5) hname = (string)"highgainqe_fit_burnin_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn] = new TF1(hname.c_str(), "[0]+[1]*(x-1.0)", 0.85, 1.05);
			afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->SetLineColor(kRed);
		}

		hname = (string)"burnin_corrected_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_burnin_corrected_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_burnin_corrected_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.85,1.05);
		a_burnin_corrected_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in");
		hname = (string)"fit_burnin_corrected_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_burnin_corrected_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-1.0)", 0.75, 1.15);
		afit_burnin_corrected_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

		hname = (string)"gainqe_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_gainqe_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_gainqe_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_gainqe_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.7,1.2);
		a_gainqe_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_gainqe_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("q.e. #cdot gain");
		hname = (string)"fit_gainqe_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_gainqe_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-0.0)", 1., 4.);
		afit_gainqe_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

		hname = (string)"gain_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_gain_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_gain_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_gain_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.7,1.2);
		a_gain_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_gain_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("gain");
		hname = (string)"fit_gain_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_gain_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-0.0)", 6., 16.);
		afit_gain_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

		hname = (string)"qe_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		a_qe_vs_VPTPNdiff[n][nn] = new TGraphErrors();
		a_qe_vs_VPTPNdiff[n][nn]->SetName(hname.c_str());
		a_qe_vs_VPTPNdiff[n][nn]->GetXaxis()->SetRangeUser(0.7,1.2);
		a_qe_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
		a_qe_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}"); a_burnin_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("q.e.");
		hname = (string)"fit_qe_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
		afit_qe_vs_VPTPNdiff[n][nn] = new TF1(hname.c_str(), "[0]+[1]*(x-0.0)", 0.0, 0.4);
		afit_qe_vs_VPTPNdiff[n][nn]->SetLineColor(kRed);

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

		//get burnin
		double burnin = 0; double mustd = 0; double muSIC = 0; int prod = -1;
		double qe = 0; double gain = 0; double gainqe = 0;
		int bin; float eta;
		//get mu_std and eta for each crystals
		if(k==0){//choose xx_EEm_yy histos
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
			if(rotatesc){
				for(unsigned int wz = 0; wz<whichscvec.size(); ++wz){
					int whichsc = whichscvec[wz];
					bool rotatesc_p90 = rotatesc_p90vec[wz];
					bool rotatesc_m90 = rotatesc_m90vec[wz];
					bool rotatesc_180 = rotatesc_180vec[wz];
					//by this we rotate virtually the sc to get correct crystal information connected to laser data.
					int isc(-1), jsc(-1);//center of sc
					for(int wx = i-2; wx<=i+2; ++wx){
					for(int wy = j-2; wy<=j+2; ++wy){
						int wb = SCcenter_EEm->FindBin(wx,wy); int wc = SCcenter_EEm->GetBinContent(wb);
						if(wc==0) continue;
						if(whichsc==wc){ isc = wx; jsc = wy; break; }
					}
					if(isc!=-1&&jsc!=-1) break;
					}
					if(isc>0&&jsc>0 && i>=(isc-2)&&i<=(isc+2) && j>=(jsc-2)&&j<=(jsc+2) ){
						int ii,jj;//intermediate i,j
						if(rotatesc_p90){
							//anticlockwise rotation of pi/2
							jj = jsc + (i-isc);// n = (i-isc) --> j = jsc + n
							ii = isc - (j-jsc);// n = (j-jsc) --> i = isc + n
						} else if(rotatesc_m90){
							//clockwise rotation of pi/2
							jj = jsc - (i-isc);// n = (i-isc) --> j = jsc - n
							ii = isc + (j-jsc);// n = (j-jsc) --> i = isc + n
						} else if(rotatesc_180){
							//rotation by pi
							ii = isc - (i-isc);
							jj = jsc - (j-jsc);
						} else {
							ii = i; jj = j;
						}
						//cout << "i:j:k " << i << ":"<<j<<":"<<k<< " new ii:jj " << ii<<":"<<jj<< " eta " << eta << " SC " << whichsc << " oldbin " << bin;
						bin = EEm_producer->FindBin(ii,jj);
						//cout << " newbin " << bin;
						prod = EEm_producer->GetBinContent(bin);
						if(prod==0 || prod>=3) continue;
						burnin = burnin_EEm    ->GetBinContent(bin);
						mustd  = mu_std_EEm    ->GetBinContent(bin);
						muSIC  = mu_SIC_EEm    ->GetBinContent(bin);
						if(useMuSICinsteafofMuStd) mustd = muSIC;
						qe     = VPT_qe_EEm    ->GetBinContent(bin);
						gainqe = VPT_gainqe_EEm->GetBinContent(bin);
						gain   = VPT_gain_EEm  ->GetBinContent(bin);
						//cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " SC " << whichsc << " and producer " << prod << "(2==chinese,1==russian) has mustd " << mustd << " and burnin " << burnin << endl;
						break;//is rotated crystal, don't check next
					}
					else bin = EEm_producer->FindBin(i,j);
				}
			} else {
			bin = EEm_producer->FindBin(i,j);
			}
			prod = EEm_producer->GetBinContent(bin);
			if(prod==0 || prod>=3) continue;
			burnin = burnin_EEm    ->GetBinContent(bin);
			mustd  = mu_std_EEm    ->GetBinContent(bin);
			muSIC  = mu_SIC_EEm    ->GetBinContent(bin);
			if(useMuSICinsteafofMuStd) mustd = muSIC;
			qe     = VPT_qe_EEm    ->GetBinContent(bin);
			gainqe = VPT_gainqe_EEm->GetBinContent(bin);
			gain   = VPT_gain_EEm  ->GetBinContent(bin);
		}
		else{//EE+
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
			if(rotatesc){
				for(unsigned int wz = 0; wz<whichscvec.size(); ++wz){
					int whichsc = whichscvec[wz];
					bool rotatesc_p90 = rotatesc_p90vec[wz];
					bool rotatesc_m90 = rotatesc_m90vec[wz];
					bool rotatesc_180 = rotatesc_180vec[wz];
					//by this we rotate virtually the sc to get correct crystal information connected to laser data.
					int isc(-1), jsc(-1);//center of sc
					for(int wx = i-2; wx<=i+2; ++wx){
					for(int wy = j-2; wy<=j+2; ++wy){
						int wb = SCcenter_EEp->FindBin(wx,wy); int wc = SCcenter_EEp->GetBinContent(wb);
						if(whichsc==wc){ isc = wx; jsc = wy; break; }
					if(isc!=-1&&jsc!=-1) break;
					}
					if(isc!=-1&&jsc!=-1) break;
					}
					if(isc>0&&jsc>0 && i>=(isc-2)&&i<=(isc+2) && j>=(jsc-2)&&j<=(jsc+2) ){
						int ii,jj;//intermediate i,j
						if(rotatesc_p90){
							//anticlockwise rotation of pi/2
							jj = jsc + (i-isc);// n = (i-isc) --> j = jsc + n
							ii = isc - (j-jsc);// n = (j-jsc) --> i = isc + n
						} else if(rotatesc_m90){
							//clockwise rotation of pi/2
							jj = jsc - (i-isc);// n = (i-isc) --> j = jsc - n
							ii = isc + (j-jsc);// n = (j-jsc) --> i = isc + n
						} else if(rotatesc_180){
							//rotation by pi
							ii = isc - (i-isc);
							jj = jsc - (j-jsc);
						} else {
							ii = i; jj = j;
						}
						//cout << "i:j:k " << i << ":"<<j<<":"<<k<< " new ii:jj " << ii<<":"<<jj<< " eta " << eta << " SC " << whichsc << " oldbin " << bin;
						bin = EEp_producer->FindBin(ii,jj);
						//cout << " newbin " << bin << endl;
						break;
					}
					else bin = EEp_producer->FindBin(i,j);

				}
			} else {
			bin = EEp_producer->FindBin(i,j);
			}
			prod = EEp_producer->GetBinContent(bin);
			if(prod==0 || prod>=3) continue;
			burnin = burnin_EEp    ->GetBinContent(bin);
			mustd  = mu_std_EEp    ->GetBinContent(bin);
			muSIC  = mu_SIC_EEp    ->GetBinContent(bin);
			if(useMuSICinsteafofMuStd) mustd = muSIC;
			qe     = VPT_qe_EEp    ->GetBinContent(bin);
			gainqe = VPT_gainqe_EEp->GetBinContent(bin);
			gain   = VPT_gain_EEp  ->GetBinContent(bin);
		}
		
		double vptdiff = 1.-y;
		int etaindex = -1;
		int highqe(0), highgainqe(0), highgain(0);
		if(fabs(eta)<1.6&&fabs(eta)>=1.4) etaindex = 0;
		else if(fabs(eta)<( 1.8))         etaindex = 1;
		else if(fabs(eta)<( 2.0))         etaindex = 2;
		else if(fabs(eta)<( 2.2))         etaindex = 3;
		else if(fabs(eta)<( 2.4))         etaindex = 4;
		else if(fabs(eta)<( 2.6))         etaindex = 5;
		else if(fabs(eta)<( 2.8))         etaindex = 6;
		else if(fabs(eta)<=(3.0))         etaindex = 7;
		if(etaindex==0)             {if(qe>0.25 ) highqe=1; if(gain>11.5) highgain=1; if(gainqe>2.7 ) highgainqe=1; }
		if(etaindex==1)             {if(qe>0.235) highqe=1; if(gain>10.5) highgain=1; if(gainqe>2.3 ) highgainqe=1; }
		if(etaindex==2)             {if(qe>0.22 ) highqe=1; if(gain>9.0 ) highgain=1; if(gainqe>2.0 ) highgainqe=1; }
		if(etaindex==3)             {if(qe>0.203) highqe=1; if(gain>9.0 ) highgain=1; if(gainqe>1.85) highgainqe=1; }
		if(etaindex==4||etaindex==5){if(qe>0.195) highqe=1; if(gain>8.5 ) highgain=1; if(gainqe>1.7 ) highgainqe=1; }
		if(etaindex==6||etaindex==7){if(qe>0.195) highqe=1; if(gain>8.5 ) highgain=1; if(gainqe>1.6 ) highgainqe=1; }
		if(etaindex>=0 && etaindex <=7 && burnin!=0){
			if(!(prod==2&&cleanviaSIC&&cleanviaSICo&&muSIC==0)){
			//if(i>95&&k==1) cout<<"ijk "<<i<<" "<<j<<" "<<k<<" has eta "<<etaindex<<" prod-1 "<<prod-1<<" burnin "<<burnin<<" vptdiff "<<vptdiff<<endl;
				int graphnbins = a_burnin_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_burnin_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, burnin, vptdiff);
				if(!fitwithouterror) a_burnin_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
				int nnn;
				if(highqe    ==1) nnn = 1; else nnn = 0;
				graphnbins = a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->GetN();;
				a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPoint(graphnbins, burnin, vptdiff);
				if(!fitwithouterror) a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPointError(graphnbins, 0., ey);
				if(highgain  ==1) nnn = 3; else nnn = 2;
				graphnbins = a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->GetN();;
				a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPoint(graphnbins, burnin, vptdiff);
				if(!fitwithouterror) a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPointError(graphnbins, 0., ey);
				if(highgainqe==1) nnn = 5; else nnn = 4;
				graphnbins = a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->GetN();;
				a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPoint(graphnbins, burnin, vptdiff);
				if(!fitwithouterror) a_splittedburnin_vs_VPTPNdiff[etaindex][prod-1][nnn]->SetPointError(graphnbins, 0., ey);
			}
		}
		if(etaindex>=0 && etaindex <=7 && mustd!=0){
			if(!(prod==2&&cleanviaSIC&&muSIC==0)){
				int graphnbins = a_mustd_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_mustd_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, mustd, vptdiff);
				if(!fitwithouterror) a_mustd_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
			}
		}
		if(etaindex>=0 && etaindex <=7 && gainqe!=0){
			if(!(prod==2&&cleanviaSIC&&cleanviaSICo&&muSIC==0)){
				int graphnbins = a_gainqe_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_gainqe_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, gainqe, vptdiff);
				if(!fitwithouterror) a_gainqe_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
			}
		}
		if(etaindex>=0 && etaindex <=7 && gain!=0){
			if(!(prod==2&&cleanviaSIC&&cleanviaSICo&&muSIC==0)){
				int graphnbins = a_gain_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_gain_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, gain, vptdiff);
				if(!fitwithouterror) a_gain_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
			}
		}
		if(etaindex>=0 && etaindex <=7 && gainqe!=0){
			if(!(prod==2&&cleanviaSIC&&cleanviaSICo&&muSIC==0)){
				int graphnbins = a_qe_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_qe_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, qe, vptdiff);
				if(!fitwithouterror) a_qe_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
			}
		}

	}}}//ijk
	//now run over all endcaps for one day
	//now make fit
	
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		double xmin(99.),xmax(-99.);
		double correff(0.), correfferr(0.);
		double xav(0.), yav(0.);
		double cov(0.), xerr(0.), yerr(0.);
		//fit now //fill correlation coefficient
		for(int nnn = 0; nnn<a_burnin_vs_VPTPNdiff[n][nn]->GetN(); ++nnn){
			double xn, yn;
			a_burnin_vs_VPTPNdiff[n][nn]->GetPoint(nnn,xn,yn);
			p_burnin_vs_VPTPNdiff[n][nn]->Fill(xn, yn);
			xav += xn; yav += yn;
			if(xmin>xn) xmin = xn;
			if(xmax<xn) xmax = xn;
		}
		xav = xav / TMath::Max(a_burnin_vs_VPTPNdiff[n][nn]->GetN(),1);
		yav = yav / TMath::Max(a_burnin_vs_VPTPNdiff[n][nn]->GetN(),1);
		if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3){
			for(int nnn = 0; nnn<a_burnin_vs_VPTPNdiff[n][nn]->GetN(); ++nnn){
				double xn, yn;
				a_burnin_vs_VPTPNdiff[n][nn]->GetPoint(nnn,xn,yn);
				cov += (xn-xav)*(yn-yav);
				xerr += (xn-xav)*(xn-xav);
				yerr += (yn-yav)*(yn-yav);
			}
			if(xerr*yerr!=0) {
				correff = cov/sqrt(xerr*yerr);
				correfferr = (1.-correff*correff)/sqrt(a_burnin_vs_VPTPNdiff[n][nn]->GetN()-3.);
				int corrcoeffbin = g_burnin_corrcoeff[n][nn]->GetN();
				g_burnin_corrcoeff[n][nn]->SetPoint(     corrcoeffbin, d, correff   );
				if(!dontploterrer) g_burnin_corrcoeff[n][nn]->SetPointError(corrcoeffbin, 0, correfferr);
			}
		}
		//continue //fit now
		if(fitprofile) {xmin = 99.; xmax = -99.;}
		int pc_burnin_vs_VPTPNdifffilled = 0;
		for(int nnn = 1; nnn<=p_burnin_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//bin no. starts from 1
			if(p_burnin_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn)<4) continue;//xxyyzz, 2 was 4
			if(p_burnin_vs_VPTPNdiff[n][nn]->GetBinContent(nnn)==0) continue;
			double x = p_burnin_vs_VPTPNdiff[n][nn]->GetBinCenter(nnn);
			double xw = p_burnin_vs_VPTPNdiff[n][nn]->GetBinWidth(nnn);
			if(fitprofile){ if(xmin>(x-xw/2.)) xmin = x - xw/2.; if(xmax<(x+xw/2.)) xmax = x + xw/2.; }
	//		pc_burnin_vs_VPTPNdiff[n][nn]->SetBinEntries(nnn,p_burnin_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn));
			pc_burnin_vs_VPTPNdiff[n][nn]->SetBinContent(nnn,p_burnin_vs_VPTPNdiff[n][nn]->GetBinContent(nnn));
			pc_burnin_vs_VPTPNdiff[n][nn]->SetBinError(  nnn,p_burnin_vs_VPTPNdiff[n][nn]->GetBinError(  nnn));
			++pc_burnin_vs_VPTPNdifffilled;
		}
	/*	for(int nnn = 1; nnn<=pc_burnin_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//give empty bins crazy uncertainty
			if(pc_burnin_vs_VPTPNdiff[n][nn]->GetBinContent(nnn)==0||pc_burnin_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn)<3){
			pc_burnin_vs_VPTPNdiff[n][nn]->SetBinEntries(nnn,1);
			pc_burnin_vs_VPTPNdiff[n][nn]->SetBinContent(nnn,-1.);
			pc_burnin_vs_VPTPNdiff[n][nn]->SetBinError(  nnn,999E6);
			}
		}
	*/
	//	if(int((d-1293836400)/(86400))>50&&int((d-1293836400)/(86400))<75) cout << "pc_burnin_vs_VPTPNdifffilled " << pc_burnin_vs_VPTPNdifffilled << endl;
		if(xmin!=99.&&xmax!=-99.) afit_burnin_vs_VPTPNdiff[n][nn]->SetRange(xmin-0.01,xmax+0.01);
		afit_burnin_vs_VPTPNdiff[n][nn]->SetRange(xmin-0.01,xmax+0.01);
		afit_burnin_vs_VPTPNdiff[n][nn]->SetRange(xmin-0.01,xmax+0.01);
		if(fitprofile){
			if(pc_burnin_vs_VPTPNdifffilled>3){//xxyyzz, 1 was 3
			//afit_burnin_vs_VPTPNdiff[n][nn]->SetParameter(0,a_burnin_vs_VPTPNdiff[n][nn]->GetMean(2));
			(pc_burnin_vs_VPTPNdiff[n][nn])->Fit(afit_burnin_vs_VPTPNdiff[n][nn], "RQ");//quiet
				/*if(int((d-1293836400)/(86400))>50&&int((d-1293836400)/(86400))<75){
					for(int nnn = 1; nnn<=pc_burnin_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//bin no. starts from 1
					cout << "cbin " << nnn << " entries " << pc_burnin_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn) << "  content " << pc_burnin_vs_VPTPNdiff[n][nn]->GetBinContent(nnn) << " +/- " << pc_burnin_vs_VPTPNdiff[n][nn]->GetBinError(nnn) << endl;
					cout << "obin " << nnn << " entries " << p_burnin_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn) << "  content " << p_burnin_vs_VPTPNdiff[n][nn]->GetBinContent(nnn) << " +/- " << p_burnin_vs_VPTPNdiff[n][nn]->GetBinError(nnn) << endl;
					}
				}*/
			}
		} else {
			if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>6){
			(a_burnin_vs_VPTPNdiff[n][nn])->Fit(afit_burnin_vs_VPTPNdiff[n][nn], "RQ");//quiet
			}
		}
		/*
		cout << "Fit " << a_burnin_vs_VPTPNdiff[n][nn]->GetName() << endl;
		cout << "Fit constant = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0) << endl;
		cout << "Fit slope    = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
		cout << "Goodness of fit: chi^2/NDF = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetChisquare() << "/" << afit_burnin_vs_VPTPNdiff[n][nn]->GetNDF() << " = " << (afit_burnin_vs_VPTPNdiff[n][nn]->GetChisquare())/(afit_burnin_vs_VPTPNdiff[n][nn]->GetNDF()) << endl;
		cout << "Fit probability: " << afit_burnin_vs_VPTPNdiff[n][nn]->GetProb() << endl;
		*/
		double burnin_slope, burnin_slopeerr, burnin_incep, burnin_inceperr;
		if((fitprofile&&pc_burnin_vs_VPTPNdifffilled>3)||(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>6)){//xxyyzz 1 was 3
			burnin_slope     = afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1);
			burnin_slopeerr  = afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1);
			burnin_incep     = afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0);
			burnin_inceperr  = afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0);
			int burnin_slopebin     = g_burnin_slopes[n][nn]    ->GetN();
			int burnin_incepbin     = g_burnin_onevalues[n][nn] ->GetN();
			//check
			if(burnin_slopebin!=burnin_incepbin) cout << __LINE__ << "  slopebin " << burnin_slopebin << "  incepbin " << burnin_incepbin << endl;
			g_burnin_slopes[n][nn]    ->SetPoint(     burnin_slopebin, d, burnin_slope   );
			if(!dontploterrer) g_burnin_slopes[n][nn]    ->SetPointError(burnin_slopebin, 0, burnin_slopeerr);
			g_burnin_onevalues[n][nn] ->SetPoint(     burnin_incepbin, d, burnin_incep   );
			if(!dontploterrer) g_burnin_onevalues[n][nn] ->SetPointError(burnin_incepbin, 0, burnin_inceperr);
			if(printtime){
			//	cout << __LINE__ << " nn " << nn << " n " << n << endl;
			//	cout << "par[0] = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0) << ", par[1] = " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1) << " +/- " << afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1) << endl;
				h_burnin_eq1[nn]  ->SetBinContent(n+1, burnin_incep   );
				h_burnin_eq1[nn]  ->SetBinError(  n+1, burnin_inceperr);
				h_burnin_slope[nn]->SetBinContent(n+1, burnin_slope   );
				h_burnin_slope[nn]->SetBinError(  n+1, burnin_slopeerr);
			}
			for(int nnn = 0; nnn<6; ++nnn){
				if(a_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetN()<=4) continue;
				(a_splittedburnin_vs_VPTPNdiff[n][nn][nnn])->Fit(afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn], "RQ");//quiet
				double splittedburnin_slope     = afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetParameter(1);
				double splittedburnin_slopeerr  = afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetParError(1);
				double splittedburnin_incep     = afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetParameter(0);
				double splittedburnin_inceperr  = afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn]->GetParError(0);
				int splittedburnin_slopebin     = g_splittedburnin_slopes[n][nn][nnn]    ->GetN();
				int splittedburnin_incepbin     = g_splittedburnin_onevalues[n][nn][nnn] ->GetN();
				g_splittedburnin_slopes[n][nn][nnn]    ->SetPoint(     splittedburnin_slopebin, d, splittedburnin_slope   );
				if(!dontploterrer) g_splittedburnin_slopes[n][nn][nnn]    ->SetPointError(splittedburnin_slopebin, 0, splittedburnin_slopeerr);
				g_splittedburnin_onevalues[n][nn][nnn] ->SetPoint(     splittedburnin_incepbin, d, splittedburnin_incep   );
				if(!dontploterrer) g_splittedburnin_onevalues[n][nn][nnn] ->SetPointError(splittedburnin_incepbin, 0, splittedburnin_inceperr);
			}
		}

		xmin=99.; xmax = -99.;
		correff=0.; correfferr=0.;
		xav=0.; yav=0.;
		cov=0.; xerr=0.; yerr=0.;
		for(int nnn = 0; nnn<a_mustd_vs_VPTPNdiff[n][nn]->GetN(); ++nnn){
			double xn, yn;
			a_mustd_vs_VPTPNdiff[n][nn]->GetPoint(nnn,xn,yn);
			p_mustd_vs_VPTPNdiff[n][nn]->Fill(xn, yn);
			xav += xn; yav += yn;
			if(xmin>xn) xmin = xn;
			if(xmax<xn) xmax = xn;
		}
		xav = xav / TMath::Max(a_mustd_vs_VPTPNdiff[n][nn]->GetN(),1);
		yav = yav / TMath::Max(a_mustd_vs_VPTPNdiff[n][nn]->GetN(),1);
		if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3){
			for(int nnn = 0; nnn<a_mustd_vs_VPTPNdiff[n][nn]->GetN(); ++nnn){
				double xn, yn;
				a_mustd_vs_VPTPNdiff[n][nn]->GetPoint(nnn,xn,yn);
				cov += (xn-xav)*(yn-yav);
				xerr += (xn-xav)*(xn-xav);
				yerr += (yn-yav)*(yn-yav);
			}
			if(xerr*yerr!=0) {
				correff = cov/sqrt(xerr*yerr);
				correfferr = (1.-correff*correff)/sqrt(a_burnin_vs_VPTPNdiff[n][nn]->GetN()-3.);
				int corrcoeffbin = g_mustd_corrcoeff[n][nn]->GetN();
				g_mustd_corrcoeff[n][nn]->SetPoint(     corrcoeffbin, d, correff   );
				if(!dontploterrer) g_mustd_corrcoeff[n][nn]->SetPointError(corrcoeffbin, 0, correfferr);
			}
		}
		//continue //fit now
		if(fitprofile) {xmin = 99.; xmax = -99.;}
		int pc_mustd_vs_VPTPNdifffilled = 0;
		for(int nnn = 1; nnn<=p_mustd_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//bin no. starts from 1
			if(p_mustd_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn)<3) continue;//xxyyzz 2 was 3
			if(p_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(nnn)==0) continue;
			double x = p_mustd_vs_VPTPNdiff[n][nn]->GetBinCenter(nnn);
			double xw = p_mustd_vs_VPTPNdiff[n][nn]->GetBinWidth(nnn);
			if(fitprofile){ if(xmin>(x-xw/2.)) xmin = x - xw/2.; if(xmax<(x+xw/2.)) xmax = x + xw/2.; }
	//		pc_mustd_vs_VPTPNdiff[n][nn]->SetBinEntries(nnn,p_mustd_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn));
			pc_mustd_vs_VPTPNdiff[n][nn]->SetBinContent(nnn,p_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(nnn));
			pc_mustd_vs_VPTPNdiff[n][nn]->SetBinError(  nnn,p_mustd_vs_VPTPNdiff[n][nn]->GetBinError(  nnn));
			++pc_mustd_vs_VPTPNdifffilled;
		}
	/*	for(int nnn = 1; nnn<=pc_mustd_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//give empty bins crazy uncertainties	
			if(pc_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(nnn)==0||pc_mustd_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn)<3){
			pc_mustd_vs_VPTPNdiff[n][nn]->SetBinEntries(nnn,1);
			pc_mustd_vs_VPTPNdiff[n][nn]->SetBinContent(nnn,-1.);
			pc_mustd_vs_VPTPNdiff[n][nn]->SetBinError(  nnn,999E6);
			}
		}
	*/
	//	if(int((d-1293836400)/(86400))>50&&int((d-1293836400)/(86400))<75) cout << "pc_mustd_vs_VPTPNdifffilled " << pc_mustd_vs_VPTPNdifffilled << endl;
		if(xmax!=-99.&&xmin!=99.) afit_mustd_vs_VPTPNdiff[n][nn]->SetRange(xmin-0.01,xmax+0.01);
		else if(xmin!=99.)        afit_mustd_vs_VPTPNdiff[n][nn]->SetRange(xmin-0.01,2.0);
		else if(xmax!=-99.)       afit_mustd_vs_VPTPNdiff[n][nn]->SetRange(0.,xmax+0.01);
		else                      afit_mustd_vs_VPTPNdiff[n][nn]->SetRange(0.,2.0);
		if(fitprofile){
			if(pc_mustd_vs_VPTPNdifffilled>3){//xxyyzz 1 was 3
		//	afit_mustd_vs_VPTPNdiff[n][nn]->SetParameter(0,a_mustd_vs_VPTPNdiff[n][nn]->GetMean(2));
			(pc_mustd_vs_VPTPNdiff[n][nn])->Fit(afit_mustd_vs_VPTPNdiff[n][nn], "RQ");//quiet
			/*	if(int((d-1293836400)/(86400))>50&&int((d-1293836400)/(86400))<75){
					for(int nnn = 1; nnn<=pc_mustd_vs_VPTPNdiff[n][nn]->GetNbinsX(); ++nnn){//bin no. starts from 1
					cout << "cbin " << nnn << " entries " << pc_mustd_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn) << "  content " << pc_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(nnn) << " +/- " << pc_mustd_vs_VPTPNdiff[n][nn]->GetBinError(nnn) << endl;
					cout << "obin " << nnn << " entries " << p_mustd_vs_VPTPNdiff[n][nn]->GetBinEntries(nnn) << "  content " << p_mustd_vs_VPTPNdiff[n][nn]->GetBinContent(nnn) << " +/- " << p_mustd_vs_VPTPNdiff[n][nn]->GetBinError(nnn) << endl;
					}
				}*/
			}
		} else {
			if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>4){
			(a_mustd_vs_VPTPNdiff[n][nn])->Fit(afit_mustd_vs_VPTPNdiff[n][nn], "RQ");//quiet
			}
		}

		double mustd_slope, mustd_slopeerr, mustd_incep, mustd_inceperr;
		if((fitprofile&&pc_mustd_vs_VPTPNdifffilled>3)||(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>4)){//xxyyzz 1 was 3
			mustd_slope     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1);
			mustd_slopeerr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1);
			mustd_incep     = afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0);
			mustd_inceperr  = afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0);
			int mustd_slopebin     = g_mustd_slopes[n][nn]    ->GetN();
			int mustd_incepbin     = g_mustd_intercepts[n][nn] ->GetN();
			g_mustd_slopes[n][nn]    ->SetPoint(     mustd_slopebin, d, mustd_slope   );
			if(!dontploterrer) g_mustd_slopes[n][nn]    ->SetPointError(mustd_slopebin, 0, mustd_slopeerr);
			g_mustd_intercepts[n][nn]->SetPoint(      mustd_incepbin, d, mustd_incep   );
			if(!dontploterrer) g_mustd_intercepts[n][nn]->SetPointError( mustd_incepbin, 0, mustd_inceperr);
			if(printtime){
				h_mustd_intercept[nn]->SetBinContent(n+1, mustd_incep   );
				h_mustd_intercept[nn]->SetBinError(  n+1, mustd_inceperr);
				h_mustd_slope[nn]    ->SetBinContent(n+1, mustd_slope   );
				h_mustd_slope[nn]    ->SetBinError(  n+1, mustd_slopeerr);
			}
			if((fitprofile&&pc_burnin_vs_VPTPNdifffilled>3)||(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>6)){
				int burninslope_mustdintercept_bin     = g_burnin_slope_vs_mustd_intercept[n][nn]    ->GetN();
				g_burnin_slope_vs_mustd_intercept[n][nn]->SetPoint(     burninslope_mustdintercept_bin, burnin_slope,    mustd_incep   );
				if(!dontploterrer) g_burnin_slope_vs_mustd_intercept[n][nn]->SetPointError(burninslope_mustdintercept_bin, burnin_slopeerr, mustd_inceperr);
				g_burnin_onevalue_vs_mustd_slope[n][nn] ->SetPoint(     burninslope_mustdintercept_bin, burnin_incep,    mustd_slope   );
				if(!dontploterrer) g_burnin_onevalue_vs_mustd_slope[n][nn] ->SetPointError(burninslope_mustdintercept_bin, burnin_inceperr, mustd_slopeerr);
			}
		}
	//	if(int((d-1293836400)/(86400))>50&&int((d-1293836400)/(86400))<75) cout << "after everything " << endl;
		if(a_gainqe_vs_VPTPNdiff[n][nn]->GetN()>4){
		(a_gainqe_vs_VPTPNdiff[n][nn])->Fit(afit_gainqe_vs_VPTPNdiff[n][nn], "RQ");//quiet
		double gainqe_slope     = afit_gainqe_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double gainqe_slopeerr  = afit_gainqe_vs_VPTPNdiff[n][nn]->GetParError(1);
		double gainqe_incep     = afit_gainqe_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double gainqe_inceperr  = afit_gainqe_vs_VPTPNdiff[n][nn]->GetParError(0);
		int gainqe_slopebin     = g_gainqe_slopes[n][nn]    ->GetN();
		int gainqe_incepbin     = g_gainqe_intercepts[n][nn] ->GetN();
		g_gainqe_slopes[n][nn]    ->SetPoint(     gainqe_slopebin, d, gainqe_slope   );
		if(!dontploterrer) g_gainqe_slopes[n][nn]    ->SetPointError(gainqe_slopebin, 0, gainqe_slopeerr);
		g_gainqe_intercepts[n][nn]->SetPoint(     gainqe_incepbin, d, gainqe_incep   );
		if(!dontploterrer) g_gainqe_intercepts[n][nn]->SetPointError(gainqe_incepbin, 0, gainqe_inceperr);
		}
		if(a_gain_vs_VPTPNdiff[n][nn]->GetN()>4){
		(a_gain_vs_VPTPNdiff[n][nn])->Fit(afit_gain_vs_VPTPNdiff[n][nn], "RQ");//quiet
		double gain_slope     = afit_gain_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double gain_slopeerr  = afit_gain_vs_VPTPNdiff[n][nn]->GetParError(1);
		double gain_incep     = afit_gain_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double gain_inceperr  = afit_gain_vs_VPTPNdiff[n][nn]->GetParError(0);
		int gain_slopebin     = g_gain_slopes[n][nn]    ->GetN();
		int gain_incepbin     = g_gain_intercepts[n][nn] ->GetN();
		g_gain_slopes[n][nn]    ->SetPoint(     gain_slopebin, d, gain_slope   );
		if(!dontploterrer) g_gain_slopes[n][nn]    ->SetPointError(gain_slopebin, 0, gain_slopeerr);
		g_gain_intercepts[n][nn]->SetPoint(     gain_incepbin, d, gain_incep   );
		if(!dontploterrer) g_gain_intercepts[n][nn]->SetPointError(gain_incepbin, 0, gain_inceperr);
		}
		if(a_qe_vs_VPTPNdiff[n][nn]->GetN()>4){
		(a_qe_vs_VPTPNdiff[n][nn])->Fit(afit_qe_vs_VPTPNdiff[n][nn], "RQ");//quiet
		double qe_slope     = afit_qe_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double qe_slopeerr  = afit_qe_vs_VPTPNdiff[n][nn]->GetParError(1);
		double qe_incep     = afit_qe_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double qe_inceperr  = afit_qe_vs_VPTPNdiff[n][nn]->GetParError(0);
		int qe_slopebin     = g_qe_slopes[n][nn]    ->GetN();
		int qe_incepbin     = g_qe_intercepts[n][nn] ->GetN();
		g_qe_slopes[n][nn]    ->SetPoint(     qe_slopebin, d, qe_slope   );
		if(!dontploterrer) g_qe_slopes[n][nn]    ->SetPointError(qe_slopebin, 0, qe_slopeerr);
		g_qe_intercepts[n][nn]->SetPoint(     qe_incepbin, d, qe_incep   );
		if(!dontploterrer) g_qe_intercepts[n][nn]->SetPointError(qe_incepbin, 0, qe_inceperr);
		}
		if(printtime){
			TString outputfile;
			TString outdirtemp;
			string hname;
			//TH1F *haxis;
			if(makeburninplots){
				if(fitlinepluserr){
					fitup   = new TF1("fitup",   "[0]+[1]*(x-1.0)", 0.8, 1.1);
					fitdown = new TF1("fitdown", "[0]+[1]*(x-1.0)", 0.8, 1.1);
					fitup  ->SetLineColor(kRed); fitup  ->SetLineStyle(3); fitup  ->SetLineWidth(2);
					fitdown->SetLineColor(kRed); fitdown->SetLineStyle(3); fitdown->SetLineWidth(2);
					fitup  ->SetParameter(0, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0)+afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0));
					fitdown->SetParameter(0, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0)-afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(0));
					fitup  ->SetParameter(1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1)-afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1));
					fitdown->SetParameter(1, afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1)+afit_burnin_vs_VPTPNdiff[n][nn]->GetParError(1));
				}
				TH1F *haxis = new TH1F("haxis","",165,0.85,1.05);
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum(0.7);
				haxis->GetXaxis()->SetTitle("burn-in");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetDirectory(0);                haxis->SetStats(0);
				haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
				haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
				haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
				haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
			/*	a_burnin_vs_VPTPNdiff[n][nn]->Draw("AP");
				haxis = a_burnin_vs_VPTPNdiff[n][nn]->GetHistogram();
				haxis->SetName(a_burnin_vs_VPTPNdiff[n][nn]->GetName());
				haxis->GetXaxis()->SetRangeUser(0.85,1.05);
				haxis->GetXaxis()->SetLimits(0.85,1.05);
				haxis->GetXaxis()->SetTitle("burn-in");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetXTitle("burn-in");
				haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum( 0.7);
			*/	haxis->Draw("hist");
				a_burnin_vs_VPTPNdiff[n][nn]->Draw("P");
				afit_burnin_vs_VPTPNdiff[n][nn]->Draw("same");
				fburnintemp->SetParameter(0,afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0));
				fburnintemp->SetParameter(1,afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1));
				fburnintemp->Draw("same");
				if(fitlinepluserr){
					fitup  ->Draw("same");
					fitdown->Draw("same");
				}
				outputfile = a_burnin_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirburnin+daydir;
				//col->Update();
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				col->cd();
				haxis->Draw("hist");
				pc_burnin_vs_VPTPNdiff[n][nn]->Draw("same");
				afit_burnin_vs_VPTPNdiff[n][nn]->Draw("same");
				fburnintemp->SetParameter(0,afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(0));
				fburnintemp->SetParameter(1,afit_burnin_vs_VPTPNdiff[n][nn]->GetParameter(1));
				fburnintemp->Draw("same");
				if(fitlinepluserr){
					fitup  ->Draw("same");
					fitdown->Draw("same");
				}
				outputfile = p_burnin_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirburnin+daydir;
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_burnin_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				delete haxis;
				if(fitlinepluserr){
					delete fitup;
					delete fitdown;
				}
			}
			if(makemustdplots){
				if(fitlinepluserr){
					fitup   = new TF1("fitup",   "[0]+[1]*(x-0.0)", 0.0, 2.0);
					fitdown = new TF1("fitdown", "[0]+[1]*(x-0.0)", 0.0, 2.0);
					fitup  ->SetLineColor(kRed); fitup  ->SetLineStyle(3); fitup  ->SetLineWidth(2);
					fitdown->SetLineColor(kRed); fitdown->SetLineStyle(3); fitdown->SetLineWidth(2);
					fitup  ->SetParameter(0, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)+afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
					fitdown->SetParameter(0, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0)-afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(0));
					fitup  ->SetParameter(1, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)+afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
					fitdown->SetParameter(1, afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1)-afit_mustd_vs_VPTPNdiff[n][nn]->GetParError(1));
				}
				TH1F *haxis = new TH1F("haxis","",100,0.,2.);
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum(1.1);
				haxis->GetXaxis()->SetTitle("#mu_{std}");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetDirectory(0);                haxis->SetStats(0);
				haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
				haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
				haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
				haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
			/*	haxis = a_mustd_vs_VPTPNdiff[n][nn]->GetHistogram();
				haxis->SetName(a_mustd_vs_VPTPNdiff[n][nn]->GetName());
				haxis->GetXaxis()->SetRangeUser(0.,2.0);
				haxis->GetXaxis()->SetLimits(0.,2.0);
				haxis->GetXaxis()->SetTitle("#mu_{std}");
				haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetXTitle("#mu_{std}");
				haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum( 1.1);
			*/	haxis->Draw("hist");
				a_mustd_vs_VPTPNdiff[n][nn]->Draw("P");
				afit_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
				fmustdtemp->SetParameter(0,afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0));
				fmustdtemp->SetParameter(1,afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1));
				fmustdtemp->Draw("same");
				if(fitlinepluserr){
					fitup  ->Draw("same");
					fitdown->Draw("same");
				}
				outputfile = a_mustd_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirmustd+daydir;
				//col->Update();
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				col->cd();
				haxis->Draw("hist");
				pc_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
				afit_mustd_vs_VPTPNdiff[n][nn]->Draw("same");
				fmustdtemp->SetParameter(0,afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(0));
				fmustdtemp->SetParameter(1,afit_mustd_vs_VPTPNdiff[n][nn]->GetParameter(1));
				fmustdtemp->Draw("same");
				if(fitlinepluserr){
					fitup  ->Draw("same");
					fitdown->Draw("same");
				}
				outputfile = p_mustd_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirmustd+daydir;
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_mustd_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				delete haxis;
				if(fitlinepluserr){
					delete fitup;
					delete fitdown;
				}
			}
			if(makegainqeplots){
				TH1F *haxis = new TH1F("haxis","",150,1.0,4.0);
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum(1.1);
				haxis->GetXaxis()->SetTitle("q.e. #cdot gain");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetDirectory(0);                haxis->SetStats(0);
				haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
				haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
				haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
				haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
			/*	haxis = a_gainqe_vs_VPTPNdiff[n][nn]->GetHistogram();
				haxis->SetName(a_gainqe_vs_VPTPNdiff[n][nn]->GetName());
				haxis->GetXaxis()->SetRangeUser(1.,4.0);
				haxis->GetXaxis()->SetLimits(1.,4.0);
				haxis->GetXaxis()->SetTitle("q.e. * gain");
				haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetXTitle("q.e. * gain");
				haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum( 1.1);
			*/	haxis->Draw("hist");
				a_gainqe_vs_VPTPNdiff[n][nn]->Draw("P");
				outputfile = a_gainqe_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirgainqe+daydir;
				//col->Update();
				if(a_gainqe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_gainqe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_gainqe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				delete haxis;
			}
			if(makegainplots){
				TH1F *haxis = new TH1F("haxis","",100,6.0,16.0);
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum(1.1);
				haxis->GetXaxis()->SetTitle("gain");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetDirectory(0);                haxis->SetStats(0);
				haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
				haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
				haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
				haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
			/*	haxis = a_gain_vs_VPTPNdiff[n][nn]->GetHistogram();
				haxis->SetName(a_gain_vs_VPTPNdiff[n][nn]->GetName());
				haxis->GetXaxis()->SetRangeUser(6.,16.0);
				haxis->GetXaxis()->SetLimits(6.,16.0);
				haxis->GetXaxis()->SetTitle("gain");
				haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetXTitle("gain");
				haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum( 1.1);
			*/	haxis->Draw("hist");
				a_gain_vs_VPTPNdiff[n][nn]->Draw("P");
				outputfile = a_gain_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirgain+daydir;
				//col->Update();
				if(a_gain_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_gain_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_gain_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				delete haxis;
			}
			if(makeqeplots){
				TH1F *haxis = new TH1F("haxis","",100,0.0,0.4);
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum(1.1);
				haxis->GetXaxis()->SetTitle("q.e");
				haxis->GetYaxis()->SetTitle("#Delta#frac{VPT}{PN}");
				haxis->SetDirectory(0);                haxis->SetStats(0);
				haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
				haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
				haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
				haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
				haxis = a_qe_vs_VPTPNdiff[n][nn]->GetHistogram();
			/*	haxis->SetName(a_qe_vs_VPTPNdiff[n][nn]->GetName());
				haxis->GetXaxis()->SetRangeUser(0.,0.4);
				haxis->GetXaxis()->SetLimits(0.,0.4);
				haxis->GetXaxis()->SetTitle("q.e.");
				haxis->GetYaxis()->SetTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetXTitle("q.e.");
				haxis->SetYTitle("#frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
				haxis->SetMinimum(-0.1);
				haxis->SetMaximum( 1.1);
			*/	haxis->Draw("hist");
				a_qe_vs_VPTPNdiff[n][nn]->Draw("P");
				outputfile = a_qe_vs_VPTPNdiff[n][nn]->GetName();
				outdirtemp = outdirqe+daydir;
				//col->Update();
				if(a_qe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
				if(a_qe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
				if(a_qe_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
				col->Clear();
				delete haxis;
			}
		}

	}}
//start BURNIN corrected
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
				ex = g_vptpn_las[i][j][k]->GetErrorX(n);	ey = g_vptpn_las[i][j][k]->GetErrorY(n);
				x = x0; y = y0; runfurther = true; }
			if(runfurther) break;
		} if(!runfurther) continue;

		double burnin = 0; double mustd = 0; double muSIC = 0; int prod = -1;
		int bin; float eta;
		if(k==0){//choose xx_EEm_yy histos
			bin = eta_EEm->FindBin(i,j);
			eta = eta_EEm->GetBinContent(bin);
			bin = EEm_producer->FindBin(i,j);
			prod = EEm_producer->GetBinContent(bin);
			if(prod==0 || prod>=3) continue;
			burnin = burnin_EEm->GetBinContent(bin);
			mustd  = mu_std_EEm ->GetBinContent(bin);
			muSIC  = mu_SIC_EEm ->GetBinContent(bin);
		}
		else{
			bin = eta_EEp->FindBin(i,j);
			eta = eta_EEp->GetBinContent(bin);
			bin = EEp_producer->FindBin(i,j);
			prod = EEp_producer->GetBinContent(bin);
			if(prod==0 || prod>=3) continue;
			burnin = burnin_EEp->GetBinContent(bin);
			mustd  = mu_std_EEp ->GetBinContent(bin);
			muSIC  = mu_SIC_EEp ->GetBinContent(bin);
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
		if(etaindex>=0 && etaindex <=7 && burnin!=0 && mustd!=0){
			if(!(prod==2&&cleanviaSIC&&muSIC==0)){
				double mustdintercept = afit_mustd_vs_VPTPNdiff[etaindex][prod-1]->GetParameter(0);
				double mustdvalue = afit_mustd_vs_VPTPNdiff[etaindex][prod-1]->Eval(mustd);
				//if(d>1315013137 && d<1315513137) cout << "VPTdiff " << vptdiff << "  VPTmustd_fitval(mustd=" << mustd<<") " << mustdvalue << " with intercept " << mustdintercept << endl;
				mustdvalue = mustdvalue-mustdintercept;//intercept subtracted --> only slope
				int graphnbins = a_burnin_corrected_vs_VPTPNdiff[etaindex][prod-1]->GetN();;
				a_burnin_corrected_vs_VPTPNdiff[etaindex][prod-1]->SetPoint(graphnbins, burnin, vptdiff-mustdvalue);
				a_burnin_corrected_vs_VPTPNdiff[etaindex][prod-1]->SetPointError(graphnbins, 0., ey);
			}
		}
	}}}//ijk
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()<=4) continue;
		(a_burnin_corrected_vs_VPTPNdiff[n][nn])->Fit(afit_burnin_corrected_vs_VPTPNdiff[n][nn], "RQ");//quiet
		double burnin_slope     = afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(1);
		double burnin_slopeerr  = afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(1);
		double burnin_incep     = afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(0);
		double burnin_inceperr  = afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(0);
		int burnin_slopebin     = g_burnin_corrected_slopes[n][nn]    ->GetN();
		int burnin_incepbin     = g_burnin_corrected_onevalues[n][nn] ->GetN();
		g_burnin_corrected_slopes[n][nn]    ->SetPoint(     burnin_slopebin, d, burnin_slope   );
		if(!dontploterrer) g_burnin_corrected_slopes[n][nn]    ->SetPointError(burnin_slopebin, 0, burnin_slopeerr);
		g_burnin_corrected_onevalues[n][nn] ->SetPoint(     burnin_incepbin, d, burnin_incep   );
		if(!dontploterrer) g_burnin_corrected_onevalues[n][nn] ->SetPointError(burnin_incepbin, 0, burnin_inceperr);
		if(printtime && makeburnincorrplots){
			string hname = (string)"profile_burnin_corrected_vs_VPTPNdiff_Eta_" + string_feta[n] + string_prod[nn];
			p_burnin_corrected_vs_VPTPNdiff[n][nn] = new TProfile(hname.c_str(),hname.c_str(), 50, 0.8, 1.1);
			p_burnin_corrected_vs_VPTPNdiff[n][nn]->SetMarkerStyle(4); //a_burnin_vs_VPTPNdiff[n][nn]->SetMarkerSize(4);
			p_burnin_corrected_vs_VPTPNdiff[n][nn]->GetYaxis()->SetTitle("VPT/PN_{max}-VPT/PN_{min}"); a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetXaxis()->SetTitle("burn-in ratio");
			for(int nnn = 0; nnn<a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN(); ++nnn){
				double xn, yn;
				a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetPoint(nnn,xn,yn);
				p_burnin_corrected_vs_VPTPNdiff[n][nn]->Fill(xn, yn);
			}
			if(fitlinepluserr){
				fitup   = new TF1("fitup",   "[0]+[1]*(x-1.0)", 0.8, 1.1);
				fitdown = new TF1("fitdown", "[0]+[1]*(x-1.0)", 0.8, 1.1);
				fitup  ->SetLineColor(kRed); fitup  ->SetLineStyle(3); fitup  ->SetLineWidth(2);
				fitdown->SetLineColor(kRed); fitdown->SetLineStyle(3); fitdown->SetLineWidth(2);
				fitup  ->SetParameter(0, afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(0)+afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(0));
				fitdown->SetParameter(0, afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(0)-afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(0));
				fitup  ->SetParameter(1, afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(1)-afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(1));
				fitdown->SetParameter(1, afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParameter(1)+afit_burnin_corrected_vs_VPTPNdiff[n][nn]->GetParError(1));
			}
			TH1F *haxis = new TH1F("haxis","",165,0.85,1.05);
			haxis->SetMinimum(-0.1);
			haxis->SetMaximum(0.7);
			haxis->GetXaxis()->SetTitle("burn-in");
			haxis->GetYaxis()->SetTitle("#mu-corrected #Delta#frac{VPT}{PN}");
			haxis->SetDirectory(0);                haxis->SetStats(0);
			haxis->GetXaxis()->SetLabelFont(42);   haxis->GetXaxis()->SetLabelSize(0.055); 
			haxis->GetXaxis()->SetTitleSize(0.06); haxis->GetXaxis()->SetTitleOffset(1.1); haxis->GetXaxis()->SetTitleFont(42);
			haxis->GetYaxis()->SetLabelFont(42);   haxis->GetYaxis()->SetLabelSize(0.055); 
			haxis->GetYaxis()->SetTitleSize(0.06); haxis->GetYaxis()->SetTitleOffset(1.1); haxis->GetYaxis()->SetTitleFont(42);
		/*	a_burnin_corrected_vs_VPTPNdiff[n][nn]->Draw("AP");
			TH1F *haxis = a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetHistogram();
			haxis->SetNameTitle(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetName(), a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetTitle());
			haxis->GetXaxis()->SetRangeUser(0.85,1.05);
			haxis->GetXaxis()->SetLimits(0.85,1.05);
			haxis->GetXaxis()->SetTitle("burn-in");
			haxis->GetYaxis()->SetTitle("corrected #frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetXTitle("burn-in");
			haxis->SetYTitle("corrected #frac{VPT}{PN} _{max} - #frac{VPT}{PN} _{min}");
			haxis->SetMinimum(-0.5);
			haxis->SetMaximum( 0.7);
		*/	haxis->Draw("hist");
			a_burnin_corrected_vs_VPTPNdiff[n][nn]->Draw("P");
			afit_burnin_corrected_vs_VPTPNdiff[n][nn]->Draw("same");
			if(fitlinepluserr){
				fitup  ->Draw("same");
				fitdown->Draw("same");
			}
			TString outputfile = a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetName();
			TString outdirtemp = outdirburninmustdcorr+daydir;
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
			col->Clear();
			col->cd();
			haxis->Draw("hist");
			p_burnin_corrected_vs_VPTPNdiff[n][nn]->Draw("same");
			afit_burnin_corrected_vs_VPTPNdiff[n][nn]->Draw("same");
				if(fitlinepluserr){
					fitup  ->Draw("same");
					fitdown->Draw("same");
				}
			outputfile = p_burnin_corrected_vs_VPTPNdiff[n][nn]->GetName();
			outdirtemp = outdirburninmustdcorr+daydir;
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_png) Util::PrintNoEPS(col, outputfile, outdirtemp);
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_eps) Util::PrintEPS(col, outputfile, outdirtemp);
			if(a_burnin_corrected_vs_VPTPNdiff[n][nn]->GetN()>3 && p_c  ) col->SaveAs(outdirtemp + "/" + outputfile + ".C");
			col->Clear();
			delete p_burnin_corrected_vs_VPTPNdiff[n][nn];
			delete haxis;
			if(fitlinepluserr){
				delete fitup;
				delete fitdown;
			}
		}
	}}
//endburnincorrected

	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		//delete now
		delete a_burnin_vs_VPTPNdiff[n][nn];
		delete afit_burnin_vs_VPTPNdiff[n][nn];
		delete a_mustd_vs_VPTPNdiff[n][nn];
		delete afit_mustd_vs_VPTPNdiff[n][nn];
		delete a_burnin_corrected_vs_VPTPNdiff[n][nn];
		delete afit_burnin_corrected_vs_VPTPNdiff[n][nn];
		delete a_gainqe_vs_VPTPNdiff[n][nn];
		delete afit_gainqe_vs_VPTPNdiff[n][nn];
		delete a_gain_vs_VPTPNdiff[n][nn];
		delete afit_gain_vs_VPTPNdiff[n][nn];
		delete a_qe_vs_VPTPNdiff[n][nn];
		delete afit_qe_vs_VPTPNdiff[n][nn];
		delete p_burnin_vs_VPTPNdiff[n][nn];
		delete pc_burnin_vs_VPTPNdiff[n][nn];
		delete p_mustd_vs_VPTPNdiff[n][nn];
		delete pc_mustd_vs_VPTPNdiff[n][nn];
		for(int nnn=0; nnn<6;++nnn){
			delete a_splittedburnin_vs_VPTPNdiff[n][nn][nnn];
			delete afit_splittedburnin_vs_VPTPNdiff[n][nn][nnn];
		}
	}}
	
	}// int d

	if(plotoverview){
	//now do plotting and saving
	TString outname;
	double slopemax(-99), incepmax(-99), slopemin(999), incepmin(999), corrcoeffmax(-99), corrcoeffmin(999);
	double timemin(10e16), timemax(-1);
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		for(int nnn=0; nnn<g_burnin_slopes[n][nn]->GetN();++nnn){
			double dx,dy;
			g_burnin_slopes[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>slopemax) slopemax = dy;
			if(dy<slopemin) slopemin = dy;
			if(dx>timemax ) timemax  = dx;
			if(dx<timemin ) timemin  = dx;
			//cout << "slope " << n << " " << nn << " point " << nnn << " has time " << (int)dx <<  " and slope " << dy << endl;
		}
		for(int nnn=0; nnn<g_burnin_onevalues[n][nn]->GetN();++nnn){
			double dx,dy;
			g_burnin_onevalues[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>incepmax) incepmax = dy;
			if(dy<incepmin) incepmin = dy;
			if(dx>timemax ) timemax  = dx;
			if(dx<timemin ) timemin  = dx;
		}
		for(int nnn=0; nnn<g_burnin_corrcoeff[n][nn]->GetN();++nnn){
			double dx,dy;
			g_burnin_corrcoeff[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>corrcoeffmax) corrcoeffmax = dy;
			if(dy<corrcoeffmin) corrcoeffmin = dy;
			//cout << "incep " << n << " " << nn << " point " << nnn << " has time " << (int)dx << " and onevalue " << dy << endl;
		}
	}}
	if(corrcoeffmax>0) corrcoeffmax = 1.2 * corrcoeffmax; else corrcoeffmax = 0.8 * corrcoeffmax;
	if(corrcoeffmin<0) corrcoeffmin = 1.2 * corrcoeffmin; else corrcoeffmin = 0.8 * corrcoeffmin;
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


	//legend
	TH1F *haxis = new TH1F("haxis","",25,1.297984e+09,1.369044e+09);
	//haxis->SetMinimum(-0.5); haxis->SetMaximum(0.3); //set separately
	//haxis->GetYaxis()->SetTitle("#rho"); //set separately
	haxis->SetStats(0);
	haxis->GetXaxis()->SetTitle("time");
	haxis->GetXaxis()->SetTimeDisplay(1);   haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}%F1970-01-1 00:00:00");
        haxis->GetXaxis()->SetNdivisions(509);
	haxis->GetXaxis()->SetLabelFont(42);    haxis->GetXaxis()->SetLabelOffset(0.03); haxis->GetXaxis()->SetLabelSize(0.06);
	haxis->GetXaxis()->SetTitleSize(0.06);  haxis->GetXaxis()->SetTitleOffset(1.44); haxis->GetXaxis()->SetTitleFont(42);
	haxis->GetYaxis()->SetLabelFont(42);    haxis->GetYaxis()->SetLabelOffset(0.01); haxis->GetYaxis()->SetLabelSize(0.06);
	haxis->GetYaxis()->SetTitleSize(0.06);  haxis->GetYaxis()->SetTitleOffset(1.03); haxis->GetYaxis()->SetTitleFont(42);
	haxis->GetZaxis()->SetLabelFont(42);    haxis->GetZaxis()->SetLabelSize(0.035);
	haxis->GetZaxis()->SetTitleSize(0.035); haxis->GetZaxis()->SetTitleFont(42);
/*
	TH1F *haxis = new TH1F("haxis", "haxis", 25, timemin, timemax);
	haxis->GetXaxis()->SetTitle("time"); haxis->GetXaxis()->SetTitleOffset(1.25); haxis->GetXaxis()->SetLabelSize(0.027);
	haxis->GetXaxis()->SetLabelOffset(0.02);haxis->GetYaxis()->SetTitleOffset(1.25); haxis->GetYaxis()->SetLabelSize(0.027);
	haxis->GetXaxis()->SetTimeDisplay(1); haxis->GetXaxis()->SetTimeFormat("%y/%m/%d%F1970-01-1 00:00:00");
	haxis->GetXaxis()->SetTimeFormat("#splitline{%Y}{%d/%m}"); haxis->GetXaxis()->SetNdivisions(510, true);
*/
	//burnin
	if(makeburninplots){
	slopemin = -2.;
	slopemax = 0.5;
	incepmin = -0.1;
	incepmax = 1.1;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("burn-in slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_burnin" );
//	col->SetTitle("slopes_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_slopes[n][0]->GetN()>0) { g_burnin_slopes[n][0]->Sort();
			g_burnin_slopes[n][0]->Draw("3"); g_burnin_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_slopes[n][0]->GetN()>0) g_burnin_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
	if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
	col->Clear();
	if(makesplittedburninplots){
		for(int nnn=0; nnn<6;++nnn){
			col->cd();
			if(nnn==0) col->SetName( "lowqe_slopes_vs_time_daily_russian_burnin" );
			if(nnn==1) col->SetName( "highqe_slopes_vs_time_daily_russian_burnin" );
			if(nnn==2) col->SetName( "lowgain_slopes_vs_time_daily_russian_burnin" );
			if(nnn==3) col->SetName( "highgain_slopes_vs_time_daily_russian_burnin" );
			if(nnn==4) col->SetName( "lowgainqe_slopes_vs_time_daily_russian_burnin" );
			if(nnn==5) col->SetName( "highgainqe_slopes_vs_time_daily_russian_burnin" );
			haxis->Draw("hist");
			if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_splittedburnin_slopes[n][0][nnn]->GetN()>0) { g_splittedburnin_slopes[n][0][nnn]->Sort();
					g_splittedburnin_slopes[n][0][nnn]->Draw("3"); g_splittedburnin_slopes[n][0][nnn]->Draw("pX"); }
					}
			} else { for(int n = 0; n<8; ++n) {if(g_splittedburnin_slopes[n][0][nnn]->GetN()>0) g_splittedburnin_slopes[n][0][nnn]->Draw("P");} }
			Legend1->Draw("same");
			col->Update();
			outname = col->GetName();
			if(p_png) Util::PrintNoEPS(col, outname, outdirburninsplitgainqe);
			if(p_eps) Util::PrintEPS(  col, outname, outdirburninsplitgainqe);
			if(p_c  ) col->SaveAs(outdirburninsplitgainqe + "/" + outname + ".C");
			col->Clear();
		}
	}
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_burnin" );
//	col->SetTitle("slopes_vs_time_daily_chinese_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_slopes[n][1]->GetN()>0) { g_burnin_slopes[n][1]->Sort();
			g_burnin_slopes[n][1]->Draw("3"); g_burnin_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_slopes[n][1]->GetN()>0) g_burnin_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
	if(p_eps) Util::PrintEPS(col, outname, outdirburnin);
	if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
	col->Clear();
	if(makesplittedburninplots){
		for(int nnn=0; nnn<6;++nnn){
			col->cd();
			if(nnn==0) col->SetName( "lowqe_slopes_vs_time_daily_chinese_burnin" );
			if(nnn==1) col->SetName( "highqe_slopes_vs_time_daily_chinese_burnin" );
			if(nnn==2) col->SetName( "lowgain_slopes_vs_time_daily_chinese_burnin" );
			if(nnn==3) col->SetName( "highgain_slopes_vs_time_daily_chinese_burnin" );
			if(nnn==4) col->SetName( "lowgainqe_slopes_vs_time_daily_chinese_burnin" );
			if(nnn==5) col->SetName( "highgainqe_slopes_vs_time_daily_chinese_burnin" );
			haxis->Draw("hist");
			if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_splittedburnin_slopes[n][1][nnn]->GetN()>0) { g_splittedburnin_slopes[n][1][nnn]->Sort();
					g_splittedburnin_slopes[n][1][nnn]->Draw("3"); g_splittedburnin_slopes[n][1][nnn]->Draw("pX"); }
					}
			} else { for(int n = 0; n<8; ++n) {if(g_splittedburnin_slopes[n][1][nnn]->GetN()>0) g_splittedburnin_slopes[n][1][nnn]->Draw("P");} }
			Legend1->Draw("same");
			col->Update();
			outname = col->GetName();
			if(p_png) Util::PrintNoEPS(col, outname, outdirburninsplitgainqe);
			if(p_eps) Util::PrintEPS(col, outname, outdirburninsplitgainqe);
			if(p_c  ) col->SaveAs(outdirburninsplitgainqe + "/" + outname + ".C");
			col->Clear();
		}
	}
	haxis->SetName("onevalues_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("burn-in = 1"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "onevalues_vs_time_daily_russian_burnin" );
//	col->SetTitle("onevalues_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_onevalues[n][0]->GetN()>0) { g_burnin_onevalues[n][0]->Sort();
			g_burnin_onevalues[n][0]->Draw("3"); g_burnin_onevalues[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_onevalues[n][0]->GetN()>0) g_burnin_onevalues[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
	if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
	col->Clear();
	if(makesplittedburninplots){
		for(int nnn=0; nnn<6;++nnn){
			col->cd();
			if(nnn==0) col->SetName( "lowqe_onevalues_vs_time_daily_russian_burnin" );
			if(nnn==1) col->SetName( "highqe_onevalues_vs_time_daily_russian_burnin" );
			if(nnn==2) col->SetName( "lowgain_onevalues_vs_time_daily_russian_burnin" );
			if(nnn==3) col->SetName( "highgain_onevalues_vs_time_daily_russian_burnin" );
			if(nnn==4) col->SetName( "lowgainqe_onevalues_vs_time_daily_russian_burnin" );
			if(nnn==5) col->SetName( "highgainqe_onevalues_vs_time_daily_russian_burnin" );
			haxis->Draw("hist");
			if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_splittedburnin_onevalues[n][0][nnn]->GetN()>0) { g_splittedburnin_onevalues[n][0][nnn]->Sort();
					g_splittedburnin_onevalues[n][0][nnn]->Draw("3"); g_splittedburnin_onevalues[n][0][nnn]->Draw("pX"); }
					}
			} else { for(int n = 0; n<8; ++n) {if(g_splittedburnin_onevalues[n][0][nnn]->GetN()>0) g_splittedburnin_onevalues[n][0][nnn]->Draw("P");} }
			Legend1->Draw("same");
			col->Update();
			outname = col->GetName();
			if(p_png) Util::PrintNoEPS(col, outname, outdirburninsplitgainqe);
			if(p_eps) Util::PrintEPS(  col, outname, outdirburninsplitgainqe);
			if(p_c  ) col->SaveAs(outdirburninsplitgainqe + "/" + outname + ".C");
			col->Clear();
		}
	}
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "onevalues_vs_time_daily_chinese_burnin" );
//	col->SetTitle("onevalues_vs_time_daily_chinese_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_onevalues[n][1]->GetN()>0) { g_burnin_onevalues[n][1]->Sort();
			g_burnin_onevalues[n][1]->Draw("3"); g_burnin_onevalues[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_onevalues[n][1]->GetN()>0) g_burnin_onevalues[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);	
	if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");	
	col->Clear();
	if(makesplittedburninplots){
		for(int nnn=0; nnn<6;++nnn){
			col->cd();
			if(nnn==0) col->SetName( "lowqe_onevalues_vs_time_daily_chinese_burnin" );
			if(nnn==1) col->SetName( "highqe_onevalues_vs_time_daily_chinese_burnin" );
			if(nnn==2) col->SetName( "lowgain_onevalues_vs_time_daily_chinese_burnin" );
			if(nnn==3) col->SetName( "highgain_onevalues_vs_time_daily_chinese_burnin" );
			if(nnn==4) col->SetName( "lowgainqe_onevalues_vs_time_daily_chinese_burnin" );
			if(nnn==5) col->SetName( "highgainqe_onevalues_vs_time_daily_chinese_burnin" );
			haxis->Draw("hist");
			if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_splittedburnin_onevalues[n][1][nnn]->GetN()>0) { g_splittedburnin_onevalues[n][1][nnn]->Sort();
					g_splittedburnin_onevalues[n][1][nnn]->Draw("3"); g_splittedburnin_onevalues[n][1][nnn]->Draw("pX"); }
					}
			} else { for(int n = 0; n<8; ++n) {if(g_splittedburnin_onevalues[n][1][nnn]->GetN()>0) g_splittedburnin_onevalues[n][1][nnn]->Draw("P");} }
			Legend1->Draw("same");
			col->Update();
			outname = col->GetName();
			if(p_png) Util::PrintNoEPS(col, outname, outdirburninsplitgainqe);
			if(p_eps) Util::PrintEPS(  col, outname, outdirburninsplitgainqe);
			if(p_c  ) col->SaveAs(outdirburninsplitgainqe + "/" + outname + ".C");
			col->Clear();
		}
	}
	if(makecorrcoeffplots){
		haxis->SetName("corrcoeff_vs_time");
		haxis->SetMinimum(corrcoeffmin);
		haxis->SetMaximum(corrcoeffmax);
	//	haxis->GetYaxis()->SetTitle("burn-in correlation coefficient"); 
		haxis->GetYaxis()->SetTitle("#rho"); 
		col->Clear();
		col->cd();
		Legend1->SetHeader("BTCP crystals");
		col->SetName( "corrcoeff_vs_time_daily_russian_burnin" );
		haxis->Draw("hist");
		if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrcoeff[n][0]->GetN()>0) { g_burnin_corrcoeff[n][0]->Sort();
				g_burnin_corrcoeff[n][0]->Draw("3"); g_burnin_corrcoeff[n][0]->Draw("pX"); }
				}
		} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrcoeff[n][0]->GetN()>0) g_burnin_corrcoeff[n][0]->Draw("P");} }
		Legend1->Draw("same");
		col->Update();
		outname = col->GetName();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
		col->Clear();
		col->cd();
		Legend1->SetHeader("SIC crystals");
		col->SetName("corrcoeff_vs_time_daily_chinese_burnin" );
		haxis->Draw("hist");
		if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrcoeff[n][1]->GetN()>0) {g_burnin_corrcoeff[n][1]->Sort();
				g_burnin_corrcoeff[n][1]->Draw("3"); g_burnin_corrcoeff[n][1]->Draw("pX"); }
				}
		} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrcoeff[n][1]->GetN()>0) g_burnin_corrcoeff[n][1]->Draw("P");} }
		Legend1->Draw("same");
		col->Update();
		outname = col->GetName();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);	
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");	
		col->Clear();
	}
	if(fprinttime){
		col->cd();
		outname = h_burnin_eq1[0]->GetName();
		h_burnin_eq1[0]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_burnin_eq1[1]->GetName();
		h_burnin_eq1[1]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_burnin_slope[0]->GetName();
		h_burnin_slope[0]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_burnin_slope[1]->GetName();
		h_burnin_slope[1]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirburnin);
		if(p_eps) Util::PrintEPS(  col, outname, outdirburnin);
		if(p_c  ) col->SaveAs(outdirburnin + "/" + outname + ".C");
		col->Clear();
	}
	}//burnin!!

	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		for(int nnn=0; nnn<g_mustd_corrcoeff[n][nn]->GetN();++nnn){
			double dx,dy;
			g_mustd_corrcoeff[n][nn]->GetPoint(nnn,dx,dy);
			if(dy>corrcoeffmax) corrcoeffmax = dy;
			if(dy<corrcoeffmin) corrcoeffmin = dy;
			//cout << "incep " << n << " " << nn << " point " << nnn << " has time " << (int)dx << " and onevalue " << dy << endl;
		}
	}}
	if(corrcoeffmax>0) corrcoeffmax = 1.2 * corrcoeffmax; else corrcoeffmax = 0.8 * corrcoeffmax;
	if(corrcoeffmin<0) corrcoeffmin = 1.2 * corrcoeffmin; else corrcoeffmin = 0.8 * corrcoeffmin;

	//mustd
	if(makemustdplots){
	slopemin = -0.15;
	slopemax = 0.6;
	incepmin = -0.15;
	incepmax = 0.6;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("#mu-slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_mustd" );
//	col->SetTitle("slopes_vs_time_daily_russian_mustd" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_slopes[n][0]->GetN()>0) { g_mustd_slopes[n][0]->Sort();
			g_mustd_slopes[n][0]->Draw("3"); g_mustd_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_mustd_slopes[n][0]->GetN()>0) g_mustd_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
	if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_mustd" );
//	col->SetTitle("slopes_vs_time_daily_chinese_mustd" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_slopes[n][1]->GetN()>0) { g_mustd_slopes[n][1]->Sort();
			g_mustd_slopes[n][1]->Draw("3"); g_mustd_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_mustd_slopes[n][1]->GetN()>0) g_mustd_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
	if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
	haxis->SetName("intercepts_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("#mu-intercept"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "intercepts_vs_time_daily_russian_mustd" );
//	col->SetTitle("intercepts_vs_time_daily_russian_mustd" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_intercepts[n][0]->GetN()>0) { g_mustd_intercepts[n][0]->Sort();
			g_mustd_intercepts[n][0]->Draw("3"); g_mustd_intercepts[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_mustd_intercepts[n][0]->GetN()>0) g_mustd_intercepts[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
	if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "intercepts_vs_time_daily_chinese_mustd" );
//	col->SetTitle("intercepts_vs_time_daily_chinese_mustd" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_intercepts[n][1]->GetN()>0) { g_mustd_intercepts[n][1]->Sort();
			g_mustd_intercepts[n][1]->Draw("3"); g_mustd_intercepts[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_mustd_intercepts[n][1]->GetN()>0) g_mustd_intercepts[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
	if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
	if(makecorrcoeffplots){
		haxis->SetName("corrcoeff_vs_time");
		haxis->SetMinimum(corrcoeffmin);
		haxis->SetMaximum(corrcoeffmax);
	//	haxis->GetYaxis()->SetTitle("#mu_{std} correlation coefficient"); 
		haxis->GetYaxis()->SetTitle("#rho"); 
		col->Clear();
		col->cd();
		Legend1->SetHeader("BTCP crystals");
		col->SetName( "corrcoeff_vs_time_daily_russian_mustd" );
		haxis->Draw("hist");
		if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_corrcoeff[n][0]->GetN()>0) { g_mustd_corrcoeff[n][0]->Sort();
				g_mustd_corrcoeff[n][0]->Draw("3"); g_mustd_corrcoeff[n][0]->Draw("pX"); }
				}
		} else { for(int n = 0; n<8; ++n) {if(g_mustd_corrcoeff[n][0]->GetN()>0) g_mustd_corrcoeff[n][0]->Draw("P");} }
		Legend1->Draw("same");
		col->Update();
		outname = col->GetName();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
		col->Clear();
		col->cd();
		Legend1->SetHeader("SIC crystals");
		col->SetName("corrcoeff_vs_time_daily_chinese_mustd" );
		haxis->Draw("hist");
		if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_mustd_corrcoeff[n][1]->GetN()>0) { g_mustd_corrcoeff[n][1]->Sort();
				g_mustd_corrcoeff[n][1]->Draw("3"); g_mustd_corrcoeff[n][1]->Draw("pX"); }
				}
		} else { for(int n = 0; n<8; ++n) {if(g_mustd_corrcoeff[n][1]->GetN()>0) g_mustd_corrcoeff[n][1]->Draw("P");} }
		Legend1->Draw("same");
		col->Update();
		outname = col->GetName();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);	
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");	
		col->Clear();
	}
	if(fprinttime){
		col->cd();
		outname = h_mustd_intercept[0]->GetName();
		h_mustd_intercept[0]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_mustd_intercept[1]->GetName();
		h_mustd_intercept[1]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_mustd_slope[0]->GetName();
		h_mustd_slope[0]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
		col->Clear();
		col->cd();
		outname = h_mustd_slope[1]->GetName();
		h_mustd_slope[1]->Draw();
		col->Update();
		if(p_png) Util::PrintNoEPS(col, outname, outdirmustd);
		if(p_eps) Util::PrintEPS(  col, outname, outdirmustd);
		if(p_c  ) col->SaveAs(outdirmustd + "/" + outname + ".C");
		col->Clear();
	}
	}//mustd

	if(makeburnincorrplots){
	//burnin corrected
	slopemin = -2.;
	slopemax = 0.5;
	incepmin = -0.1;
	incepmax = 1.1;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("#mu-corrected burn-in slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_burnin" );
//	col->SetTitle("slopes_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrected_slopes[n][0]->GetN()>0) { g_burnin_corrected_slopes[n][0]->Sort();
			g_burnin_corrected_slopes[n][0]->Draw("3"); g_burnin_corrected_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrected_slopes[n][0]->GetN()>0) g_burnin_corrected_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburninmustdcorr);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburninmustdcorr);
	if(p_c  ) col->SaveAs(outdirburninmustdcorr + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_burnin" );
//	col->SetTitle("slopes_vs_time_daily_chinese_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrected_slopes[n][1]->GetN()>0) { g_burnin_corrected_slopes[n][1]->Sort();
			g_burnin_corrected_slopes[n][1]->Draw("3"); g_burnin_corrected_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrected_slopes[n][1]->GetN()>0) g_burnin_corrected_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburninmustdcorr);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburninmustdcorr);
	if(p_c  ) col->SaveAs(outdirburninmustdcorr + "/" + outname + ".C");
	haxis->SetName("onevalues_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("#mu-corrected burn-in = 1"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "onevalues_vs_time_daily_russian_burnin" );
//	col->SetTitle("onevalues_vs_time_daily_russian_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrected_onevalues[n][0]->GetN()>0) { g_burnin_corrected_onevalues[n][0]->Sort();
			g_burnin_corrected_onevalues[n][0]->Draw("3"); g_burnin_corrected_onevalues[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrected_onevalues[n][0]->GetN()>0) g_burnin_corrected_onevalues[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburninmustdcorr);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburninmustdcorr);
	if(p_c  ) col->SaveAs(outdirburninmustdcorr + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "onevalues_vs_time_daily_chinese_burnin" );
//	col->SetTitle("onevalues_vs_time_daily_chinese_burnin" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_corrected_onevalues[n][1]->GetN()>0) { g_burnin_corrected_onevalues[n][1]->Sort();
			g_burnin_corrected_onevalues[n][1]->Draw("3"); g_burnin_corrected_onevalues[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_corrected_onevalues[n][1]->GetN()>0) g_burnin_corrected_onevalues[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirburninmustdcorr);
	if(p_eps) Util::PrintEPS(  col, outname, outdirburninmustdcorr);	
	if(p_c  ) col->SaveAs(outdirburninmustdcorr + "/" + outname + ".C");	
	}

	if(makeburninvsmustdplots){
	//correlate mustd fits with burninfits
	TH1F *haxis1 = new TH1F("haxis1", "", 25, -1.2, 0.6);//burnin slopes vs. mustd intercepts
	haxis1->SetMinimum(-0.1); haxis1->SetMaximum(0.7);
	haxis1->GetXaxis()->SetLabelFont(42);   haxis1->GetXaxis()->SetLabelSize(0.06);
	haxis1->GetXaxis()->SetTitleSize(0.06); haxis1->GetXaxis()->SetTitleOffset(1.15); haxis1->GetXaxis()->SetTitleFont(42);
	haxis1->GetYaxis()->SetLabelFont(42);   haxis1->GetYaxis()->SetLabelSize(0.05);
	haxis1->GetYaxis()->SetTitleSize(0.06); haxis1->GetYaxis()->SetTitleOffset(1.15); haxis1->GetYaxis()->SetTitleFont(42);
	TH1F *haxis2 = new TH1F("haxis2", "", 25, -0.1, 1.5);//burnin = 1    vs. mustd slopes
	haxis2->SetMinimum(-0.1); haxis2->SetMaximum(0.4);
	haxis2->GetXaxis()->SetLabelFont(42);   haxis2->GetXaxis()->SetLabelSize(0.06);
	haxis2->GetXaxis()->SetTitleSize(0.06); haxis2->GetXaxis()->SetTitleOffset(1.15); haxis2->GetXaxis()->SetTitleFont(42);
	haxis2->GetYaxis()->SetLabelFont(42);   haxis2->GetYaxis()->SetLabelSize(0.05);
	haxis2->GetYaxis()->SetTitleSize(0.06); haxis2->GetYaxis()->SetTitleOffset(1.15); haxis2->GetYaxis()->SetTitleFont(42);
//	haxis1->SetNameTitle("burninslopes_vs_mustdintercepts", "burnin slopes vs. #mu_{std} intercepts");
	haxis1->SetName("burninslopes_vs_mustdintercepts");
//	haxis1->SetMinimum(slopemin);
//	haxis1->SetMaximum(slopemax);
	haxis1->GetXaxis()->SetTitle("burn-in slope"); 
	haxis1->GetYaxis()->SetTitle("#mu-intercept"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "burninslopes_vs_mustdintercepts_daily_russian" );
//	col->SetTitle("burninslopes_vs_mustdintercepts_daily_russian" );
	haxis1->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_slope_vs_mustd_intercept[n][0]->GetN()>0) { g_burnin_slope_vs_mustd_intercept[n][0]->Sort();
			g_burnin_slope_vs_mustd_intercept[n][0]->Draw("3"); g_burnin_slope_vs_mustd_intercept[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_slope_vs_mustd_intercept[n][0]->GetN()>0) g_burnin_slope_vs_mustd_intercept[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustdcorrburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustdcorrburnin);
	if(p_c  ) col->SaveAs(outdirmustdcorrburnin + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "burninslopes_vs_mustdintercepts_daily_chinese" );
//	col->SetTitle("burninslopes_vs_mustdintercepts_daily_chinese" );
	haxis1->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_slope_vs_mustd_intercept[n][1]->GetN()>0) { g_burnin_slope_vs_mustd_intercept[n][1]->Sort();
			g_burnin_slope_vs_mustd_intercept[n][1]->Draw("3"); g_burnin_slope_vs_mustd_intercept[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_slope_vs_mustd_intercept[n][1]->GetN()>0) g_burnin_slope_vs_mustd_intercept[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustdcorrburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustdcorrburnin);
	if(p_c  ) col->SaveAs(outdirmustdcorrburnin + "/" + outname + ".C");
//	haxis2->SetNameTitle("burninonevalues_vs_mustdslopes", "burnin=1 equivalents vs. #mu_{std} slopes");
	haxis2->SetName("burninonevalues_vs_mustdslopes");
//	haxis2->SetMinimum(incepmin);
//	haxis2->SetMaximum(incepmax);
	haxis2->GetXaxis()->SetTitle("burn-in = 1"); 
	haxis2->GetYaxis()->SetTitle("#mu-slope"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "burninonevalues_vs_mustdslopes_daily_russian" );
//	col->SetTitle("burninonevalues_vs_mustdslopes_daily_russian" );
	haxis2->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_onevalue_vs_mustd_slope[n][0]->GetN()>0) { g_burnin_onevalue_vs_mustd_slope[n][0]->Sort();
			g_burnin_onevalue_vs_mustd_slope[n][0]->Draw("3"); g_burnin_onevalue_vs_mustd_slope[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_onevalue_vs_mustd_slope[n][0]->GetN()>0) g_burnin_onevalue_vs_mustd_slope[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustdcorrburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustdcorrburnin);
	if(p_c  ) col->SaveAs(outdirmustdcorrburnin + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "burninonevalues_vs_mustdslopes_daily_chinese" );
//	col->SetTitle("burninonevalues_vs_mustdslopes_daily_chinese" );
	haxis2->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_burnin_onevalue_vs_mustd_slope[n][1]->GetN()>0) { g_burnin_onevalue_vs_mustd_slope[n][1]->Sort();
			g_burnin_onevalue_vs_mustd_slope[n][1]->Draw("3"); g_burnin_onevalue_vs_mustd_slope[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_burnin_onevalue_vs_mustd_slope[n][1]->GetN()>0) g_burnin_onevalue_vs_mustd_slope[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirmustdcorrburnin);
	if(p_eps) Util::PrintEPS(  col, outname, outdirmustdcorrburnin);	
	if(p_c  ) col->SaveAs(outdirmustdcorrburnin + "/" + outname + ".C");
	delete haxis1;
	delete haxis2;
	}

	if(makegainqeplots){
	//qe*gain
	slopemin = -0.4;
	slopemax = 0.3;
	incepmin = -0.1;
	incepmax = 1.0;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("q.e.#cdot gain slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_gainqe" );
//	col->SetTitle("slopes_vs_time_daily_russian_gainqe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gainqe_slopes[n][0]->GetN()>0) { g_gainqe_slopes[n][0]->Sort();
			g_gainqe_slopes[n][0]->Draw("3"); g_gainqe_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gainqe_slopes[n][0]->GetN()>0) g_gainqe_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgainqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgainqe);
	if(p_c  ) col->SaveAs(outdirgainqe + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_gainqe" );
//	col->SetTitle("slopes_vs_time_daily_chinese_gainqe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gainqe_slopes[n][1]->GetN()>0) { g_gainqe_slopes[n][1]->Sort();
			g_gainqe_slopes[n][1]->Draw("3"); g_gainqe_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gainqe_slopes[n][1]->GetN()>0) g_gainqe_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgainqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgainqe);
	if(p_c  ) col->SaveAs(outdirgainqe + "/" + outname + ".C");
	haxis->SetName("intercepts_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("q.e.#cdot gain intercept"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "intercepts_vs_time_daily_russian_gainqe" );
//	col->SetTitle("intercepts_vs_time_daily_russian_gainqe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gainqe_intercepts[n][0]->GetN()>0) { g_gainqe_intercepts[n][0]->Sort();
			g_gainqe_intercepts[n][0]->Draw("3"); g_gainqe_intercepts[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gainqe_intercepts[n][0]->GetN()>0) g_gainqe_intercepts[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgainqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgainqe);
	if(p_c  ) col->SaveAs(outdirgainqe + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "intercepts_vs_time_daily_chinese_gainqe" );
//	col->SetTitle("intercepts_vs_time_daily_chinese_gainqe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gainqe_intercepts[n][1]->GetN()>0) { g_gainqe_intercepts[n][1]->Sort();
			g_gainqe_intercepts[n][1]->Draw("3"); g_gainqe_intercepts[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gainqe_intercepts[n][1]->GetN()>0) g_gainqe_intercepts[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgainqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgainqe);
	if(p_c  ) col->SaveAs(outdirgainqe + "/" + outname + ".C");
	}

	if(makegainplots){
	//gain
	slopemin = -0.2;
	slopemax = 0.2;
	incepmin = -0.1;
	incepmax = 1.0;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("gain slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_gain" );
//	col->SetTitle("slopes_vs_time_daily_russian_gain" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gain_slopes[n][0]->GetN()>0) { g_gain_slopes[n][0]->Sort();
			g_gain_slopes[n][0]->Draw("3"); g_gain_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gain_slopes[n][0]->GetN()>0) g_gain_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgain);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgain);
	if(p_c  ) col->SaveAs(outdirgain + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_gain" );
//	col->SetTitle("slopes_vs_time_daily_chinese_gain" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gain_slopes[n][1]->GetN()>0) { g_gain_slopes[n][1]->Sort();
			g_gain_slopes[n][1]->Draw("3"); g_gain_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gain_slopes[n][1]->GetN()>0) g_gain_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgain);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgain);
	if(p_c  ) col->SaveAs(outdirgain + "/" + outname + ".C");
	haxis->SetName("intercepts_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("gain intercept"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "intercepts_vs_time_daily_russian_gain" );
//	col->SetTitle("intercepts_vs_time_daily_russian_gain" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gain_intercepts[n][0]->GetN()>0) { g_gain_intercepts[n][0]->Sort();
			g_gain_intercepts[n][0]->Draw("3"); g_gain_intercepts[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gain_intercepts[n][0]->GetN()>0) g_gain_intercepts[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgain);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgain);
	if(p_c  ) col->SaveAs(outdirgain + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "intercepts_vs_time_daily_chinese_gain" );
//	col->SetTitle("intercepts_vs_time_daily_chinese_gain" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_gain_intercepts[n][1]->GetN()>0) { g_gain_intercepts[n][1]->Sort();
			g_gain_intercepts[n][1]->Draw("3"); g_gain_intercepts[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_gain_intercepts[n][1]->GetN()>0) g_gain_intercepts[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirgain);
	if(p_eps) Util::PrintEPS(  col, outname, outdirgain);
	if(p_c  ) col->SaveAs(outdirgain + "/" + outname + ".C");
	}

	if(makeqeplots){
	//qe
	slopemin = -1.2;
	slopemax = 0.4;
	incepmin = -1.;
	incepmax = 1.;
	haxis->SetName("slopes_vs_time");
	haxis->SetMinimum(slopemin);
	haxis->SetMaximum(slopemax);
	haxis->GetYaxis()->SetTitle("qe slope"); 
	col->cd();
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "slopes_vs_time_daily_russian_qe" );
//	col->SetTitle("slopes_vs_time_daily_russian_qe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_qe_slopes[n][0]->GetN()>0) { g_qe_slopes[n][0]->Sort();
			g_qe_slopes[n][0]->Draw("3"); g_qe_slopes[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_qe_slopes[n][0]->GetN()>0) g_qe_slopes[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirqe);
	if(p_c  ) col->SaveAs(outdirqe + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("SIC crystals");
	col->SetName( "slopes_vs_time_daily_chinese_qe" );
//	col->SetTitle("slopes_vs_time_daily_chinese_qe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_qe_slopes[n][1]->GetN()>0) { g_qe_slopes[n][1]->Sort();
			g_qe_slopes[n][1]->Draw("3"); g_qe_slopes[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_qe_slopes[n][1]->GetN()>0) g_qe_slopes[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirqe);
	if(p_c  ) col->SaveAs(outdirqe + "/" + outname + ".C");
	haxis->SetName("intercepts_vs_time");
	haxis->SetMinimum(incepmin);
	haxis->SetMaximum(incepmax);
	haxis->GetYaxis()->SetTitle("qe intercept"); 
	col->Clear();
	col->cd();
	Legend1->SetHeader("BTCP crystals");
	col->SetName( "intercepts_vs_time_daily_russian_qe" );
//	col->SetTitle("intercepts_vs_time_daily_russian_qe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_qe_intercepts[n][0]->GetN()>0) { g_qe_intercepts[n][0]->Sort();
			g_qe_intercepts[n][0]->Draw("3"); g_qe_intercepts[n][0]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_qe_intercepts[n][0]->GetN()>0) g_qe_intercepts[n][0]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirqe);
	if(p_c  ) col->SaveAs(outdirqe + "/" + outname + ".C");
	col->Clear();
	col->cd();
	Legend1->SetHeader("chinese crystals");
	col->SetName( "intercepts_vs_time_daily_chinese_qe" );
//	col->SetTitle("intercepts_vs_time_daily_chinese_qe" );
	haxis->Draw("hist");
	if(ploterrshaded) { for(int n = 0; n<8; ++n) { if(g_qe_intercepts[n][1]->GetN()>0) { g_qe_intercepts[n][1]->Sort();
			g_qe_intercepts[n][1]->Draw("3"); g_qe_intercepts[n][1]->Draw("pX"); }
			}
	} else { for(int n = 0; n<8; ++n) {if(g_qe_intercepts[n][1]->GetN()>0) g_qe_intercepts[n][1]->Draw("P");} }
	Legend1->Draw("same");
	col->Update();
	outname = col->GetName();
	if(p_png) Util::PrintNoEPS(col, outname, outdirqe);
	if(p_eps) Util::PrintEPS(  col, outname, outdirqe);
	if(p_c  ) col->SaveAs(outdirqe + "/" + outname + ".C");
	}
	delete haxis;
	}//plotoverview

	if(savefinalplots){
	TFile *newfile;
	if(fitprofile&&fitprofile2) newfile = new TFile ("FileOutputs/20130402_everything_vs_VPTPNdiff_VectorNTuples_profilefits2_daily.root", "RECREATE");
	else if(fitprofile) newfile = new TFile ("FileOutputs/20130402_everything_vs_VPTPNdiff_VectorNTuples_profilefits_daily.root", "RECREATE");
	else           newfile = new TFile ("FileOutputs/20130402_everything_vs_VPTPNdiff_VectorNTuples_fits_daily.root", "RECREATE");
	if(useMuSICinsteafofMuStd){
		if(fitprofile&&fitprofile2) newfile = new TFile ("FileOutputs/20130402_everythingSIC_vs_VPTPNdiff_VectorNTuples_profilefits2_daily.root", "RECREATE");
		else if(fitprofile) newfile = new TFile ("FileOutputs/20130402_everythingSIC_vs_VPTPNdiff_VectorNTuples_profilefits_daily.root", "RECREATE");
		else           newfile = new TFile ("FileOutputs/20130402_everythingSIC_vs_VPTPNdiff_VectorNTuples_fits_daily.root", "RECREATE");
	}
	newfile->cd();
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		g_mustd_corrcoeff[n][nn]->Write();
		g_burnin_slopes[n][nn]->Write();
		g_burnin_onevalues[n][nn]->Write();
		g_burnin_corrcoeff[n][nn]->Write();
		g_mustd_slopes[n][nn]->Write();
		g_mustd_intercepts[n][nn]->Write();
		g_burnin_corrected_slopes[n][nn]->Write();
		g_burnin_corrected_onevalues[n][nn]->Write();
		g_burnin_slope_vs_mustd_intercept[n][nn]->Write();
		g_burnin_onevalue_vs_mustd_slope[n][nn] ->Write();
		g_gainqe_slopes[n][nn]->Write();
		g_gainqe_intercepts[n][nn]->Write();
		g_gain_slopes[n][nn]->Write();
		g_gain_intercepts[n][nn]->Write();
		g_qe_slopes[n][nn]->Write();
		g_qe_intercepts[n][nn]->Write();
		for(int nnn=0; nnn<6; ++nnn){
			g_splittedburnin_slopes[n][nn][nnn]->Write();
			g_splittedburnin_onevalues[n][nn][nnn]->Write();
		}
	}}
	newfile->Close();
	cout << "file saved as " << newfile->GetName() << endl;
	}
	
}