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
#include "LaserAnalysisVectorTrees.C"
#include "RootMacros/Utilities.hh"


using namespace std;


//overlay slopes/intercept and potentially different other variables
// for various eta bins

void EEstudies_overlay_orange_laser(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	bool startdiff         = true;//only true if normalized == false, always diff w.r.t. 2011 beginning
	bool normalized        = false;
	const int period       = 7;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11
				//2011: 5 periods: 2011+2012: 7 periods
	bool samelumi          = true;//if same lumi = true, use only last three periods with first 1.1fb-1
	int  whatoverlay       = 1; // 1: different lumi periods (separate eta plots) 2: different eta bins (separate lumi period plots)
	bool only2011          = false;

	bool plotscaled        = false;//save also plots where VPT/PN has been scaled

	TFile *oldFile[period+1];
	TFile *oldledFile[period+1];
	TH1D *intercepts[period+1][2];
	TH1D *orangeleds[period+1][2];
	TGraphErrors *orange_vs_laser[8][period+1][2];

	string string_prod[2] = {"__russ", "__chin"};
	string string_aeta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_acap[1] = {"_"};

	//period == 0 is technical stop
	for(int i = 0; i<=period; ++i){
		if(i==1||i==2||i==6||i==7) continue;//temporary
		TString periodTS, periods; string pe; string histname;
		if(i==1) {periodTS = "/period110315_110327"; periods = "_period110315_110327"; pe = "_p1";}
		if(i==2) {periodTS = "/period110412_110505"; periods = "_period110412_110505"; pe = "_p2";}
		if(i==3) {periodTS = "/period110513_110701"; periods = "_period110513_110701"; pe = "_p3";}
		if(i==4) {periodTS = "/period110715_110825"; periods = "_period110715_110825"; pe = "_p4";}
		if(i==5) {periodTS = "/period110904_111102"; periods = "_period110904_111102"; pe = "_p5";}
		if(i==6) {periodTS = "/period120406_120420"; periods = "_period120406_120420"; pe = "_p6";}
		if(i==7) {periodTS = "/period120505_120618"; periods = "_period120505_120618"; pe = "_p7";}
		if(i==4&&samelumi&& only2011){periodTS = "/period110715_110808_1p1fb-1"; periods = "_period110715_110808_1p1fb-1"; pe = "_p41";}
		if(i==5&&samelumi&& only2011){periodTS = "/period110904_110923_1p1fb-1"; periods = "_period110904_110923_1p1fb-1"; pe = "_p51";}
		if(i==3&&samelumi&&!only2011){periodTS = "/period110513_110624_950pb-1"; periods = "_period110513_110624_950pb-1"; pe = "_p32";}
		if(i==4&&samelumi&&!only2011){periodTS = "/period110715_110807_950pb-1"; periods = "_period110715_110807_950pb-1"; pe = "_p42";}
		if(i==5&&samelumi&&!only2011){periodTS = "/period110904_110922_950pb-1"; periods = "_period110904_110922_950pb-1"; pe = "_p52";}
		if(i==7&&samelumi&&!only2011){periodTS = "/period120505_110517_950pb-1"; periods = "_period120505_110517_950pb-1"; pe = "_p72";}
		if(i==0){periodTS = "/period110822_110907"; periods = "_period110822_110907"; pe = "_p0";}//technical stop Sept.2011
		if(startdiff) oldledFile[i] = TFile::Open("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_diffwrt2011begin_"+periods+".root");
		if(normalized) oldledFile[i] = TFile::Open("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_"+periods+".root");
		else if(!startdiff) oldledFile[i] = TFile::Open("~/ECAL/Laser/Analysis/FileOutputs/20121023_orange_VPTPNmaxMinusVPTPNmin_nonormalization_"+periods+".root");
		oldFile[i] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_diffwrt2011begin"+periods+".root");
		if(normalized) oldFile[i] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted"+periods+".root");
		else if(!startdiff) oldFile[i] = TFile::Open("FileOutputs/20121023_mustd_vs_VPTPNmaxMinusVPTPNmin_VectorNTuples_fitted_nonormalization"+periods+".root");
		histname = "intercept_vs_AbsEta___russ";
		intercepts[i][0] = (TH1D*)oldFile[i]->Get(histname.c_str());
		histname = histname + pe;
		intercepts[i][0]->SetName((histname).c_str());
		histname = "intercept_vs_AbsEta___chin";
		intercepts[i][1] = (TH1D*)oldFile[i]->Get(histname.c_str());
		histname = histname + pe;
		intercepts[i][1]->SetName((histname).c_str());
		histname = "h_VPTPNdiff_vs_eta_russ";
		orangeleds[i][0] = (TH1D*)oldledFile[i]->Get(histname.c_str());
		histname = histname + pe;
		orangeleds[i][0]->SetName((histname).c_str());
		histname = "h_VPTPNdiff_vs_eta_chin";
		orangeleds[i][1] = (TH1D*)oldledFile[i]->Get(histname.c_str());
		histname = histname + pe;
		orangeleds[i][1]->SetName((histname).c_str());
		for(int j = 0; j<8; ++j){
			histname = "orange_vs_laser_russ";
			orange_vs_laser[j][i][0] = new TGraphErrors();
			histname = histname+pe+string_aeta[j];
			orange_vs_laser[j][i][0]->SetName((histname).c_str());
			orange_vs_laser[j][i][0]->SetLineWidth(2); orange_vs_laser[j][i][0]->SetMarkerStyle(20);
			histname = "orange_vs_laser_chin";
			orange_vs_laser[j][i][1] = new TGraphErrors();
			histname = histname+pe+string_aeta[j];
			orange_vs_laser[j][i][1]->SetName((histname).c_str());
			orange_vs_laser[j][i][1]->SetLineWidth(2); orange_vs_laser[j][i][1]->SetMarkerStyle(20);
			int bin = j+1;
			double interceptvalue = intercepts[i][0]->GetBinContent(bin);
			double intercepterror = intercepts[i][0]->GetBinError(bin);
			double orangeledvalue = orangeleds[i][0]->GetBinContent(bin);
			double orengelederror = orangeleds[i][0]->GetBinError(bin);
			int n = orange_vs_laser[j][i][0]->GetN();//dummy - should be always 0
			orange_vs_laser[j][i][0]->SetPoint(     n, orangeledvalue, interceptvalue);
			orange_vs_laser[j][i][0]->SetPointError(n, orengelederror, intercepterror);
			       interceptvalue = intercepts[i][1]->GetBinContent(bin);
			       intercepterror = intercepts[i][1]->GetBinError(bin);
			       orangeledvalue = orangeleds[i][1]->GetBinContent(bin);
			       orengelederror = orangeleds[i][1]->GetBinError(bin);
			     n = orange_vs_laser[j][i][1]->GetN();//dummy - should be always 0
			orange_vs_laser[j][i][1]->SetPoint(     n, orangeledvalue, interceptvalue);
			orange_vs_laser[j][i][1]->SetPointError(n, orengelederror, intercepterror);
		}
		//eta colors
		orange_vs_laser[0][i][0]->SetMarkerColor(kViolet);   orange_vs_laser[0][i][0]->SetLineColor(kViolet);
		orange_vs_laser[0][i][0]->SetFillColor(kViolet); 
		orange_vs_laser[0][i][1]->SetMarkerColor(kViolet);   orange_vs_laser[0][i][1]->SetLineColor(kViolet); 
		orange_vs_laser[0][i][1]->SetFillColor(kViolet); 
		orange_vs_laser[1][i][0]->SetMarkerColor(kBlue);     orange_vs_laser[1][i][0]->SetLineColor(kBlue); 
		orange_vs_laser[1][i][0]->SetFillColor(kBlue); 
		orange_vs_laser[1][i][1]->SetMarkerColor(kBlue);     orange_vs_laser[1][i][1]->SetLineColor(kBlue); 
		orange_vs_laser[1][i][1]->SetFillColor(kBlue); 
		orange_vs_laser[2][i][0]->SetMarkerColor(kCyan+2);   orange_vs_laser[2][i][0]->SetLineColor(kCyan+2); 
		orange_vs_laser[2][i][0]->SetFillColor(kCyan+2); 
		orange_vs_laser[2][i][1]->SetMarkerColor(kCyan+2);   orange_vs_laser[2][i][1]->SetLineColor(kCyan+2); 
		orange_vs_laser[2][i][1]->SetFillColor(kCyan+2); 
		orange_vs_laser[3][i][0]->SetMarkerColor(kGreen+2);  orange_vs_laser[3][i][0]->SetLineColor(kGreen+2); 
		orange_vs_laser[3][i][0]->SetFillColor(kGreen+2); 
		orange_vs_laser[3][i][1]->SetMarkerColor(kGreen+2);  orange_vs_laser[3][i][1]->SetLineColor(kGreen+2); 
		orange_vs_laser[3][i][1]->SetFillColor(kGreen+2); 
		orange_vs_laser[4][i][0]->SetMarkerColor(kYellow+3); orange_vs_laser[4][i][0]->SetLineColor(kYellow+3); 
		orange_vs_laser[4][i][0]->SetFillColor(kYellow+3); 
		orange_vs_laser[4][i][1]->SetMarkerColor(kYellow+3); orange_vs_laser[4][i][1]->SetLineColor(kYellow+3); 
		orange_vs_laser[4][i][1]->SetFillColor(kYellow+3); 
		orange_vs_laser[5][i][0]->SetMarkerColor(kOrange+5); orange_vs_laser[5][i][0]->SetLineColor(kOrange+5); 
		orange_vs_laser[5][i][0]->SetFillColor(kOrange+5); 
		orange_vs_laser[5][i][1]->SetMarkerColor(kOrange+5); orange_vs_laser[5][i][1]->SetLineColor(kOrange+5); 
		orange_vs_laser[5][i][1]->SetFillColor(kOrange+5); 
		orange_vs_laser[6][i][0]->SetMarkerColor(kOrange);   orange_vs_laser[6][i][0]->SetLineColor(kOrange); 
		orange_vs_laser[6][i][0]->SetFillColor(kOrange); 
		orange_vs_laser[6][i][1]->SetMarkerColor(kOrange);   orange_vs_laser[6][i][1]->SetLineColor(kOrange); 
		orange_vs_laser[6][i][1]->SetFillColor(kOrange); 
		orange_vs_laser[7][i][0]->SetMarkerColor(kRed);      orange_vs_laser[7][i][0]->SetLineColor(kRed); 
		orange_vs_laser[7][i][0]->SetFillColor(kRed); 
		orange_vs_laser[7][i][1]->SetMarkerColor(kRed);      orange_vs_laser[7][i][1]->SetLineColor(kRed); 
		orange_vs_laser[7][i][1]->SetFillColor(kRed); 
	}
	TString outputdir;
	outputdir = "Plots/20121023/EEstudies/continued/Orange_vs_BlueLaser/diffwrt2011begin/";
	if(normalized) outputdir = "Plots/20121023/EEstudies/continued/Orange_vs_BlueLaser/normalized/";
	else if(!startdiff) outputdir = "Plots/20121023/EEstudies/continued/Orange_vs_BlueLaser/notnormalized/";



	Util::MakeOutputDir(outputdir);

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);


	string string_eta_leg[8] = {"1.4 #leq |#eta| < 1.6", "1.6 #leq |#eta| < 1.8", "1.8 #leq |#eta| < 2.0", "2.0 #leq |#eta| < 2.2", "2.2 #leq |#eta| < 2.4", "2.4 #leq |#eta| < 2.6", "2.6 #leq |#eta| < 2.8", "2.8 #leq |#eta|"};
	string string_per_leg[period+1] = {"22.8.-7.9.2011", "15.3.-27.3.2011", "12.4.-5.5.2011", "13.5.-1.7.2011", "15.7.-22.8.2011", "7.9.-11.11.2011", "6.4.-20.4.2012", "5.5.-18.6.2012"};
	if(samelumi &&  only2011){ string_per_leg[4] = "15.7.-8.8.2011"; string_per_leg[5] = "7.9.-23.9.2011"; }
	if(samelumi && !only2011){ string_per_leg[3] = "13.5.-24.6.2011"; string_per_leg[4] = "15.7.-7.8.2011"; string_per_leg[5] = "7.9.-22.9.2011"; string_per_leg[7] = "5.5.-17.5.2012"; }

	string aname;
	TH1F *haxis;
	if(startdiff) haxis = new TH1F("haxis", "", 13, -0.05, 0.6);
	else          haxis = new TH1F("haxis", "", 13, -0.05, 0.4);
	if(only2011){
		if(normalized)     haxis->SetMaximum(0.25);
		else if(startdiff) haxis->SetMaximum(0.4);
		else               haxis->SetMaximum(0.185);
	} else{
		if(normalized)     haxis->SetMaximum(0.27);
		else if(startdiff) haxis->SetMaximum(/*0.55*/0.35);
		else               haxis->SetMaximum(/*0.22*/0.16);
	}
	haxis->SetMinimum(0.0);
	haxis->GetXaxis()->SetTitle("VPTdiff (orange led)");
	haxis->GetYaxis()->SetTitle("intercept");
	TLegend *Legend1 = new TLegend(.68,.5,.88,.88);
	TLegend *Legend2 = new TLegend(.68,.5,.88,.88);
	Legend1 -> SetFillColor(0);
	Legend1 -> SetBorderSize(0);
	Legend2 -> SetFillColor(0);
	Legend2 -> SetBorderSize(0);
 
	for(int e = 0; e<8; ++e){
		Legend1->AddEntry(orange_vs_laser[e][0][0], (string_eta_leg[e]).c_str(), "f");
	}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	for(int i = 0; i<=period; ++i){
		for(int nn = 0; nn<2; ++nn){
		if(i==1||i==2||i==6||i==7) continue;//temporary
			for(int e = 0; e<8; ++e){
				orange_vs_laser[e][i][nn]->Draw("AP");
				orange_vs_laser[e][i][nn]->SetHistogram(haxis);
			}
		}
	}
	//overlay VPTdiff_oled vs intercept of different eta in lumiperiod bins
	for(int i = 0; i<=period; ++i){
		if(i==1||i==2||i==6||i==7) continue;//temporary
		TString periods;
		if(i==1) periods = "_period110315_110327";	if(i==2) periods = "_period110412_110505";	if(i==3) periods = "_period110513_110701";
		if(i==4) periods = "_period110715_110825";	if(i==5) periods = "_period110904_111102";	if(i==6) periods = "_period120406_120420";
		if(i==7) periods = "_period120505_120618";	if(i==4&&samelumi&& only2011) periods = "_period110715_110808_1p1fb-1";
		if(i==5&&samelumi&& only2011) periods = "_period110904_110923_1p1fb-1";	if(i==3&&samelumi&&!only2011) periods = "_period110513_110624_950pb-1";	if(i==4&&samelumi&&!only2011) periods = "_period110715_110807_950pb-1";	if(i==5&&samelumi&&!only2011) periods = "_period110904_110922_950pb-1";	if(i==7&&samelumi&&!only2011) periods = "_period120505_110517_950pb-1";
		if(i==0) periods = "_period110822_110907";

		for(int nn = 0; nn<2; ++nn){
			TString outputfile, outputname;
			     if(nn==0) { Legend1->SetHeader("russian crystals");   outputfile  = "_Overlay_Russian"; }
			else if(nn==1) { Legend1->SetHeader("chinese crystals");   outputfile  = "_Overlay_Chinese"; }
			col->cd();

			outputname = "Intercept_vs_VPTdiffOrangeLED" + outputfile + periods;
			haxis->SetTitle(outputname);
			haxis->Draw();
			for(int e = 0; e<8; ++e){
				if(orange_vs_laser[e][i][nn]->GetN()>0) orange_vs_laser[e][i][nn]->Draw("P");
			}
			Legend1->Draw("same");
			col->Update();
			Util::PrintNoEPS(col, outputname, outputdir);
			Util::PrintEPS(col, outputname, outputdir);
			col->Clear();
		}
	}
	//make color scheme for period overlay
		//period colors
	for(int e = 0; e<8; ++e){
		//if(i==1||i==2||i==6||i==7) continue;//temporary
		orange_vs_laser[e][0][0]->SetMarkerColor(kViolet);   orange_vs_laser[e][0][0]->SetLineColor(kViolet); 
		orange_vs_laser[e][0][0]->SetFillColor(kViolet);
		orange_vs_laser[e][0][1]->SetMarkerColor(kViolet);   orange_vs_laser[e][0][1]->SetLineColor(kViolet); 
		orange_vs_laser[e][0][1]->SetFillColor(kViolet);
	//	orange_vs_laser[e][1][0]->SetMarkerColor(kBlue);     orange_vs_laser[e][1][0]->SetLineColor(kBlue); 
	//	orange_vs_laser[e][1][0]->SetFillColor(kBlue);
	//	orange_vs_laser[e][1][1]->SetMarkerColor(kBlue);     orange_vs_laser[e][1][1]->SetLineColor(kBlue); 
	//	orange_vs_laser[e][1][1]->SetFillColor(kBlue);
	//	orange_vs_laser[e][2][0]->SetMarkerColor(kCyan+2);   orange_vs_laser[e][2][0]->SetLineColor(kCyan+2); 
	//	orange_vs_laser[e][2][0]->SetFillColor(kCyan+2);
	//	orange_vs_laser[e][2][1]->SetMarkerColor(kCyan+2);   orange_vs_laser[e][2][1]->SetLineColor(kCyan+2); 
	//	orange_vs_laser[e][2][1]->SetFillColor(kCyan+2);
		orange_vs_laser[e][3][0]->SetMarkerColor(kGreen+2);  orange_vs_laser[e][3][0]->SetLineColor(kGreen+2); 
		orange_vs_laser[e][3][0]->SetFillColor(kGreen+2);
		orange_vs_laser[e][3][1]->SetMarkerColor(kGreen+2);  orange_vs_laser[e][3][1]->SetLineColor(kGreen+2); 
		orange_vs_laser[e][3][1]->SetFillColor(kGreen+2);
		orange_vs_laser[e][4][0]->SetMarkerColor(kYellow+3); orange_vs_laser[e][4][0]->SetLineColor(kYellow+3); 
		orange_vs_laser[e][4][0]->SetFillColor(kYellow+3);
		orange_vs_laser[e][4][1]->SetMarkerColor(kYellow+3); orange_vs_laser[e][4][1]->SetLineColor(kYellow+3); 
		orange_vs_laser[e][4][1]->SetFillColor(kYellow+3);
		orange_vs_laser[e][5][0]->SetMarkerColor(kOrange+5); orange_vs_laser[e][5][0]->SetLineColor(kOrange+5); 
		orange_vs_laser[e][5][0]->SetFillColor(kOrange+5);
		orange_vs_laser[e][5][1]->SetMarkerColor(kOrange+5); orange_vs_laser[e][5][1]->SetLineColor(kOrange+5); 
		orange_vs_laser[e][5][1]->SetFillColor(kOrange+5);
	//	orange_vs_laser[e][6][0]->SetMarkerColor(kOrange);   orange_vs_laser[e][6][0]->SetLineColor(kOrange); 
	//	orange_vs_laser[e][6][0]->SetFillColor(kOrange);
	//	orange_vs_laser[e][6][1]->SetMarkerColor(kOrange);   orange_vs_laser[e][6][1]->SetLineColor(kOrange); 
	//	orange_vs_laser[e][6][1]->SetFillColor(kOrange);
	//	orange_vs_laser[e][7][0]->SetMarkerColor(kRed);      orange_vs_laser[e][7][0]->SetLineColor(kRed); 
	//	orange_vs_laser[e][7][0]->SetFillColor(kRed);
	//	orange_vs_laser[e][7][1]->SetMarkerColor(kRed);      orange_vs_laser[e][7][1]->SetLineColor(kRed); 
	//	orange_vs_laser[e][7][1]->SetFillColor(kRed);
	}
	for(int p = 0; p<=period; ++p){
		if(p==1||p==2||p==6||p==7) continue;//temporary
		Legend2->AddEntry(orange_vs_laser[0][p][0], (string_per_leg[p]).c_str(), "f");
	}
	for(int e = 0; e<8; ++e){
		for(int nn = 0; nn<2; ++nn){
			TString outputfile, outputname;
			     if(nn==0) { Legend2->SetHeader("russian crystals");   outputfile  = "_Overlay_Russian"; }
			else if(nn==1) { Legend2->SetHeader("chinese crystals");   outputfile  = "_Overlay_Chinese"; }
			col->cd();

			outputname = "Intercept_vs_VPTdiffOrangeLED" + outputfile + string_aeta[e];
			haxis->SetTitle(outputname);
			haxis->Draw();
			for(int p = 0; p<=period; ++p){
				if(p==1||p==2||p==6||p==7) continue;//temporary
				if(orange_vs_laser[e][p][nn]->GetN()>0) orange_vs_laser[e][p][nn]->Draw("P");
			}
			Legend2->Draw("same");
			col->Update();
			Util::PrintNoEPS(col, outputname, outputdir);
			Util::PrintEPS(col, outputname, outputdir);
			col->Clear();
		}
	}
}
