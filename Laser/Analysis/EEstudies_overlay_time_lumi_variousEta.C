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


//overlay slopes/intercept and potentially different other variables
// for various eta bins

void EEstudies_overlay_time_lumi_variousEta(){

	bool startdiff         = false;//only true if normalized == false, always diff w.r.t. 2011 beginning
	bool normalized        = false;
	const int period       = 7;//0: original, 1: 15/03-27/03, 2: 12/04-05/05, 3: 13/05-01/07, 4: 15/07-25/08, 5: 04/09-02/11
				//2011: 5 periods: 2011+2012: 7 periods
	bool samelumi          = true;//if same lumi = true, use only last three periods with first 1.1fb-1
	int  whatoverlay       = 1; // 1: integrated Lumi 2: time
	bool only2011          = false;


	bool plotscaled        = false;//save also plots where VPT/PN has been scaled

	TFile *oldFile;
	if(only2011){
	if(whatoverlay==1){
	if(samelumi){
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin_only2011_sameLumi1p1fb-1.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized_only2011_sameLumi1p1fb-1.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_sameLumi1p1fb-1.root");
	} else{
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_diffwrt2011begin.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011_normalized.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_only2011.root");
	}
	} else if(whatoverlay==2){
	if(samelumi){
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_diffwrt2011begin_only2011_sameLumi1p1fb-1.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_normalized_only2011_sameLumi1p1fb-1.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_only2011_sameLumi1p1fb-1.root");
	}
	else{
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_only2011_diffwrt2011begin.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_only2011_normalized.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_only2011.root");
	}
	}
	} else{
	if(whatoverlay==1){
	if(samelumi){
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin_sameLumi950pb-1.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized_sameLumi950pb-1.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_sameLumi950pb-1.root");
	} else{
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_diffwrt2011begin.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples_normalized.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_integrLumi_VectorNTuples.root");
	}
	} else if(whatoverlay==2){
	if(samelumi){
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_diffwrt2011begin_sameLumi950pb-1.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_normalized_sameLumi950pb-1.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_sameLumi950pb-1.root");
	}
	else{
	if(startdiff)  oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_diffwrt2011begin.root");
	if(normalized) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples_normalized.root");
	else if(!startdiff) oldFile = TFile::Open("FileOutputs/20121023_intercept_slopes_VS_time_VectorNTuples.root");
	}
	}
	}

	TString outputdir;
	outputdir = "Plots/20121023/EEstudies/continued/diffwrt2011begin/InterceptSlopesVSintegrLumi";
	if(normalized) outputdir = "Plots/20121023/EEstudies/continued/normalized/InterceptSlopesVSintegrLumi";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
	else if(!startdiff) outputdir = "Plots/20121023/EEstudies/continued/notnormalized/InterceptSlopesVSintegrLumi";
	if(whatoverlay==2){
		outputdir = "Plots/20121023/EEstudies/continued/diffwrt2011begin/InterceptSlopesVStime";
		if(normalized) outputdir = "Plots/20121023/EEstudies/continued/normalized/InterceptSlopesVStime";//normalized if normalize VPT/PN // continued/normalized <-> continued/SIC/normalized
		else if(!startdiff) outputdir = "Plots/20121023/EEstudies/continued/notnormalized/InterceptSlopesVStime";
	}
	if(samelumi &&  only2011) outputdir = outputdir + (TString)"_sameLumi1p1fb-1";
	if(samelumi && !only2011) outputdir = outputdir + (TString)"_sameLumi950pb-1";

	Util::MakeOutputDir(outputdir);

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	string string_prod[2] = {"__russ", "__chin"};
	string string_aeta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_acap[1] = {"_"};

	string string_eta_leg[8] = {"1.4 #leq |#eta| < 1.6", "1.6 #leq |#eta| < 1.8", "1.8 #leq |#eta| < 2.0", "2.0 #leq |#eta| < 2.2", "2.2 #leq |#eta| < 2.4", "2.4 #leq |#eta| < 2.6", "2.6 #leq |#eta| < 2.8", "2.8 #leq |#eta|"};


	string aname;

//	const int Ntimebins = 11;
//	double timebins[Ntimebins+1] = {0., 45., 77., 78., 337., 338., 1444., 1445., 2924., 2925., 6137., 6138.};
//	if(samelumi){ timebins[8] = 2493.; timebins[10] = 3986.; }
//	if(whatoverlay==2){
//		timebins[0] = 1293840000.; timebins[1] = 1300200000.; timebins[2] = 1301210000.; timebins[3] = 1302590000.; timebins[4] = 1304620000.; timebins[5] = 1305260000.; timebins[6] = 1309540000.; timebins[7] = 1310710000.; timebins[8] = 1314270000.; timebins[9] = 1315120000.; timebins[10] = 1320250000.; timebins[11] = 1325376000.; 
//		if(samelumi){
//			timebins[8] = 1312800923.; timebins[10] = 1316813137.;
//		}
//	}
	const int Ntimebins = 15;
	double timebins[Ntimebins+1] = {0., 45., 77., 78., 337., 338., 1444., 1445., 2924., 2925., 6137., 6177.,7080.,7178.,12777.,15000.};
	if(only2011){//last ones are dummys for 2011
		timebins[11] = 6138.; timebins[12] = 6139.; timebins[13] = 6140.; timebins[14] = 6141.; timebins[15] = 7500.; }
	if( only2011 && samelumi){ timebins[8] = 2538.; timebins[10] = 3986.; }
	if(!only2011 && samelumi){ timebins[6] = 1300.; timebins[8] = 2400.; timebins[10] = 3875.; timebins[14] = 8130.;}
	if(whatoverlay==2){
		timebins[0] = 1293840000.; timebins[1] = 1300200000.; timebins[2] = 1301210000.; timebins[3] = 1302590000.; timebins[4] = 1304620000.; timebins[5] = 1305260000.; timebins[6] = 1309540000.; timebins[7] = 1310710000.; timebins[8] = 1314270000.; timebins[9] = 1315120000.; timebins[10] = 1320250000.; timebins[11] = 1333740000.; timebins[12] = 1334910000.; timebins[13] = 1336230000.; timebins[14] = 1340000000.; timebins[15] = 1345000000.;
		if(only2011){
		timebins[11] = 1321000000.; timebins[12] = 1322000000.; timebins[13] = 1323000000.; timebins[14] = 1324000000.; timebins[15] = 1325376000.;
		}
		if(samelumi &&  only2011){ timebins[8] = 1312800923.; timebins[10] = 1316813137.; }
		if(samelumi && !only2011){ timebins[6] = 1308900000.; timebins[8] = 1312720000.; timebins[10] = 1316695000.; timebins[14] = 1337260000.; }
	}

	TH1D *haxis = new TH1D("haxis", "", Ntimebins, timebins);
	haxis->SetMinimum(0.);
	if(     whatoverlay==1) haxis->GetXaxis()->SetTitle("L_{int} [pb^{-1}]");
	else if(whatoverlay==2) haxis->GetXaxis()->SetTitle("time");
	TLegend *Legend1 = new TLegend(.68,.5,.88,.88);
	Legend1 -> SetFillColor(0);
	Legend1 -> SetBorderSize(0);
   	TH1F* h[8]; 

	TH1D * h_intercept_vs_time[8][2];
	TH1D * h_intercept_vs_time_scaled[8][2];
	TH1D * h_slope_vs_time[8][2];
	TH1D * h_slope_vs_time_scaled[8][2];
	for(int e = 0; e<8; ++e){
	   h[e] = new TH1F("", (string_eta_leg[e]).c_str(), 1, 0, 1); 
	   for(int nn = 0; nn<2; ++nn){
		aname = (string)"intercept_vs_integrLumi_" + string_aeta[e] + string_prod[nn];
		if(whatoverlay==2) aname = (string)"intercept_vs_time_" + string_aeta[e] + string_prod[nn];
		h_intercept_vs_time[e][nn] = (TH1D*)oldFile->Get(aname.c_str());
		aname = (string)"intercept_vs_integrLumi_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
		if(whatoverlay==2) aname = (string)"intercept_vs_time_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
		h_intercept_vs_time_scaled[e][nn] = (TH1D*)oldFile->Get(aname.c_str());
		aname = (string)"slope_vs_integrLumi_" + string_aeta[e] + string_prod[nn];
		if(whatoverlay==2) aname = (string)"slope_vs_time_" + string_aeta[e] + string_prod[nn];
		h_slope_vs_time[e][nn] = (TH1D*)oldFile->Get(aname.c_str());
		aname = (string)"slope_vs_integrLumi_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
		if(whatoverlay==2) aname = (string)"slope_vs_time_" + string_aeta[e] + (string)"_scaled" + string_prod[nn];
		h_slope_vs_time_scaled[e][nn] = (TH1D*)oldFile->Get(aname.c_str());
		if(e==0) { h_intercept_vs_time[e][nn]       ->SetLineColor(kViolet);   h_intercept_vs_time[e][nn]       ->SetMarkerColor(kViolet);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kViolet);   h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kViolet);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kViolet);   h_slope_vs_time[e][nn]           ->SetMarkerColor(kViolet);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kViolet);   h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kViolet);
			   h[e]                             ->SetFillColor(kViolet);   }
		if(e==1) { h_intercept_vs_time[e][nn]       ->SetLineColor(kBlue);     h_intercept_vs_time[e][nn]       ->SetMarkerColor(kBlue);  
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kBlue);     h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kBlue);  
			   h_slope_vs_time[e][nn]           ->SetLineColor(kBlue);     h_slope_vs_time[e][nn]           ->SetMarkerColor(kBlue);  
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kBlue);     h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kBlue);  
			   h[e]                             ->SetFillColor(kBlue);     }
		if(e==2) { h_intercept_vs_time[e][nn]       ->SetLineColor(kCyan+2);   h_intercept_vs_time[e][nn]       ->SetMarkerColor(kCyan+2);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kCyan+2);   h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kCyan+2);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kCyan+2);   h_slope_vs_time[e][nn]           ->SetMarkerColor(kCyan+2);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kCyan+2);   h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kCyan+2);
			   h[e]                             ->SetFillColor(kCyan+2);   }
		if(e==3) { h_intercept_vs_time[e][nn]       ->SetLineColor(kGreen+2);  h_intercept_vs_time[e][nn]       ->SetMarkerColor(kGreen+2);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kGreen+2);  h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kGreen+2);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kGreen+2);  h_slope_vs_time[e][nn]           ->SetMarkerColor(kGreen+2);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kGreen+2);  h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kGreen+2);
			   h[e]                             ->SetFillColor(kGreen+2);  }
		if(e==4) { h_intercept_vs_time[e][nn]       ->SetLineColor(kYellow+3); h_intercept_vs_time[e][nn]       ->SetMarkerColor(kYellow+3);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kYellow+3); h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kYellow+3);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kYellow+3); h_slope_vs_time[e][nn]           ->SetMarkerColor(kYellow+3);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kYellow+3); h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kYellow+3);
			   h[e]                             ->SetFillColor(kYellow+3); }
		if(e==5) { h_intercept_vs_time[e][nn]       ->SetLineColor(kOrange+5); h_intercept_vs_time[e][nn]       ->SetMarkerColor(kOrange+5);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kOrange+5); h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kOrange+5);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kOrange+5); h_slope_vs_time[e][nn]           ->SetMarkerColor(kOrange+5);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kOrange+5); h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kOrange+5);
			   h[e]                             ->SetFillColor(kOrange+5); }
		if(e==6) { h_intercept_vs_time[e][nn]       ->SetLineColor(kOrange);   h_intercept_vs_time[e][nn]       ->SetMarkerColor(kOrange);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kOrange);   h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kOrange);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kOrange);   h_slope_vs_time[e][nn]           ->SetMarkerColor(kOrange);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kOrange);   h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kOrange);
			   h[e]                             ->SetFillColor(kOrange);   }
		if(e==7) { h_intercept_vs_time[e][nn]       ->SetLineColor(kRed);      h_intercept_vs_time[e][nn]       ->SetMarkerColor(kRed);
			   h_intercept_vs_time_scaled[e][nn]->SetLineColor(kRed);      h_intercept_vs_time_scaled[e][nn]->SetMarkerColor(kRed);
			   h_slope_vs_time[e][nn]           ->SetLineColor(kRed);      h_slope_vs_time[e][nn]           ->SetMarkerColor(kRed);
			   h_slope_vs_time_scaled[e][nn]    ->SetLineColor(kRed);      h_slope_vs_time_scaled[e][nn]    ->SetMarkerColor(kRed);
			   h[e]                             ->SetFillColor(kRed);      }
	   }
	}

	for(int e = 0; e<8; ++e){
		Legend1->AddEntry(h[e], h[e]->GetTitle(), "f");
	}

	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
	for(int nn = 0; nn<2; ++nn){
		TString outputfile, outputname;
		     if(nn==0) { Legend1->SetHeader("russian crystals"); outputfile  = "_Overlay_Russian"; }
		else if(nn==1) { Legend1->SetHeader("chinese crystals");   outputfile  = "_Overlay_Chinese"; }
		col->cd();
		if(only2011){
			if(normalized)     haxis->SetMaximum(0.25);
			else if(startdiff) haxis->SetMaximum(0.4);
			else               haxis->SetMaximum(0.185);
		} else{
			if(normalized)     haxis->SetMaximum(0.27);
			else if(startdiff) haxis->SetMaximum(0.55);
			else               haxis->SetMaximum(0.22);
		}
		outputname = "Intercept" + outputfile;
		haxis->SetTitle(outputname);
		haxis->GetYaxis()->SetTitle("intercept");
		haxis->Draw();
		for(int e = 0; e<8; ++e){
			h_intercept_vs_time[e][nn]->Draw("same");
		}
		Legend1->Draw("same");
		col->Update();
		Util::PrintNoEPS(col, outputname, outputdir);
		Util::PrintEPS(col, outputname, outputdir);
		col->Clear();
		if(plotscaled){
		col->cd();
		if(only2011){
			if(normalized)     haxis->SetMaximum(6.);
			else if(startdiff) haxis->SetMaximum(5.);
			else               haxis->SetMaximum(7.);
		} else{
			if(normalized)     haxis->SetMaximum(6.);
			else if(startdiff) haxis->SetMaximum(5.5);
			else               haxis->SetMaximum(7.);
		}
		outputname = "InterceptScaled" + outputfile;
		haxis->SetTitle(outputname);
		haxis->GetYaxis()->SetTitle("intercept (scaled)");
		haxis->Draw();
		for(int e = 0; e<8; ++e){
			h_intercept_vs_time_scaled[e][nn]->Draw("same");
		}
		Legend1->Draw("same");
		col->Update();
		Util::PrintNoEPS(col, outputname, outputdir);
		Util::PrintEPS(col, outputname, outputdir);
		col->Clear();
		}
		col->cd();
		if(only2011){
			if(normalized)     haxis->SetMaximum(0.16);
			else if(startdiff) haxis->SetMaximum(0.2);
			else               haxis->SetMaximum(0.1);
		} else{
			if(normalized)     haxis->SetMaximum(0.17);
			else if(startdiff) haxis->SetMaximum(0.24);
			else               haxis->SetMaximum(0.1);
		}
		outputname = "Slope" + outputfile;
		haxis->SetTitle(outputname);
		haxis->GetYaxis()->SetTitle("slope");
		haxis->Draw();
		for(int e = 0; e<8; ++e){
			h_slope_vs_time[e][nn]->Draw("same");
		}
		Legend1->Draw("same");
		col->Update();
		Util::PrintNoEPS(col, outputname, outputdir);
		Util::PrintEPS(col, outputname, outputdir);
		col->Clear();
		col->cd();
		if(plotscaled){
		if(only2011){
			if(normalized)     haxis->SetMaximum(10.);
			else if(startdiff) haxis->SetMaximum(15.);
			else               haxis->SetMaximum(10.);
		} else{
			if(normalized)     haxis->SetMaximum(10.);
			else if(startdiff) haxis->SetMaximum(15.);
			else               haxis->SetMaximum(10.);
		}
		outputname = "SlopeScaled" + outputfile;
		haxis->SetTitle(outputname);
		haxis->GetYaxis()->SetTitle("slope (scaled)");
		haxis->Draw();
		for(int e = 0; e<8; ++e){
			h_slope_vs_time_scaled[e][nn]->Draw("same");
		}
		Legend1->Draw("same");
		col->Update();
		Util::PrintNoEPS(col, outputname, outputdir);
		Util::PrintEPS(col, outputname, outputdir);
		col->Clear();
		}
	}

}