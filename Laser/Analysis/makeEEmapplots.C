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
//#include "LaserAnalysisReducedTrees.C"
#include "RootMacros/Utilities.hh"

using namespace std;


void makeEEmapplots(){

	bool startdiff         = false;//only true if normalized == false
	bool normalized        = false;

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

	TFile *SCfile          = TFile::Open("/shome/haweber/ECAL/DataLaser/20130118_EESC.root");
	TH2D *scpos_EEp = (TH2D*)SCfile->Get("scpos_EEp");
	TH2D *scpos_EEm = (TH2D*)SCfile->Get("scpos_EEm");
	TH2D *sc_EEp = (TH2D*)SCfile->Get("sc_EEp");
	TH2D *sc_EEm = (TH2D*)SCfile->Get("sc_EEm");
	TH2D *dee_EEp = (TH2D*)SCfile->Get("dee_EEp");
	TH2D *dee_EEm = (TH2D*)SCfile->Get("dee_EEm");
	TH2D *scbarcode_EEp = (TH2D*)SCfile->Get("scbarcode_EEp");
	TH2D *scbarcode_EEm = (TH2D*)SCfile->Get("scbarcode_EEm");


	TH2D *mu_std_EEplus        = new TH2D("mu_std_EEplus ",       "mu_std_EEplus ",       101, 0, 101, 101, 0, 101);
	TH2D *mu_std_EEminus       = new TH2D("mu_std_EEminus",       "mu_std_EEminus",       101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEplus        = new TH2D("mu_SIC_EEplus ",       "mu_SIC_EEplus ",       101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEminus       = new TH2D("mu_SIC_EEminus",       "mu_SIC_EEminus",       101, 0, 101, 101, 0, 101);
	TH2D *VPT_gainqe_EEplus    = new TH2D("VPT_gainqe_EEplus ",   "VPT_gainqe_EEplus ",   101, 0, 101, 101, 0, 101);
	TH2D *VPT_gainqe_EEminus   = new TH2D("VPT_gainqe_EEminus",   "VPT_gainqe_EEminus",   101, 0, 101, 101, 0, 101);
	TH2D *VPT_gain_EEplus      = new TH2D("VPT_gain_EEplus ",     "VPT_gain_EEplus ",     101, 0, 101, 101, 0, 101);
	TH2D *VPT_gain_EEminus     = new TH2D("VPT_gain_EEminus",     "VPT_gain_EEminus",     101, 0, 101, 101, 0, 101);
	TH2D *VPT_qe_EEplus        = new TH2D("VPT_qe_EEplus ",       "VPT_qe_EEplus ",       101, 0, 101, 101, 0, 101);
	TH2D *VPT_qe_EEminus       = new TH2D("VPT_qe_EEminus",       "VPT_qe_EEminus",       101, 0, 101, 101, 0, 101);
	TH2D *EEplus_producer     = new TH2D("EEplus_producer",     "EEplus_producer",     101, 0, 101, 101, 0, 101);
	TH2D *EEminus_producer     = new TH2D("EEminus_producer",     "EEminus_producer",     101, 0, 101, 101, 0, 101);

	TH2D *burnin_EEminus       = new TH2D("burnin_EEminus",       "burnin_EEminus",       101, 0, 101, 101, 0, 101);
	TH2D *burnin_EEplus        = new TH2D("burnin_EEplus ",       "burnin_EEplus ",       101, 0, 101, 101, 0, 101);
	TH2D *eta_EEminus       = new TH2D("eta_EEminus",       "eta_EEminus",       101, 0, 101, 101, 0, 101);
	TH2D *eta_EEplus        = new TH2D("eta_EEplus ",       "eta_EEplus ",       101, 0, 101, 101, 0, 101);

	TH2D *scpos_EEplus  = new TH2D("scpos_EEplus ", "scpos_EEplus ", 101, 0, 101, 101, 0, 101);
	TH2D *scpos_EEminus = new TH2D("scpos_EEminus", "scpos_EEminus", 101, 0, 101, 101, 0, 101);
	TH2D *sc_EEplus  = new TH2D("sc_EEplus ", "sc_EEplus ", 101, 0, 101, 101, 0, 101);
	TH2D *sc_EEminus = new TH2D("sc_EEminus", "sc_EEminus", 101, 0, 101, 101, 0, 101);
	TH2D *dee_EEplus  = new TH2D("dee_EEplus ", "dee_EEplus ", 101, 0, 101, 101, 0, 101);
	TH2D *dee_EEminus = new TH2D("dee_EEminus", "dee_EEminus", 101, 0, 101, 101, 0, 101);
	TH2D *scbarcode_EEplus  = new TH2D("scbarcode_EEplus ", "scbarcode_EEplus ", 101, 0, 101, 101, 0, 101);
	TH2D *scbarcode_EEminus = new TH2D("scbarcode_EEminus", "scbarcode_EEminus", 101, 0, 101, 101, 0, 101);
	TH2D *scbarcode_EEplus2  = new TH2D("scbarcode_EEplus2 ", "scbarcode_EEplus2 ", 101, 0, 101, 101, 0, 101);
	TH2D *scbarcode_EEminus2 = new TH2D("scbarcode_EEminus2", "scbarcode_EEminus2", 101, 0, 101, 101, 0, 101);

	double mustdmax(-10e16), musicmax(-10e16), burninmax(-10e16), deemax(4), scmax(-10e16), scposmax(-10e16), scbarmax(-10e16);
	double mustdmin( 10e16), musicmin( 10e16), burninmin( 10e16), deemin(1), scmin( 10e16), scposmin( 10e16), scbarmin( 10e16);
	for(int i=0; i<=101; ++i){
	for(int j=0; j<=101; ++j){
		int bin = eta_EEm->FindBin(i,j);
		double etap = eta_EEp->GetBinContent(bin);
		double etam = eta_EEm->GetBinContent(bin);
		int prodp = EEp_producer->GetBinContent(bin);
		int prodm = EEm_producer->GetBinContent(bin);
		if(prodp==2 && fabs(etap)>=1.8 && fabs(etap)<=2.0){
			cout << "crystal ix:iy:iz " << i<<":"<<j<<":"<<1<< " at eta = " << etap << " has following values:" << endl;
			double mustd = mu_std_EEp->GetBinContent(bin);
			double music = mu_SIC_EEp->GetBinContent(bin);
			double burnin= burnin_EEp->GetBinContent(bin);
			int dee      = dee_EEp->GetBinContent(bin);
			double sc    = sc_EEp->GetBinContent(bin);
			double scpos = scpos_EEp->GetBinContent(bin);
			double scbar = scbarcode_EEp->GetBinContent(bin)-33201000016000;
			cout << "mustd " << mustd << " music " << music << " burnin " << burnin << " dee " << dee << endl;
			cout << "sc " << sc << "("<<(Long64_t)sc<<") sc_pos " << scpos <<"("<<(Long64_t)scpos<<") sc_bar " << scbar << "("<<(Long64_t)scbar<<")"<<endl;
			if(mustd>mustdmax) mustdmax = mustd; if(mustd<mustdmin) mustdmin = mustd;
			if(music>musicmax) musicmax = music; if(music<musicmin) musicmin = music;
			if(burnin>burninmax) burninmax = burnin; if(burnin<burninmin) burninmin = burnin;
			if(sc>scmax) scmax = sc; if(sc<scmin) scmin = sc;
			if(scpos>scposmax) scposmax = scpos; if(scpos<scposmin) scposmin = scpos;
			if(scbar>scbarmax) scbarmax = scbar; if(scbar<scbarmin) scbarmin = scbar;
			if(mustd>0) mu_std_EEplus->SetBinContent(i,j,mustd);
			mu_SIC_EEplus->SetBinContent(i,j,music);
			burnin_EEplus->SetBinContent(i,j,burnin);
			dee_EEplus->SetBinContent(i,j,dee);
			sc_EEplus->SetBinContent(i,j,sc);
			scpos_EEplus->SetBinContent(i,j,scpos);
			scbarcode_EEplus->SetBinContent(i,j,scbar);
			if(mustd>0) scbarcode_EEplus2->SetBinContent(i,j,scbar);
		}
		if(prodm==2 && fabs(etam)>=1.8 && fabs(etam)<=2.0){
			cout << "crystal ix:iy:iz " << i<<":"<<j<<":"<<-1<< " at eta = " << etap << " has following values:" << endl;
			double mustd = mu_std_EEm->GetBinContent(bin);
			double music = mu_SIC_EEm->GetBinContent(bin);
			double burnin= burnin_EEm->GetBinContent(bin);
			int dee      = dee_EEm->GetBinContent(bin);
			double sc    = sc_EEm->GetBinContent(bin);
			double scpos = scpos_EEm->GetBinContent(bin);
			double scbar = scbarcode_EEm->GetBinContent(bin)-33201000016000;
			cout << "mustd " << mustd << " music " << music << " burnin " << burnin << " dee " << dee << endl;
			cout << "sc " << sc << "("<<(Long64_t)sc<<") sc_pos " << scpos <<"("<<(Long64_t)scpos<<") sc_bar " << scbar << "("<<(Long64_t)scbar<<")"<<endl;
			if(mustd>mustdmax) mustdmax = mustd; if(mustd<mustdmin) mustdmin = mustd;
			if(music>musicmax) musicmax = music; if(music<musicmin) musicmin = music;
			if(burnin>burninmax) burninmax = burnin; if(burnin<burninmin) burninmin = burnin;
			if(sc>scmax) scmax = sc; if(sc<scmin) scmin = sc;
			if(scpos>scposmax) scposmax = scpos; if(scpos<scposmin) scposmin = scpos;
			if(scbar>scbarmax) scbarmax = scbar; if(scbar<scbarmin) scbarmin = scbar;
			if(mustd>0) mu_std_EEminus->SetBinContent(i,j,mustd);
			mu_SIC_EEminus->SetBinContent(i,j,music);
			burnin_EEminus->SetBinContent(i,j,burnin);
			dee_EEminus->SetBinContent(i,j,dee);
			sc_EEminus->SetBinContent(i,j,sc);
			scpos_EEminus->SetBinContent(i,j,scpos);
			scbarcode_EEminus->SetBinContent(i,j,scbar);
			if(mustd>0) scbarcode_EEminus2->SetBinContent(i,j,scbar);
		}
	}}
	

	mu_std_EEplus        ->SetMaximum(mustdmax );
	mu_std_EEminus       ->SetMaximum(mustdmax );
	mu_SIC_EEplus        ->SetMaximum(musicmax );
	mu_SIC_EEminus       ->SetMaximum(musicmax );
	burnin_EEminus       ->SetMaximum(burninmax);
	burnin_EEplus        ->SetMaximum(burninmax);
	scpos_EEplus         ->SetMaximum(scposmax );
	scpos_EEminus        ->SetMaximum(scposmax );
	sc_EEplus            ->SetMaximum(scmax    );
	sc_EEminus           ->SetMaximum(scmax    );
	dee_EEplus           ->SetMaximum(deemax   );
	dee_EEminus          ->SetMaximum(deemax   );
	scbarcode_EEplus     ->SetMaximum(scbarmax );
	scbarcode_EEminus    ->SetMaximum(scbarmax );
	scbarcode_EEplus2    ->SetMaximum(scbarmax );
	scbarcode_EEminus2   ->SetMaximum(scbarmax );

	mu_std_EEplus        ->SetMinimum(mustdmin );
	mu_std_EEminus       ->SetMinimum(mustdmin );
	mu_SIC_EEplus        ->SetMinimum(musicmin );
	mu_SIC_EEminus       ->SetMinimum(musicmin );
	burnin_EEminus       ->SetMinimum(burninmin);
	burnin_EEplus        ->SetMinimum(burninmin);
	scpos_EEplus         ->SetMinimum(scposmin );
	scpos_EEminus        ->SetMinimum(scposmin );
	sc_EEplus            ->SetMinimum(scmin    );
	sc_EEminus           ->SetMinimum(scmin    );
	dee_EEplus           ->SetMinimum(deemin   );
	dee_EEminus          ->SetMinimum(deemin   );
	scbarcode_EEplus     ->SetMinimum(scbarmin );
	scbarcode_EEminus    ->SetMinimum(scbarmin );
	scbarcode_EEplus2    ->SetMinimum(scbarmin );
	scbarcode_EEminus2   ->SetMinimum(scbarmin );

	TCanvas *c1 = new TCanvas(); c1->cd();
	mu_std_EEminus->Draw("COLZ");
	TCanvas *c2 = new TCanvas(); c2->cd();
	mu_SIC_EEminus->Draw("COLZ");
	TCanvas *c3 = new TCanvas(); c3->cd();
	burnin_EEminus->Draw("COLZ");
	TCanvas *c4 = new TCanvas(); c4->cd();
	scpos_EEminus->Draw("COLZ");
	TCanvas *c5 = new TCanvas(); c5->cd();
	sc_EEminus->Draw("COLZ");
	TCanvas *c6 = new TCanvas(); c6->cd();
	dee_EEminus->Draw("COLZ");
	TCanvas *c7 = new TCanvas(); c7->cd();
	scbarcode_EEminus->Draw("COLZ");
	TCanvas *c8 = new TCanvas(); c8->cd();
	scbarcode_EEminus2->Draw("COLZ");


//	TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
/*
	mu_std_EEplus        
	mu_std_EEminus       
	mu_SIC_EEplus        
	mu_SIC_EEminus       
	burnin_EEminus       
	burnin_EEplus        
	scpos_EEplus         
	scpos_EEminus        
	sc_EEplus            
	sc_EEminus           
	dee_EEplus           
	dee_EEminus          
	scbarcode_EEplus     
	scbarcode_EEminus    
*/
}


