#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "RootMacros/Utilities.hh"

using namespace std;

void PlotEtaProfile();
void Plot_EE_mu_maps_2();

void Plot_EE_mu_maps_2(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

//	TFile *outf = new TFile("/shome/haweber/ECAL/Laser/Analysis/Plots/20121126_EEmustd.root", "RECREATE");
	TFile *outf = new TFile("/shome/haweber/ECAL/DataLaser/20121116_EEmustd.root", "RECREATE");

	TFile *file = TFile::Open("/shome/haweber/ECAL/DataLaser/ana_rad_ecal_v3.root");
	TTree *tree = (TTree*)file->Get("EE");

	TH2D *mu_ECAL_EEp_russ = new TH2D("mu_ECAL_EEp_russ", "mu_ECAL_EEp_russ", 101, 0, 101, 101, 0, 101);
	TH2D *mu_ECAL_EEp_chin = new TH2D("mu_ECAL_EEp_chin", "mu_ECAL_EEp_chin", 101, 0, 101, 101, 0, 101);
	TH2D *mu_ECAL_EEm_russ = new TH2D("mu_ECAL_EEm_russ", "mu_ECAL_EEm_russ", 101, 0, 101, 101, 0, 101);
	TH2D *mu_ECAL_EEm_chin = new TH2D("mu_ECAL_EEm_chin", "mu_ECAL_EEm_chin", 101, 0, 101, 101, 0, 101);

	TH2D *mu_SIC_EEp_russ  = new TH2D("mu_SIC_EEp_russ",  "mu_SIC_EEp_russ",  101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEp_chin  = new TH2D("mu_SIC_EEp_chin",  "mu_SIC_EEp_chin",  101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEm_russ  = new TH2D("mu_SIC_EEm_russ",  "mu_SIC_EEm_russ",  101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEm_chin  = new TH2D("mu_SIC_EEm_chin",  "mu_SIC_EEm_chin",  101, 0, 101, 101, 0, 101);


	TH2D *mu_std_EEp       = new TH2D("mu_std_EEp",       "mu_std_EEp",       101, 0, 101, 101, 0, 101);
	TH2D *mu_std_EEm       = new TH2D("mu_std_EEm",       "mu_std_EEm",       101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEp       = new TH2D("mu_SIC_EEp",       "mu_SIC_EEp",       101, 0, 101, 101, 0, 101);
	TH2D *mu_SIC_EEm       = new TH2D("mu_SIC_EEm",       "mu_SIC_EEm",       101, 0, 101, 101, 0, 101);

	TH2D *VPT_gainqe_EEp   = new TH2D("VPT_gainqe_EEp",   "VPT_gainqe_EEp",   101, 0, 101, 101, 0, 101);
	TH2D *VPT_gainqe_EEm   = new TH2D("VPT_gainqe_EEm",   "VPT_gainqe_EEm",   101, 0, 101, 101, 0, 101);
	TH2D *VPT_gain_EEp     = new TH2D("VPT_gain_EEp",     "VPT_gain_EEp",     101, 0, 101, 101, 0, 101);
	TH2D *VPT_gain_EEm     = new TH2D("VPT_gain_EEm",     "VPT_gain_EEm",     101, 0, 101, 101, 0, 101);
	TH2D *VPT_qe_EEp       = new TH2D("VPT_qe_EEp",       "VPT_qe_EEp",       101, 0, 101, 101, 0, 101);
	TH2D *VPT_qe_EEm       = new TH2D("VPT_qe_EEm",       "VPT_qe_EEm",       101, 0, 101, 101, 0, 101);

	TH2D *EEp_producer     = new TH2D("EEp_producer",     "EEp_producer",     101, 0, 101, 101, 0, 101);
	TH2D *EEm_producer     = new TH2D("EEm_producer",     "EEm_producer",     101, 0, 101, 101, 0, 101);

	int nev;
	TString variable, var;
	TString selection, sel;

   	Long64_t nentries = tree->GetEntriesFast();
	cout << "Tree contains "<< nentries << " entries" << endl;

	int ix,iy,iz;
	double mu_ECAL, mu_SIC;
	float vpt_gain, vpt_qe, vpt_pg;
	int producer;
	TBranch *xx     = tree->GetBranch("ix");
	TBranch *yy     = tree->GetBranch("iy");
	TBranch *zz     = tree->GetBranch("iz");
	TBranch *muECAL = tree->GetBranch("mu_ECAL");
	TBranch *muSIC  = tree->GetBranch("mu_SIC");
	TBranch *prod   = tree->GetBranch("producer");
	TBranch *gain   = tree->GetBranch("vpt_gain");
	TBranch *qe     = tree->GetBranch("vpt_qe");
	TBranch *pg     = tree->GetBranch("vpt_pg");
	xx    ->SetAddress(&ix);
	yy    ->SetAddress(&iy);
	zz    ->SetAddress(&iz);
	muECAL->SetAddress(&mu_ECAL);
	muSIC ->SetAddress(&mu_SIC);
	prod  ->SetAddress(&producer);
	gain  ->SetAddress(&vpt_gain);
	qe    ->SetAddress(&vpt_qe);
	pg    ->SetAddress(&vpt_pg);

   	Long64_t nbytes = 0, nb = 0;
   	    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      		//Long64_t ientry = tree->LoadTree(jentry);
      		//if (ientry < 0) break;
		//if(jentry%(nentries/25)==0) cout << "Processing event " << jentry << endl;
		nb = tree->GetEntry(jentry);   nbytes += nb;

	//	int iX = ix;
	//	int iY = iy;
	//	int iZ = iz;
	//	if(mu_ECAL>0){
	//	cout << "ix " << iX << " iy " << iY << " iz " << iZ << endl;
	//	cout << "mustd " << mu_ECAL << " muSIC " << mu_SIC << " producer " << producer << endl;
	//	}
		int bin = mu_ECAL_EEp_russ->FindBin(ix,iy);

		if(mu_ECAL>=0 && iz==1 && producer==1){//EE+, russian
			if(mu_ECAL_EEp_russ->GetBinContent(bin)>0) cout << "mu_ECAL: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_ECAL_EEp_russ->GetBinContent(bin) << " and now is filled with " << mu_ECAL << endl;
			mu_ECAL_EEp_russ->SetBinContent(bin, mu_ECAL);
		}
		if(mu_ECAL>=0 && iz==-1 && producer==1){//EE-, russian
			if(mu_ECAL_EEm_russ->GetBinContent(bin)>0) cout << "mu_ECAL: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_ECAL_EEm_russ->GetBinContent(bin) << " and now is filled with " << mu_ECAL << endl;
			mu_ECAL_EEm_russ->SetBinContent(bin, mu_ECAL);
		}
		if(mu_ECAL>=0 && iz==1 && producer==2){//EE+, chinese
			if(mu_ECAL_EEp_chin->GetBinContent(bin)>0) cout << "mu_ECAL: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_ECAL_EEp_chin->GetBinContent(bin) << " and now is filled with " << mu_ECAL << endl;
			mu_ECAL_EEp_chin->SetBinContent(bin, mu_ECAL);
		}
		if(mu_ECAL>=0 && iz==-1 && producer==2){//EE-, chinese
			if(mu_ECAL_EEm_chin->GetBinContent(bin)>0) cout << "mu_ECAL: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_ECAL_EEm_chin->GetBinContent(bin) << " and now is filled with " << mu_ECAL << endl;
			mu_ECAL_EEm_chin->SetBinContent(bin, mu_ECAL);
		}
		if(mu_SIC>=0 && iz==1 && producer==1){//EE+, russian
			if(mu_SIC_EEp_russ->GetBinContent(bin)>0) cout << "mu_SIC: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_SIC_EEp_russ->GetBinContent(bin) << " and now is filled with " << mu_SIC << endl;
			mu_SIC_EEp_russ->SetBinContent(bin, mu_SIC);
		}
		if(mu_SIC>=0 && iz==-1 && producer==1){//EE-, russian
			if(mu_SIC_EEm_russ->GetBinContent(bin)>0) cout << "mu_SIC: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_SIC_EEm_russ->GetBinContent(bin) << " and now is filled with " << mu_SIC << endl;
			mu_SIC_EEm_russ->SetBinContent(bin, mu_SIC);
		}
		if(mu_SIC>=0 && iz==1 && producer==2){//EE+, chinese
			if(mu_SIC_EEp_chin->GetBinContent(bin)>0) cout << "mu_SIC: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_SIC_EEp_chin->GetBinContent(bin) << " and now is filled with " << mu_SIC << endl;
			mu_SIC_EEp_chin->SetBinContent(bin, mu_SIC);
		}
		if(mu_SIC>=0 && iz==-1 && producer==2){//EE-, chinese
			if(mu_SIC_EEm_chin->GetBinContent(bin)>0) cout << "mu_SIC: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << mu_SIC_EEm_chin->GetBinContent(bin) << " and now is filled with " << mu_SIC << endl;
			mu_SIC_EEm_chin->SetBinContent(bin, mu_SIC);
		}
		if(iz==1 && (producer==1||producer==2)){
			if(EEp_producer->GetBinContent(bin)>0) cout << "producer: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << EEp_producer->GetBinContent(bin) << " and now is filled with " << producer << endl;
			EEp_producer->SetBinContent(bin, producer);
		}
		if(iz==-1 && (producer==1||producer==2)){
			if(EEm_producer->GetBinContent(bin)>0) cout << "producer: bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << EEm_producer->GetBinContent(bin) << " and now is filled with " << producer << endl;
			EEm_producer->SetBinContent(bin, producer);
		}

		if(vpt_pg>=0 && iz==1){
			if(VPT_gainqe_EEp->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_gainqe_EEp->GetBinContent(bin) << " and now is filled with " << vpt_pg << endl;
			VPT_gainqe_EEp->SetBinContent(bin, vpt_pg);
		}
		if(vpt_pg>=0 && iz==-1){
			if(VPT_gainqe_EEm->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_gainqe_EEm->GetBinContent(bin) << " and now is filled with " << vpt_pg << endl;
			VPT_gainqe_EEm->SetBinContent(bin, vpt_pg);
		}
		if(vpt_qe>=0 && iz==1){
			if(VPT_qe_EEp->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_qe_EEp->GetBinContent(bin) << " and now is filled with " << vpt_qe << endl;
			VPT_qe_EEp->SetBinContent(bin, vpt_qe);
		}
		if(vpt_qe>=0 && iz==-1){
			if(VPT_qe_EEm->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_qe_EEm->GetBinContent(bin) << " and now is filled with " << vpt_qe << endl;
			VPT_qe_EEm->SetBinContent(bin, vpt_qe);
		}
		if(vpt_gain>=0 && iz==1){
			if(VPT_gain_EEp->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_gain_EEp->GetBinContent(bin) << " and now is filled with " << vpt_gain << endl;
			VPT_gain_EEp->SetBinContent(bin, vpt_gain);
		}
		if(vpt_gain>=0 && iz==-1){
			if(VPT_gain_EEm->GetBinContent(bin)>0) cout << __LINE__ << ": bin " << bin << " (ix,iy,iz=" << ix <<","<<iy<<","<<iz<<") was filled with " << VPT_gain_EEm->GetBinContent(bin) << " and now is filled with " << vpt_gain << endl;
			VPT_gain_EEm->SetBinContent(bin, vpt_gain);
		}
	   }

	mu_std_EEp->Add(mu_ECAL_EEp_russ,1);
	mu_std_EEp->Add(mu_ECAL_EEp_chin,1);
	mu_std_EEm->Add(mu_ECAL_EEm_russ,1);
	mu_std_EEm->Add(mu_ECAL_EEm_chin,1);
	mu_SIC_EEp->Add(mu_SIC_EEp_russ, 1);
	mu_SIC_EEp->Add(mu_SIC_EEp_chin, 1);
	mu_SIC_EEm->Add(mu_SIC_EEm_russ, 1);
	mu_SIC_EEm->Add(mu_SIC_EEm_chin, 1);

	mu_std_EEp->SetMaximum(2.); mu_std_EEp->GetZaxis()->SetRangeUser(0.,2.);
	mu_std_EEm->SetMaximum(2.); mu_std_EEm->GetZaxis()->SetRangeUser(0.,2.);
	mu_SIC_EEp->SetMaximum(2.); mu_SIC_EEp->GetZaxis()->SetRangeUser(0.,2.);
	mu_SIC_EEm->SetMaximum(2.); mu_SIC_EEm->GetZaxis()->SetRangeUser(0.,2.);

	TString outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots";

	if(mu_ECAL_EEp_russ->GetEntries()>0){
	TCanvas *can1 = new TCanvas("can1", "mu_ECAL_EEp_russ", 0, 0, 900, 700);
	can1->cd();
	mu_ECAL_EEp_russ->Draw("COLZ");
	//Util::PrintNoEPS(can1, can1->GetTitle(), outputdir);
	Util::PrintEPS(can1, can1->GetTitle(), outputdir);
	} if(mu_ECAL_EEp_chin->GetEntries()>0){
	TCanvas *can2 = new TCanvas("can2", "mu_ECAL_EEp_chin", 0, 0, 900, 700);
	can2->cd();
	mu_ECAL_EEp_chin->Draw("COLZ");
	//Util::PrintNoEPS(can2, can2->GetTitle(), outputdir);
	Util::PrintEPS(can2, can2->GetTitle(), outputdir);
	} if(mu_ECAL_EEm_russ->GetEntries()>0){
	TCanvas *can3 = new TCanvas("can3", "mu_ECAL_EEm_russ", 0, 0, 900, 700);
	can3->cd();
	mu_ECAL_EEm_russ->Draw("COLZ");
	//Util::PrintNoEPS(can3, can3->GetTitle(), outputdir);
	Util::PrintEPS(can3, can3->GetTitle(), outputdir);
	} if(mu_ECAL_EEm_chin->GetEntries()>0){
	TCanvas *can4 = new TCanvas("can4", "mu_ECAL_EEm_chin", 0, 0, 900, 700);
	can4->cd();
	mu_ECAL_EEm_chin->Draw("COLZ");
	//Util::PrintNoEPS(can4, can4->GetTitle(), outputdir);
	Util::PrintEPS(can4, can4->GetTitle(), outputdir);
	} if(mu_SIC_EEp_russ->GetEntries()>0){
	TCanvas *can5 = new TCanvas("can5", "mu_SIC_EEp_russ", 0, 0, 900, 700);
	can5->cd();
	mu_SIC_EEp_russ->Draw("COLZ");
	//Util::PrintNoEPS(can5, can5->GetTitle(), outputdir);
	Util::PrintEPS(can5, can5->GetTitle(), outputdir);
	} if(mu_SIC_EEp_chin->GetEntries()>0){
	TCanvas *can6 = new TCanvas("can6", "mu_SIC_EEp_chin", 0, 0, 900, 700);
	can6->cd();
	mu_SIC_EEp_chin->Draw("COLZ");
	//Util::PrintNoEPS(can6, can6->GetTitle(), outputdir);
	Util::PrintEPS(can6, can6->GetTitle(), outputdir);	
	} if(mu_SIC_EEm_russ->GetEntries()>0){
	TCanvas *can7 = new TCanvas("can7", "mu_SIC_EEm_russ", 0, 0, 900, 700);
	can7->cd();
	mu_SIC_EEm_russ->Draw("COLZ");
	//Util::PrintNoEPS(can7, can7->GetTitle(), outputdir);
	Util::PrintEPS(can7, can7->GetTitle(), outputdir);
	} if(mu_SIC_EEm_chin->GetEntries()>0){
	TCanvas *can8 = new TCanvas("can8", "mu_SIC_EEm_chin", 0, 0, 900, 700);
	can8->cd();
	mu_SIC_EEm_chin->Draw("COLZ");
	//Util::PrintNoEPS(can8, can8->GetTitle(), outputdir);
	Util::PrintEPS(can8, can8->GetTitle(), outputdir);
	}

	TCanvas *c1 = new TCanvas("c1", "mu_SIC_EEm_chin", 0, 0, 900, 700);
	c1->cd();
	mu_std_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c1, mu_std_EEp->GetTitle(), outputdir);
	Util::PrintEPS(c1, mu_std_EEp->GetTitle(), outputdir);
//	c1->Clear();
	TCanvas *c2 = new TCanvas("c2", "mu_SIC_EEm_chin", 0, 0, 900, 700);
	c2->cd();
	mu_std_EEm->Draw("COLZ");
	//Util::PrintNoEPS(c2, mu_std_EEm->GetTitle(), outputdir);
	Util::PrintEPS(c2, mu_std_EEm->GetTitle(), outputdir);
//	c2->Clear();
	TCanvas *c3 = new TCanvas("c3", "mu_SIC_EEm_chin", 0, 0, 900, 700);
	c3->cd();
	mu_SIC_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c3, mu_SIC_EEp->GetTitle(), outputdir);
	Util::PrintEPS(c3, mu_SIC_EEp->GetTitle(), outputdir);
//	c3->Clear();
	TCanvas *c4 = new TCanvas("c4", "mu_SIC_EEm_chin", 0, 0, 900, 700);
	c4->cd();
	mu_SIC_EEm->Draw("COLZ");
	//Util::PrintNoEPS(c4, mu_SIC_EEm->GetTitle(), outputdir);
	Util::PrintEPS(c4, mu_SIC_EEm->GetTitle(), outputdir);

	TCanvas *c5 = new TCanvas("c5", "VPT_gainqe_EEp", 0, 0, 900, 700);
	c5->cd();
	VPT_gainqe_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c5, VPT_gainqe_EEp->GetTitle(), outputdir);
	Util::PrintEPS(c5, VPT_gainqe_EEp->GetTitle(), outputdir);
	TCanvas *c6 = new TCanvas("c6", "VPT_gainqe_EEm", 0, 0, 900, 700);
	c6->cd();
	VPT_gainqe_EEm->Draw("COLZ");
	//Util::PrintNoEPS(c6, VPT_gainqe_EEm->GetTitle(), outputdir);
	Util::PrintEPS(c6, VPT_gainqe_EEm->GetTitle(), outputdir);
	TCanvas *c7 = new TCanvas("c7", "VPT_gain_EEp", 0, 0, 900, 700);
	c7->cd();
	VPT_gain_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c7, VPT_gain_EEp->GetTitle(), outputdir);
	Util::PrintEPS(c7, VPT_gain_EEp->GetTitle(), outputdir);
	TCanvas *c8 = new TCanvas("c8", "VPT_gain_EEm", 0, 0, 900, 700);
	c8->cd();
	VPT_gain_EEm->Draw("COLZ");
	//Util::PrintNoEPS(c8, VPT_gain_EEm->GetTitle(), outputdir);
	Util::PrintEPS(c8, VPT_gain_EEm->GetTitle(), outputdir);
	TCanvas *c9 = new TCanvas("c9", "VPT_qe_EEp", 0, 0, 900, 700);
	c9->cd();
	VPT_qe_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c9, VPT_qe_EEp->GetTitle(), outputdir);
	Util::PrintEPS(c9, VPT_qe_EEp->GetTitle(), outputdir);
	TCanvas *c0 = new TCanvas("c0", "VPT_qe_EEm", 0, 0, 900, 700);
	c0->cd();
	VPT_qe_EEp->Draw("COLZ");
	//Util::PrintNoEPS(c0, VPT_qe_EEm->GetTitle(), outputdir);
	Util::PrintEPS(c0, VPT_qe_EEm->GetTitle(), outputdir);
//	c4->Clear();
	
	outf			->cd();
	mu_ECAL_EEp_russ	->Write();
	mu_ECAL_EEp_chin	->Write();
	mu_ECAL_EEm_russ	->Write();
	mu_ECAL_EEm_chin	->Write();
	mu_SIC_EEp_russ		->Write();
	mu_SIC_EEp_chin		->Write();
	mu_SIC_EEm_russ		->Write();
	mu_SIC_EEm_chin		->Write();
	EEp_producer		->Write();
	EEm_producer		->Write();
	mu_std_EEp		->Write();
	mu_std_EEm		->Write();
	mu_SIC_EEp		->Write();
	mu_SIC_EEm		->Write();
	VPT_gainqe_EEp		->Write();
	VPT_gainqe_EEm		->Write();
	VPT_gain_EEp		->Write();
	VPT_gain_EEm		->Write();
	VPT_qe_EEp		->Write();
	VPT_qe_EEm		->Write();
	cout << "saved histograms in: " << outf->GetName() << endl;
	outf->Close();

	PlotEtaProfile();

} 

void PlotEtaProfile(){

	TFile *outf = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEmustd.root");

	TH2D *mu_std_EEp       = (TH2D*)outf->Get("mu_std_EEp");
	TH2D *mu_std_EEm       = (TH2D*)outf->Get("mu_std_EEm");
	TH2D *mu_SIC_EEp       = (TH2D*)outf->Get("mu_SIC_EEp");
	TH2D *mu_SIC_EEm       = (TH2D*)outf->Get("mu_SIC_EEm");

	TH2D *VPT_gainqe_EEp   = (TH2D*)outf->Get("VPT_gainqe_EEp");
	TH2D *VPT_gainqe_EEm   = (TH2D*)outf->Get("VPT_gainqe_EEm");
	TH2D *VPT_gain_EEp     = (TH2D*)outf->Get("VPT_gain_EEp");
	TH2D *VPT_gain_EEm     = (TH2D*)outf->Get("VPT_gain_EEm");
	TH2D *VPT_qe_EEp       = (TH2D*)outf->Get("VPT_qe_EEp");
	TH2D *VPT_qe_EEm       = (TH2D*)outf->Get("VPT_qe_EEm");

	TH2D *EEp_producer     = (TH2D*)outf->Get("EEp_producer");
	TH2D *EEm_producer     = (TH2D*)outf->Get("EEm_producer");

	TFile *etaFile         = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root");
	TH2F *eta_EEp          = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm          = (TH2F*)etaFile->Get("eta_EEm");

	TFile *burninfile      = TFile::Open("/shome/haweber/ECAL/DataLaser/20121116_EEburnin.root");
	TH2D *burnin_EEm       = (TH2D*)burninfile->Get("burnin_EEm");
	TH2D *burnin_EEp       = (TH2D*)burninfile->Get("burnin_EEp");

	TH1D *mu_std[8][2];
	TH1D *eta_mu_std[8][2];
	TH1D *mu_SIC[8][2];
	TH1D *burn_in[8][2];
	TH1D *gain_qe[8][2];
	TH1D *gain_[8][2];
	TH1D *qe_[8][2];

	string string_feta[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string string_prod[2] = {"__russ", "__chin"};
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
		string hname;
		hname = "mu_std_eta"+string_feta[n]+string_prod[nn];
		mu_std[n][nn] = new TH1D(hname.c_str(), hname.c_str(),10,0.0,2.0);
		hname = "eta_mu_std_eta"+string_feta[n]+string_prod[nn];
		eta_mu_std[n][nn] = new TH1D(hname.c_str(), hname.c_str(),4,0.0,0.2);
		hname = "mu_SIC_eta"+string_feta[n]+string_prod[nn];
		mu_SIC[n][nn] = new TH1D(hname.c_str(), hname.c_str(),10,0.0,2.0);
		hname = "burnin_eta"+string_feta[n]+string_prod[nn];
		burn_in[n][nn] = new TH1D(hname.c_str(), hname.c_str(),10,0.9,1.0);
		hname = "gainqe_eta"+string_feta[n]+string_prod[nn];
		gain_qe[n][nn] = new TH1D(hname.c_str(), hname.c_str(),30,1.,4.);
		hname = "gain_eta"+string_feta[n]+string_prod[nn];
		gain_[n][nn] = new TH1D(hname.c_str(), hname.c_str(),28,6.,20.);
		hname = "qe_eta"+string_feta[n]+string_prod[nn];
		qe_[n][nn] = new TH1D(hname.c_str(), hname.c_str(),30,0.15,0.35);
	}}
	//EE+
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		int bin = EEp_producer->GetBin(i,j);
		int    prod   = EEp_producer  ->GetBinContent(bin);
		double eta    = fabs(eta_EEp  ->GetBinContent(bin));
		double mustd  = mu_std_EEp    ->GetBinContent(bin);
		double music  = mu_SIC_EEp    ->GetBinContent(bin);
		double burnin = burnin_EEp    ->GetBinContent(bin);
		double gain   = VPT_gain_EEp  ->GetBinContent(bin);
		double qe     = VPT_qe_EEp    ->GetBinContent(bin);
		double gainqe = VPT_gainqe_EEp->GetBinContent(bin);
		int p;
		if(prod==1) p = 0;
		else if(prod==2) p = 1;
		else continue;
		int e;
		double el;
		//cout << "gain " << gain << "  qe " << qe << "  gainqe " << gainqe << endl;
		if(eta<1.4) continue;
		else if(eta<1.6) {e = 0; el = 1.4; }
		else if(eta<1.8) {e = 1; el = 1.6; }
		else if(eta<2.0) {e = 2; el = 1.8; }
		else if(eta<2.2) {e = 3; el = 2.0; }
		else if(eta<2.4) {e = 4; el = 2.2; }
		else if(eta<2.6) {e = 5; el = 2.4; }
		else if(eta<2.8) {e = 6; el = 2.6; }
		else             {e = 7; el = 2.8; }
		if(mustd>0.001&&mustd<=2.) mu_std[e][p]->Fill(mustd);
		if(music>0.001&&music<=2.) mu_SIC[e][p]->Fill(music);
		burn_in[e][p]->Fill(burnin);
		gain_qe[e][p]->Fill(gainqe);
		gain_[e][p]->Fill(gain);
		qe_[e][p]->Fill(qe);
		if(mustd>0.001&&mustd<=2.) eta_mu_std[e][p]->Fill(eta-el);
	   }
	  }
	//EE-
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		int bin = EEm_producer->GetBin(i,j);
		int    prod   = EEm_producer  ->GetBinContent(bin);
		double eta    = fabs(eta_EEm  ->GetBinContent(bin));
		double mustd  = mu_std_EEm    ->GetBinContent(bin);
		double music  = mu_SIC_EEm    ->GetBinContent(bin);
		double burnin = burnin_EEm    ->GetBinContent(bin);
		double gain   = VPT_gain_EEm  ->GetBinContent(bin);
		double qe     = VPT_qe_EEm    ->GetBinContent(bin);
		double gainqe = VPT_gainqe_EEm->GetBinContent(bin);
		int p;
		if(prod==1) p = 0;
		else if(prod==2) p = 1;
		else continue;
		int e;
		double el;
		if(eta<1.4) continue;
		else if(eta<1.6) {e = 0; el = 1.4; }
		else if(eta<1.8) {e = 1; el = 1.6; }
		else if(eta<2.0) {e = 2; el = 1.8; }
		else if(eta<2.2) {e = 3; el = 2.0; }
		else if(eta<2.4) {e = 4; el = 2.2; }
		else if(eta<2.6) {e = 5; el = 2.4; }
		else if(eta<2.8) {e = 6; el = 2.6; }
		else             {e = 7; el = 2.8; }
		if(mustd>0.001&&mustd<=2.) mu_std[e][p]->Fill(mustd);
		if(music>0.001&&music<=2.) mu_SIC[e][p]->Fill(music);
		burn_in[e][p]->Fill(burnin);
		gain_qe[e][p]->Fill(gainqe);
		gain_[e][p]->Fill(gain);
		qe_[e][p]->Fill(qe);
		if(mustd>0.001&&mustd<=2.) eta_mu_std[e][p]->Fill(eta-el);
	   }
	  }

	TFile *outfile = new TFile("/shome/haweber/ECAL/DataLaser/testqualityparameter.root", "RECREATE");
	outfile->cd();
	for(int n = 0; n<8; ++n){//eta bin
	for(int nn = 0; nn<2; ++ nn){//producer
	mu_std[n][nn]->Write();
	eta_mu_std[n][nn]->Write();
	mu_SIC[n][nn]->Write();
	burn_in[n][nn]->Write();
	gain_qe[n][nn]->Write();
	gain_[n][nn]->Write();
	qe_[n][nn]->Write();
	}}
}
