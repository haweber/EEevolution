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
#include <TAxis.h>
#include <TSystem.h>
#include <TMap.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>

#include <map>
#include <string>
#include <iomanip>
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

float avg(vector<float> x){

  float sum = 0;
  for (unsigned int i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}

float sigma(vector<float> x, float mean){

  float sum = 0;
  for (unsigned int i = 0; i<x.size(); ++i){
    sum += (x[i]-mean)*(x[i]-mean);
  }
  if (x.size()!=0)  return sqrt(sum/(float)x.size());
  return 0;
}

float rms(vector<float> x){

  float sum = 0;
  for (unsigned int i = 0; i<x.size(); ++i){
    sum += x[i]*x[i];
  }
  if  (x.size()!=0)   return sqrt(sum/(float)x.size());
  return 0;
}

void EEstudies_continued_part2(){

	Bool_t plotting = true;

	TFile *oldFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/test.root");

	TH2F *mu_VPTPN_EEp           = (TH2F*)oldFile->Get("mu_VPTPN_EEp");
	TH2F *mu_VPTPN_EEm           = (TH2F*)oldFile->Get("mu_VPTPN_EEm");
	TH2F *VPTPN_ratio_EEp        = (TH2F*)oldFile->Get("VPTPN_ratio_EEp");
	TH2F *VPTPN_ratio_EEm        = (TH2F*)oldFile->Get("VPTPN_ratio_EEm");

	TH2F *mu_VPTPNnorm_EEp_chin  = (TH2F*)oldFile->Get("mu_VPTPNnorm_EEp_chin");
	TH2F *mu_VPTPNnorm_EEp_russ  = (TH2F*)oldFile->Get("mu_VPTPNnorm_EEp_russ");
	TH2F *mu_VPTPNnorm_EEm_chin  = (TH2F*)oldFile->Get("mu_VPTPNnorm_EEm_chin");
	TH2F *mu_VPTPNnorm_EEm_russ  = (TH2F*)oldFile->Get("mu_VPTPNnorm_EEm_russ");

	TFile *etaFile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_eta_in_ix_iy.root");
	TH2F *eta_EEp                = (TH2F*)etaFile->Get("eta_EEp");
	TH2F *eta_EEm                = (TH2F*)etaFile->Get("eta_EEm");

	TFile *stdmapsfile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20111219_EEplots.root");
	TH2D *mu_ECAL_EEp_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_russ");
	TH2D *mu_ECAL_EEp_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_chin");
	TH2D *mu_ECAL_EEm_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_russ");
	TH2D *mu_ECAL_EEm_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_chin");

	TH2D *EEp_producer     = (TH2D*)stdmapsfile->Get("EEp_producer");
	TH2D *EEm_producer     = (TH2D*)stdmapsfile->Get("EEm_producer");


	map<string, TH1F*> histos;
	map<string, TH2F*> histos2;

    string cap[2]  = {"EEm", "EEp"};
    string prod[3] = {"chin", "russ", "both"};
    string etap[9] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8", "all"};
    string etam[9] = {"3p0", "2p8", "2p6", "2p4", "2p2", "2p0", "1p8", "1p6", "all"};
    string etaname[9] = {"mu_std_vs_eta", "mu_ind_vs_eta", "mu_indSelected_vs_eta", "mu_ind_div_mu_std_vs_eta", "RMSmu_ind_div_RMSmu_std_vs_eta", "RMSmu_ind_vs_eta", "RMSmu_indSelected_vs_eta", "RMSmu_std_vs_eta", "RMS_mu_ind_div_mu_std_vs_eta"};
    string restname[] = {"mu_ind", "mu_indSelected", "mu_ind_vs_mu_std", "ratio_vs_mu_std", "mu_std", "mu_ind_div_mu_std",  "RMSmu_ind_vs_RMSmu_std", "mu_std_vs_eta_2d", "mu_ind_vs_eta_2d", "mu_indSelected_vs_eta_2d", "mu_ind_div_mu_std_vs_eta_2d"};
	const int Nmuind = 543;
	double muind[Nmuind + 1];
	for(int i = 0; i<28; ++i)    muind[i] = (-1.6+(double)i*0.05);//up to -0.205
	for(int i = 28; i<528; ++i)  muind[i] = (-0.2+(double)(i-28)*0.0005);//up to 0.0495
	for(int i = 528; i<544; ++i) muind[i] = (0.05+(double)(i-528)*0.05);//up to 0.8
	const int Nnorm = 870;
	double norm[Nnorm+1];
	for(int i = 0; i<50; ++i)    norm[i] = (-0.8+(double)i*0.01);//up to -0.3
	for(int i = 50; i<850; ++i)  norm[i] = (-0.3+(double)(i-50)*0.0005);//up to 0.1
	for(int i = 850; i<871; ++i) norm[i] = (0.1+(double)(i-850)*0.05);//up to 0.3
	const int Nrms = 280;
	double rms[Nrms + 1];
	for(int i = 0; i<40; ++i)    rms[i] = (     (double)i*0.02);//up to 0.8
	for(int i = 40; i<240; ++i)  rms[i] = (0.8+(double)(i-40)*0.0005);//up to 1.2 ->799 entries
	for(int i = 240; i<281; ++i) rms[i] = (1.2+(double)(i-240)*0.02);//up to 2
	const int Nrat = 103;
	double rat[Nrat + 1];
	for(int i = 0; i<25; ++i)    rat[i] = (0.6 +(double)i*0.01);//up to 0.85
	for(int i = 25; i<71; ++i)   rat[i] = (0.85+(double)(i-25)*0.005);//up to 1.08
	for(int i = 71; i<104; ++i)  rat[i] = (1.08+(double)(i-71)*0.01);//up to 1.4
	//const double muind[Nmuind +1 ] = hmuind;
	//const double norm[Nnorm +1 ]   = hnorm;
	//const double rms[Nrms +1 ]     = hrms;
	//const double rat[Nrat +1 ]     = hrat;
	const int Neta = 8;
	const double peta[Neta+1]      ={ 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
	const double neta[Neta+1]      ={-3.0,-2.8,-2.6,-2.4,-2.2,-2.0,-1.8,-1.6,-1.4};
	string hs, hs1;

	for(int c = 0; c<2; ++c){
	for(int p = 0; p<3; ++p){
	double etal, etau;
	if(c==0){ etal = -3.0; etau = -1.4; }
	else    { etal =  1.4; etau =  3.0; }
	hs1 = string("-")+cap[c]+string("-")+prod[p];
	hs = string("RMSmu_ind_vs_RMSmu_std")+hs1;
       	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 50, 0., 0.5, 80, -4., 0.);
	histos2[hs]->GetXaxis()->SetRangeUser(0.0, 0.5); histos2[hs]->GetYaxis()->SetRangeUser(-4.0, 0.0);
	histos2[hs]->SetXTitle("RMS/av mu_std"); histos2[hs]->SetYTitle("RMS/av mu_ind"); 
	hs = string("mu_ind_vs_mu_std") + hs1;
       	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 500, 0., 2.0, Nmuind, muind);
	histos2[hs]->GetXaxis()->SetRangeUser(0.5, 1.5); histos2[hs]->GetYaxis()->SetRangeUser(-0.7, 0.5);
	histos2[hs]->SetXTitle("mu_std"); histos2[hs]->SetYTitle("mu_ind"); 
	hs = string("ratio_vs_mu_std") + hs1;
       	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 500, 0., 2.0, Nrat, rat);
	histos2[hs]->GetXaxis()->SetRangeUser(0.5, 1.5); histos2[hs]->GetYaxis()->SetRangeUser(0.95, 1.2);
	histos2[hs]->SetXTitle("mu_std"); histos2[hs]->SetYTitle("ratio VPT/PN"); 
	hs = string("mu_ind_vs_eta_2d") + hs1;
	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 32, etal, etau, Nmuind, muind);
	histos2[hs]->GetYaxis()->SetRangeUser(-0.7, 0.5);
	histos2[hs]->SetXTitle("eta"); histos2[hs]->SetYTitle("mu_ind"); 
	hs = string("mu_indSelected_vs_eta_2d") + hs1;
	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 32, etal, etau, Nmuind, muind);
	histos2[hs]->GetYaxis()->SetRangeUser(-0.7, 0.5);
	histos2[hs]->SetXTitle("eta"); histos2[hs]->SetYTitle("mu_ind (where mu_std!=0)"); 
	hs = string("mu_ind_div_mu_std_vs_eta_2d")+hs1;
	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 32, etal, etau, Nnorm, norm);
	histos2[hs]->GetYaxis()->SetRangeUser(-0.7, 0.5);
	histos2[hs]->SetXTitle("eta"); histos2[hs]->SetYTitle("mu_ind/mu_std"); 
	hs = string("mu_std_vs_eta_2d")+hs1;
	if(histos2.count(hs) == 0 ) histos2[hs] = new TH2F(hs.c_str(), hs.c_str(), 32, etal, etau, 500, 0., 2.0);
	histos2[hs]->GetYaxis()->SetRangeUser(0.5, 1.5);
	histos2[hs]->SetXTitle("eta"); histos2[hs]->SetYTitle("mu_std"); 
	hs = string("mu_std_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(0.5, 1.5);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("av mu_std"); 
	hs = string("mu_ind_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-0.25, 0.1);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("av mu_ind"); 
	hs = string("mu_indSelected_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-0.25, 0.1);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("av mu_ind (where mu_std!=0)"); 
	hs = string("mu_ind_div_mu_std_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-0.25, 0.1);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("av mu_ind/mu_std"); 
	hs = string("RMSmu_ind_div_RMSmu_std_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-6.0, 1.0);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("(RMS/av mu_ind)/(RMS/av mu_std)"); 
	hs = string("RMSmu_ind_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-6.0, 1.0);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("RMS/av mu_ind"); 
	hs = string("RMSmu_indSelected_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-6.0, 1.0);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("RMS/av mu_ind (where mu_std!=0)"); 
	hs = string("RMSmu_std_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(0.0, 0.5);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("RMS/av mu_std"); 
	hs = string("RMS_mu_ind_div_mu_std_vs_eta")+hs1;
        if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 8, etal, etau);
	histos[hs]->GetYaxis()->SetRangeUser(-6.0, 1.0);
	histos[hs]->SetXTitle("eta"); histos[hs]->SetYTitle("RMS/av mu_ind/mu_std"); 
	for(int e = 0; e<9; ++e){
	if(c==0) hs1 = string("-")+cap[c]+string("-")+prod[p]+string("-")+etam[e];
	else     hs1 = string("-")+cap[c]+string("-")+prod[p]+string("-")+etap[e];
	hs = string("mu_std")+hs1;
       	if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), 500, 0., 2.0);
	if(etap[e]!="1p4") {histos[hs]->GetXaxis()->SetRangeUser(0.0, 2.0); histos[hs]->GetYaxis()->SetRangeUser(0., 15.); }
	else                histos[hs]->GetXaxis()->SetRangeUser(0.0, 2.0); histos[hs]->GetYaxis()->SetRangeUser(0., 25.);
	histos[hs]->SetXTitle("mu_std"); histos[hs]->SetYTitle("entries"); 
	hs = string("mu_ind")+hs1;
       	if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), Nmuind, muind);
	if((c==0 && e>=5) || (c==1 && e<=3)){ histos[hs]->GetXaxis()->SetRangeUser(-0.7, 0.5); histos[hs]->GetYaxis()->SetRangeUser(0, 85); }//eta<=2.0
	else                                  histos[hs]->GetXaxis()->SetRangeUser(-0.7, 0.5); histos[hs]->GetYaxis()->SetRangeUser(0, 25);//eta> 2.0
	histos[hs]->SetXTitle("mu_ind"); histos[hs]->SetYTitle("entries"); 
	hs = string("mu_indSelected")+hs1;
       	if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), Nmuind, muind);
	histos[hs]->GetXaxis()->SetRangeUser(-0.7, 0.5); histos[hs]->GetYaxis()->SetRangeUser(0, 20);
	histos[hs]->SetXTitle("mu_ind (where mu_std!=0)"); histos[hs]->SetYTitle("entries"); 
	hs = string("mu_ind_div_mu_std")+hs1;
       	if(histos.count(hs) == 0 ) histos[hs] = new TH1F(hs.c_str(), hs.c_str(), Nnorm, norm);
	histos[hs]->GetXaxis()->SetRangeUser(-0.7, 0.5); histos[hs]->GetYaxis()->SetRangeUser(0, 20);
	histos[hs]->SetXTitle("mu_ind/mu_std"); histos[hs]->SetYTitle("entries"); 
	}//e
	}//p
	}//c

	for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
		h->second->Sumw2();}
	for(map<string,TH2F*>::iterator h=histos2.begin(); h!=histos2.end();++h){
		h->second->Sumw2();}

	for(int x = 0; x<101; ++x){
		for(int y = 0; y<101; ++y){
			//EE+
			float value_std, value_rat, value_ind, value_nor;
			int bin; float eta; string seta; int producer;
			bin = eta_EEp->FindBin(x,y);
			eta = eta_EEp->GetBinContent(bin);
			//if(eta==0) continue;//no crystal
			producer = EEp_producer->GetBinContent(bin);
			seta = "all";
			if(eta<1.4 && eta!=0) { cout << "what eta " << eta << endl; seta = "1p4"; }
			else if(eta>=1.4 && eta<1.6) seta = "1p4";
			else if(eta>=1.6 && eta<1.8) seta = "1p6";
			else if(eta>=1.8 && eta<2.0) seta = "1p8";
			else if(eta>=2.0 && eta<2.2) seta = "2p0";
			else if(eta>=2.2 && eta<2.4) seta = "2p2";
			else if(eta>=2.4 && eta<2.6) seta = "2p4";
			else if(eta>=2.6 && eta<2.8) seta = "2p6";
			else if(eta>=2.8 && eta<3.0) seta = "2p8";
			else if(eta!=0) { cout << "what eta " << eta << endl; seta = "2p8"; }
			//russ
			value_std = mu_ECAL_EEp_russ->GetBinContent(bin);
			value_rat = VPTPN_ratio_EEp->GetBinContent(bin);
			value_ind = mu_VPTPN_EEp->GetBinContent(bin);
			value_nor = mu_VPTPNnorm_EEp_russ->GetBinContent(bin);
			if(value_std!=0){
				histos[(string)"mu_std-EEp-russ-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEp-russ-all"]->Fill(value_std);
				histos[(string)"mu_std-EEp-both-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEp-both-all"]->Fill(value_std);
				histos[(string)"mu_indSelected-EEp-russ-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-russ-all"]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-both-all"]->Fill(value_ind);
				histos[(string)"mu_ind_div_mu_std-EEp-russ-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-russ-all"]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-both-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-both-all"]->Fill(value_nor);
				histos2[(string)"mu_ind_vs_mu_std-EEp-russ"]->Fill(value_std, value_ind);
				histos2[(string)"mu_ind_vs_mu_std-EEp-both"]->Fill(value_std, value_ind);
				histos2[(string)"ratio_vs_mu_std-EEp-russ"]->Fill(value_std, value_rat);
				histos2[(string)"ratio_vs_mu_std-EEp-both"]->Fill(value_std, value_rat);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEp-russ"]->Fill(eta, value_ind);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEp-both"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEp-russ"]->Fill(eta, value_nor);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEp-both"]->Fill(eta, value_nor);
				histos2[(string)"mu_std_vs_eta_2d-EEp-russ"]->Fill(eta, value_std);
				histos2[(string)"mu_std_vs_eta_2d-EEp-both"]->Fill(eta, value_std);
			}
			if(value_ind!=0){
				if(producer==1) histos[(string)"mu_ind-EEp-russ-"+seta]->Fill(value_ind);
				if(producer==1) histos[(string)"mu_ind-EEp-russ-all"]->Fill(value_ind);
				histos[(string)"mu_ind-EEp-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_ind-EEp-both-all"]->Fill(value_ind);
				if(producer==1) histos2[(string)"mu_ind_vs_eta_2d-EEp-russ"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_vs_eta_2d-EEp-both"]->Fill(eta, value_ind);
			}
			//chin
			value_std = mu_ECAL_EEp_chin->GetBinContent(bin);
			value_rat = VPTPN_ratio_EEp->GetBinContent(bin);
			value_ind = mu_VPTPN_EEp->GetBinContent(bin);
			value_nor = mu_VPTPNnorm_EEp_chin->GetBinContent(bin);
			if(value_std!=0){
				histos[(string)"mu_std-EEp-chin-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEp-chin-all"]->Fill(value_std);
				histos[(string)"mu_std-EEp-both-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEp-both-all"]->Fill(value_std);
				histos[(string)"mu_indSelected-EEp-chin-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-chin-all"]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEp-both-all"]->Fill(value_ind);
				histos[(string)"mu_ind_div_mu_std-EEp-chin-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-chin-all"]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-both-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEp-both-all"]->Fill(value_nor);
				histos2[(string)"mu_ind_vs_mu_std-EEp-chin"]->Fill(value_std, value_ind);
				histos2[(string)"mu_ind_vs_mu_std-EEp-both"]->Fill(value_std, value_ind);
				histos2[(string)"ratio_vs_mu_std-EEp-chin"]->Fill(value_std, value_rat);
				histos2[(string)"ratio_vs_mu_std-EEp-both"]->Fill(value_std, value_rat);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEp-chin"]->Fill(eta, value_ind);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEp-both"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEp-chin"]->Fill(eta, value_nor);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEp-both"]->Fill(eta, value_nor);
				histos2[(string)"mu_std_vs_eta_2d-EEp-chin"]->Fill(eta, value_std);
				histos2[(string)"mu_std_vs_eta_2d-EEp-both"]->Fill(eta, value_std);
			}
			if(value_ind!=0){
				if(producer==2) histos[(string)"mu_ind-EEp-chin-"+seta]->Fill(value_ind);
				if(producer==2) histos[(string)"mu_ind-EEp-chin-all"]->Fill(value_ind);
				if(producer==2) histos2[(string)"mu_ind_vs_eta_2d-EEp-chin"]->Fill(eta, value_ind);
			}
			//EE-
			bin = eta_EEm->FindBin(x,y);
			eta = eta_EEm->GetBinContent(bin);
			//if(eta==0) continue;//no crystal
			producer = EEm_producer->GetBinContent(bin);
			seta = "all";
			if(eta<-3.0) { cout << "what eta " << eta << endl; seta = "3p0"; }
			else if(eta>=-3.0 && eta<-2.8) seta = "3p0";
			else if(eta>=-2.8 && eta<-2.6) seta = "2p8";
			else if(eta>=-2.6 && eta<-2.4) seta = "2p4";
			else if(eta>=-2.4 && eta<-2.2) seta = "2p6";
			else if(eta>=-2.2 && eta<-2.0) seta = "2p2";
			else if(eta>=-2.0 && eta<-1.8) seta = "2p0";
			else if(eta>=-1.8 && eta<-1.6) seta = "1p8";
			else if(eta>=-1.6 && eta<-1.4) seta = "1p6";
			else if(eta!=0) { cout << "what eta " << eta << endl; seta = "1p6"; }
			//russ
			value_std = mu_ECAL_EEm_russ->GetBinContent(bin);
			value_rat = VPTPN_ratio_EEm->GetBinContent(bin);
			value_ind = mu_VPTPN_EEm->GetBinContent(bin);
			value_nor = mu_VPTPNnorm_EEm_russ->GetBinContent(bin);
			if(value_std!=0){
				histos[(string)"mu_std-EEm-russ-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEm-russ-all"]->Fill(value_std);
				histos[(string)"mu_std-EEm-both-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEm-both-all"]->Fill(value_std);
				histos[(string)"mu_indSelected-EEm-russ-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-russ-all"]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-both-all"]->Fill(value_ind);
				histos[(string)"mu_ind_div_mu_std-EEm-russ-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-russ-all"]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-both-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-both-all"]->Fill(value_nor);
				histos2[(string)"mu_ind_vs_mu_std-EEm-russ"]->Fill(value_std, value_ind);
				histos2[(string)"mu_ind_vs_mu_std-EEm-both"]->Fill(value_std, value_ind);
				histos2[(string)"ratio_vs_mu_std-EEm-russ"]->Fill(value_std, value_rat);
				histos2[(string)"ratio_vs_mu_std-EEm-both"]->Fill(value_std, value_rat);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEm-russ"]->Fill(eta, value_ind);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEm-both"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEm-russ"]->Fill(eta, value_nor);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEm-both"]->Fill(eta, value_nor);
				histos2[(string)"mu_std_vs_eta_2d-EEm-russ"]->Fill(eta, value_std);
				histos2[(string)"mu_std_vs_eta_2d-EEm-both"]->Fill(eta, value_std);
			}
			if(value_ind!=0){
				if(producer==1) histos[(string)"mu_ind-EEm-russ-"+seta]->Fill(value_ind);
				if(producer==1) histos[(string)"mu_ind-EEm-russ-all"]->Fill(value_ind);
				histos[(string)"mu_ind-EEm-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_ind-EEm-both-all"]->Fill(value_ind);
				if(producer==1) histos2[(string)"mu_ind_vs_eta_2d-EEm-russ"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_vs_eta_2d-EEm-both"]->Fill(eta, value_ind);
			}
			//chin
			value_std = mu_ECAL_EEm_chin->GetBinContent(bin);
			value_rat = VPTPN_ratio_EEm->GetBinContent(bin);
			value_ind = mu_VPTPN_EEm->GetBinContent(bin);
			value_nor = mu_VPTPNnorm_EEm_chin->GetBinContent(bin);
			if(value_std!=0){
				histos[(string)"mu_std-EEm-chin-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEm-chin-all"]->Fill(value_std);
				histos[(string)"mu_std-EEm-both-"+seta]->Fill(value_std);
				histos[(string)"mu_std-EEm-both-all"]->Fill(value_std);
				histos[(string)"mu_indSelected-EEm-chin-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-chin-all"]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-both-"+seta]->Fill(value_ind);
				histos[(string)"mu_indSelected-EEm-both-all"]->Fill(value_ind);
				histos[(string)"mu_ind_div_mu_std-EEm-chin-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-chin-all"]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-both-"+seta]->Fill(value_nor);
				histos[(string)"mu_ind_div_mu_std-EEm-both-all"]->Fill(value_nor);
				histos2[(string)"mu_ind_vs_mu_std-EEm-chin"]->Fill(value_std, value_ind);
				histos2[(string)"mu_ind_vs_mu_std-EEm-both"]->Fill(value_std, value_ind);
				histos2[(string)"ratio_vs_mu_std-EEm-chin"]->Fill(value_std, value_rat);
				histos2[(string)"ratio_vs_mu_std-EEm-both"]->Fill(value_std, value_rat);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEm-chin"]->Fill(eta, value_ind);
				histos2[(string)"mu_indSelected_vs_eta_2d-EEm-both"]->Fill(eta, value_ind);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEm-chin"]->Fill(eta, value_nor);
				histos2[(string)"mu_ind_div_mu_std_vs_eta_2d-EEm-both"]->Fill(eta, value_nor);
				histos2[(string)"mu_std_vs_eta_2d-EEm-chin"]->Fill(eta, value_std);
				histos2[(string)"mu_std_vs_eta_2d-EEm-both"]->Fill(eta, value_std);
			}
			if(value_ind!=0){
				if(producer==2) histos[(string)"mu_ind-EEm-chin-"+seta]->Fill(value_ind);
				if(producer==2) histos[(string)"mu_ind-EEm-chin-all"]->Fill(value_ind);
				if(producer==2) histos2[(string)"mu_ind_vs_eta_2d-EEm-chin"]->Fill(eta, value_ind);
			}
		}//y
	}//x
	//EE+
	for(int c = 0; c<2; ++c){//EE+,EE-
	for(int p = 0; p<3; ++p){//chin, russ, both
	for(int e = 0; e<8; ++e){
	int etabin;
	if(c==0)      etabin = histos["mu_std_vs_eta-EEm-both"]->FindBin(-2.9+ e *0.2);
	else if(c==1) etabin = histos["mu_std_vs_eta-EEp-both"]->FindBin(1.5+ e *0.2);
	float av_std, sig_std, av_ind, sig_ind, av_indS, sig_indS, av_norm, sig_norm;

    	string cc[2]  = {"EEm", "EEp"};
    	string pp[3] = {"chin", "russ", "both"};
//	string ee[8] = (c==0) ? {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"} : {"3p0", "2p8", "2p6", "2p4", "2p2", "2p0", "1p8", "1p6"};
    	string eep[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
    	string eem[8] = {"3p0", "2p8", "2p6", "2p4", "2p2", "2p0", "1p8", "1p6"};
//    	string ee[8] = {"1p4", "1p6", "1p8", "2p0", "2p2", "2p4", "2p6", "2p8"};
	string s = "-";
	//both
	if(c==0){
		av_ind  = histos[(string)"mu_ind"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetMean();
		sig_ind = histos[(string)"mu_ind"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetRMS();
		av_std  = histos[(string)"mu_std"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetMean();
		sig_std = histos[(string)"mu_std"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetRMS();
		av_norm = histos[(string)"mu_ind_div_mu_std"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetMean();
		sig_norm= histos[(string)"mu_ind_div_mu_std"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetRMS();
		av_indS = histos[(string)"mu_indSelected"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetMean();
		sig_indS= histos[(string)"mu_indSelected"+s+cc[c]+s+pp[p]+s+eem[e] ]->GetRMS();
	}
	else if(c==1){
		av_ind  = histos[(string)"mu_ind"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetMean();
		sig_ind = histos[(string)"mu_ind"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetRMS();
		av_std  = histos[(string)"mu_std"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetMean();
		sig_std = histos[(string)"mu_std"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetRMS();
		av_norm = histos[(string)"mu_ind_div_mu_std"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetMean();
		sig_norm= histos[(string)"mu_ind_div_mu_std"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetRMS();
		av_indS = histos[(string)"mu_indSelected"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetMean();
		sig_indS= histos[(string)"mu_indSelected"+s+cc[c]+s+pp[p]+s+eep[e] ]->GetRMS();
	}
	histos[(string)"mu_indSelected_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, av_indS);
	histos[(string)"mu_indSelected_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinError(etabin, sig_indS);
	histos[(string)"mu_ind_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, av_ind);
	histos[(string)"mu_ind_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinError(etabin, sig_ind);
	histos[(string)"mu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, av_std);
	histos[(string)"mu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinError(etabin, sig_std);
	histos[(string)"mu_ind_div_mu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, av_norm);
	histos[(string)"mu_ind_div_mu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinError(etabin, sig_norm);
	if(av_std!=0) histos[(string)"RMSmu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, sig_std/av_std);
	if(av_ind!=0) histos[(string)"RMSmu_ind_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, sig_ind/av_ind);
	if(av_indS!=0) histos[(string)"RMSmu_indSelected_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, sig_indS/av_indS);
	if(av_norm!=0) histos[(string)"RMS_mu_ind_div_mu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, sig_norm/av_norm);
	if(sig_std!=0 && av_indS!=0 && av_std!=0 ) histos[(string)"RMSmu_ind_div_RMSmu_std_vs_eta"+s+cc[c]+s+pp[p] ]->SetBinContent(etabin, (sig_indS/av_indS)/(sig_std/av_std) );
	if(av_std!=0 && av_ind!=0) histos2[(string)"RMSmu_ind_vs_RMSmu_std"+s+cc[c]+s+pp[p] ]->Fill((sig_std/av_std), (sig_ind/av_ind) );
	else cout << "why :(" << "     " << "av_std " << av_std << " av_ind " << av_ind << endl;
	}//eta
	}//producer
	}//endcap

   	 TFile *fsavefile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120224_EEstudies_continued_part2.root","RECREATE");
	fsavefile->cd();
	for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
		if(h->second->GetEntries()>0) h->second->Write();
	}
	for(map<string,TH2F*>::iterator h=histos2.begin(); h!=histos2.end();++h){
		if(h->second->GetEntries()>0) h->second->Write();
	}
	fsavefile->Close();
	cout << "Saved histograms in " << "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/" << "20120224_EEstudies_continued_part2.root" << endl;



	//plotting
	if(plotting){
		gROOT->SetStyle("Plain");
		gStyle->SetOptStat(0);
		//gStyle->SetPalette(1);

/*	//old plotting: plot by plot
		TString outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/part2";
		Util::MakeOutputDir(outputdir);
		for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
			TString outputfile = (TString)h->first;
			TCanvas *col = new TCanvas("col", outputfile, 0, 0, 1600, 975);
			h->second->Draw("");
			gPad->RedrawAxis();
			col ->Update();
		//	Util::PrintNoEPS(col, outputfile, outputdir, false);
		//	Util::PrintEPS(col, outputfile, outputdir);
			delete col;
			getchar();
		}
		for(map<string,TH2F*>::iterator h=histos2.begin(); h!=histos2.end();++h){
			TString outputfile = (TString)h->first;
			TCanvas *col = new TCanvas("col", outputfile, 0, 0, 900, 700);
			h->second->Draw("");
			gPad->RedrawAxis();
			col ->Update();
		//	Util::PrintNoEPS(col, outputfile, outputdir, false);
		//	Util::PrintEPS(col, outputfile, outputdir);
			delete col;
		}
*/
	//new plotting, 5x4 plots per figure

		TString outputdir = "/shome/haweber/ECAL/Laser/Analysis/Plots/20120224/EEstudies/continued/part2/packed";
		Util::MakeOutputDir(outputdir);
			TCanvas *col = new TCanvas("col", "col", 60, 0, 1575, 980);
			col->cd();
			TPad *pad[5][4];
			for(int xx = 0; xx<5; ++xx){
				for(int yy = 0; yy<4; ++yy){
					TString padname = TString::Format("pad_x%i_y%i",xx,yy);
					pad[xx][yy] = new TPad(padname, padname, (float)xx*0.2, 0.75-(float)yy*0.25, (float)xx*0.2+0.2, 1.-(float)yy*0.25, 0,0);
					//cout << "xlow " << (float)xx*0.2 << " xup " << (float)xx*0.2+0.2 << " ylow " << 0.75-(float)yy*0.25 << " yup " << 1.-(float)yy*0.25 << endl;
					pad[xx][yy]->SetBottomMargin(0.1);
					pad[xx][yy]->SetTopMargin(0.1);
					pad[xx][yy]->Draw();
				}
			}
			int x = 0; int y = 0; int count = 1;
			for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
				pad[x][y]->Clear();
				pad[x][y]->cd();
				h->second->Draw("");
				if(x!=4) {++x;}
				else if(x==4 && y!=3) {++y; x=0; }
				else { x=0; y =0;
					TString outputfile = TString::Format("EEstudies_plot_no%i",count);
					col->Update();
					Util::PrintNoEPS(col, outputfile, outputdir, false);
					Util::PrintEPS(col, outputfile, outputdir);
					++count;
				}
			}	
			for(map<string,TH2F*>::iterator h=histos2.begin(); h!=histos2.end();++h){
				pad[x][y]->Clear();//delete all previous entries
				pad[x][y]->cd();
				h->second->Draw("");
				if(x!=4) {++x;}
				else if(x==4 && y!=3) {++y; x=0; }
				else { x=0; y =0;//reset values
					TString outputfile = TString::Format("EEstudies_plot_no%i",count);
					col->Update();
					Util::PrintNoEPS(col, outputfile, outputdir, false);
					Util::PrintEPS(col, outputfile, outputdir);
					++count;
				}
			}
			for(int xx = 0; xx<5; ++xx){
				for(int yy = 0; yy<4; ++yy){
					delete pad[xx][yy];
				}
			}
			delete col;	

	}


}


/*
	TH1I *mu_VPTPN_EEp_entries = new TH1I("mu_VPTPN_EEp_entries", "mu_VPTPN_EEp_entries", 8, 1.4, 3.0);
	TH1I *mu_VPTPN_EEm_entries = new TH1I("mu_VPTPN_EEm_entries", "mu_VPTPN_EEm_entries", 8, -3.0, -1.4);
	TH1F *mu_VPTPN_EEp_average = new TH1F("mu_VPTPN_EEp_average", "mu_VPTPN_EEp_average", 8, 1.4, 3.0);
	TH1F *mu_VPTPN_EEm_average = new TH1F("mu_VPTPN_EEm_average", "mu_VPTPN_EEm_average", 8, -3.0, -1.4);

	TH1I *mu_VPTPNnorm_EEp_chin_entries = new TH1I("mu_VPTPNnorm_EEp_chin_entries", "mu_VPTPNnorm_EEp_chin_entries", 8, 1.4, 3.0);
	TH1I *mu_VPTPNnorm_EEm_chin_entries = new TH1I("mu_VPTPNnorm_EEm_chin_entries", "mu_VPTPNnorm_EEm_chin_entries", 8, -3.0, -1.4);
	TH1F *mu_VPTPNnorm_EEp_chin_average = new TH1F("mu_VPTPNnorm_EEp_chin_average", "mu_VPTPNnorm_EEp_chin_average", 8, 1.4, 3.0);
	TH1F *mu_VPTPNnorm_EEm_chin_average = new TH1F("mu_VPTPNnorm_EEm_chin_average", "mu_VPTPNnorm_EEm_chin_average", 8, -3.0, -1.4);
	TH1I *mu_VPTPNnorm_EEp_russ_entries = new TH1I("mu_VPTPNnorm_EEp_russ_entries", "mu_VPTPNnorm_EEp_russ_entries", 8, 1.4, 3.0);
	TH1I *mu_VPTPNnorm_EEm_russ_entries = new TH1I("mu_VPTPNnorm_EEm_russ_entries", "mu_VPTPNnorm_EEm_russ_entries", 8, -3.0, -1.4);
	TH1F *mu_VPTPNnorm_EEp_russ_average = new TH1F("mu_VPTPNnorm_EEp_russ_average", "mu_VPTPNnorm_EEp_russ_average", 8, 1.4, 3.0);
	TH1F *mu_VPTPNnorm_EEm_russ_average = new TH1F("mu_VPTPNnorm_EEm_russ_average", "mu_VPTPNnorm_EEm_russ_average", 8, -3.0, -1.4);

	mu_VPTPN_EEp_entries->Sumw2();
	mu_VPTPN_EEm_entries->Sumw2();
	mu_VPTPN_EEp_average->Sumw2();
	mu_VPTPN_EEm_average->Sumw2();
	mu_VPTPNnorm_EEp_chin_entries->Sumw2();
	mu_VPTPNnorm_EEm_chin_entries->Sumw2();
	mu_VPTPNnorm_EEp_chin_average->Sumw2();
	mu_VPTPNnorm_EEm_chin_average->Sumw2();
	mu_VPTPNnorm_EEp_russ_entries->Sumw2();
	mu_VPTPNnorm_EEm_russ_entries->Sumw2();
	mu_VPTPNnorm_EEp_russ_average->Sumw2();
	mu_VPTPNnorm_EEm_russ_average->Sumw2();


	vector<float> mu_VPTPN_EEp_eta13, mu_VPTPN_EEp_eta14, mu_VPTPN_EEp_eta15, mu_VPTPN_EEp_eta16, mu_VPTPN_EEp_eta17, mu_VPTPN_EEp_eta18, mu_VPTPN_EEp_eta19, mu_VPTPN_EEp_eta20, mu_VPTPN_EEp_eta21, mu_VPTPN_EEp_eta22, mu_VPTPN_EEp_eta23, mu_VPTPN_EEp_eta24, mu_VPTPN_EEp_eta25, mu_VPTPN_EEp_eta26, mu_VPTPN_EEp_eta27, mu_VPTPN_EEp_eta28, mu_VPTPN_EEp_eta29;
	vector<float> mu_VPTPN_EEm_eta15, mu_VPTPN_EEm_eta16, mu_VPTPN_EEm_eta17, mu_VPTPN_EEm_eta18, mu_VPTPN_EEm_eta19, mu_VPTPN_EEm_eta20, mu_VPTPN_EEm_eta21, mu_VPTPN_EEm_eta22, mu_VPTPN_EEm_eta23, mu_VPTPN_EEm_eta24, mu_VPTPN_EEm_eta25, mu_VPTPN_EEm_eta26, mu_VPTPN_EEm_eta27, mu_VPTPN_EEm_eta28, mu_VPTPN_EEm_eta29, mu_VPTPN_EEm_eta30, mu_VPTPN_EEm_eta31;
	vector<float> mu_VPTPN_EEp_chin_eta13, mu_VPTPN_EEp_chin_eta14, mu_VPTPN_EEp_chin_eta15, mu_VPTPN_EEp_chin_eta16, mu_VPTPN_EEp_chin_eta17, mu_VPTPN_EEp_chin_eta18, mu_VPTPN_EEp_chin_eta19, mu_VPTPN_EEp_chin_eta20, mu_VPTPN_EEp_chin_eta21, mu_VPTPN_EEp_chin_eta22, mu_VPTPN_EEp_chin_eta23, mu_VPTPN_EEp_chin_eta24, mu_VPTPN_EEp_chin_eta25, mu_VPTPN_EEp_chin_eta26, mu_VPTPN_EEp_chin_eta27, mu_VPTPN_EEp_chin_eta28, mu_VPTPN_EEp_chin_eta29;
	vector<float> mu_VPTPN_EEm_chin_eta15, mu_VPTPN_EEm_chin_eta16, mu_VPTPN_EEm_chin_eta17, mu_VPTPN_EEm_chin_eta18, mu_VPTPN_EEm_chin_eta19, mu_VPTPN_EEm_chin_eta20, mu_VPTPN_EEm_chin_eta21, mu_VPTPN_EEm_chin_eta22, mu_VPTPN_EEm_chin_eta23, mu_VPTPN_EEm_chin_eta24, mu_VPTPN_EEm_chin_eta25, mu_VPTPN_EEm_chin_eta26, mu_VPTPN_EEm_chin_eta27, mu_VPTPN_EEm_chin_eta28, mu_VPTPN_EEm_chin_eta29, mu_VPTPN_EEm_chin_eta30, mu_VPTPN_EEm_chin_eta31;
	vector<float> mu_VPTPN_EEp_russ_eta13, mu_VPTPN_EEp_russ_eta14, mu_VPTPN_EEp_russ_eta15, mu_VPTPN_EEp_russ_eta16, mu_VPTPN_EEp_russ_eta17, mu_VPTPN_EEp_russ_eta18, mu_VPTPN_EEp_russ_eta19, mu_VPTPN_EEp_russ_eta20, mu_VPTPN_EEp_russ_eta21, mu_VPTPN_EEp_russ_eta22, mu_VPTPN_EEp_russ_eta23, mu_VPTPN_EEp_russ_eta24, mu_VPTPN_EEp_russ_eta25, mu_VPTPN_EEp_russ_eta26, mu_VPTPN_EEp_russ_eta27, mu_VPTPN_EEp_russ_eta28, mu_VPTPN_EEp_russ_eta29;
	vector<float> mu_VPTPN_EEm_russ_eta15, mu_VPTPN_EEm_russ_eta16, mu_VPTPN_EEm_russ_eta17, mu_VPTPN_EEm_russ_eta18, mu_VPTPN_EEm_russ_eta19, mu_VPTPN_EEm_russ_eta20, mu_VPTPN_EEm_russ_eta21, mu_VPTPN_EEm_russ_eta22, mu_VPTPN_EEm_russ_eta23, mu_VPTPN_EEm_russ_eta24, mu_VPTPN_EEm_russ_eta25, mu_VPTPN_EEm_russ_eta26, mu_VPTPN_EEm_russ_eta27, mu_VPTPN_EEm_russ_eta28, mu_VPTPN_EEm_russ_eta29, mu_VPTPN_EEm_russ_eta30, mu_VPTPN_EEm_russ_eta31;
	for(int x = 0; x<101; ++x){
		for(int y = 0; y<101; ++y){
			int bin = eta_EEp->FindBin(x,y);
			//EEp
			float eta = eta_EEp->GetBinContent(bin);
			if(eta<1.4){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta13.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta13.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta13.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<1.5){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta13.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta14.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta13.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta14.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta13.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta14.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<1.6){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta14.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta15.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta14.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta15.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta14.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta15.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<1.7){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta15.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta16.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta15.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta16.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta15.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta16.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<1.8){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta16.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta17.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta16.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta17.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta16.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta17.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<1.9){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta17.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta18.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta17.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta18.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta17.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta18.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.0){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta18.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta19.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta18.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta19.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta18.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta19.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.1){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta19.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta20.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta19.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta20.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta19.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta20.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.2){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta20.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta21.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta20.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta21.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta20.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta21.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.3){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta21.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta22.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta21.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta22.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta21.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta22.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.4){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta22.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta23.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta22.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta23.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta22.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta23.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.5){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta23.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta24.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta23.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta24.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta23.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta24.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.6){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta24.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta25.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta24.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta25.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta24.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta25.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.7){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta25.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta26.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta25.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta26.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta25.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta26.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.8){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta26.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta27.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta26.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta27.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta26.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta27.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<2.9){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta27.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta28.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta27.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta28.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta27.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta28.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<3.0){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta28.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta29.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta28.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta29.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta28.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta29.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta<3.1){
				if(mu_VPTPN_EEp->GetBinContent(bin)!=0) mu_VPTPN_EEp_eta29.push_back(mu_VPTPN_EEp->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_chin->GetBinContent(bin)!=0) mu_VPTPN_EEp_chin_eta29.push_back(mu_VPTPNnorm_EEp_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEp_russ->GetBinContent(bin)!=0) mu_VPTPN_EEp_russ_eta29.push_back(mu_VPTPNnorm_EEp_russ->GetBinContent(bin));
			}
			else if(eta!=0) cout << "eta " << eta << endl;
			//EE-
			bin = eta_EEm->FindBin(x,y);
			eta = eta_EEm->GetBinContent(bin);
			if(eta<(-3.0)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta31.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta31.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta31.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.9)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta31.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta30.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta31.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta30.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta31.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta30.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.8)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta30.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta29.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta30.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta29.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta30.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta29.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.7)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta29.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta28.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta29.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta28.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta29.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta28.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.6)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta28.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta27.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta28.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta27.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta28.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta27.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.5)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta27.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta26.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta27.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta26.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta27.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta26.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.4)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta26.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta25.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta26.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta25.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta26.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta25.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.3)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta25.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta24.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta25.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta24.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta25.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta24.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.2)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta24.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta23.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta24.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta23.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta24.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta23.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.1)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta23.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta22.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta23.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta22.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta23.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta22.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-2.0)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta22.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta21.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta22.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta21.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta22.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta21.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.9)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta21.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta20.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta21.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta20.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta21.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta20.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.8)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta20.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta19.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta20.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta19.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta20.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta19.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.7)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta19.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta18.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta19.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta18.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta19.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta18.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.6)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta18.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta17.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta18.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta17.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta18.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta17.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.5)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta17.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta16.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta17.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta16.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta17.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta16.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.4)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta16.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta15.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta16.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta15.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta16.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta15.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta<(-1.3)){
				if(mu_VPTPN_EEm->GetBinContent(bin)!=0) mu_VPTPN_EEm_eta15.push_back(mu_VPTPN_EEm->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_chin->GetBinContent(bin)!=0) mu_VPTPN_EEm_chin_eta15.push_back(mu_VPTPNnorm_EEm_chin->GetBinContent(bin));
				if(mu_VPTPNnorm_EEm_russ->GetBinContent(bin)!=0) mu_VPTPN_EEm_russ_eta15.push_back(mu_VPTPNnorm_EEm_russ->GetBinContent(bin));
			}
			else if(eta!=0) cout << "eta = " << eta << endl;
		}//y
	}//x
	float av = 0;
	float sig= 0;
	int   n  = 0;
	int   bin= 0;
	
	n = mu_VPTPN_EEp_eta13.size(); av = avg(mu_VPTPN_EEp_eta13); sig = sigma(mu_VPTPN_EEp_eta13, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(1.31); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta13.size(); av = avg(mu_VPTPN_EEp_chin_eta13); sig = sigma(mu_VPTPN_EEp_chin_eta13, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(1.31); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta13.size(); av = avg(mu_VPTPN_EEp_russ_eta13); sig = sigma(mu_VPTPN_EEp_russ_eta13, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(1.31); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta15.size(); av = avg(mu_VPTPN_EEp_eta15); sig = sigma(mu_VPTPN_EEp_eta15, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(1.51); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta15.size(); av = avg(mu_VPTPN_EEp_chin_eta15); sig = sigma(mu_VPTPN_EEp_chin_eta15, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(1.51); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta15.size(); av = avg(mu_VPTPN_EEp_russ_eta15); sig = sigma(mu_VPTPN_EEp_russ_eta15, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(1.51); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta17.size(); av = avg(mu_VPTPN_EEp_eta17); sig = sigma(mu_VPTPN_EEp_eta17, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(1.71); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta17.size(); av = avg(mu_VPTPN_EEp_chin_eta17); sig = sigma(mu_VPTPN_EEp_chin_eta17, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(1.71); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta17.size(); av = avg(mu_VPTPN_EEp_russ_eta17); sig = sigma(mu_VPTPN_EEp_russ_eta17, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(1.71); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta19.size(); av = avg(mu_VPTPN_EEp_eta19); sig = sigma(mu_VPTPN_EEp_eta19, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(1.91); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta19.size(); av = avg(mu_VPTPN_EEp_chin_eta19); sig = sigma(mu_VPTPN_EEp_chin_eta19, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(1.91); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta19.size(); av = avg(mu_VPTPN_EEp_russ_eta19); sig = sigma(mu_VPTPN_EEp_russ_eta19, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(1.91); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta21.size(); av = avg(mu_VPTPN_EEp_eta21); sig = sigma(mu_VPTPN_EEp_eta21, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(2.11); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta21.size(); av = avg(mu_VPTPN_EEp_chin_eta21); sig = sigma(mu_VPTPN_EEp_chin_eta21, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(2.11); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta21.size(); av = avg(mu_VPTPN_EEp_russ_eta21); sig = sigma(mu_VPTPN_EEp_russ_eta21, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(2.11); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta23.size(); av = avg(mu_VPTPN_EEp_eta23); sig = sigma(mu_VPTPN_EEp_eta23, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(2.31); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta23.size(); av = avg(mu_VPTPN_EEp_chin_eta23); sig = sigma(mu_VPTPN_EEp_chin_eta23, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(2.31); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta23.size(); av = avg(mu_VPTPN_EEp_russ_eta23); sig = sigma(mu_VPTPN_EEp_russ_eta23, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(2.31); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta25.size(); av = avg(mu_VPTPN_EEp_eta25); sig = sigma(mu_VPTPN_EEp_eta25, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(2.51); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta25.size(); av = avg(mu_VPTPN_EEp_chin_eta25); sig = sigma(mu_VPTPN_EEp_chin_eta25, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(2.51); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta25.size(); av = avg(mu_VPTPN_EEp_russ_eta25); sig = sigma(mu_VPTPN_EEp_russ_eta25, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(2.51); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta27.size(); av = avg(mu_VPTPN_EEp_eta27); sig = sigma(mu_VPTPN_EEp_eta27, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(2.71); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta27.size(); av = avg(mu_VPTPN_EEp_chin_eta27); sig = sigma(mu_VPTPN_EEp_chin_eta27, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(2.71); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta27.size(); av = avg(mu_VPTPN_EEp_russ_eta27); sig = sigma(mu_VPTPN_EEp_russ_eta27, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(2.71); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta29.size(); av = avg(mu_VPTPN_EEp_eta29); sig = sigma(mu_VPTPN_EEp_eta29, av);
	bin = mu_VPTPN_EEp_entries_2->FindBin(2.91); mu_VPTPN_EEp_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_2->SetBinContent(bin, av); mu_VPTPN_EEp_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta29.size(); av = avg(mu_VPTPN_EEp_chin_eta29); sig = sigma(mu_VPTPN_EEp_chin_eta29, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_2->FindBin(2.91); mu_VPTPNnorm_EEp_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta29.size(); av = avg(mu_VPTPN_EEp_russ_eta29); sig = sigma(mu_VPTPN_EEp_russ_eta29, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_2->FindBin(2.91); mu_VPTPNnorm_EEp_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_2->SetBinError(bin, sig);


	n = mu_VPTPN_EEp_eta14.size(); av = avg(mu_VPTPN_EEp_eta14); sig = sigma(mu_VPTPN_EEp_eta14, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(1.41); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta14.size(); av = avg(mu_VPTPN_EEp_chin_eta14); sig = sigma(mu_VPTPN_EEp_chin_eta14, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(1.41); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta14.size(); av = avg(mu_VPTPN_EEp_russ_eta14); sig = sigma(mu_VPTPN_EEp_russ_eta14, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(1.41); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta16.size(); av = avg(mu_VPTPN_EEp_eta16); sig = sigma(mu_VPTPN_EEp_eta16, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(1.61); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta16.size(); av = avg(mu_VPTPN_EEp_chin_eta16); sig = sigma(mu_VPTPN_EEp_chin_eta16, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(1.61); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta16.size(); av = avg(mu_VPTPN_EEp_russ_eta16); sig = sigma(mu_VPTPN_EEp_russ_eta16, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(1.61); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta18.size(); av = avg(mu_VPTPN_EEp_eta18); sig = sigma(mu_VPTPN_EEp_eta18, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(1.81); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta18.size(); av = avg(mu_VPTPN_EEp_chin_eta18); sig = sigma(mu_VPTPN_EEp_chin_eta18, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(1.81); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta18.size(); av = avg(mu_VPTPN_EEp_russ_eta18); sig = sigma(mu_VPTPN_EEp_russ_eta18, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(1.81); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta20.size(); av = avg(mu_VPTPN_EEp_eta20); sig = sigma(mu_VPTPN_EEp_eta20, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(2.01); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta20.size(); av = avg(mu_VPTPN_EEp_chin_eta20); sig = sigma(mu_VPTPN_EEp_chin_eta20, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(2.01); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta20.size(); av = avg(mu_VPTPN_EEp_russ_eta20); sig = sigma(mu_VPTPN_EEp_russ_eta20, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(2.01); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta22.size(); av = avg(mu_VPTPN_EEp_eta22); sig = sigma(mu_VPTPN_EEp_eta22, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(2.21); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta22.size(); av = avg(mu_VPTPN_EEp_chin_eta22); sig = sigma(mu_VPTPN_EEp_chin_eta22, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(2.21); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta22.size(); av = avg(mu_VPTPN_EEp_russ_eta22); sig = sigma(mu_VPTPN_EEp_russ_eta22, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(2.21); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta24.size(); av = avg(mu_VPTPN_EEp_eta24); sig = sigma(mu_VPTPN_EEp_eta24, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(2.41); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta24.size(); av = avg(mu_VPTPN_EEp_chin_eta24); sig = sigma(mu_VPTPN_EEp_chin_eta24, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(2.41); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta24.size(); av = avg(mu_VPTPN_EEp_russ_eta24); sig = sigma(mu_VPTPN_EEp_russ_eta24, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(2.41); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta26.size(); av = avg(mu_VPTPN_EEp_eta26); sig = sigma(mu_VPTPN_EEp_eta26, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(2.61); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta26.size(); av = avg(mu_VPTPN_EEp_chin_eta26); sig = sigma(mu_VPTPN_EEp_chin_eta26, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(2.61); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta26.size(); av = avg(mu_VPTPN_EEp_russ_eta26); sig = sigma(mu_VPTPN_EEp_russ_eta26, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(2.61); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEp_eta28.size(); av = avg(mu_VPTPN_EEp_eta28); sig = sigma(mu_VPTPN_EEp_eta28, av);
	bin = mu_VPTPN_EEp_entries_1->FindBin(2.81); mu_VPTPN_EEp_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEp_average_1->SetBinContent(bin, av); mu_VPTPN_EEp_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_chin_eta28.size(); av = avg(mu_VPTPN_EEp_chin_eta28); sig = sigma(mu_VPTPN_EEp_chin_eta28, av);
	bin = mu_VPTPNnorm_EEp_chin_entries_1->FindBin(2.81); mu_VPTPNnorm_EEp_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEp_russ_eta28.size(); av = avg(mu_VPTPN_EEp_russ_eta28); sig = sigma(mu_VPTPN_EEp_russ_eta28, av);
	bin = mu_VPTPNnorm_EEp_russ_entries_1->FindBin(2.81); mu_VPTPNnorm_EEp_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEp_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEp_russ_average_1->SetBinError(bin, sig);



	n = mu_VPTPN_EEm_eta31.size(); av = avg(mu_VPTPN_EEm_eta31); sig = sigma(mu_VPTPN_EEm_eta31, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-3.09); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta31.size(); av = avg(mu_VPTPN_EEm_chin_eta31); sig = sigma(mu_VPTPN_EEm_chin_eta31, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-3.09); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta31.size(); av = avg(mu_VPTPN_EEm_russ_eta31); sig = sigma(mu_VPTPN_EEm_russ_eta31, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-3.09); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta29.size(); av = avg(mu_VPTPN_EEm_eta29); sig = sigma(mu_VPTPN_EEm_eta29, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-2.89); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta29.size(); av = avg(mu_VPTPN_EEm_chin_eta29); sig = sigma(mu_VPTPN_EEm_chin_eta29, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-2.89); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta29.size(); av = avg(mu_VPTPN_EEm_russ_eta29); sig = sigma(mu_VPTPN_EEm_russ_eta29, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-2.89); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta27.size(); av = avg(mu_VPTPN_EEm_eta27); sig = sigma(mu_VPTPN_EEm_eta27, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-2.69); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta27.size(); av = avg(mu_VPTPN_EEm_chin_eta27); sig = sigma(mu_VPTPN_EEm_chin_eta27, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-2.69); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta27.size(); av = avg(mu_VPTPN_EEm_russ_eta27); sig = sigma(mu_VPTPN_EEm_russ_eta27, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-2.69); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta25.size(); av = avg(mu_VPTPN_EEm_eta25); sig = sigma(mu_VPTPN_EEm_eta25, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-2.49); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta25.size(); av = avg(mu_VPTPN_EEm_chin_eta25); sig = sigma(mu_VPTPN_EEm_chin_eta25, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-2.49); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta25.size(); av = avg(mu_VPTPN_EEm_russ_eta25); sig = sigma(mu_VPTPN_EEm_russ_eta25, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-2.49); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta23.size(); av = avg(mu_VPTPN_EEm_eta23); sig = sigma(mu_VPTPN_EEm_eta23, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-2.29); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta23.size(); av = avg(mu_VPTPN_EEm_chin_eta23); sig = sigma(mu_VPTPN_EEm_chin_eta23, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-2.29); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta23.size(); av = avg(mu_VPTPN_EEm_russ_eta23); sig = sigma(mu_VPTPN_EEm_russ_eta23, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-2.29); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta21.size(); av = avg(mu_VPTPN_EEm_eta21); sig = sigma(mu_VPTPN_EEm_eta21, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-2.09); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta21.size(); av = avg(mu_VPTPN_EEm_chin_eta21); sig = sigma(mu_VPTPN_EEm_chin_eta21, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-2.09); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta21.size(); av = avg(mu_VPTPN_EEm_russ_eta21); sig = sigma(mu_VPTPN_EEm_russ_eta21, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-2.09); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta19.size(); av = avg(mu_VPTPN_EEm_eta19); sig = sigma(mu_VPTPN_EEm_eta19, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-1.89); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta19.size(); av = avg(mu_VPTPN_EEm_chin_eta19); sig = sigma(mu_VPTPN_EEm_chin_eta19, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-1.89); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta19.size(); av = avg(mu_VPTPN_EEm_russ_eta19); sig = sigma(mu_VPTPN_EEm_russ_eta19, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-1.89); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta17.size(); av = avg(mu_VPTPN_EEm_eta17); sig = sigma(mu_VPTPN_EEm_eta17, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-1.69); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta17.size(); av = avg(mu_VPTPN_EEm_chin_eta17); sig = sigma(mu_VPTPN_EEm_chin_eta17, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-1.69); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta17.size(); av = avg(mu_VPTPN_EEm_russ_eta17); sig = sigma(mu_VPTPN_EEm_russ_eta17, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-1.69); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta15.size(); av = avg(mu_VPTPN_EEm_eta15); sig = sigma(mu_VPTPN_EEm_eta15, av);
	bin = mu_VPTPN_EEm_entries_2->FindBin(-1.49); mu_VPTPN_EEm_entries_2->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_2->SetBinContent(bin, av); mu_VPTPN_EEm_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta15.size(); av = avg(mu_VPTPN_EEm_chin_eta15); sig = sigma(mu_VPTPN_EEm_chin_eta15, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_2->FindBin(-1.49); mu_VPTPNnorm_EEm_chin_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_2->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta15.size(); av = avg(mu_VPTPN_EEm_russ_eta15); sig = sigma(mu_VPTPN_EEm_russ_eta15, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_2->FindBin(-1.49); mu_VPTPNnorm_EEm_russ_entries_2->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_2->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_2->SetBinError(bin, sig);


	n = mu_VPTPN_EEm_eta30.size(); av = avg(mu_VPTPN_EEm_eta30); sig = sigma(mu_VPTPN_EEm_eta30, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-2.99); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta30.size(); av = avg(mu_VPTPN_EEm_chin_eta30); sig = sigma(mu_VPTPN_EEm_chin_eta30, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-2.99); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta30.size(); av = avg(mu_VPTPN_EEm_russ_eta30); sig = sigma(mu_VPTPN_EEm_russ_eta30, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-2.99); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta28.size(); av = avg(mu_VPTPN_EEm_eta28); sig = sigma(mu_VPTPN_EEm_eta28, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-2.79); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta28.size(); av = avg(mu_VPTPN_EEm_chin_eta28); sig = sigma(mu_VPTPN_EEm_chin_eta28, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-2.79); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta28.size(); av = avg(mu_VPTPN_EEm_russ_eta28); sig = sigma(mu_VPTPN_EEm_russ_eta28, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-2.79); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta26.size(); av = avg(mu_VPTPN_EEm_eta26); sig = sigma(mu_VPTPN_EEm_eta26, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-2.59); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta26.size(); av = avg(mu_VPTPN_EEm_chin_eta26); sig = sigma(mu_VPTPN_EEm_chin_eta26, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-2.59); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta26.size(); av = avg(mu_VPTPN_EEm_russ_eta26); sig = sigma(mu_VPTPN_EEm_russ_eta26, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-2.59); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta24.size(); av = avg(mu_VPTPN_EEm_eta24); sig = sigma(mu_VPTPN_EEm_eta24, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-2.39); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta24.size(); av = avg(mu_VPTPN_EEm_chin_eta24); sig = sigma(mu_VPTPN_EEm_chin_eta24, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-2.39); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta24.size(); av = avg(mu_VPTPN_EEm_russ_eta24); sig = sigma(mu_VPTPN_EEm_russ_eta24, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-2.39); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta22.size(); av = avg(mu_VPTPN_EEm_eta22); sig = sigma(mu_VPTPN_EEm_eta22, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-2.19); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta22.size(); av = avg(mu_VPTPN_EEm_chin_eta22); sig = sigma(mu_VPTPN_EEm_chin_eta22, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-2.19); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta22.size(); av = avg(mu_VPTPN_EEm_russ_eta22); sig = sigma(mu_VPTPN_EEm_russ_eta22, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-2.19); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta20.size(); av = avg(mu_VPTPN_EEm_eta20); sig = sigma(mu_VPTPN_EEm_eta20, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-1.99); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta20.size(); av = avg(mu_VPTPN_EEm_chin_eta20); sig = sigma(mu_VPTPN_EEm_chin_eta20, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-1.99); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta20.size(); av = avg(mu_VPTPN_EEm_russ_eta20); sig = sigma(mu_VPTPN_EEm_russ_eta20, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-1.99); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta18.size(); av = avg(mu_VPTPN_EEm_eta18); sig = sigma(mu_VPTPN_EEm_eta18, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-1.79); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta18.size(); av = avg(mu_VPTPN_EEm_chin_eta18); sig = sigma(mu_VPTPN_EEm_chin_eta18, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-1.79); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta18.size(); av = avg(mu_VPTPN_EEm_russ_eta18); sig = sigma(mu_VPTPN_EEm_russ_eta18, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-1.79); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);

	n = mu_VPTPN_EEm_eta16.size(); av = avg(mu_VPTPN_EEm_eta16); sig = sigma(mu_VPTPN_EEm_eta16, av);
	bin = mu_VPTPN_EEm_entries_1->FindBin(-1.59); mu_VPTPN_EEm_entries_1->SetBinContent(bin, n);
	mu_VPTPN_EEm_average_1->SetBinContent(bin, av); mu_VPTPN_EEm_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_chin_eta16.size(); av = avg(mu_VPTPN_EEm_chin_eta16); sig = sigma(mu_VPTPN_EEm_chin_eta16, av);
	bin = mu_VPTPNnorm_EEm_chin_entries_1->FindBin(-1.59); mu_VPTPNnorm_EEm_chin_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_chin_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_chin_average_1->SetBinError(bin, sig);
	n = mu_VPTPN_EEm_russ_eta16.size(); av = avg(mu_VPTPN_EEm_russ_eta16); sig = sigma(mu_VPTPN_EEm_russ_eta16, av);
	bin = mu_VPTPNnorm_EEm_russ_entries_1->FindBin(-1.59); mu_VPTPNnorm_EEm_russ_entries_1->SetBinContent(bin, n);
	mu_VPTPNnorm_EEm_russ_average_1->SetBinContent(bin, av); mu_VPTPNnorm_EEm_russ_average_1->SetBinError(bin, sig);


	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120224_EE_mu_eta_averages.root", "RECREATE");
	newFile->cd();
	mu_VPTPN_EEp_entries_1->Write();
	mu_VPTPN_EEp_entries_2->Write();
	mu_VPTPN_EEm_entries_1->Write();
	mu_VPTPN_EEm_entries_2->Write();
	mu_VPTPN_EEp_average_1->Write();
	mu_VPTPN_EEp_average_2->Write();
	mu_VPTPN_EEm_average_1->Write();
	mu_VPTPN_EEm_average_2->Write();
	mu_VPTPNnorm_EEp_chin_entries_1->Write();
	mu_VPTPNnorm_EEp_chin_entries_2->Write();
	mu_VPTPNnorm_EEm_chin_entries_1->Write();
	mu_VPTPNnorm_EEm_chin_entries_2->Write();
	mu_VPTPNnorm_EEp_chin_average_1->Write();
	mu_VPTPNnorm_EEp_chin_average_2->Write();
	mu_VPTPNnorm_EEm_chin_average_1->Write();
	mu_VPTPNnorm_EEm_chin_average_2->Write();
	mu_VPTPNnorm_EEp_russ_entries_1->Write();
	mu_VPTPNnorm_EEp_russ_entries_2->Write();
	mu_VPTPNnorm_EEm_russ_entries_1->Write();
	mu_VPTPNnorm_EEm_russ_entries_2->Write();
	mu_VPTPNnorm_EEp_russ_average_1->Write();
	mu_VPTPNnorm_EEp_russ_average_2->Write();
	mu_VPTPNnorm_EEm_russ_average_1->Write();
	mu_VPTPNnorm_EEm_russ_average_2->Write();

}

//}


//}

*/