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
#include <TH2I.h>
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
#include <TBranch.h>

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

using namespace std;

void load();

//vector<TFile *> files;
vector<string> files;

double avg(vector<double> x){

  double sum = 0;
  for (Int_t i = 0; i<x.size(); ++i){
    sum += x[i];
  }
  if (x.size()!=0)  return sum/(float)x.size();
  return 0;
}
double sigma(vector<double> x, double mean){

  double sum = 0;
  for (Int_t i = 0; i<x.size(); ++i){
    sum += (x[i]-mean)*(x[i]-mean);
  }
  if (x.size()!=0)  return sqrt(sum/(float)x.size());
  return 0;
}

struct EEgeo{//tagged jets

	//Default constructor will be needed if inserting into vector
	EEgeo() : ix(), iy(), iz(), eta(), detId(), fed(), elecId(), harness() {}
	EEgeo( int x, int y, int z, float e, int d, int f, int eI, int h ) : ix(x), iy(y), iz(z), eta(e), detId(d), fed(f), elecId(eI), harness(h) {}
	//sorting needs operator < to compile
	bool operator<( EEgeo const& rhs ) const
	{ return eta < rhs.eta; }//ascending sorting by discr_value

	int ix;
	int iy;
	int iz;
	float eta;
	float R(){ 
		return sqrt(ix*ix + iy*iy);
	}
	int detId;
	int fed;
	int elecId;
	int harness;

};

struct EBgeo{//tagged jets

	//Default constructor will be needed if inserting into vector
	EBgeo() : ieta(), iphi(),eta(), detId(), fed(), elecId(), harness() {}
	EBgeo( int x, int y, float e, int d, int f, int eI, int h ) : ieta(x), iphi(y), eta(e), detId(d), fed(f), elecId(eI), harness(h) {}
	//sorting needs operator < to compile
	bool operator<( EBgeo const& rhs ) const
	{ return eta < rhs.eta; }//ascending sorting by discr_value

	int ieta;
	int iphi;
	float eta;
	int detId;
	int fed;
	int elecId;
	int harness;

};

/*struct EEgeo{
	int ix;
	int iy;
	int iz;
	float eta;
};*/

bool EEetasort (EEgeo i,EEgeo j) { return (i.eta<j.eta); }//sort by lowest eta
bool EBetasort (EBgeo i,EBgeo j) { return (i.eta<j.eta); }//sort by lowest eta


void EEstudies(){

	vector<EEgeo> EEgeoVec;
	vector<EBgeo> EBgeoVec;
	//vector<int> ixxx;
	//vector<int> iyyy;
	//vector<int> izzz;
	//vector<float> eeta;


	//load files
	load();
/*
//	TFile *stdmapsfile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20111219_EEplots.root");
	TFile *stdmapsfile = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/Plots/20120403_EEmustdplots.root");
	TH2D *mu_ECAL_EEp_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_russ");
	TH2D *mu_ECAL_EEp_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEp_chin");
	TH2D *mu_ECAL_EEm_russ = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_russ");
	TH2D *mu_ECAL_EEm_chin = (TH2D*)stdmapsfile->Get("mu_ECAL_EEm_chin");

	int steptime = 1314807448;
	int lasttime = 1315398150;
	vector<float > firstaverage;
	vector<float > lastaverage;
	float average1;
	float average2;

	TH2F *VPTPN_firstaverage_EEp = new TH2F("VPTPN_firstaverage_EEp", "VPTPN_firstaverage_EEp", 101, 0, 101, 101, 0, 101);
	TH2F *VPTPN_firstaverage_EEm = new TH2F("VPTPN_firstaverage_EEm", "VPTPN_firstaverage_EEm", 101, 0, 101, 101, 0, 101);
	TH2F *VPTPN_lastaverage_EEp  = new TH2F("VPTPN_lastaverage_EEp" , "VPTPN_lastaverage_EEp" , 101, 0, 101, 101, 0, 101);
	TH2F *VPTPN_lastaverage_EEm  = new TH2F("VPTPN_lastaverage_EEm" , "VPTPN_lastaverage_EEm" , 101, 0, 101, 101, 0, 101);

	TH2F *mu_VPTPN_EEp           = new TH2F("mu_VPTPN_EEp", "mu_VPTPN_EEp", 101, 0, 101, 101, 0, 101);
	TH2F *mu_VPTPN_EEm           = new TH2F("mu_VPTPN_EEm", "mu_VPTPN_EEm", 101, 0, 101, 101, 0, 101);

	TH2F *mu_VPTPNnorm_EEp_chin  = new TH2F("mu_VPTPNnorm_EEp_chin", "mu_VPTPNnorm_EEp_chin", 101, 0, 101, 101, 0, 101);
	TH2F *mu_VPTPNnorm_EEp_russ  = new TH2F("mu_VPTPNnorm_EEp_russ", "mu_VPTPNnorm_EEp_russ", 101, 0, 101, 101, 0, 101);
	TH2F *mu_VPTPNnorm_EEm_chin  = new TH2F("mu_VPTPNnorm_EEm_chin", "mu_VPTPNnorm_EEm_chin", 101, 0, 101, 101, 0, 101);
	TH2F *mu_VPTPNnorm_EEm_russ  = new TH2F("mu_VPTPNnorm_EEm_russ", "mu_VPTPNnorm_EEm_russ", 101, 0, 101, 101, 0, 101);
*/
/* //detId etc. connections
	TH2F *eta_EEp                = new TH2F("eta_EEp", "eta_EEp", 101, 0, 101, 101, 0, 101);
	TH2F *eta_EEm                = new TH2F("eta_EEm", "eta_EEm", 101, 0, 101, 101, 0, 101);
	TH2I *detId_EEp              = new TH2I("detId_EEp", "detId_EEp", 101, 0, 101, 101, 0, 101);
	TH2I *detId_EEm              = new TH2I("detId_EEm", "detId_EEm", 101, 0, 101, 101, 0, 101);
	TH2I *fed_EEp                = new TH2I("fed_EEp", "fed_EEp", 101, 0, 101, 101, 0, 101);
	TH2I *fed_EEm                = new TH2I("fed_EEm", "fed_EEm", 101, 0, 101, 101, 0, 101);
	TH2I *elecId_EEp             = new TH2I("elecId_EEp", "elecId_EEp", 101, 0, 101, 101, 0, 101);
	TH2I *elecId_EEm             = new TH2I("elecId_EEm", "elecId_EEm", 101, 0, 101, 101, 0, 101);
	TH2I *harness_EEp            = new TH2I("harness_EEp", "harness_EEp", 101, 0, 101, 101, 0, 101);
	TH2I *harness_EEm            = new TH2I("harness_EEm", "harness_EEm", 101, 0, 101, 101, 0, 101);

	TH1I *detId_ieta_EB          = new TH1I("detId_ieta_EB", "detId_ieta_EB", 108917, 838861306, 838970223);
	TH1I *detId_iphi_EB          = new TH1I("detId_iphi_EB", "detId_iphi_EB", 108917, 838861306, 838970223);
	TH1I *detId_ix_EEm           = new TH1I("detId_ix_EEm", "detId_ix_EEm",   12694,  872415400, 872428094);
	TH1I *detId_iy_EEm           = new TH1I("detId_iy_EEm", "detId_iy_EEm",   12694,  872415400, 872428094);
	TH1I *detId_ix_EEp           = new TH1I("detId_ix_EEp", "detId_ix_EEp",   12694,  872431783, 872444477);
	TH1I *detId_iy_EEp           = new TH1I("detId_iy_EEp", "detId_iy_EEp",   12694,  872431783, 872444477);

	TH2F *eta_EB                 = new TH2F("eta_EB", "eta_EB", 171, -85, 86, 361, 0, 361);
	TH2I *detId_EB               = new TH2I("detId_EB", "detId_EB", 171, -85, 86, 361, 0, 361);
	TH2I *fed_EB                 = new TH2I("fed_EB", "fed_EB", 171, -86, 85, 361, 0, 361);
	TH2I *elecId_EB              = new TH2I("elecId_EB", "elecId_EB", 171, -85, 86, 361, 0, 361);
	TH2I *harness_EB             = new TH2I("harness_EB", "harness_EB", 171, -85, 86, 361, 0, 361);
*/

//	TGraph* g_vptpn_las[101][101][2];
//	TGraph* g_vptpn_led[101][101][2];
	TGraph* g_vptpn_oled[101][101][2];
//	TGraph* g_vptpn_rat[101][101][2];
//	TGraph* g_vptpn_orat[101][101][2];
//	TGraph* g_EE_ampl[101][101][2];
//	TGraph* g_EE_fwhm[101][101][2];
/*
	TGraph* g_apdpn_las[171][361];
	TGraph* g_apdpn_led[171][361];
	TGraph* g_apdpn_rat[171][361];
	TGraph* g_EB_ampl[171][361];
	TGraph* g_EB_fwhm[171][361];
*/

//	TGraphErrors* g_vptpn_las_day[101][101][2];
	TGraphErrors* g_vptpn_oled_day[101][101][2];

	char gname[101];
	bool errors = true;//this needs also change of g_vptpn_las_day[101][101][2];
	bool saveAllChannels=true;
	bool saveeverypoint = false;
	bool saveoneday     = true;//save one (avg) point per day
	int oneday = 86400;
	vector<double> (value[101][101][2]);
	int daystart[101][101][2];
	
//	TGraph *lumi_vs_time;

	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
	      (value[i][j][k]).clear();
	      daystart[i][j][k] = -1;
	      if (saveAllChannels){
	      if(saveeverypoint){
/*		if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEp", i, j);
		g_vptpn_las[i][j][k] = new TGraph();
		g_vptpn_las[i][j][k] -> SetTitle(gname);
		g_vptpn_las[i][j][k] -> SetName(gname);
		g_vptpn_las[i][j][k] -> SetMarkerColor(kBlue);
		g_vptpn_las[i][j][k] -> SetMarkerStyle(22);
		
		if (k==0) sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEp", i, j);
		g_vptpn_led[i][j][k] = new TGraph();
		g_vptpn_led[i][j][k] -> SetTitle(gname);
		g_vptpn_led[i][j][k] -> SetName(gname);
		g_vptpn_led[i][j][k] -> SetMarkerColor(kCyan);
		g_vptpn_led[i][j][k] -> SetMarkerStyle(23);
		//g_vptpn_led[i][j][k] -> SetMarkerSize(0.5);
*/

		if (k==0) sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEp", i, j);
		g_vptpn_oled[i][j][k] = new TGraph();
		g_vptpn_oled[i][j][k] -> SetTitle(gname);
		g_vptpn_oled[i][j][k] -> SetName(gname);
		g_vptpn_oled[i][j][k] -> SetMarkerColor(kOrange);
		g_vptpn_oled[i][j][k] -> SetMarkerStyle(7);
/*
		if (k==0) sprintf(gname,"g_vptpn_rat_ix%d_iy%d_EEm", i, j);//ratio las/led
		else sprintf(gname,"g_vptpn_rat_ix%d_iy%d_EEp", i, j);
		g_vptpn_orat[i][j][k] = new TGraph();
		g_vptpn_orat[i][j][k] -> SetTitle(gname);
		g_vptpn_orat[i][j][k] -> SetName(gname);
		g_vptpn_orat[i][j][k] -> SetMarkerColor(kBlack);
		g_vptpn_orat[i][j][k] -> SetMarkerStyle(20);

		if (k==0) sprintf(gname,"g_vptpn_orat_ix%d_iy%d_EEm", i, j);//ratio las/oled
		else sprintf(gname,"g_vptpn_orat_ix%d_iy%d_EEp", i, j);
		g_vptpn_rat[i][j][k] = new TGraph();
		g_vptpn_rat[i][j][k] -> SetTitle(gname);
		g_vptpn_rat[i][j][k] -> SetName(gname);
		g_vptpn_rat[i][j][k] -> SetMarkerColor(kBlack);
		g_vptpn_rat[i][j][k] -> SetMarkerStyle(20);

		if (k==0) sprintf(gname,"g_EE_ampl_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_EE_ampl_ix%d_iy%d_EEp", i, j);
		g_EE_ampl[i][j][k] = new TGraph();
		g_EE_ampl[i][j][k] -> SetTitle(gname);
		g_EE_ampl[i][j][k] -> SetName(gname);
		g_EE_ampl[i][j][k] -> SetMarkerColor(kGreen+3);
		g_EE_ampl[i][j][k] -> SetMarkerStyle(22);

		if (k==0) sprintf(gname,"g_EE_fwhm_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_EE_fwhm_ix%d_iy%d_EEp", i, j);
		g_EE_fwhm[i][j][k] = new TGraph();
		g_EE_fwhm[i][j][k] -> SetTitle(gname);
		g_EE_fwhm[i][j][k] -> SetName(gname);
		g_EE_fwhm[i][j][k] -> SetMarkerColor(kMagenta+2);
		g_EE_fwhm[i][j][k] -> SetMarkerStyle(23);*/
	      }
	      if(saveoneday){
		if (k==0) sprintf(gname,"tg_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"tg_vptpn_oled_ix%d_iy%d_EEp", i, j);
		g_vptpn_oled_day[i][j][k] = new TGraphErrors();
		g_vptpn_oled_day[i][j][k] -> SetTitle(gname);
		g_vptpn_oled_day[i][j][k] -> SetName(gname);
		g_vptpn_oled_day[i][j][k] -> SetMarkerColor(kOrange);
		g_vptpn_oled_day[i][j][k] -> SetMarkerStyle(7);
              }
	      }
	    }
	  }
	  }
		//EB
	/*  for (int i = 0; i < 171; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 361; ++j){//as iX goes from 1 to 100
	      if (saveAllChannels){
		sprintf(gname,"g_apdpn_las_ieta%d_iphi%d", i, j);
		g_apdpn_las[i][j] = new TGraph();
		g_apdpn_las[i][j] -> SetTitle(gname);
		g_apdpn_las[i][j] -> SetName(gname);
		g_apdpn_las[i][j] -> SetMarkerColor(kBlue);
		g_apdpn_las[i][j] -> SetMarkerStyle(22);
		
		sprintf(gname,"g_apdpn_led_ieta%d_iphi%d", i, j);
		g_apdpn_led[i][j] = new TGraph();
		g_apdpn_led[i][j] -> SetTitle(gname);
		g_apdpn_led[i][j] -> SetName(gname);
		g_apdpn_led[i][j] -> SetMarkerColor(kRed);
		g_apdpn_led[i][j] -> SetMarkerStyle(23);

		sprintf(gname,"g_apdpn_rat_ieta%d_iphi%d", i, j);
		g_apdpn_rat[i][j] = new TGraph();
		g_apdpn_rat[i][j] -> SetTitle(gname);
		g_apdpn_rat[i][j] -> SetName(gname);
		g_apdpn_rat[i][j] -> SetMarkerColor(kBlack);
		g_apdpn_rat[i][j] -> SetMarkerStyle(20);

		sprintf(gname,"g_EB_ampl_ieta%d_iphi%d", i, j);
		g_EB_ampl[i][j] = new TGraph();
		g_EB_ampl[i][j] -> SetTitle(gname);
		g_EB_ampl[i][j] -> SetName(gname);
		g_EB_ampl[i][j] -> SetMarkerColor(kGreen+3);
		g_EB_ampl[i][j] -> SetMarkerStyle(22);

		sprintf(gname,"g_EB_fwhm_ieta%d_iphi%d", i, j);
		g_EB_fwhm[i][j] = new TGraph();
		g_EB_fwhm[i][j] -> SetTitle(gname);
		g_EB_fwhm[i][j] -> SetName(gname);
		g_EB_fwhm[i][j] -> SetMarkerColor(kMagenta+2);
		g_EB_fwhm[i][j] -> SetMarkerStyle(23);

	      }
	    }
	  }*/

	


//	TFile *f = TFile::Open("/shome/haweber/ECAL/DataLaser/small_ntu_data_00173378-00176049.root");
//	f->Add((TFile*)TFile::Open("/shome/haweber/ECAL/DataLaser/small_ntu_data_00173378-00176049.root"), );
//	TTree *tr = (TTree*)f->Get("x");

	for(int fileindex=0; fileindex<files.size(); ++fileindex){
	cout << endl << "*************************************************" << endl;
	cout << "Processing file " << fileindex + 1 << "/" << files.size() << endl;

	TFile *file = TFile::Open((files[fileindex]).c_str());

	//TTree *tr = (TTree*)files[fileindex]->Get("x");
	TTree *tr = (TTree*)file->Get("x");


	
	LaserAnalysisReducedTrees x((TTree*)tr);

   	Long64_t nentries = tr->GetEntriesFast();
	cout << "Tree contains "<< nentries << " entries" << endl;

   	Long64_t nbytes = 0, nb = 0;
   	    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      		Long64_t ientry = x.LoadTree(jentry);
      		if (ientry < 0) break;
		if(jentry%(nentries/25)==0) cout << "Processing event " << jentry << endl;
		nb = x.GetEntry(jentry);   nbytes += nb;

		int iX = x.ix;
		int iY = x.iy;
		int iZ = x.iz;
/*		if(iZ!=0){
			iX = iX + 50 ;
			if(x.ix < 0) iX = iX + 1 ;
			iY = iY * iZ;
		}
		if(iZ==0){
			iY = iY + 10;
			if(x.ix < 0) iY = iY + 1;
		}

		if( ((iX<1 || iX>100 || iY<1 || iY>100)&&iZ!=0) || (iZ<-1 || iZ>1) || (iZ==0&&(iY<0 || iY>360 || iX<-85 || iX>85)) ) 
			cout << "iX " << iX << " iY " << iY << " iZ " << iZ << " ix " << x.ix << " iy " << x.iy << " iz " << x.iz << endl;
*/
/* //detID etc connections
		if(iZ==0) {
			bool add2 = true;
			for(int l = 0; l<EBgeoVec.size(); ++l){
				if((EBgeoVec[l]).ieta==iX && (EBgeoVec[l]).iphi==iY) {
					add2 = false;
					break;
				}
			}
			if(add2){
				EBgeo help2;
				help2.ieta = iX;
				help2.iphi = iY;
				help2.eta = x.eta;
				help2.detId = x.detId;
				help2.fed = x.fed;
				help2.elecId = x.elecId;
				help2.harness = x.harness;
				EBgeoVec.push_back(help2);
			}
			continue;//don't run over barrel
		}

		bool add = true;
		for(int l = 0; l<EEgeoVec.size(); ++l){
			if((EEgeoVec[l]).ix==iX && (EEgeoVec[l]).iy==iY && (EEgeoVec[l]).iz==iZ) {
				add = false;
				break;
			}
		}
		if(add){
			EEgeo help;
			help.ix = iX;
			help.iy = iY;
			help.iz = iZ;
			help.eta = x.eta;
			help.detId = x.detId;
			help.fed = x.fed;
			help.elecId = x.elecId;
			help.harness = x.harness;
			//ixxx.push_back(iX);
			//iyyy.push_back(iY);
			//izzz.push_back(iZ);
			//eeta.push_back(x.eta);
			EEgeoVec.push_back(help);
		}
*/
		float apdpnlas  = x.apdpnAB[0];
		float apdpnled  = x.apdpnAB[1];
		float apdpnoled = x.apdpnAB[2];
		float apdpnrat  = -99.9;
		float apdpnorat = -99.9;
		float ampl      = x.l_ampli[0];
		float fwhm      = x.l_fwhm[0];
		if(apdpnled>0) apdpnrat = apdpnlas/apdpnled;
		if(apdpnoled>0)apdpnorat= apdpnlas/apdpnled;
		int tlas       = x.time[0];
		int tled       = x.time[1];
		int toled      = x.time[2];
			if (saveAllChannels){
				if(saveeverypoint){
				if(iZ!=0){
				int k = iZ;
				if(iZ<0) k = 0;
		//		if(apdpnlas>0)  g_vptpn_las[iX][iY][k]  -> SetPoint(g_vptpn_las[iX][iY][k]->GetN(),  tlas , apdpnlas);
		//		if(apdpnled>0)  g_vptpn_led[iX][iY][k]  -> SetPoint(g_vptpn_led[iX][iY][k]->GetN(),  tled , apdpnled);
				if(apdpnoled>0) g_vptpn_oled[iX][iY][k] -> SetPoint(g_vptpn_oled[iX][iY][k]->GetN(), toled, apdpnoled);
			//	if(apdpnrat>0)  g_vptpn_rat[iX][iY][k]  -> SetPoint(g_vptpn_rat[iX][iY][k]->GetN(),  tlas , apdpnrat);
			//	if(apdpnorat>0) g_vptpn_rat[iX][iY][k]  -> SetPoint(g_vptpn_orat[iX][iY][k]->GetN(), tlas , apdpnorat);
			//	if(ampl>0)      g_EE_ampl[iX][iY][k]    -> SetPoint(g_EE_ampl[iX][iY][k]->GetN(),    tlas , ampl);
			//	if(fwhm>0)      g_EE_fwhm[iX][iY][k]    -> SetPoint(g_EE_fwhm[iX][iY][k]->GetN(),    tlas , fwhm);
				}
				//EB
				/*else{
				int iXX = iX+85;//as ix=ieta starts at -85
				if(apdpnlas>0) g_apdpn_las[iXX][iY] -> SetPoint(g_apdpn_las[iXX][iY]->GetN(), tlas , apdpnlas);
				}*/
				}
				if(saveoneday){
				if(iZ!=0){
					int k = iZ;
					if(iZ<0) k = 0;
					if(daystart[iX][iY][k]<=0){
					daystart[iX][iY][k]           = x.time[0];
					if(iZ!=0) daystart[iX][iY][k] = x.time[73];
					}
					if((tlas-daystart[iX][iY][k])<0) cout << "Error, tlas-daystart " << (int)tlas << "-" << (int)daystart[iX][iY][k] << " = " << int(tlas-daystart[iX][iY][k]) << endl;
				//	if((tlas-daystart[iX][iY][k])<oneday){ (value[iX][iY][k]).push_back(apdpnlas); }
					if((toled-daystart[iX][iY][k])<oneday){  if(apdpnoled>0) (value[iX][iY][k]).push_back(apdpnoled); }
						//if(iX==3&&iY==44) cout << "good " << tlas << " " << daystart[iX][iY][k] << endl;}
					else {
						//if(iX==3&&iY==44) cout << "bad  " << tlas << " " << daystart[iX][iY][k] << " fullday " << oneday <<  endl;
						double average = avg(value[iX][iY][k]);
						double sigmaav = sigma(value[iX][iY][k], average);
					//	int NNN = g_vptpn_las_day[iX][iY][k]->GetN();
					//	g_vptpn_las_day[iX][iY][k]->SetPoint(NNN, daystart[iX][iY][k]+oneday/2, average);
					//	if(errors) g_vptpn_las_day[iX][iY][k]->SetPointError(NNN, oneday/2, sigmaav);
						int NNN = g_vptpn_oled_day[iX][iY][k]->GetN();
						if((daystart[iX][iY][k]+oneday/2)>1294000000) g_vptpn_oled_day[iX][iY][k]->SetPoint(NNN, daystart[iX][iY][k]+oneday/2, average);
						if(errors && (daystart[iX][iY][k]+oneday/2)>1294000000) g_vptpn_oled_day[iX][iY][k]->SetPointError(NNN, oneday/2, sigmaav);
						//start new day
						value[iX][iY][k].clear();
					//	(value[iX][iY][k]).push_back(apdpnlas);
					//	daystart[iX][iY][k] = tlas;
						if(apdpnoled>0) (value[iX][iY][k]).push_back(apdpnoled);
						daystart[iX][iY][k] = toled;
						if(NNN%1000==0 && NNN!=0) cout << "ix iy z " << iX << " "<<iY<<" "<<k <<" has " << NNN << " entries" << endl;
					}
				}}
			}

   	    }//event loop
	file->Close();

  	}//file index

//	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120202_EE_VPT_over_PN.root", "RECREATE");
//	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN_all.root", "RECREATE");
//	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121102_EE_VPT_over_PN_WG1_oled.root", "RECREATE");
	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20121130_EE_VPT_over_PN_WG1_oled_daily.root", "RECREATE");
	newFile->cd();
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){
	    for (int j = 0; j < 101; ++j){
	      if (saveAllChannels){
		//if(g_vptpn_las[i][j][k]->GetN()>0){
//			g_vptpn_las[i][j][k]->Sort();
//			g_vptpn_las[i][j][k]->Write();
		//}
		//if(g_vptpn_led[i][j][k]->GetN()>0){ 
//			g_vptpn_led[i][j][k]->Sort();
//			g_vptpn_led[i][j][k]->Write();
		//}
		//	g_vptpn_oled[i][j][k]->Sort();
		//	g_vptpn_oled[i][j][k]->Write();
/*
			g_vptpn_rat[i][j][k]->Sort();
			g_vptpn_rat[i][j][k]->Write();

			g_vptpn_orat[i][j][k]->Sort();
			g_vptpn_orat[i][j][k]->Write();

			g_EE_ampl[i][j][k]->Sort();
			g_EE_ampl[i][j][k]->Write();

			g_EE_fwhm[i][j][k]->Sort();
			g_EE_fwhm[i][j][k]->Write();
*/

			g_vptpn_oled_day[i][j][k]->Sort();
			g_vptpn_oled_day[i][j][k]->Write();
	      }
	    }
	  }
	}
	cout << "saved histograms in: " << newFile->GetName() << endl;

/*	//EB
	TFile *newFileEB = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EB_APD_over_PN_file3.root", "RECREATE");
	newFileEB->cd();
	//for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 171; ++i){
	    for (int j = 0; j < 361; ++j){
	      if (saveAllChannels){
		//if(g_vptpn_las[i][j][k]->GetN()>0){
			g_apdpn_las[i][j]->Sort();
			g_apdpn_las[i][j]->Write();
		//}
		//if(g_vptpn_led[i][j][k]->GetN()>0){ 
			g_apdpn_led[i][j]->Sort();
			g_apdpn_led[i][j]->Write();
		//}

			g_apdpn_rat[i][j]->Sort();
			g_apdpn_rat[i][j]->Write();

			g_EB_ampl[i][j]->Sort();
			g_EB_ampl[i][j]->Write();

			g_EB_fwhm[i][j]->Sort();
			g_EB_fwhm[i][j]->Write();
	      }
	    }
	  }
	//}
	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EB_APD_over_PN.root" << endl;
*/
/*
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){
	    for (int j = 0; j < 101; ++j){
		if(g_vptpn_las[i][j][k]->GetN()>0){
			firstaverage.clear();
			lastaverage.clear();
			average1=0;
			average2=0;
			for(int n = 0; n<g_vptpn_las[i][j][k]->GetN(); ++n){//start from beginning
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(x<1314807448) continue;//add points only after jump
				firstaverage.push_back(y);
				if(firstaverage.size()>=10) break;//do not average more than 10 points
			}
			for(int n = g_vptpn_las[i][j][k]->GetN()-1; n>=0; --n){//start from end
				double x,y;
				g_vptpn_las[i][j][k]->GetPoint(n,x,y);
				if(x>1315398150) continue;//add points only before data taking
				lastaverage.push_back(y);
				if(lastaverage.size()>=10) break;//do not average more than 10 points
			}
			if(firstaverage.size()==0 || lastaverage.size() == 0) {
				cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " (iz=0->EE-) has no entries between 1314807448 and 1315398150" << endl;
				continue;
			}
			if(firstaverage.size()<10 || lastaverage.size() <10) {
				cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " (iz=0->EE-) has only " << firstaverage.size() << "(1st av.) / " << lastaverage.size() << "(last av.) entries between 1314807448 and 1315398150" << endl;
			}
			for(int n = 0; n<firstaverage.size(); ++n) average1 += (firstaverage[n])/(firstaverage.size());
			for(int n = 0; n<lastaverage.size(); ++n)  average2 += (lastaverage[n] )/(lastaverage.size() );
			if(average1==average2){
				cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " (iz=0->EE-) has same (1st av.) as last av. between 1314807448 and 1315398150" << endl;
				cout << "firstaverage = " << average1 << " (out of " << firstaverage.size() << "entries [";
				for(int n = 0; n<firstaverage.size(); ++n) cout << firstaverage[n] << " ";
				cout << "])" << endl;
				cout << "lasstaverage  = " << average2 << " (out of " << lastaverage.size() << "entries [";
				for(int n = 0; n<lastaverage.size(); ++n) cout << lastaverage[n] << " ";
				cout << "])" << endl;
			}
			//int bin = VPTPN_firstaverage_EEp->FindBin(i,j);
			//int binx(-1), biny(-1), binz(-1);
			//VPTPN_firstaverage_EEp->GetBinXYZ(bin, binx, biny, binz);
			//cout << "crystal " << i << ":" << j << ":" << k << " has bins " << binx << ":" << biny << ":" << binz << " (global " << bin << ")" << endl;
			if(k==0){//EE-
				VPTPN_firstaverage_EEm->SetBinContent(i+1, j+1, average1);
				VPTPN_lastaverage_EEm ->SetBinContent(i+1, j+1, average2);
				if(average1!=0 && average2!=0) {
					if(-(1./0.22)*log(average2/average1)==0){
						cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " (iz=0->EE-) has same (1st av.) as last av. between 1314807448 and 1315398150 - x" << endl;
						cout << "firstaverage = " << average1 << " (out of " << firstaverage.size() << "entries [";
						for(int n = 0; n<firstaverage.size(); ++n) cout << firstaverage[n] << " ";
						cout << "])" << endl;
						cout << "lasstaverage  = " << average2 << " (out of " << lastaverage.size() << "entries [";
						for(int n = 0; n<lastaverage.size(); ++n) cout << lastaverage[n] << " ";
						cout << "])" << endl;
					}
					mu_VPTPN_EEm->SetBinContent(i+1, j+1, -(1./0.22)*log(average2/average1) );
					if((mu_ECAL_EEm_russ->GetBinContent(i+1, j+1))!=0) mu_VPTPNnorm_EEm_russ->SetBinContent(i+1, j+1, (-(1./0.22)*log(average2/average1))/(mu_ECAL_EEm_russ->GetBinContent(i+1, j+1)) );
					if((mu_ECAL_EEm_chin->GetBinContent(i+1, j+1))!=0) mu_VPTPNnorm_EEm_chin->SetBinContent(i+1, j+1, (-(1./0.22)*log(average2/average1))/(mu_ECAL_EEm_chin->GetBinContent(i+1, j+1)) );
				}
			}
			else{//EE+
				VPTPN_firstaverage_EEp->SetBinContent(i+1, j+1, average1);
				VPTPN_lastaverage_EEp ->SetBinContent(i+1, j+1, average2);
				if(average1!=0 && average2!=0) {
					if(-(1./0.22)*log(average2/average1)==0){
						cout << "crystal ix:iy:iz = " << i << ":" << j << ":" << k << " (iz=0->EE-) has same (1st av.) as last av. between 1314807448 and 1315398150 - y" << endl;
						cout << "firstaverage = " << average1 << " (out of " << firstaverage.size() << "entries [";
						for(int n = 0; n<firstaverage.size(); ++n) cout << firstaverage[n] << " ";
						cout << "])" << endl;
						cout << "lasstaverage  = " << average2 << " (out of " << lastaverage.size() << "entries [";
						for(int n = 0; n<lastaverage.size(); ++n) cout << lastaverage[n] << " ";
						cout << "])" << endl;
					}
					mu_VPTPN_EEp->SetBinContent(i+1, j+1, -(1./0.22)*log(average2/average1) );
					if((mu_ECAL_EEp_russ->GetBinContent(i+1, j+1))!=0) mu_VPTPNnorm_EEp_russ->SetBinContent(i+1, j+1, (-(1./0.22)*log(average2/average1))/(mu_ECAL_EEp_russ->GetBinContent(i+1, j+1)) );
					if((mu_ECAL_EEp_chin->GetBinContent(i+1, j+1))!=0) mu_VPTPNnorm_EEp_chin->SetBinContent(i+1, j+1, (-(1./0.22)*log(average2/average1))/(mu_ECAL_EEp_chin->GetBinContent(i+1, j+1)) );
				}
			}
		}//graph N>0
	    }//i
	  }//j
	}//k

	TFile *newFile2 = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120202_EE_mu_normalized.root", "RECREATE");
	newFile2->cd();
	VPTPN_firstaverage_EEp ->Write();
	VPTPN_firstaverage_EEm ->Write();
	VPTPN_lastaverage_EEp  ->Write();
	VPTPN_lastaverage_EEm  ->Write();
	mu_VPTPN_EEp           ->Write();
	mu_VPTPN_EEm           ->Write();
	mu_VPTPNnorm_EEp_chin  ->Write();
	mu_VPTPNnorm_EEp_russ  ->Write();
	mu_VPTPNnorm_EEm_chin  ->Write();
	mu_VPTPNnorm_EEm_russ  ->Write();
	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120202_EE_mu_normalized.root" << endl;
*/
/* //detId etc. connections
	sort(EEgeoVec.begin(), EEgeoVec.end(), EEetasort);
	sort(EBgeoVec.begin(), EBgeoVec.end(), EBetasort);
	for(int l = 0; l<EEgeoVec.size(); ++l){
		cout << "ix:iy:iz = " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << " (R = sqrt(ix^2+iy^2) = " << (EEgeoVec[l]).R() << ") --> eta = " << (EEgeoVec[l]).eta << endl;
		int bin = eta_EEm->FindBin((EEgeoVec[l]).ix,(EEgeoVec[l]).iy);
		if((EEgeoVec[l]).iz==-1) { 
			if(eta_EEm->GetBinContent(bin)!=0) 
				cout << "content " << eta_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).eta << endl;
			eta_EEm->SetBinContent(bin, (EEgeoVec[l]).eta);
			if(detId_EEm->GetBinContent(bin)!=0) 
				cout << "content " << detId_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).detId << endl;
			detId_EEm->SetBinContent(bin, (EEgeoVec[l]).detId);
			if(fed_EEm->GetBinContent(bin)!=0) 
				cout << "content " << fed_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).fed << endl;
			fed_EEm->SetBinContent(bin, (EEgeoVec[l]).fed);
			if(elecId_EEm->GetBinContent(bin)!=0) 
				cout << "content " << elecId_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).elecId << endl;
			elecId_EEm->SetBinContent(bin, (EEgeoVec[l]).elecId);
			if(harness_EEm->GetBinContent(bin)!=0) 
				cout << "content " << harness_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).harness << endl;
			harness_EEm->SetBinContent(bin, (EEgeoVec[l]).harness);
		}
		else if((EEgeoVec[l]).iz==+1) { 
			if(eta_EEp->GetBinContent(bin)!=0) 
				cout << "content " << eta_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).eta << endl;
			eta_EEp->SetBinContent(bin, (EEgeoVec[l]).eta);
			if(detId_EEp->GetBinContent(bin)!=0) 
				cout << "content " << detId_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).detId << endl;
			detId_EEp->SetBinContent(bin, (EEgeoVec[l]).detId);
			if(elecId_EEp->GetBinContent(bin)!=0) 
				cout << "content " << elecId_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).elecId << endl;
			elecId_EEp->SetBinContent(bin, (EEgeoVec[l]).elecId);
			if(fed_EEp->GetBinContent(bin)!=0) 
				cout << "content " << fed_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).fed << endl;
			fed_EEp->SetBinContent(bin, (EEgeoVec[l]).fed);
			if(harness_EEp->GetBinContent(bin)!=0) 
				cout << "content " << harness_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).harness << endl;
			harness_EEp->SetBinContent(bin, (EEgeoVec[l]).harness);
		}
		else cout << " should not happen: ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << endl;
		if((EEgeoVec[l]).iz==-1) { 
			bin = detId_ix_EEm->FindBin((EEgeoVec[l]).detId);
			if(detId_ix_EEm->GetBinContent(bin)!=0) 
				cout << "content " << detId_ix_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).ix << endl;
			detId_ix_EEm->SetBinContent(bin, (EEgeoVec[l]).ix);
			if(detId_iy_EEm->GetBinContent(bin)!=0) 
				cout << "content " << detId_iy_EEm->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).iy << endl;
			detId_iy_EEm->SetBinContent(bin, (EEgeoVec[l]).iy);
		}
		else if((EEgeoVec[l]).iz==+1) { 
			bin = detId_ix_EEp->FindBin((EEgeoVec[l]).detId);
			if(detId_ix_EEp->GetBinContent(bin)!=0) 
				cout << "content " << detId_ix_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).ix << endl;
			detId_ix_EEp->SetBinContent(bin, (EEgeoVec[l]).ix);
			if(detId_iy_EEp->GetBinContent(bin)!=0) 
				cout << "content " << detId_iy_EEp->GetBinContent(bin) << "for ix:iy:iz " << (EEgeoVec[l]).ix << ":" << (EEgeoVec[l]).iy << ":" << (EEgeoVec[l]).iz << ", now will fill with " << (EEgeoVec[l]).iy << endl;
			detId_iy_EEp->SetBinContent(bin, (EEgeoVec[l]).iy);
		}
	}
	for(int l = 0; l<EBgeoVec.size(); ++l){
		cout << "ieta:iphi = " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << " --> eta = " << (EBgeoVec[l]).eta << endl;
		int bin = eta_EB->FindBin((EBgeoVec[l]).ieta,(EBgeoVec[l]).iphi);
		if(eta_EB->GetBinContent(bin)!=0) 
			cout << "content " << eta_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).eta << endl;
		eta_EB->SetBinContent(bin, (EBgeoVec[l]).eta);
		if(detId_EB->GetBinContent(bin)!=0) 
			cout << "content " << detId_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).detId << endl;
		detId_EB->SetBinContent(bin, (EBgeoVec[l]).detId);
		if(fed_EB->GetBinContent(bin)!=0) 
			cout << "content " << fed_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).fed << endl;
		fed_EB->SetBinContent(bin, (EBgeoVec[l]).fed);
		if(elecId_EB->GetBinContent(bin)!=0) 
			cout << "content " << elecId_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).elecId << endl;
		elecId_EB->SetBinContent(bin, (EBgeoVec[l]).elecId);
		if(harness_EB->GetBinContent(bin)!=0) 
			cout << "content " << harness_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).harness << endl;
		harness_EB->SetBinContent(bin, (EBgeoVec[l]).harness);
		bin = detId_ieta_EB->FindBin((EBgeoVec[l]).detId);
		if(detId_ieta_EB->GetBinContent(bin)!=0) 
			cout << "content " << detId_ieta_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).ieta << endl;
		detId_ieta_EB->SetBinContent(bin, (EBgeoVec[l]).ieta);
		if(detId_iphi_EB->GetBinContent(bin)!=0) 
			cout << "content " << detId_iphi_EB->GetBinContent(bin) << "for ieta:iphi " << (EBgeoVec[l]).ieta << ":" << (EBgeoVec[l]).iphi << ", now will fill with " << (EBgeoVec[l]).iphi << endl;
		detId_iphi_EB->SetBinContent(bin, (EBgeoVec[l]).iphi);
	}
	TFile *newFile3 = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/EE_EB_IDconnection_histograms.root", "RECREATE");
	newFile3->cd();
	eta_EEm->Write();
	eta_EEp->Write();
	detId_EEp->Write();
	detId_EEm->Write();
	fed_EEp->Write();
	fed_EEm->Write();
	elecId_EEp->Write();
	elecId_EEm->Write();
	harness_EEp->Write();
	harness_EEm->Write();
	detId_ieta_EB->Write();
	detId_iphi_EB->Write();
	detId_ix_EEm->Write();
	detId_iy_EEm->Write();
	detId_ix_EEp->Write();
	detId_iy_EEp->Write();
	eta_EB->Write();
	detId_EB->Write();
	fed_EB->Write();
	elecId_EB->Write();
	harness_EB->Write();
*/
}



void load(){

//	files.push_back(TFile::Open("/shome/haweber/ECAL/DataLaser/small_ntu_data_00173378-00176049.root"));//deleted
//	files.push_back(TFile::Open("/shome/haweber/ECAL/DataLaser/small_ntu_data_00173X-00175X.root"));
//	files.push_back(TFile::Open("/shome/haweber/ECAL/DataLaser/small_ntu_data_00176X-00178X.root"));
//	files.push_back(TFile::Open("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00179X-00185X.root"));

//	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00173X-00175X.root"));
//	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00176X-00178X.root"));
//	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00179X-00185X.root"));

	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_0015X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00160X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00161X-00162X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00163X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00164X-00165X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00166X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00167X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00168X-00170X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00171X-00172X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00173X-00175X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00176X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00177X-00178X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00179X-00180X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00181X.root"));
	files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00182X-00185X.root"));


	//files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00176X.root"));
	//files.push_back(string("dcap://t3se01.psi.ch:22125/pnfs/psi.ch/cms/trivcat/store/user/haweber/ECAL/small_ntu_data_00177X-00178X.root"));



}