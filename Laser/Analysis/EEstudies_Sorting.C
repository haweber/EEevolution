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

using namespace std;

void load();

TFile* file1;
TFile* file2;
//do this with hadd newfilename.root oldfile1.root oldfile2.root ...

//This is sorting only
void EEstudies_Sorting(){

	load();

	TGraph* g_vptpn_las[101][101][2];
	TGraph* g_vptpn_led[101][101][2];
	TGraph* g_vptpn_oled[101][101][2];
	TGraph* g_vptpn_rat[101][101][2];
	TGraph* g_EE_ampl[101][101][2];
	TGraph* g_EE_fwhm[101][101][2];

	TGraph* g_apdpn_las[171][361];
	TGraph* g_apdpn_led[171][361];
	TGraph* g_apdpn_rat[171][361];
	TGraph* g_EB_ampl[171][361];
	TGraph* g_EB_fwhm[171][361];

	char gname[101];
	bool saveAllChannels=true;
	
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
	      if (saveAllChannels){
	/*	if (k==0) sprintf(gname,"gg_vptpn_las_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_vptpn_las_ix%d_iy%d_EEp", i, j);
		g_vptpn_las[i][j][k] = new TGraph();
		g_vptpn_las[i][j][k] -> SetTitle(gname);
		g_vptpn_las[i][j][k] -> SetName(gname);
		g_vptpn_las[i][j][k] -> SetMarkerColor(kBlue);
		g_vptpn_las[i][j][k] -> SetMarkerStyle(22);
	*/ /*	
		if (k==0) sprintf(gname,"gg_vptpn_led_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_vptpn_led_ix%d_iy%d_EEp", i, j);
		g_vptpn_led[i][j][k] = new TGraph();
		g_vptpn_led[i][j][k] -> SetTitle(gname);
		g_vptpn_led[i][j][k] -> SetName(gname);
		g_vptpn_led[i][j][k] -> SetMarkerColor(kCyan);
		g_vptpn_led[i][j][k] -> SetMarkerStyle(23);
		//g_vptpn_led[i][j][k] -> SetMarkerSize(0.5);

		if (k==0) sprintf(gname,"gg_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_vptpn_oled_ix%d_iy%d_EEp", i, j);
		g_vptpn_oled[i][j][k] = new TGraph();
		g_vptpn_oled[i][j][k] -> SetTitle(gname);
		g_vptpn_oled[i][j][k] -> SetName(gname);
		g_vptpn_oled[i][j][k] -> SetMarkerColor(kOrange);
		g_vptpn_oled[i][j][k] -> SetMarkerStyle(7);

		if (k==0) sprintf(gname,"gg_vptpn_rat_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_vptpn_rat_ix%d_iy%d_EEp", i, j);
		g_vptpn_rat[i][j][k] = new TGraph();
		g_vptpn_rat[i][j][k] -> SetTitle(gname);
		g_vptpn_rat[i][j][k] -> SetName(gname);
		g_vptpn_rat[i][j][k] -> SetMarkerColor(kBlack);
		g_vptpn_rat[i][j][k] -> SetMarkerStyle(20);

		if (k==0) sprintf(gname,"gg_EE_ampl_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_EE_ampl_ix%d_iy%d_EEp", i, j);
		g_EE_ampl[i][j][k] = new TGraph();
		g_EE_ampl[i][j][k] -> SetTitle(gname);
		g_EE_ampl[i][j][k] -> SetName(gname);
		g_EE_ampl[i][j][k] -> SetMarkerColor(kGreen+3);
		g_EE_ampl[i][j][k] -> SetMarkerStyle(22);

		if (k==0) sprintf(gname,"gg_EE_fwhm_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"gg_EE_fwhm_ix%d_iy%d_EEp", i, j);
		g_EE_fwhm[i][j][k] = new TGraph();
		g_EE_fwhm[i][j][k] -> SetTitle(gname);
		g_EE_fwhm[i][j][k] -> SetName(gname);
		g_EE_fwhm[i][j][k] -> SetMarkerColor(kMagenta+2);
		g_EE_fwhm[i][j][k] -> SetMarkerStyle(23);
	*/      }
	    }
	  }
	  }
	  for (int i = 0; i < 171; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 361; ++j){//as iX goes from 1 to 100
	      if (saveAllChannels){
/*		sprintf(gname,"gg_apdpn_las_ieta%d_iphi%d", i, j);
		g_apdpn_las[i][j] = new TGraph();
		g_apdpn_las[i][j] -> SetTitle(gname);
		g_apdpn_las[i][j] -> SetName(gname);
		g_apdpn_las[i][j] -> SetMarkerColor(kBlue);
		g_apdpn_las[i][j] -> SetMarkerStyle(22);
		
		sprintf(gname,"gg_apdpn_led_ieta%d_iphi%d", i, j);
		g_apdpn_led[i][j] = new TGraph();
		g_apdpn_led[i][j] -> SetTitle(gname);
		g_apdpn_led[i][j] -> SetName(gname);
		g_apdpn_led[i][j] -> SetMarkerColor(kRed);
		g_apdpn_led[i][j] -> SetMarkerStyle(23);

		sprintf(gname,"gg_apdpn_rat_ieta%d_iphi%d", i, j);
		g_apdpn_rat[i][j] = new TGraph();
		g_apdpn_rat[i][j] -> SetTitle(gname);
		g_apdpn_rat[i][j] -> SetName(gname);
		g_apdpn_rat[i][j] -> SetMarkerColor(kBlack);
		g_apdpn_rat[i][j] -> SetMarkerStyle(20);

		sprintf(gname,"gg_EB_ampl_ix%d_iy%d", i, j);
		g_EB_ampl[i][j] = new TGraph();
		g_EB_ampl[i][j] -> SetTitle(gname);
		g_EB_ampl[i][j] -> SetName(gname);
		g_EB_ampl[i][j] -> SetMarkerColor(kGreen+3);
		g_EB_ampl[i][j] -> SetMarkerStyle(22);

		sprintf(gname,"gg_EB_fwhm_ix%d_iy%d", i, j);
		g_EB_fwhm[i][j] = new TGraph();
		g_EB_fwhm[i][j] -> SetTitle(gname);
		g_EB_fwhm[i][j] -> SetName(gname);
		g_EB_fwhm[i][j] -> SetMarkerColor(kMagenta+2);
		g_EB_fwhm[i][j] -> SetMarkerStyle(23);
*/
	      }
	    }
	  }

	

//now open Tgraphs from each single file

//	for(int fileindex = 0; fileindex<files.size(); ++fileindex){
//	TGraph* vptpn_las;
	/*TGraph* vptpn_led;
	TGraph* vptpn_oled;
	TGraph* vptpn_rat;
	TGraph* EE_ampl;
	TGraph* EE_fwhm;

	TGraph* apdpn_las;
	TGraph* apdpn_led;
	TGraph* apdpn_rat;
	TGraph* EB_ampl;*/
//	TGraph* EB_fwhm;

//	double x,y;
//	cout << "file #" << fileindex << "/"<< files.size() << endl;
	cout << "EE..." << endl;
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){//as iX goes from 1 to 100
	    for (int j = 0; j < 101; ++j){//as iX goes from 1 to 100
		if (k==0) sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_oled_ix%d_iy%d_EEp", i, j);
		g_vptpn_oled[i][j][k] = (TGraph*)file1->Get(gname);
		if (k==0) sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_las_ix%d_iy%d_EEp", i, j);
		g_vptpn_las[i][j][k] = (TGraph*)file1->Get(gname);
		if (k==0) sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_led_ix%d_iy%d_EEp", i, j);
		g_vptpn_led[i][j][k] = (TGraph*)file1->Get(gname);
		if (k==0) sprintf(gname,"g_vptpn_rat_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_vptpn_rat_ix%d_iy%d_EEp", i, j);
		g_vptpn_rat[i][j][k] = (TGraph*)file1->Get(gname);
		if (k==0) sprintf(gname,"g_EE_ampl_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_EE_ampl_ix%d_iy%d_EEp", i, j);
		g_EE_ampl[i][j][k] = (TGraph*)file1->Get(gname);
		if (k==0) sprintf(gname,"g_EE_fwhm_ix%d_iy%d_EEm", i, j);
		else sprintf(gname,"g_EE_fwhm_ix%d_iy%d_EEp", i, j);
		g_EE_fwhm[i][j][k] = (TGraph*)file1->Get(gname);


		g_vptpn_las[i][j][k]->Sort();
		g_vptpn_led[i][j][k]->Sort();
		g_vptpn_oled[i][j][k]->Sort();
		g_vptpn_rat[i][j][k]->Sort();
		g_EE_ampl[i][j][k]->Sort();
		g_EE_fwhm[i][j][k]->Sort();


	     }//i
	    }//j
	   }//k
	cout << " done. EB..." << endl;
	  for (int i = 0; i < 171; ++i){
	    for (int j = 0; j < 361; ++j){
		/*sprintf(gname,"g_apdpn_las_ieta%d_iphi%d", i, j);
		apdpn_las = (TGraph*)files2[fileindex]->Get(gname);
		sprintf(gname,"g_apdpn_led_ieta%d_iphi%d", i, j);
		apdpn_led = (TGraph*)files2[fileindex]->Get(gname);
		sprintf(gname,"g_apdpn_rat_ieta%d_iphi%d", i, j);
		apdpn_rat = (TGraph*)files2[fileindex]->Get(gname);
		sprintf(gname,"g_EB_ampl_ix%d_iy%d", i, j);
		EB_ampl = (TGraph*)files2[fileindex]->Get(gname);
		sprintf(gname,"g_EB_fwhm_ix%d_iy%d", i, j);
		EB_fwhm = (TGraph*)files2[fileindex]->Get(gname);

		for(int n = 0; n<apdpn_led->GetN(); ++n){
			apdpn_led->GetPoint(n,x,y);
			if(x!=0 && y!=0){
				g_apdpn_led[i][j]->SetPoint(g_apdpn_led[i][j]->GetN(), x,y);
			}
		}
		for(int n = 0; n<apdpn_las->GetN(); ++n){
			apdpn_las->GetPoint(n,x,y);
			if(x!=0 && y!=0){
				g_apdpn_las[i][j]->SetPoint(g_apdpn_las[i][j]->GetN(), x,y);
			}
		}
		for(int n = 0; n<apdpn_rat->GetN(); ++n){
			apdpn_rat->GetPoint(n,x,y);
			if(x!=0 && y!=0){
				g_apdpn_rat[i][j]->SetPoint(g_apdpn_rat[i][j]->GetN(), x,y);
			}
		}
		for(int n = 0; n<EB_ampl->GetN(); ++n){
			EB_ampl->GetPoint(n,x,y);
			if(x!=0 && y!=0){
				g_EB_ampl[i][j]->SetPoint(g_EB_ampl[i][j]->GetN(), x,y);
			}
		}
		for(int n = 0; n<EB_fwhm->GetN(); ++n){
			EB_fwhm->GetPoint(n,x,y);
			if(x!=0 && y!=0){
				g_EB_fwhm[i][j]->SetPoint(g_EB_fwhm[i][j]->GetN(), x,y);
			}
		}*/
	     }//i
	    }//j
	cout << " done." << endl;
//	}//fileindex


	TFile *newFile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN_new2sorted.root", "RECREATE");
	newFile->cd();
	for (int k = 0; k < 2; ++k){
	  for (int i = 0; i < 101; ++i){
	    for (int j = 0; j < 101; ++j){
	      if (saveAllChannels){
			g_vptpn_las[i][j][k]->Write();
			g_vptpn_led[i][j][k]->Write();
			g_vptpn_oled[i][j][k]->Write();
			g_vptpn_rat[i][j][k]->Write();
			g_EE_ampl[i][j][k]->Write();
			g_EE_fwhm[i][j][k]->Write();
	      }
	    }
	  }
	}
	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN_new2sorted.root" << endl;

//	TFile *newFileEB = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EB_APD_over_PN_new2sorted.root", "RECREATE");
//	newFileEB->cd();
	//for (int k = 0; k < 2; ++k){
//	  for (int i = 0; i < 171; ++i){
//	    for (int j = 0; j < 361; ++j){
//	      if (saveAllChannels){
			/*g_apdpn_las[i][j]->Write();
			g_apdpn_led[i][j]->Write();
			g_apdpn_rat[i][j]->Write();
			g_EB_ampl[i][j]->Write();
			g_EB_fwhm[i][j]->Write();*/
//	      }
//	    }
//	  }
	//}
//	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EB_APD_over_PN_new2sorted.root" << endl;

}


void load(){

	file1 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EE_VPT_over_PN_new2.root");

//	file2 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20120403_EB_APD_over_PN_new2.root");
}