#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include <time.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TString.h"

using namespace std;

inline TString MakeOutputDir(TString dir){
	if(!dir.EndsWith("/")) dir += "/";
	// Create directory if needed
	//  >> NOTE: This function needs to be called before the booking functions!
	char cmd[100];
	sprintf(cmd,"mkdir -p %s", dir.Data());
	system(cmd);
	return dir;
}
//DOES NOT WORK!!!!!!!!!!!!!!!!!!!!!!!!

const Int_t    nfile                 = 11;  // max file = 11
const Int_t    feds                  = 72; //max 72 
const Int_t    lasers                = 3;//max 3
const Int_t    variables             = 3;//max 3

Bool_t         veto                  = true;//true: veto inf entries, false: do not veto

TString        newfilename           = "StepCalibration.root";
TString        newfiledir            = "/shome/haweber/ECAL/Laser/Analysis/FileOutputs/20111026/";


void FileAdding(){

	TString laser[3];
	TString variable[3];
	TString name[72][3][3];


	laser[0] = "bluelaser";
	laser[1] = "blueLED";
	laser[2] = "bluelaserVsLED";
	variable[0] = "apdpnAB";
	variable[1] = "ampl";
	variable[2] = "fwhm";

	for(int l = 0; l<3; ++l){
		for(int v = 0; v<3; ++v){

			name[0][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed601_harness9";
			name[1][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed601_harness10";
			name[2][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed601_harness11";
			name[3][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed601_harness12";
			name[4][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed602_harness5";
			name[5][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed602_harness6";
			name[6][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed602_harness7";
			name[7][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed602_harness8";
			name[8][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed603_harness1";
			name[9][v][l]  = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed603_harness2";
			name[10][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed603_harness3";
			name[11][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed603_harness4";
			name[12][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed604_harness1";
			name[13][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed604_harness2";
			name[14][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed604_harness3";
			name[15][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed604_harness4";
			name[16][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed605_harness5";
			name[17][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed605_harness6";
			name[18][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed605_harness7";
			name[19][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed605_harness8";
			name[20][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed606_harness9";
			name[21][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed606_harness10";
			name[22][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed606_harness11";
			name[23][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed606_harness12";
			name[24][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed607_harness14";
			name[25][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed607_harness15";
			name[26][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed607_harness16";
			name[27][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed607_harness17";
			name[28][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed608_harness18";
			name[29][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed608_harness19";
			name[30][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed601_harness13";
			name[31][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed646_harness13";
			name[32][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed609_harness14";
			name[33][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed609_harness15";
			name[34][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed609_harness16";
			name[35][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed609_harness17";
			name[36][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed646_harness9";
			name[37][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed646_harness10";
			name[38][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed646_harness11";
			name[39][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed646_harness12";
			name[40][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed647_harness5";
			name[41][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed647_harness6";
			name[42][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed647_harness7";
			name[43][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed647_harness8";
			name[44][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed648_harness1";
			name[45][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed648_harness2";
			name[46][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed648_harness3";
			name[47][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed648_harness4";
			name[48][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed649_harness1";
			name[49][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed649_harness2";
			name[50][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed649_harness3";
			name[51][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed649_harness4";
			name[52][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed650_harness5";
			name[53][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed650_harness6";
			name[54][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed650_harness7";
			name[55][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed650_harness8";
			name[56][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed651_harness9";
			name[57][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed651_harness10";
			name[58][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed651_harness11";
			name[59][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed651_harness12";
			name[60][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed652_harness14";
			name[61][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed652_harness15";
			name[62][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed652_harness16";
			name[63][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed652_harness17";
			name[64][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed653_harness18";
			name[65][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed653_harness19";
			name[66][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed606_harness13";
			name[67][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed651_harness13";
			name[68][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed654_harness14";
			name[69][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed654_harness15";
			name[70][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed654_harness16";
			name[71][v][l] = "h_"+variable[v]+"_vs_time_"+laser[l]+"_fed654_harness17";
		}
	}

	TString fname[11];
	fname[0]  ="../FileOutputs/20111026/StepCalibration_runs158851_160450.root";
	fname[1]  ="../FileOutputs/20111026/StepCalibration_runs160454_161897.root";
	fname[2]  ="../FileOutputs/20111026/StepCalibration_runs161991_163392.root";
	fname[3]  ="../FileOutputs/20111026/StepCalibration_runs163397_165154.root";
	fname[4]  ="../FileOutputs/20111026/StepCalibration_runs165158_166246.root";
	fname[5]  ="../FileOutputs/20111026/StepCalibration_runs166429_167541.root";
	fname[6]  ="../FileOutputs/20111026/StepCalibration_runs167543_170876.root";
	fname[7]  ="../FileOutputs/20111026/StepCalibration_runs170896_173377.root";
	fname[8]  ="../FileOutputs/20111026/StepCalibration_runs173378_176049.root";
	fname[9]  ="../FileOutputs/20111026/StepCalibration_runs176051_177493.root";
	fname[10] ="../FileOutputs/20111026/StepCalibration_runs177497_178888.root";

	TGraph *mygraph_new[feds][variables][lasers];


	for(int ifile = 0; ifile < nfile; ++ifile){

		TFile *f = TFile::Open(fname[ifile]);
		TGraph *mygraph[feds][variables][lasers];

		std::cout << "Read out file " << fname[ifile] << std::endl;

		for (Int_t i= 0 ; i<feds; ++i){
			for(Int_t v= 0 ; v<variables; ++v){
				for(Int_t l= 0 ; l<lasers; ++l){
					if(l!=0 && v>0) continue;
cout << "came here i " << i << " v " << v << " l " << l << endl;
					mygraph[i][v][l]= (TGraph*)f->Get(name[i][v][l]);
					vector<double> xv,yv;
					int N = mygraph[i][v][l]->GetN();
					int nn = mygraph_new[i][v][l]->GetN();
cout << "came here2" << endl;
					for(int p = 0; p<N; ++p){
						Double_t x,y;
						mygraph[i][v][l]->GetPoint(p,x,y);
						if(veto && fabs(y)>999999.) continue;
						xv.push_back(x);
						yv.push_back(y);
					}
cout << "came here2a" << endl;
cout << "came here2b" << endl;

					for(int p = 0; p<xv.size(); ++p){
						mygraph_new[i][v][l]->SetPoint(nn+p,xv[p],yv[p]);
					}
cout << "came here2c" << endl;

					xv.clear(); yv.clear();
					

cout << "came here3" << endl;
				}
			}
		}
		std::cout << "Readed file " << fname[ifile] << std::endl;
		//f->Close();
	}
	MakeOutputDir(newfiledir);
	TFile *newfile = new TFile(newfiledir+newfilename, "RECREATE");
	newfile->cd();
	for (Int_t i= 0 ; i<feds; ++i){
		for(Int_t v= 0 ; v<variables; ++v){
			for(Int_t l= 0 ; l<lasers; ++l){
				if(mygraph_new[i][v][l]->GetN()>0) mygraph_new[i][v][l]->Write();
			}
		}
	}
	newfile->Close();
}
