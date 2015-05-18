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
#include "TH2I.h"

using namespace std;

void getdetId(){

	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	TFile *outf = new TFile("/shome/haweber/ECAL/Laser/Analysis/detideb.root", "RECREATE");
//	TFile *outf = new TFile("/shome/haweber/ECAL/Laser/Analysis/detidee.root", "RECREATE");

	TFile *file = TFile::Open("/shome/haweber/ECAL/DataLaser/ana_rad_ecal_v3.root");
	TTree *tree = (TTree*)file->Get("EB");
//	TTree *tree = (TTree*)file->Get("EE");


	//EB
	TH2I *detId_eb = new TH2I("detId_eb", "detId_eb", 172, -86, 86, 361, 0, 361);
	TH3I *detId_eb3 = new TH3I("detId_eb3", "detId_eb3", 172, -86, 86, 361, 0, 361, 128925, 838851300, 838980225);


	int nev;
	TString variable, var;
	TString selection, sel;

	var = "rawId:iphi:ieta";//z-axis, y-axis, x-axis
	sel = "rawId>0";//EE+, russian
	variable  = TString::Format("%s>>%s",var.Data(),detId_eb3->GetName());
	selection = sel.Data();
	cout << " drawing: " << variable << ", sel: "<< selection << endl;
	nev = tree->Draw(variable,selection,"goff");
	cout << "Events found: " << detId_eb3->Integral() << endl;


	for(int ix = 1; ix<=172; ++ix){
		for(int iy = 1; iy<=361; ++iy){
				int counter = 0;
			for(int iz = 1; iz<=108925; ++iz){
				int bin2d = detId_eb->GetBin(ix,iy);
				int bin3d = detId_eb3->GetBin(ix,iy,iz);
				double content;
				int contentbin;
				contentbin = detId_eb3->GetBinContent(ix,iy,iz);
				content    = detId_eb3->GetZaxis()->GetBinLowEdge(iz);
				if(contentbin>0){
					if(counter>0) cout << "Already " << counter << " entries; value " << detId_eb->GetBinContent(ix,iy) << " adding " << content << " and bincontent " << contentbin << endl;
					if( detId_eb->GetBinContent(ix,iy) >0 ) cout << "value " << detId_eb->GetBinContent(ix,iy) << " adding " << content << endl;
					counter += 1;
					if(content!=int(content)) cout << content << endl;
					detId_eb->SetBinContent(ix, iy,  (int)content);
					//cout << "ix " << ix << " iy " << iy << " content " << content << endl;
				}
			}
		}
	}
	
	outf			->cd();
	if(detId_eb   ->GetEntries()>0) detId_eb	->Write();

	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/detideb.root" << endl;


/* //EE
	TH2I *detId_eem = new TH2I("detId_eem", "detId_eem", 101, 0, 101, 101, 0, 101);
	TH2I *detId_eep = new TH2I("detId_eep", "detId_eep", 101, 0, 101, 101, 0, 101);
	TH3I *detId_eem3  = new TH3I("detId_eem3",  "detId_eem3",  101, 0, 101, 101, 0, 101, 29080, 872415400, 872444480);
	TH3I *detId_eep3  = new TH3I("detId_eep3",  "detId_eep3",  101, 0, 101, 101, 0, 101, 29080, 872415400, 872444480);


	int nev;
	TString variable, var;
	TString selection, sel;

	var = "rawId:iy:ix";//z-axis, y-axis, x-axis
	sel = "iz==1";//EE+, russian
	variable  = TString::Format("%s>>%s",var.Data(),detId_eep3->GetName());
	selection = sel.Data();
	cout << " drawing: " << variable << ", sel: "<< selection << endl;
	nev = tree->Draw(variable,selection,"goff");
	cout << "Events found: " << detId_eep3->Integral() << endl;

	var = "rawId:iy:ix";//z-axis, y-axis, x-axis
	sel = "iz==-1";//EE+, russian
	variable  = TString::Format("%s>>%s",var.Data(),detId_eem3->GetName());
	selection = sel.Data();
	cout << " drawing: " << variable << ", sel: "<< selection << endl;
	nev = tree->Draw(variable,selection,"goff");
	cout << "Events found: " << detId_eem3->Integral() << endl;

	for(int ix = 1; ix<=101; ++ix){
		for(int iy = 1; iy<=101; ++iy){
				int counter[2] = {0,0};
			for(int iz = 1; iz<=29080; ++iz){
				int bin2d = detId_eem->GetBin(ix,iy);
				int bin3d = detId_eem3->GetBin(ix,iy,iz);
				double content;
				int contentbin;
				contentbin = detId_eem3->GetBinContent(ix,iy,iz);
				content    = detId_eem3->GetZaxis()->GetBinLowEdge(iz);
				if(contentbin>0){
					if(counter[0]>0) cout << "Already " << counter[0] << " entries; value " << detId_eem->GetBinContent(ix,iy) << " adding " << content << " and bincontent " << contentbin << endl;
					if( detId_eem->GetBinContent(ix,iy) >0 ) cout << "value " << detId_eem->GetBinContent(ix,iy) << " adding " << content << endl;
					counter[0] += 1;
					if(content!=int(content)) cout << content << endl;
					detId_eem->SetBinContent(ix, iy,  (int)content);
					//cout << "ix " << ix << " iy " << iy << " content " << content << endl;
				}
	
				contentbin = detId_eep3->GetBinContent(ix,iy,iz);
				content    = detId_eep3->GetZaxis()->GetBinLowEdge(iz);
				if(contentbin>0){
					if(counter[1]>0) cout << "Already " << counter[1] << " entries; value " << detId_eep->GetBinContent(ix,iy) << " adding " << content << " and bincontent " << contentbin << endl;
					if( detId_eep->GetBinContent(ix,iy) >0 ) cout << "value " << detId_eep->GetBinContent(ix,iy) << " adding " << content << endl;
					counter[1] += 1;
					if(content!=int(content)) cout << content << endl;
					detId_eep->SetBinContent(ix, iy,  (int)content);
					//cout << "ix " << ix << " iy " << iy << " content " << content << endl;
				}
			}
		}
	}
	
	outf			->cd();
	if(detId_eem   ->GetEntries()>0) detId_eem	->Write();
	if(detId_eep   ->GetEntries()>0) detId_eep	->Write();

	cout << "saved histograms in: " << "/shome/haweber/ECAL/Laser/Analysis/detidee.root" << endl;
*/

//	Long64_t nentries =  tree->GetEntries();
//	Long64_t nbytes = 0, nb = 0;
//	int nev =0;

//    	for (Long64_t jentry=0; jentry<nentries;jentry++) {
//      		nb = Sample.tree->GetEntry(jentry);   nbytes += nb;

		

}