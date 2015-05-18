#include <TH1.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TFile.h>
#include <vector>
#include <iostream>

void RebinHistos(){

	//define Histos here
//	TFile *f1 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/detId_startvalues_qmax_apdpnAB_fwhm.root");
//	TFile *f2 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/detId_startvalues_apdpnA_apdpnB.root");
//	TFile *f3 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/detId_startvalues_pnA_pnB.root");
//	TFile *f4 = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/detId_startvalues_qmaxDivpnA_qmaxDivpnB.root");
	TFile *f = TFile::Open("/shome/haweber/ECAL/Laser/Analysis/UsedFiles/DetIDCrystalconnections.root");

	TH1I *old1 = (TH1I*)f->Get("h_detId_vs_ieta_ix");
	TH1I *old2 = (TH1I*)f->Get("h_detId_vs_iphi_iy");
	TH1I *old3 = (TH1I*)f->Get("h_detId_iz");

	//const int nhistos = 9;
	//TH1F *old[nhistos];

/*	TH1F *old1  = (TH1F*)f1->Get("h_detId_vs_qmax_wl0");
	TH1F *old2  = (TH1F*)f1->Get("h_detId_vs_qmax_wl1");
	TH1F *old3  = (TH1F*)f1->Get("h_detId_vs_qmax_wl2");
	TH1F *old4  = (TH1F*)f1->Get("h_detId_vs_apdpnAB_wl0");
	TH1F *old5  = (TH1F*)f1->Get("h_detId_vs_apdpnAB_wl1");
	TH1F *old6  = (TH1F*)f1->Get("h_detId_vs_apdpnAB_wl2");
	old[6]  = (TH1F*)f2->Get("h_detId_vs_apdpnA_wl0");
	old[7]  = (TH1F*)f2->Get("h_detId_vs_apdpnA_wl1");
	old[8]  = (TH1F*)f2->Get("h_detId_vs_apdpnA_wl2");
	old[9]  = (TH1F*)f2->Get("h_detId_vs_apdpnB_wl0");
	old[10] = (TH1F*)f2->Get("h_detId_vs_apdpnB_wl1");
	old[11] = (TH1F*)f2->Get("h_detId_vs_apdpnB_wl2");
	old[12] = (TH1F*)f3->Get("h_detId_vs_pnAqmax_wl0");
	old[13] = (TH1F*)f3->Get("h_detId_vs_pnAqmax_wl1");
	old[14] = (TH1F*)f3->Get("h_detId_vs_pnAqmax_wl2");
	old[15] = (TH1F*)f3->Get("h_detId_vs_pnBqmax_wl0");
	old[16] = (TH1F*)f3->Get("h_detId_vs_pnBqmax_wl1");
	old[17] = (TH1F*)f3->Get("h_detId_vs_pnBqmax_wl2");
	old[18] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnA_qmax_wl0");
	old[19] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnA_qmax_wl1");
	old[20] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnA_qmax_wl2");
	old[21] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnB_qmax_wl0");
	old[22] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnB_qmax_wl1");
	old[23] = (TH1F*)f4->Get("h_detId_vs_qmaxDivpnB_qmax_wl2");
	TH1F *old7 = (TH1F*)f1->Get("h_detId_vs_fwhm_wl0");
	TH1F *old8 = (TH1F*)f1->Get("h_detId_vs_fwhm_wl1");
	TH1F *old9 = (TH1F*)f1->Get("h_detId_vs_fwhm_wl2");*/



	//define binvectors
//	vector<double> vec[nhistos];

	vector<double> vec1;
	vector<double> vec2;
	vector<double> vec3;
	/*vector<double> vec4;
	vector<double> vec5;
	vector<double> vec6;
	vector<double> vec7;
	vector<double> vec8;
	vector<double> vec9;
	vector<double> vec0;*/ 
	//define entries
//	vector<double> vecc[nhistos];

	vector<int> vecc1;
	vector<int> vecc2;
	vector<int> vecc3;
/*	vector<double> vecc4;
	vector<double> vecc5;
	vector<double> vecc6;
	vector<double> vecc7;
	vector<double> vecc8;
	vector<double> vecc9;
	vector<double> vecc0; */

/*	int size[nhistos];

	for(int k = 0; k<nhistos; ++k){
		double lowedge;
		for(int i = 1; i<=old[k]->GetNbinsX(); ++i){
			if(old[k]->GetBinContent(i)!=0){
				lowedge = old[k]->GetBinLowEdge(i);
				if((vec[k]).size()==0) (vec[k]).push_back(lowedge-1.);
				(vec[k]).push_back(lowedge);
				(vecc[k]).push_back(old[k]->GetBinContent(i));
			}
		}
		(vec[k]).push_back((vec[k])[(vec[k]).size()-1]+1.);
		size[k] = (vec[k]).size();
	}

	TH1F *hnew[nhistos];

	for(int k = 0; k<nhistos; ++k){
		double arr[(size[k])];
		for(int i = 0; i<(size[k]); ++i) arr[i] = (vec[k])[i];
		hnew[k] = new TH1F(old[k]->GetName()+TString("_rebinned"), old[k]->GetName()+TString("_rebinned"), size[k]-1, arr);
		for(int i = 2; i<= hnew[k]->GetNbinsX(); ++i) hnew[k]->SetBinContent(i, (vecc[k])[i-2]);
	}

	TFile *newfile = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/detId_startvalues_rebinned1.root", "RECREATE");
	newfile->cd();
	for(int k = 0; k<nhistos; ++k){
		hnew[k]->Write();
	}
	newfile->Close();
*/

	//find binning -> add a first and a last bin
	//here all histos have the apriori same bin
		//BREAKS SOMEWHERE HERE
	for(int i = 1; i<=old1->GetNbinsX(); ++i){
		double lowedge;
		if(old1->GetBinContent(i)!=0){
			lowedge = old1->GetBinLowEdge(i);
			if(vec1.size()==0) vec1.push_back(lowedge-1.);
			vec1.push_back(lowedge);
			vecc1.push_back(old1->GetBinContent(i));
		}
		if(old2->GetBinContent(i)!=0){
			lowedge = old2->GetBinLowEdge(i);
			if(vec2.size()==0) vec2.push_back(lowedge-1.);
			vec2.push_back(lowedge);
			vecc2.push_back(old2->GetBinContent(i));
		}
		if(old3->GetBinContent(i)!=0){
			lowedge = old3->GetBinLowEdge(i);
			if(vec3.size()==0) vec3.push_back(lowedge-1.);
			vec3.push_back(lowedge);
			vecc3.push_back(old3->GetBinContent(i));
		}
/*		if(old4->GetBinContent(i)!=0){
			lowedge = old4->GetBinLowEdge(i);
			if(vec4.size()==0) vec4.push_back(lowedge-1.);
			vec4.push_back(lowedge);
			vecc4.push_back(old4->GetBinContent(i));
		}
		if(old5->GetBinContent(i)!=0){
			lowedge = old5->GetBinLowEdge(i);
			if(vec5.size()==0) vec5.push_back(lowedge-1.);
			vec5.push_back(lowedge);
			vecc5.push_back(old5->GetBinContent(i));
		}
		if(old6->GetBinContent(i)!=0){
			lowedge = old6->GetBinLowEdge(i);
			if(vec6.size()==0) vec6.push_back(lowedge-1.);
			vec6.push_back(lowedge);
			vecc6.push_back(old6->GetBinContent(i));
		}
		if(old7->GetBinContent(i)!=0){
			lowedge = old7->GetBinLowEdge(i);
			if(vec7.size()==0) vec7.push_back(lowedge-1.);
			vec7.push_back(lowedge);
			vecc7.push_back(old7->GetBinContent(i));
		}
		if(old8->GetBinContent(i)!=0){
			lowedge = old8->GetBinLowEdge(i);
			if(vec8.size()==0) vec8.push_back(lowedge-1.);
			vec8.push_back(lowedge);
			vecc8.push_back(old8->GetBinContent(i));
		}
		if(old9->GetBinContent(i)!=0){
			lowedge = old9->GetBinLowEdge(i);
			if(vec9.size()==0) vec9.push_back(lowedge-1.);
			vec9.push_back(lowedge);
			vecc9.push_back(old9->GetBinContent(i));
		}*/
		/*if(old0->GetBinContent(i)!=0){
			lowedge = old0->GetBinLowEdge(i);
			if(vec0.size()==0) vec0.push_back(lowedge-1.);
			vec0.push_back(lowedge);
			vecc0.push_back(old0->GetBinContent(i));
		}*/
	}
	vec1.push_back(vec1[vec1.size()-1]+1.);
	vec2.push_back(vec2[vec2.size()-1]+1.);
	vec3.push_back(vec3[vec3.size()-1]+1.);
	/*vec4.push_back(vec4[vec4.size()-1]+1.);
	vec5.push_back(vec5[vec5.size()-1]+1.);
	vec6.push_back(vec6[vec6.size()-1]+1.);
	vec7.push_back(vec7[vec7.size()-1]+1.);
	vec8.push_back(vec8[vec8.size()-1]+1.);
	vec9.push_back(vec9[vec9.size()-1]+1.);*/
	//vec0.push_back(vec0[vec0.size()-1]+1.);

	//define arrays and size
	int size1 = vec1.size();
	int size2 = vec2.size();
	int size3 = vec3.size();
	/*int size4 = vec4.size();
	int size5 = vec5.size();
	int size6 = vec6.size();
	int size7 = vec7.size();
	int size8 = vec8.size();
	int size9 = vec9.size();
	int size0 = vec0.size();*/
	double arr1[size1];
	double arr2[size2];
	double arr3[size3];
	/*double arr4[size4];
	double arr5[size5];
	double arr6[size6];
	double arr7[size7];
	double arr8[size8];
	double arr9[size9];*/
	//double arr0[size0];
	for(int i = 0; i<size1; ++i) arr1[i] = vec1[i];
	for(int i = 0; i<size2; ++i) arr2[i] = vec2[i];
	for(int i = 0; i<size3; ++i) arr3[i] = vec3[i];
	/*for(int i = 0; i<size4; ++i) arr4[i] = vec4[i];
	for(int i = 0; i<size5; ++i) arr5[i] = vec5[i];
	for(int i = 0; i<size6; ++i) arr6[i] = vec6[i];
	for(int i = 0; i<size7; ++i) arr7[i] = vec7[i];
	for(int i = 0; i<size8; ++i) arr8[i] = vec8[i];
	for(int i = 0; i<size9; ++i) arr9[i] = vec9[i];*/
	//for(int i = 0; i<size0; ++i) arr0[i] = vec0[i];


	//define new histos here
	TH1I *new1 = new TH1I(old1->GetName()+TString("_rebinned"), old1->GetName()+TString("_rebinned"), size1-1, arr1);
	TH1I *new2 = new TH1I(old2->GetName()+TString("_rebinned"), old2->GetName()+TString("_rebinned"), size2-1, arr2);
	TH1I *new3 = new TH1I(old3->GetName()+TString("_rebinned"), old3->GetName()+TString("_rebinned"), size3-1, arr3);
/*	TH1F *new1 = new TH1F(old1->GetName()+TString("_rebinned"), old1->GetName()+TString("_rebinned"), size1-1, arr1);
	TH1F *new2 = new TH1F(old2->GetName()+TString("_rebinned"), old2->GetName()+TString("_rebinned"), size2-1, arr2);
	TH1F *new3 = new TH1F(old3->GetName()+TString("_rebinned"), old3->GetName()+TString("_rebinned"), size3-1, arr3);
	TH1F *new4 = new TH1F(old4->GetName()+TString("_rebinned"), old4->GetName()+TString("_rebinned"), size4-1, arr4);
	TH1F *new5 = new TH1F(old5->GetName()+TString("_rebinned"), old5->GetName()+TString("_rebinned"), size5-1, arr5);
	TH1F *new6 = new TH1F(old6->GetName()+TString("_rebinned"), old6->GetName()+TString("_rebinned"), size6-1, arr6);
	TH1F *new7 = new TH1F(old7->GetName()+TString("_rebinned"), old7->GetName()+TString("_rebinned"), size7-1, arr7);
	TH1F *new8 = new TH1F(old8->GetName()+TString("_rebinned"), old8->GetName()+TString("_rebinned"), size8-1, arr8);
	TH1F *new9 = new TH1F(old9->GetName()+TString("_rebinned"), old9->GetName()+TString("_rebinned"), size9-1, arr9);*/
	//TH1F *new0 = new TH1F(old0->GetName()+TString("_rebinned"), old0->GetName()+TString("_rebinned"), size0-1, arr0);

	//fill the new histos
	//note: first bin should be empty
	for(int i = 2; i<= new1->GetNbinsX(); ++i) new1->SetBinContent(i, vecc1[i-2]);
	for(int i = 2; i<= new2->GetNbinsX(); ++i) new2->SetBinContent(i, vecc2[i-2]);
	for(int i = 2; i<= new3->GetNbinsX(); ++i) new3->SetBinContent(i, vecc3[i-2]);
/*	for(int i = 2; i<= new4->GetNbinsX(); ++i) new4->SetBinContent(i, vecc4[i-2]);
	for(int i = 2; i<= new5->GetNbinsX(); ++i) new5->SetBinContent(i, vecc5[i-2]);
	for(int i = 2; i<= new6->GetNbinsX(); ++i) new6->SetBinContent(i, vecc6[i-2]);
	for(int i = 2; i<= new7->GetNbinsX(); ++i) new7->SetBinContent(i, vecc7[i-2]);
	for(int i = 2; i<= new8->GetNbinsX(); ++i) new8->SetBinContent(i, vecc8[i-2]);
	for(int i = 2; i<= new9->GetNbinsX(); ++i) new9->SetBinContent(i, vecc9[i-2]);*/
	//for(int i = 2; i<= new0->GetNbinsX(); ++i) new0->SetBinContent(i, vecc0[i-2]);

//	std::cout << "sizes " << size1 << "/" << size2 << "/" << size3 << "/" << size4 << "/" << size5 << "/" << size6 << "/" << size7 << "/" << size8 << "/" << size9 << "/" << size0 << std::endl;
//	std::cout << "veccs " << vecc1.size() << "/" << vecc2.size() << "/" << vecc3.size() << "/" << vecc4.size() << "/" << vecc5.size() << "/" << vecc6.size() << "/" << vecc7.size() << "/" << vecc8.size() << "/" << vecc9.size() << "/" << vecc0.size() << std::endl;

	std::cout << "sizes " << size1 << "/" << size2 << "/" << size3 << "  veccs " << vecc1.size() << "/" << vecc2.size() << "/" << vecc3.size() << std::endl;
	//save histos in new file
	TFile *g = new TFile("/shome/haweber/ECAL/Laser/Analysis/FileOutputs/DetIDCrystalconnections_rebinned.root", "RECREATE");
	g->cd();
	new1->Write();
	new2->Write();
	new3->Write();
	/*new4->Write();
	new5->Write();
	new6->Write();
	new7->Write();
	new8->Write();
	new9->Write();*/
	//new0->Write();
	g->Close();  
}