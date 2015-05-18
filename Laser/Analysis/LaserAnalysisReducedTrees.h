//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  2 16:25:10 2012 by ROOT version 5.27/06b
// from TTree x/xtals info
// found on file: small_ntu_data_00173378-00176049.root
//////////////////////////////////////////////////////////

#ifndef LaserAnalysisReducedTrees_h
#define LaserAnalysisReducedTrees_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class LaserAnalysisReducedTrees {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           fed;
   Int_t           ix;
   Int_t           iy;
   Int_t           iz;
   Int_t           detId;
   Int_t           elecId;
   Int_t           harness;
   Float_t         eta;
   Float_t         alpha;
   Int_t           time[3];
   Float_t         lumi[3];
   Float_t         qmax[3];
   Float_t         tmax[3];
 //  Float_t         apdpnA[3];
 //  Float_t         apdpnB[3];
   Float_t         apdpnAB[3];
 //  Float_t         apdpnABS[3];
   Float_t         l_ampli[3];
   Float_t         l_rise_time[3];
   Float_t         l_fwhm[3];
   Float_t         l_prepulse[3];

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_fed;   //!
   TBranch        *b_ix;   //!
   TBranch        *b_iy;   //!
   TBranch        *b_iz;   //!
   TBranch        *b_detId;   //!
   TBranch        *b_elecId;   //!
   TBranch        *b_harness;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_alpha;   //!
   TBranch        *b_time;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_qmax;   //!
   TBranch        *b_tmax;   //!
 //  TBranch        *b_apdpnA;   //!
 //  TBranch        *b_apdpnB;   //!
   TBranch        *b_apdpnAB;   //!
 //  TBranch        *b_apdpnABS;   //!
   TBranch        *b_l_ampli;   //!
   TBranch        *b_l_rise_time;   //!
   TBranch        *b_l_fwhm;   //!
   TBranch        *b_l_prepulse;   //!

   LaserAnalysisReducedTrees(TTree *tree=0);
   virtual ~LaserAnalysisReducedTrees();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef LaserAnalysisReducedTrees_cxx
LaserAnalysisReducedTrees::LaserAnalysisReducedTrees(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("small_ntu_data_00173378-00176049.root");
      if (!f) {
         f = new TFile("small_ntu_data_00173378-00176049.root");
      }
      tree = (TTree*)gDirectory->Get("x");

   }
   Init(tree);
}

LaserAnalysisReducedTrees::~LaserAnalysisReducedTrees()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LaserAnalysisReducedTrees::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LaserAnalysisReducedTrees::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void LaserAnalysisReducedTrees::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("fed", &fed, &b_fed);
   fChain->SetBranchAddress("ix", &ix, &b_ix);
   fChain->SetBranchAddress("iy", &iy, &b_iy);
   fChain->SetBranchAddress("iz", &iz, &b_iz);
   fChain->SetBranchAddress("detId", &detId, &b_detId);
   fChain->SetBranchAddress("elecId", &elecId, &b_elecId);
   fChain->SetBranchAddress("harness", &harness, &b_harness);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("alpha", &alpha, &b_alpha);
   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("lumi", lumi, &b_lumi);
   fChain->SetBranchAddress("qmax", qmax, &b_qmax);
   fChain->SetBranchAddress("tmax", tmax, &b_tmax);
//   fChain->SetBranchAddress("apdpnA", apdpnA, &b_apdpnA);
//   fChain->SetBranchAddress("apdpnB", apdpnB, &b_apdpnB);
   fChain->SetBranchAddress("apdpnAB", apdpnAB, &b_apdpnAB);
//   fChain->SetBranchAddress("apdpnABS", apdpnABS, &b_apdpnABS);
   fChain->SetBranchAddress("l_ampli", l_ampli, &b_l_ampli);
   fChain->SetBranchAddress("l_rise_time", l_rise_time, &b_l_rise_time);
   fChain->SetBranchAddress("l_fwhm", l_fwhm, &b_l_fwhm);
   fChain->SetBranchAddress("l_prepulse", l_prepulse, &b_l_prepulse);
   Notify();
}

Bool_t LaserAnalysisReducedTrees::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LaserAnalysisReducedTrees::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LaserAnalysisReducedTrees::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LaserAnalysisReducedTrees_cxx
