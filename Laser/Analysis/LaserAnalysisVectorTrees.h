//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 20 09:52:16 2011 by ROOT version 5.27/06b
// from TTree LDB/Dump of Laser data from file
// found on file: /shome/donega/ECAL/DATA_corrected/DumpLaserDB.root
//////////////////////////////////////////////////////////

#ifndef LaserAnalysisVectorTrees_h
#define LaserAnalysisVectorTrees_h


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//Tree Reader for Vector NTuples (i.e. Francesca C.'s ntuples coming from Frederico)

class LaserAnalysisVectorTrees {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           time[92];
   Int_t           detId[75848];
   Int_t           eta[75848];
   Int_t           phi[75848];
   Int_t           x[75848];
   Int_t           y[75848];
   Int_t           z[75848];
   Float_t         cor[75848];

   // List of branches
   TBranch        *b_time;   //!
   TBranch        *b_detId;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_cor;   //!

   LaserAnalysisVectorTrees(TTree *tree=0);
  //define here all functions used in laser.C
   virtual ~LaserAnalysisVectorTrees();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef LaserAnalysisVectorTrees_cxx
LaserAnalysisVectorTrees::LaserAnalysisVectorTrees(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/shome/donega/ECAL/DATA_corrected/DumpLaserDB.root");
      if (!f) {
         f = new TFile("/shome/donega/ECAL/DATA_corrected/DumpLaserDB.root");
      }
      tree = (TTree*)gDirectory->Get("LDB");

   }
   Init(tree);
}

LaserAnalysisVectorTrees::~LaserAnalysisVectorTrees()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t LaserAnalysisVectorTrees::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t LaserAnalysisVectorTrees::LoadTree(Long64_t entry)
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

void LaserAnalysisVectorTrees::Init(TTree *tree)
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

   fChain->SetBranchAddress("time", time, &b_time);
   fChain->SetBranchAddress("detId", detId, &b_detId);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("cor", cor, &b_cor);
   Notify();
}

Bool_t LaserAnalysisVectorTrees::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void LaserAnalysisVectorTrees::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t LaserAnalysisVectorTrees::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef LaserAnalysisVectorTrees_cxx
