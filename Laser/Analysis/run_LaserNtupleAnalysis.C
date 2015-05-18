#include "TROOT.h"
#include <iostream>
#include <TChain.h>


void run_LaserNtupleAnalysis(){
  gROOT->Reset();

  std::cout << "Loading data files" << std::endl;
  TChain data("x");
//  gROOT->ProcessLine(".x ./NtupleChain_158851_160450.C");
//  gROOT->ProcessLine(".x ./NtupleChain_160454_161897.C");
//  gROOT->ProcessLine(".x ./NtupleChain_161991_163392.C");
//  gROOT->ProcessLine(".x ./NtupleChain_163397_165154.C");
//  gROOT->ProcessLine(".x ./NtupleChain_165158_166426.C");
//  gROOT->ProcessLine(".x ./NtupleChain_166429_167541.C");
//  gROOT->ProcessLine(".x ./NtupleChain_167543_170876.C");
//  gROOT->ProcessLine(".x ./NtupleChain_170896_173377.C");
//  gROOT->ProcessLine(".x ./NtupleChain_173378_176049.C");
//  gROOT->ProcessLine(".x ./NtupleChain_176051_177493.C");
  gROOT->ProcessLine(".x ./NtupleChain_177497_178888.C");

  std::cout << "Compiling" << std::endl;
  gROOT->ProcessLine(".L ./LaserNtupleAnalysis.C++");

  std::cout << "Process data" << std::endl;
  gROOT->ProcessLine("LaserNtupleAnalysis a((TTree*) data)");
  gROOT->ProcessLine("a.Loop()");
  
}