/*
 * Merge a number of .root files into one output file.
 *
 * Execute the function like this from the command line:
 *
 * echo outputfile.root input1.root input2.root... | root -b merge_root.C+
 *
 * This macro assumes that the TTree is named 'TMBTree'. Change
 * the corresponding line below if this is different for your case.
 */

#include <string>
#include <iostream>
#include <fstream>

#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

// Change this line if your TTree has a different name
const char *TreeName = "x";

void merge_root()//don't merge to many files - cpu memory overload
{
    using namespace std;

    string outfile = "/shome/haweber/ECAL/DataLaser/ntu-data_00177026-00177954.root";//ADD here the outputfilename
    const char *listname = "Filelist_173853-177954.txt";

    TList tree_list;
    std::string filename;


    Int_t total_events = 0;

    TString rootFile;
    ifstream is(listname);
    while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
	cout << "Filename: " << rootFile << endl;
	if(rootFile[0] == '#') continue;
   // while(cin >> filename) {
    //    TFile *f = new TFile(filename.c_str());
	
	TFile *f = new TFile(rootFile);
	//TFile *f = TFile::Open(rootFile);

        if(TTree *tree = (TTree *)f->Get(TreeName)) {

            cout << "Adding file: " << rootFile << endl;
            tree_list.Add(tree);

            total_events += (Int_t )tree->GetEntries();

        } else {
            cout << "File has no TTree named x" << endl;
        }
    }

    cout << "Opening output file: " << outfile << endl;
    TFile output(outfile.c_str(), "RECREATE");

    cout << "Merging trees...patience..." << endl;
    TTree::MergeTrees(&tree_list);
    output.Write();
    output.Close();

    cout << "Total Events: " << total_events << endl;
}