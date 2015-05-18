#include <iostream>
#include <cmath>
#include <vector>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <THistPainter.h>
#include <TGraphPainter.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TString.h>
#include <TMath.h>
#include <TVector.h>
#include "Utilities.hh"
#include <time.h>
#include <stdio.h>
#include <string>
#include <sstream>

using namespace std;

/*
//initial routine which has been modified below to timeconvertedintochar
void timeconverter(time_t utctime){

//    time_t     now;
    struct tm  ts;
//    char       buf[80];
	char buf[28];

    // Get current time
//    time(&now);
 //   TString time1 = TString(now);
    //char *s=ctime(&now);
    //s[strlen(s)-1]=0;

    // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
//    ts = *localtime(&now);
    ts = *localtime(&utctime);
    strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
    //strftime(buf, sizeof(buf), "%y/%m/%d" , &ts);
    printf("%s\n", buf);
//	printf("%s\n", buf.size());
    //printf("%s\n", s);


    return;
}
*/

//usage example, after program has been compiled under root
//char X[10];
//timeconverterintochar(1300000000, X);
// X has now form yyyy/mm/dd as char

void timeconvertedintochar(time_t utctime, char *buf){

//    time_t     now;//time is now function imput
    struct tm  ts;
//    char       buf[80];//buf is now function input

      
    // Get current time//not anymore
//    time(&now);
 //   TString time1 = TString(now);
    //char *s=ctime(&now);
    //s[strlen(s)-1]=0;

    // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
//    ts = *localtime(&now);//if now is wanted time
    ts = *localtime(&utctime);
    //strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", &ts);//if other display is wanted
    strftime(buf, 80, "%Y-%m-%d" , &ts);



    return;
}

//output is a TString of design "yyyy/mm/dd", input utctime
TString timeconverter(time_t utctime){

	char buf[32];

	timeconvertedintochar(utctime, buf);

	stringstream timestream;

	string timestring;
	timestream << buf;
	timestream >> timestring;
	TString output = timestring;

	//std::cout << "String " << output << " has length " << output.Length() << std::endl;

	return output;

}
