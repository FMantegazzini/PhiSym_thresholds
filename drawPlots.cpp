/*************************************************************
C++ program to read root file with histos and draw a final plot
compile with ---> c++ -o drawPlots `root-config --cflags --glibs` drawPlots.cpp
or with ---> c++ -o drawPlots.cpp drawPlots `root-config --cflags --glibs`
run with ---> ./drawPlots root_list.txt
*************************************************************/

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <map>
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"


void DrawPlot(std::map<int,TGraphErrors*> g1, std::string EBEE, std::string PM, std::string Gain, std::vector<int> runId);

int main(int argc, char** argv) {
  
  gStyle->SetOptStat(0000);
  
  char* inputLIST = argv[1];
  std::string inputList = std::string(inputLIST);
  std::cout << "inputList = " << inputList << std::endl;
  
  std::vector<std::string> inputFiles;
  char dumps[500];
  FILE *f_dumps;
  f_dumps = fopen(inputList.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading inputFile: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles.push_back(DUMPS);
  }

  //flat mean, mean_s1, mean_s2, maxEnergy;
  Double_t mean, sigma, maxEnergy;
  static const int kEndcEtaRings = 39;
  static const int kBarlRings = 170;

  for(unsigned int ii = 0; ii < inputFiles.size(); ii++){ //loop over root files

    std::cout << "Reading inputFile: " << inputFiles.at(ii) << std::endl;    
    TFile *inFile = TFile::Open(inputFiles.at(ii).c_str());

    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetGrid();   
    c1->SetFillColor(0);

    TGraph *mean_EB = new TGraph(); //energy mean vs ring for EB
    ostringstream t;    
    t << "mean_EB_" << ii+1;
    mean_EB->SetName(t.str().c_str());   

    TGraph *mean_EC = new TGraph(); //energy mean vs ring for EC
    t << "mean_EC_" << ii+1;
    mean_EC->SetName(t.str().c_str());   

    for(int i = 0; i < kBarlRings; i++) { //EB-
      TH1F *h = (TH1F*)inFile.Get(("eCut_spectrum_b_" + i+1).c_str());
      mean = h->GetMean();
      sigma = h->GetStdDev();
      mean_EB->SetPoint(i,i-85,mean);
    }

  } //loop over the root files
 
 
}
