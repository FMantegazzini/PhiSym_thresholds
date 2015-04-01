/*************************************************************
C++ program to read root file with histos and draw a final plot
compile with ---> c++ -o drawPlots `root-config --cflags --glibs` drawPlots.cpp
or with ---> c++ -o drawPlots.cpp drawPlots `root-config --cflags --glibs`
run with ---> ./drawPlots root_list.txt
*************************************************************/

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <map>
#include <vector>
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
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
    std::cout << "Filling the list with: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles.push_back(DUMPS);
  }

  //float mean, mean_s1, mean_s2, maxEnergy;
  Double_t mean, sigma, maxEnergy;
  static const int kEndcEtaRings = 39;
  static const int kBarlRings = 170;

  for(unsigned int ii = 0; ii < inputFiles.size(); ii++){ //loop over root_list

    std::cout << "Reading inputFile: " << inputFiles.at(ii) << std::endl;    
    TFile *inFile = TFile::Open(inputFiles.at(ii).c_str());
    
    TGraph *mean_EB = new TGraph(); //energy mean vs ring for EB    
    TGraph *mean_1s_EB = new TGraph(); //energy mean + 1sigma vs ring for EB
    TGraph *mean_2s_EB = new TGraph(); //energy mean + 2sigma vs ring for EB
    TGraph *maxEnergy_EB = new TGraph(); //max energy vs ring for EB

    TGraph *mean_EE = new TGraph(); //energy mean vs ring for EE
    TGraph *mean_1s_EE = new TGraph(); //energy mean + 1sigma vs ring for EE
    TGraph *mean_2s_EE = new TGraph(); //energy mean + 2sigma vs ring for EE
    TGraph *maxEnergy_EE = new TGraph(); //max energy vs ring for EE

    std::cout << "TGraphs created" << std::endl;
    
    ostringstream t;

    for(int i = 0; i < kBarlRings; i++) { //EB-
      if (i < 85) { //EB-
	t << "eCut_spectrum_b_" << i-85;
	TH1F *h = (TH1F*)inFile->Get(t.str().c_str());
	mean = h->GetMean();
	sigma = h->GetStdDev();
	mean_EB->SetPoint(i,i-85,mean);
	mean_1s_EB->SetPoint(i,i-85,mean+sigma);
	mean_2s_EB->SetPoint(i,i-85,mean+(2*sigma));
      }
      else if (i >= 85) { //EB+
	t << "eCut_spectrum_b_" << i-84;
	TH1F *h = (TH1F*)inFile->Get(t.str().c_str());
	mean = h->GetMean();
	sigma = h->GetStdDev();
	mean_EB->SetPoint(i,i-84,mean);
	mean_1s_EB->SetPoint(i,i-84,mean+sigma);
	mean_2s_EB->SetPoint(i,i-84,mean+(2*sigma));
      }
    }
    std::cout << "TGraphs for EB filled" << std::endl;

    for(int i = 0; i < kEndcEtaRings; i++) { //EE-
      t << "eCut_spectrum_e_" << i+1;
      TH1F *h = (TH1F*)inFile->Get(t.str().c_str());
      mean = h->GetMean();
      sigma = h->GetStdDev();
      mean_EE->SetPoint(i,i+1,mean);
      mean_1s_EE->SetPoint(i,i+1,mean+sigma);
      mean_2s_EE->SetPoint(i,i+1,mean+(2*sigma));
    }
    std::cout << "TGraphs for EE filled" << std::endl;
    
    //EB plot    
    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetGrid();   
    c1->SetFillColor(0);
        
    mean_EB->SetMarkerColor(kRed);
    mean_EB->SetMarkerStyle(20);
    mean_EB->SetMarkerSize(1);
    mean_EB->SetTitle("EB;Ring (#);Energy (MeV)");
    mean_EB->Draw("AP");
    delete mean_EB;

    mean_1s_EB->SetMarkerColor(kBlue);
    mean_1s_EB->SetMarkerStyle(20);
    mean_1s_EB->SetMarkerSize(1);
    mean_1s_EB->Draw("P");
    delete mean_1s_EB;

    mean_2s_EB->SetMarkerColor(kBlack);
    mean_2s_EB->SetMarkerStyle(20);
    mean_2s_EB->SetMarkerSize(1);
    mean_2s_EB->Draw("P");
    delete mean_2s_EB;

    /*
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->SetFillColor(0);
    //leg->SetHeader("Legend");
    leg->AddEntry(mean_EB, "Energy mean", "p");
    leg->AddEntry(mean_1s_EB, "Energy mean + 1 sigma", "p");
    leg->AddEntry(mean_2s_EB, "Energy mean + 2 sigma", "p");
    leg->Draw();
    */
    
    t << "plot_EB_" << ii+1 << ".pdf";
    c1->Draw();
    c1->Print((t.str().c_str(),"pdf");
    delete c1;

    //EE plot
    TCanvas *c2 = new TCanvas("c2","c2"); 
    c2->SetGrid();   
    c2->SetFillColor(0);
    
    mean_EE->SetMarkerColor(kRed);
    mean_EE->SetMarkerStyle(20);
    mean_EE->SetMarkerSize(1);
    mean_EE->SetTitle("EE;Ring (#);Energy (MeV)");
    mean_EE->Draw("AP");
    delete mean_EE;
    
    mean_1s_EB->SetMarkerColor(kBlue);
    mean_1s_EB->SetMarkerStyle(20);
    mean_1s_EB->SetMarkerSize(1);
    mean_1s_EB->Draw("P");
    delete mean_1s_EE;

    mean_2s_EB->SetMarkerColor(kBlack);
    mean_2s_EB->SetMarkerStyle(20);
    mean_2s_EB->SetMarkerSize(1);
    mean_2s_EB->Draw("P");
    delete mean_2s_EE;
    
    t << "plot_EE_" << ii+1 << ".pdf";
    c2->Draw();
    c2->Print((t.str().c_str(),"pdf");
    delete c2;
    
  } //loop over root_list
 
 
}//main
