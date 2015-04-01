// C++ code for offline thresholds studies
//Output: energy distributions for every ring (EB and EE) + final plots (energy vs rings)
// to compile: c++ -o OfflineTHR_studies `root-config --cflags --glibs` OfflineTHR_studies.cpp geometryUtils.cc
// to run: ./OfflineTHR_studies APDPN_list.txt alpha_list.txt IC_list.txt ADC_list.txt ChStatus_list.txt

using namespace std;

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
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMultiGraph.h"

#include "geometryUtils.h"

string uintToString (unsigned int);

int main(int argc, char** argv)
{

  //lists of input files
  char* inputLIST_APDPN = argv[1];
  char* inputLIST_alpha = argv[2];
  char* inputLIST_IC = argv[3];
  char* inputLIST_ADC = argv[4];
  char* inputLIST_ChStat = argv[5];
  
  std::string inputList_APDPN = std::string(inputLIST_APDPN);
  std::string inputList_alpha = std::string(inputLIST_alpha);
  std::string inputList_IC = std::string(inputLIST_IC);
  std::string inputList_ADC = std::string(inputLIST_ADC);
  std::string inputList_ChStat = std::string(inputLIST_ChStat);
   
  std::cout << "inputList_APDPN = " << inputList_APDPN << std::endl;
  std::cout << "inputList_alpha = " << inputList_alpha << std::endl;
  std::cout << "inputList_IC = " << inputList_IC << std::endl;
  std::cout << "inputList_ADC = " << inputList_ADC << std::endl;
  std::cout << "inputList_ChStat = " << inputList_ChStat << std::endl << std::endl;
  
  //files vectors for every list
  std::vector<std::string> inputFiles_APDPN;
  std::vector<std::string> inputFiles_alpha;
  std::vector<std::string> inputFiles_IC;
  std::vector<std::string> inputFiles_ADC;
  std::vector<std::string> inputFiles_ChStat;
 
  char dumps[500];
  FILE *f_dumps;

  f_dumps = fopen(inputList_APDPN.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading file from APDPN_list: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_APDPN.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_alpha.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading file from alpha_list: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_alpha.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_IC.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading inputFile_IC: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_IC.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_ADC.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading file from ADC_list: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_ADC.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_ChStat.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading file from ChStatus_list: " << dumps << std::endl << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_ChStat.push_back(DUMPS);
  }  
  
  //variables
  int ieta, iphi, iz, chStatus;
  long int rawid;
  float ADCAmp_b = 8.;
  float ADCAmp_e = 12.;
  float eCut_b = 0.;
  float eCut_e = 0.;
  int cutChStatus = 3;

  static const int  kEndcEtaRings  = 39;
  static const int  kBarlRings  = 170;
  Int_t nBins_b = 50;  
  Double_t eMin_b = 0.;
  Double_t eMax_b = 500.;
  Int_t nBins_e = 75;   
  Double_t eMin_e = 0.;
  Double_t eMax_e = 1500.;
  
  TEndcapRegions* eeId = new TEndcapRegions();
  
  //make maps RawId --> ieta, iphi, iz  (ix, iy, iz)
  //make map RawId --> EB or EE (1 for EB, 0 for EE)
  std::map<long int,int> mapRawId_ieta; //ix for EE
  std::map<long int,int> mapRawId_iphi; //iy for EE
  std::map<long int,int> mapRawId_iz;
  std::map<long int,int> mapRawId_EB;
  
  char cieta[100], ciphi[100], ciz[100], crawid[100], cped12[100], crms12[100], cped6[100], crms6[100], cped1[100], crms1[100];
 
  std::ifstream infile("/afs/cern.ch/work/f/fedmante/public/PhiSym_OfflineThres/dump_EcalPedestals_hlt__since_00208206_till_00208206.dat");
  int iline = 0;
  while(infile >> cieta >> ciphi >> ciz >> cped12 >> crms12 >> cped6 >> crms6 >> cped1 >> crms1 >> crawid) {

    iline++;

    if(iline > 1 && iline <= 61201){ //EB
      ieta = atoi(cieta);
      iphi = atoi(ciphi);
      iz = atoi(ciz);
      rawid = atol(crawid);
      
      mapRawId_ieta[rawid] = ieta;
      mapRawId_iphi[rawid] = iphi;
      mapRawId_iz[rawid] = iz;
      mapRawId_EB[rawid] = 1;    
    }
    
    else if(iline > 61201){ //EE
      ieta = atoi(cieta); //ix
      iphi = atoi(ciphi); //iy
      iz = atoi(ciz);
      rawid = atol(crawid);
      
      mapRawId_ieta[rawid] = ieta; //ix
      mapRawId_iphi[rawid] = iphi; //iy
      mapRawId_iz[rawid] = iz;
      mapRawId_EB[rawid] = 0;
    }
  }
  
  std::cout << "Maps rawid --> coordinates, done" << std::endl;
  std::cout << "Map rawid --> EB, done" << std::endl << std::endl;
 
  std::cout << "LOOP OVER THE APDPN FILES LIST:" << std::endl;
  
  for (unsigned int ii = 0; ii < inputFiles_APDPN.size(); ii++) { //loop over the apdpn list (index ii)
    
    std::cout << ">>>>>> Loop number: " << ii+1 << std::endl;    
    std::cout << "Reading inputFile_APDPN: " << inputFiles_APDPN.at(ii) << std::endl;
    std::cout << "Reading inputFile_alpha: " << inputFiles_alpha.at(ii) << std::endl;
    std::cout << "Reading inputFile_ADC: " << inputFiles_ADC.at(ii) << std::endl;
    std::cout << "Reading inputFile_IC: " << inputFiles_IC.at(ii) << std::endl;
    std::cout << "Reading inputFile_ChStat: " << inputFiles_ChStat.at(ii) << std::endl;
 
    //create the output file
    std::string str = uintToString(ii+1);
    TFile *outputFile = new TFile (("ESpectra_" + str + ".root").c_str(),"RECREATE");
    std::cout << "New file created: " << "ESpectra_" << str << ".root" << std::endl; 
  
    //create spectra
    std::vector<TH1F*> eCut_spectrum_b_histos;
    std::vector<TH1F*> eCut_spectrum_e_histos;

    eCut_spectrum_b_histos.resize(kBarlRings);
    eCut_spectrum_e_histos.resize(kEndcEtaRings);


    ostringstream t;
    for (int i=0; i<kBarlRings; i++) { //EB
      if ( i < 85) { //EB-
	t << "eCut_spectrum_b_" << i-85;
	eCut_spectrum_b_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_b,eMin_b,eMax_b);
	t.str("");
      }

      if ( i >= 85) { //EB+
	t << "eCut_spectrum_b_" << i-84;
	eCut_spectrum_b_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_b,eMin_b,eMax_b); 
	t.str("");
      }
    }

    for (int i=0; i<kEndcEtaRings; i++) { //EE
      t << "eCut_spectrum_e_" << i+1;
      eCut_spectrum_e_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_e,eMin_e,eMax_e); //number of bins?
      t.str("");
    }
 
    std::cout << "Spectra created" << std::endl; 
     
    //read the Channel Status file and create ChStatus map: coord-->channel status
    std::map<int,std::map<int,std::map<int,int> > > ichStatus;
    std::ifstream infileChStatus(inputFiles_ChStat.at(ii).c_str());
    std::cout << "Reading Channel Status file" << std::endl;
    while(infileChStatus >> ieta >> iphi >> iz >> chStatus) {
      ichStatus[ieta][iphi][iz]=chStatus; //for EC ieta =ix, iphi=iy
    }
    
    //get ADCToGeV corrections (costants)
    std::ifstream infileADC(inputFiles_ADC.at(ii).c_str());
    std::string s1,s2;
    float ADCToGeV_b, ADCToGeV_e;
    std::cout << "Reading ADCToGeV file:" << std::endl; 
    while(infileADC >> s1 >> ADCToGeV_b >> s2 >> ADCToGeV_e) {
      std::cout << "ADCToGeV constants: for EB = " <<  ADCToGeV_b << ", for EE = " << ADCToGeV_e << std::endl;
    }
      
    //read the IC file and create IC map: coord-->IC
    std::map<int,std::map<int,std::map<int,float> > > ICmap;
    std::ifstream infileIC(inputFiles_IC.at(ii).c_str());  
    float IC;
    while(infileIC >> ieta >> iphi >> iz >> IC) {
      ICmap[ieta][iphi][iz]=IC; //for EC ieta =ix, iphi=iy 
    }
    
    //read the alpha file and create alpha map: coord --> alpha
    std::map<int,std::map<int,std::map<int,float> > > alphaMap; 
    std::ifstream infileAlpha(inputFiles_alpha.at(ii).c_str());
    float alpha;
    std::cout << "Reading alpha file" << std::endl;
    while(infileAlpha >> ieta >> iphi >> iz >> alpha) {
      alphaMap[ieta][iphi][iz]=alpha; //for EC ieta =ix, iphi=iy 
    }
    
    //read APDPN file and create apdpnratio map: rawid --> apdpnratio
    std::map<long int,float> apdpnMap;
    FILE *fAPDPN;
    fAPDPN = fopen (inputFiles_APDPN.at(ii).c_str(),"r");
    std::cout << "Reading APDPN file" << std::endl;
    int BUF = 5000;
    char line [ BUF ];
    int nline = 0;
    float apdpnratio;
    
    while ( fgets(line, sizeof line, fAPDPN) != NULL ) { //reading the file
    
      int n;
      int count = 0;
      char field [ BUF / 2 ], *ptr = strtok(line, "\n"); 
      int isT = 0;
           
      nline++;
     
      if (nline > 75849 ) break; 
      while ( sscanf(ptr, "%1999[^ ]%n", field, &n) == 1 ) { //reading every line splitted into fields separated by space
	count++;
	ptr += n;
	isT = 0;

	if ( *ptr != ' ' ) break;
	if ( *field == 'T' ) {
	  isT = 1;	
	  break;
	}

	if (count == 2) {
	  rawid = atol(field); //get my rawid
	}
	if (count == 4) {
	  apdpnratio = atof(field); //get the apdnpratio coefficient
	}
	++ptr;
      }//reading the line
      
      if (isT == 0) { 
	apdpnMap[rawid] = apdpnratio;
      }
      
    }//reading the file 

                 
    //loop on the apdpn map --> get rawid and apdpnratio
    std::cout << ">>>>>> Beginning iteration over the crystals" << std::endl;
    for (std::map<long int,float>::iterator it=apdpnMap.begin(); it!=apdpnMap.end(); ++it) { //apdpn map iterator

      rawid = it->first;
      apdpnratio = it->second;
      
      //get coordinates for the rawid from the RawId map
      ieta = mapRawId_ieta.find(rawid)->second;
      iphi = mapRawId_iphi.find(rawid)->second;
      iz = mapRawId_iz.find(rawid)->second;
      
      float a = alphaMap[ieta][iphi][iz]; //get alpha from the map
      float LC = pow (apdpnratio, a); //get LC coefficient
      float IC = ICmap[ieta][iphi][iz]; //get IC coefficient from the map
       
      if(ichStatus[ieta][iphi][iz] > cutChStatus) continue;
      
      if ( mapRawId_EB[rawid]==1 ) { //EB
	eCut_b = ADCAmp_b * LC * IC * ADCToGeV_b;
	
	if (ieta < 0) { //EB- histo from 0 to 84
	  eCut_spectrum_b_histos[ieta+85]->Fill(eCut_b*1000.); 
	}
	else if (ieta > 0) { //EB+ histo 85 to 169
	  eCut_spectrum_b_histos[ieta+84]->Fill(eCut_b*1000.); 
	}      
      }

      else if ( mapRawId_EB[rawid]==0 ) { //EC
	eCut_e = ADCAmp_e * LC * IC * ADCToGeV_e;
	int iring = eeId->GetEndcapRing(ieta,iphi,iz); //actually ieta is ix and iphi is iy
	eCut_spectrum_e_histos[iring]->Fill(eCut_e*1000.);  
      }
      
    }//apdpn map iterator
    
    std::cout << ">>>>>> End of the iteration over the crystals" << std::endl;
      
    //write ouput root file with energy distributions for every ring
    outputFile->cd();
    for(int i=0;i<kBarlRings;i++){
      eCut_spectrum_b_histos[i]->Write();
    }
    for(int i=0;i<kEndcEtaRings;i++){
      eCut_spectrum_e_histos[i]->Write();
    }
    
    std::cout << "Output file ESpectra written" << std::endl;

    //MAKE PLOTS
    Double_t mean, sigma, maxEnergy;

    TGraph *mean_EB = new TGraph(); //energy mean vs ring for EB
    TGraph *mean_1s_EB = new TGraph(); //energy mean + 1sigma vs ring for EB
    TGraph *mean_2s_EB = new TGraph(); //energy mean + 2sigma vs ring for EB
    TGraph *maxEnergy_EB = new TGraph(); //max energy vs ring for EB
    
    TGraph *mean_EE = new TGraph(); //energy mean vs ring for EE
    TGraph *mean_1s_EE = new TGraph(); //energy mean + 1sigma vs ring for EE
    TGraph *mean_2s_EE = new TGraph(); //energy mean + 2sigma vs ring for EE
    TGraph *maxEnergy_EE = new TGraph(); //max energy vs ring for EE

    std::cout << "TGraphs created" << std::endl;
        
    for(int i = 0; i < kBarlRings; i++) { //EB-
      if (i < 85) { //EB-
	mean = eCut_spectrum_b_histos[i]->GetMean();
	sigma = eCut_spectrum_b_histos[i]->GetStdDev();
	maxEnergy = (eCut_spectrum_b_histos[i]->FindLastBinAbove()) * (eMax_b - eMin_b) / nBins_b;
	
	mean_EB->SetPoint(i,i-85,mean);
	mean_1s_EB->SetPoint(i,i-85,mean+sigma);
	mean_2s_EB->SetPoint(i,i-85,mean+(2*sigma));
	maxEnergy_EB->SetPoint(i,i-85,maxEnergy);
      }

      else if (i >= 85) { //EB+
	mean = eCut_spectrum_b_histos[i]->GetMean();
	sigma = eCut_spectrum_b_histos[i]->GetStdDev();
	maxEnergy = (eCut_spectrum_b_histos[i]->FindLastBinAbove()) * (eMax_b - eMin_b) / nBins_b;
	
	mean_EB->SetPoint(i,i-84,mean);
	mean_1s_EB->SetPoint(i,i-84,mean+sigma);
	mean_2s_EB->SetPoint(i,i-84,mean+(2*sigma));
	maxEnergy_EB->SetPoint(i,i-84,maxEnergy);
      }
    }
    std::cout << "TGraphs for EB filled" << std::endl;
    
    for(int i = 0; i < kEndcEtaRings; i++) { //EE
      mean = eCut_spectrum_b_histos[i]->GetMean();
      sigma = eCut_spectrum_b_histos[i]->GetStdDev();
      maxEnergy = (eCut_spectrum_b_histos[i]->FindLastBinAbove()) * (eMax_e - eMin_e) / nBins_e;
      
      mean_EE->SetPoint(i,i+1,mean);
      mean_1s_EE->SetPoint(i,i+1,mean+sigma);
      mean_2s_EE->SetPoint(i,i+1,mean+(2*sigma));
      maxEnergy_EE->SetPoint(i,i+1,maxEnergy);
    }
    std::cout << "TGraphs for EE filled" << std::endl;
    
    outputFile->Close();
       
    //draw EB plot
    TCanvas *c1 = new TCanvas("c1","c1");
    c1->SetGrid();
    c1->SetFillColor(0);

    TMultiGraph *mg1 = new TMultiGraph();

    mg1->Add(mean_EB,"lp");
    mg1->Add(mean_1s_EB,"cp");
    mg1->Add(mean_2s_EB,"cp");
    mg1->Add(maxEnergy_EB,"cp");
    
    mean_EB->SetMarkerColor(kRed);
    mean_EB->SetMarkerStyle(20);
    mean_EB->SetMarkerSize(0.8);
    
    mean_1s_EB->SetMarkerColor(kBlue);
    mean_1s_EB->SetMarkerStyle(20);
    mean_1s_EB->SetMarkerSize(0.8);
    
    mean_2s_EB->SetMarkerColor(kBlack);
    mean_2s_EB->SetMarkerStyle(20);
    mean_2s_EB->SetMarkerSize(0.8);

    maxEnergy_EB->SetMarkerColor(kGreen);
    maxEnergy_EB->SetMarkerStyle(20);
    maxEnergy_EB->SetMarkerSize(0.8);
       
    mg1->SetTitle("EB;Rings (#); Energy (MeV)");
    mg1->Draw("APC");
    t << "plot_EB_" << ii+1 << ".pdf";
    c1->Draw();
    c1->Print(t.str().c_str(),"pdf");
    t.str("");

    delete mean_EB;
    delete mean_1s_EB;
    delete mean_2s_EB;    
    delete maxEnergy_EB;
    delete c1;
    delete mg1;
         
    //draw EE plot
    TCanvas *c2 = new TCanvas("c2","c2");
    c2->SetGrid();
    c2->SetFillColor(0);
    
    TMultiGraph *mg2 = new TMultiGraph();

    mg2->Add(mean_EE,"lp");
    mg2->Add(mean_1s_EE,"cp");
    mg2->Add(mean_2s_EE,"cp");
    mg2->Add(maxEnergy_EE,"cp");
    
    mean_EE->SetMarkerColor(kRed);
    mean_EE->SetMarkerStyle(20);
    mean_EE->SetMarkerSize(0.8);
    
    mean_1s_EE->SetMarkerColor(kBlue);
    mean_1s_EE->SetMarkerStyle(20);
    mean_1s_EE->SetMarkerSize(0.8);
   
    mean_2s_EE->SetMarkerColor(kBlack);
    mean_2s_EE->SetMarkerStyle(20);
    mean_2s_EE->SetMarkerSize(0.8);
    
    maxEnergy_EE->SetMarkerColor(kGreen);
    maxEnergy_EE->SetMarkerStyle(20);
    maxEnergy_EE->SetMarkerSize(0.8);
    
    mg2->SetTitle("EE;Rings (#); Energy (MeV)");
    mg2->Draw("APC");
    t << "plot_EE_" << ii+1 << ".pdf";
    c2->Draw();
    c2->Print(t.str().c_str(),"pdf");
    
    delete mean_EE;
    delete mean_1s_EE;
    delete mean_2s_EE;
    delete maxEnergy_EE;
    delete c2;
    delete mg2;
    
  }//loop over the list (index ii)
 
    
}//main 

// function to convert unsigned int into string
string uintToString(unsigned int val)
{
  char buff[500];
  sprintf(buff, "%u", val);
  string str = buff;
  return(str);
}
