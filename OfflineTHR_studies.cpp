// C++ code for offline thresholds studies
// Output: energy distributions for every ring (EB and EE) + final plots (energy vs rings) for EBM, EBP, EEM, EEP
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
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStyle.h"
#include "TMultiGraph.h"
#include "TPaveStats.h"

#include "geometryUtils.h"

string uintToString (unsigned int);
void drawGraphs(TGraph* g1,TGraph* g2,TGraph* g3, TGraph* g4, std::string Title, int xmin, int xmax, int ymin, int ymax);
void drawChannelsMaps (TH2F *h2_IC, TH2F *h2_LC, std::string run, std::string EBEE);

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
  float ADCAmp_b = 1.;
  float ADCAmp_e = 1.;
  float eCut_b = 0.;
  float eCut_e = 0.;
  int cutChStatus = 0;

  static const int  EE_Rings = 39;
  static const int  EB_Rings  = 85;

  Int_t nBins_b = 50;  
  Double_t eMin_b = 0.;
  Double_t eMax_b = 100.;
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
    std::vector<TH1F*> eCut_spectrum_BP_histos;
    std::vector<TH1F*> eCut_spectrum_BM_histos;
    std::vector<TH1F*> eCut_spectrum_EP_histos;
    std::vector<TH1F*> eCut_spectrum_EM_histos;

    eCut_spectrum_BP_histos.resize(EB_Rings);
    eCut_spectrum_BM_histos.resize(EB_Rings);
    eCut_spectrum_EP_histos.resize(EE_Rings);
    eCut_spectrum_EM_histos.resize(EE_Rings);
    
    ostringstream t;
    for (int i=0; i<EB_Rings; i++) { //EB-
     	t << "eCut_spectrum_BM_" << i-85;
	eCut_spectrum_BM_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_b,eMin_b,eMax_b);

	t.str("");
    }

    for (int i=0; i<EB_Rings; i++) { //EB+
	t << "eCut_spectrum_BP_" << i+1;
	eCut_spectrum_BP_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_b,eMin_b,eMax_b); 
	t.str("");
    }
    
    for (int i=0; i<EE_Rings; i++) { //EE-
      t << "eCut_spectrum_EM_" << i+1;
      eCut_spectrum_EM_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_e,eMin_e,eMax_e);
      t.str("");
    }

    for (int i=0; i<EE_Rings; i++) { //EE+
      t << "eCut_spectrum_EP_" << i+1;
      eCut_spectrum_EP_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",nBins_e,eMin_e,eMax_e);
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
	  rawid = atol(field); //get the rawid
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

    //TH2F for IC and LC Channel maps
    TH2F *h2_IC_EB = new TH2F("h2_IC_EB","h2_IC_EB",360,0,360,170,-85,85);
    TH2F *h2_LC_EB = new TH2F("h2_LC_EB","h2_LC_EB",360,0,360,170,-85,85);
    TH2F *h2_IC_EEM = new TH2F("h2_IC_EEM","h2_IC_EEM",100,1,101,100,1,101);
    TH2F *h2_LC_EEM = new TH2F("h2_LC_EEM","h2_LC_EEM",100,1,101,100,1,101);
    TH2F *h2_IC_EEP = new TH2F("h2_IC_EEP","h2_IC_EEP",100,1,101,100,1,101);
    TH2F *h2_LC_EEP = new TH2F("h2_LC_EEP","h2_LC_EEP",100,1,101,100,1,101);
                  
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
      float LC = pow (apdpnratio, -a); //get LC coefficient
      float IC = ICmap[ieta][iphi][iz]; //get IC coefficient from the map
       
      if(ichStatus[ieta][iphi][iz] > cutChStatus) continue;
      
      if ( mapRawId_EB[rawid]==1 ) { //EB
	eCut_b = ADCAmp_b * LC * IC * ADCToGeV_b;
	if (IC > 1.5){
	  h2_IC_EB->SetBinContent(h2_IC_EB->FindBin(iphi,ieta),IC); 
	}
	if (LC > 5.) {
	  h2_LC_EB->SetBinContent(h2_LC_EB->FindBin(iphi,ieta),LC); 
	}
	if (ieta < 0) { //EB-
	  eCut_spectrum_BM_histos[ieta+85]->Fill(eCut_b*1000.); 
	}
	else if (ieta > 0) { //EB+
	  eCut_spectrum_BP_histos[ieta-1]->Fill(eCut_b*1000.); 
	}      
      }

      else if ( mapRawId_EB[rawid]==0 ) { //EC
	eCut_e = ADCAmp_e * LC * IC * ADCToGeV_e;
	int iring = eeId->GetEndcapRing(ieta,iphi,iz); //actually ieta is ix and iphi is iy	
	if (iz < 0) { // EE-
	  eCut_spectrum_EM_histos[iring]->Fill(eCut_e*1000.);
	  if (IC > 1.5) {
	    h2_IC_EEM->SetBinContent(h2_IC_EEM->FindBin(ieta,iphi),IC); //ieta=ix, iphi=iy
	  }
	  if (LC > 5.){
	    h2_LC_EEM->SetBinContent(h2_LC_EEM->FindBin(ieta,iphi),LC); //ieta=ix, iphi=iy
	  }
	}
	else if (iz > 0) { // EE+
	  eCut_spectrum_EP_histos[iring]->Fill(eCut_e*1000.); 
	  if (IC > 1.5){
	    h2_IC_EEP->SetBinContent(h2_IC_EEP->FindBin(ieta,iphi),IC); //ieta=ix, iphi=iy
	  }
	  if (LC > 5.){
	    h2_LC_EEP->SetBinContent(h2_LC_EEP->FindBin(ieta,iphi),LC); //ieta=ix, iphi=iy
	  }
	}
      }
      
    }//apdpn map iterator
    
    std::cout << ">>>>>> End of the iteration over the crystals" << std::endl;
      
    //write ouput root file with energy distributions for every ring
    outputFile->cd();
    for(int i=0;i<EB_Rings;i++){
      eCut_spectrum_BP_histos[i]->Write();
    }

    for(int i=0;i<EB_Rings;i++){
      eCut_spectrum_BM_histos[i]->Write();
    }

    for(int i=0;i<EE_Rings;i++){
      eCut_spectrum_EP_histos[i]->Write();
    }

    for(int i=0;i<EE_Rings;i++){
      eCut_spectrum_EM_histos[i]->Write();
    }
    
    std::cout << "Output file ESpectra written" << std::endl;

    //draw TH2F for IC and LC maps
    std::string run;
    if (ii == 0)
      run = "203777";
    else if (ii == 1)
      run = "208686";

    //drawChannelsMaps (h2_IC_EB, h2_LC_EB, run, "EB");
    //drawChannelsMaps (h2_IC_EEM, h2_LC_EEM, run, "EEM");
    //drawChannelsMaps (h2_IC_EEP, h2_LC_EEP, run, "EEP");
    
    //MAKE FINAL PLOTS
    Double_t mean, sigma, maxEnergy;

    //EB
    TGraph *mean_EBP = new TGraph(); //energy mean vs ring for EB+
    TGraph *mean_EBM = new TGraph(); //energy mean vs ring for EB-
    TGraph *mean_1s_EBP = new TGraph(); //energy mean + 1sigma vs ring for EB+
    TGraph *mean_1s_EBM = new TGraph(); //energy mean + 1sigma vs ring for EB-
    TGraph *mean_2s_EBP = new TGraph(); //energy mean + 2sigma vs ring for EB+
    TGraph *mean_2s_EBM = new TGraph(); //energy mean + 2sigma vs ring for EB-
    TGraph *maxEnergy_EBP = new TGraph(); //max energy vs ring for EB+
    TGraph *maxEnergy_EBM = new TGraph(); //max energy vs ring for EB-
    
    //EE
    TGraph *mean_EEP = new TGraph(); //energy mean vs ring for EE+
    TGraph *mean_EEM = new TGraph(); //energy mean vs ring for EE-
    TGraph *mean_1s_EEP = new TGraph(); //energy mean + 1sigma vs ring for EE+
    TGraph *mean_1s_EEM = new TGraph(); //energy mean + 1sigma vs ring for EE-
    TGraph *mean_2s_EEP = new TGraph(); //energy mean + 2sigma vs ring for EE+
    TGraph *mean_2s_EEM = new TGraph(); //energy mean + 2sigma vs ring for EE-
    TGraph *maxEnergy_EEP = new TGraph(); //max energy vs ring for EE+
    TGraph *maxEnergy_EEM = new TGraph(); //max energy vs ring for EE-
        
    std::cout << "TGraphs created" << std::endl;
        
    for(int i = 0; i < EB_Rings; i++) { //EB-
      mean = eCut_spectrum_BM_histos[i]->GetMean();
      sigma = eCut_spectrum_BM_histos[i]->GetStdDev();
      maxEnergy = (eCut_spectrum_BM_histos[i]->FindLastBinAbove()) * (eMax_b - eMin_b) / nBins_b;
      
      mean_EBM->SetPoint(i,i-85,mean);
      mean_1s_EBM->SetPoint(i,i-85,mean+sigma);
      mean_2s_EBM->SetPoint(i,i-85,mean+(2*sigma));
      maxEnergy_EBM->SetPoint(i,i-85,maxEnergy);
    }

    for(int i = 0; i < EB_Rings; i++) { //EB+
      mean = eCut_spectrum_BP_histos[i]->GetMean();
      sigma = eCut_spectrum_BP_histos[i]->GetStdDev();
      maxEnergy = (eCut_spectrum_BP_histos[i]->FindLastBinAbove()) * (eMax_b - eMin_b) / nBins_b;
      
      mean_EBP->SetPoint(i,i+1,mean);
      mean_1s_EBP->SetPoint(i,i+1,mean+sigma);
      mean_2s_EBP->SetPoint(i,i+1,mean+(2*sigma));
      maxEnergy_EBP->SetPoint(i,i+1,maxEnergy);
    }
    std::cout << "TGraphs for EB filled" << std::endl;
    
    for(int i = 0; i < EE_Rings; i++) { //EE-
      mean = eCut_spectrum_EM_histos[i]->GetMean();
      sigma = eCut_spectrum_EM_histos[i]->GetStdDev();
      maxEnergy = (eCut_spectrum_EM_histos[i]->FindLastBinAbove()) * (eMax_e - eMin_e) / nBins_e;
      
      mean_EEM->SetPoint(i,i,mean);
      mean_1s_EEM->SetPoint(i,i,mean+sigma);
      mean_2s_EEM->SetPoint(i,i,mean+(2*sigma));
      maxEnergy_EEM->SetPoint(i,i,maxEnergy);
    }

    for(int i = 0; i < EE_Rings; i++) { //EE+
      mean = eCut_spectrum_EP_histos[i]->GetMean();
      sigma = eCut_spectrum_EP_histos[i]->GetStdDev();
      maxEnergy = (eCut_spectrum_EP_histos[i]->FindLastBinAbove()) * (eMax_e - eMin_e) / nBins_e;
      
      mean_EEP->SetPoint(i,i,mean);
      mean_1s_EEP->SetPoint(i,i,mean+sigma);
      mean_2s_EEP->SetPoint(i,i,mean+(2*sigma));
      maxEnergy_EEP->SetPoint(i,i,maxEnergy);
    }
    std::cout << "TGraphs for EE filled" << std::endl;
    
    outputFile->Close();

    //draw plots   
    drawGraphs(mean_EBM,mean_1s_EBM,mean_2s_EBM,maxEnergy_EBM,std::string("EBM_"+run),-87,0,20,110);   
    drawGraphs(mean_EBP,mean_1s_EBP,mean_2s_EBP,maxEnergy_EBP,std::string("EBP_"+run),-1,86,20,110);   
    drawGraphs(mean_EEM,mean_1s_EEM,mean_2s_EEM,maxEnergy_EEM,std::string("EEM_"+run),-1,40,0,1000);   
    drawGraphs(mean_EEP,mean_1s_EEP,mean_2s_EEP,maxEnergy_EEP,std::string("EEP_"+run),-1,40,0,1450);   
      
  }//loop over the list (index ii)
 
    
}//main 

// function to convert unsigned int into string
string uintToString(unsigned int val) {
  char buff[500];
  sprintf(buff, "%u", val);
  string str = buff;
  return(str);
}

void drawGraphs(TGraph* g1,TGraph* g2,TGraph* g3,TGraph* g4, std::string Title, int xmin, int xmax, int ymin, int ymax) {
    
  //fit
  //TF1 *f = new TF1("exp", "[0]+[1]*exp([2]*x)", 0, 38);
  //g1 -> Fit("exp","","", 0, 38);
  
  TF1 *f = new TF1("parabola", "pol2", 0, 38);
  g1 -> Fit("parabola","","", 0, 38);

  //TF1 *f = new TF1("f3", "pol3", 0, 38);
  //g1 -> Fit("f3","","", 0, 38);
  
  gStyle -> SetOptFit (00111);
  gStyle -> SetOptStat ("");
  gStyle -> SetStatX (.90);
  gStyle -> SetStatY (.90);
  gStyle -> SetStatW (.15);
     
  g1 -> SetTitle(Title.c_str());
  g1 -> GetXaxis() -> SetLabelSize(0.04);
  g1 -> GetYaxis() -> SetLabelSize(0.04);
  g1 -> GetXaxis() -> SetTitleSize(0.05);
  g1 -> GetYaxis() -> SetTitleSize(0.05);
  g1 -> GetYaxis() -> SetTitleOffset(1.);
  g1 -> GetYaxis() -> SetRangeUser(ymin,ymax);
  g1 -> GetXaxis() -> SetRangeUser(xmin,xmax);
 
  g1 -> GetXaxis() -> SetTitle("iRing");
  g1 -> GetYaxis() -> SetTitle("ADCtoEnergy (MeV/ADC)");
   
  g1 -> SetMarkerStyle(20);
  g1 -> SetMarkerSize(0.6);
  g1 -> SetMarkerColor(kBlack);
  g1 -> SetLineColor(kBlack);
  g1 -> SetLineWidth(1.8);

  g2 -> SetMarkerStyle(20);
  g2 -> SetMarkerSize(0.6);
  g2 -> SetMarkerColor(kBlue);
  g2 -> SetLineColor(kBlue);
  g2 -> SetLineWidth(1.8);
 
  g3 -> SetMarkerStyle(20);
  g3 -> SetMarkerSize(0.6);
  g3 -> SetMarkerColor(kRed);
  g3 -> SetLineColor(kRed);
  g3 -> SetLineWidth(1.8);

  g4 -> SetMarkerStyle(20);
  g4 -> SetMarkerSize(0.6);
  g4 -> SetMarkerColor(kGreen);
  g4 -> SetLineColor(kGreen);
  g4 -> SetLineWidth(1.8);
  

  TLegend* legend = new TLegend(.65, 0.84, 1., 0.96);
  legend -> SetFillColor(kWhite);
  legend -> SetFillStyle(1000);
  legend -> SetLineWidth(0); 
  legend -> SetLineColor(kWhite);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.04);
  legend -> AddEntry(g1,"mean energy","L");
  legend -> AddEntry(g2,"mean energy + 1 #sigma","L");
  legend -> AddEntry(g3,"mean energy + 2 #sigma","L");
  legend -> AddEntry(g4,"max energy","L");

  TCanvas* c1 = new TCanvas("c1","c1");
  c1 -> cd();

  g1 -> Draw("APL");
  //g2 -> Draw("PL");
  //g3 -> Draw("PL");
  //g4 -> Draw("PL");
  //legend -> Draw("same");
  
  c1 -> Print((Title+".png").c_str(),"png");
  c1 -> Print((Title+".pdf").c_str(),"pdf");
    
  delete c1;

}

void drawChannelsMaps (TH2F *h2_IC, TH2F *h2_LC, std::string run, std::string EBEE) {

  std::string TitleIC;
  std::string TitleLC;

  if (EBEE == "EB") {

    TitleIC = "IC_EB_"+run;
    TitleLC = "LC_EB_"+run;

    h2_IC -> SetTitle(TitleIC.c_str());
    h2_IC -> GetXaxis() -> SetTitle("iphi");
    h2_IC -> GetYaxis() -> SetTitle("ieta");
    h2_IC -> GetZaxis() -> SetTitle("IC");
    h2_IC -> GetZaxis() -> SetRangeUser(1.5,3.);
    h2_IC -> SetMarkerStyle(20);
    h2_IC -> SetMarkerColor(kRed);
    h2_IC -> SetMarkerSize(0.6);
    
    h2_LC -> SetTitle(TitleLC.c_str());
    h2_LC -> GetXaxis() -> SetTitle("iphi");
    h2_LC -> GetYaxis() -> SetTitle("ieta");
    h2_LC -> GetZaxis() -> SetTitle("LC");
    h2_LC -> GetZaxis() -> SetRangeUser(5.,12.5);
    h2_LC -> SetMarkerStyle(20);
    h2_LC -> SetMarkerColor(kRed);
    h2_LC -> SetMarkerSize(0.6);
    
  }

  if (EBEE == "EEM") {

    TitleIC = "IC_EEM_"+run;
    TitleLC = "LC_EEM_"+run;
   
    h2_IC -> SetTitle(TitleIC.c_str());
    h2_IC -> GetXaxis() -> SetTitle("ix");
    h2_IC -> GetYaxis() -> SetTitle("iy");
    h2_IC -> GetZaxis() -> SetTitle("IC");
    h2_IC -> GetZaxis() -> SetRangeUser(1.5,3.);
    h2_IC -> SetMarkerStyle(20);
    h2_IC -> SetMarkerColor(kRed);
    h2_IC -> SetMarkerSize(0.6);
    
    h2_LC -> SetTitle(TitleLC.c_str());
    h2_LC -> GetXaxis() -> SetTitle("ix");
    h2_LC -> GetYaxis() -> SetTitle("iy");
    h2_LC -> GetZaxis() -> SetTitle("LC");
    h2_LC -> GetZaxis() -> SetRangeUser(5.,12.5);
    h2_LC -> SetMarkerStyle(20);
    h2_LC -> SetMarkerColor(kRed);
    h2_LC -> SetMarkerSize(0.6);

  }

  if (EBEE == "EEP") {

    TitleIC = "IC_EEP_"+run;
    TitleLC = "LC_EEP_"+run;  

    h2_IC -> SetTitle(TitleIC.c_str());
    h2_IC -> GetXaxis() -> SetTitle("ix");
    h2_IC -> GetYaxis() -> SetTitle("iy");
    h2_IC -> GetZaxis() -> SetTitle("IC");
    h2_IC -> GetZaxis() -> SetRangeUser(1.5,3.);
    h2_IC -> SetMarkerStyle(20);
    h2_IC -> SetMarkerColor(kRed);
    h2_IC -> SetMarkerSize(0.6);

    h2_LC -> SetTitle(TitleLC.c_str());
    h2_LC -> GetXaxis() -> SetTitle("ix");
    h2_LC -> GetYaxis() -> SetTitle("iy");
    h2_LC -> GetZaxis() -> SetTitle("LC");
    h2_LC -> GetZaxis() -> SetRangeUser(5.,12.5);
    h2_LC -> SetMarkerStyle(20);
    h2_LC -> SetMarkerColor(kRed);
    h2_LC -> SetMarkerSize(0.6);
 
  }

  TCanvas* c1 = new TCanvas("c1","c1");
  c1 -> cd();
  h2_IC -> Draw("colz");
  
  c1 -> Print((TitleIC+".png").c_str(),"png");
  c1 -> Print((TitleIC+".pdf").c_str(),"pdf");
    
  delete c1;

  TCanvas* c2 = new TCanvas("c2","c2");
  c2 -> cd();
  h2_LC -> Draw("colz");
  
  c2 -> Print((TitleLC+".png").c_str(),"png");
  c2 -> Print((TitleLC+".pdf").c_str(),"pdf");

  delete c2;

}

  
