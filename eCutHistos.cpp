// C++ code to plot energy thresholds distributions for every ECAL ring (EB and EC)
// to compile: c++ -o eCutHistos `root-config --cflags --glibs` eCutHistos.cpp
// to run: ./eCutHistos.cpp APDPN_list.txt alpha_list.txt IC_list.txt ADC_list.txt ChStatus_list.txt

using namespace std;

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <map>
#include "TFile.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTree.h"
#include "TStyle.h"
#include <vector>

std::vector<std::string> splitString(std::string& inputName, char split_char);

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
  std::cout << "inputList_ChStat = " << inputList_ChStat << std::endl;

  std::vector<std::string> inputFiles_APDPN;
  std::vector<std::string> inputFiles_alpha;
  std::vector<std::string> inputFiles_IC;
  std::vector<std::string> inputFiles_ADC;
  std::vector<std::string> inputFiles_ChStat;

  char dumps[500];

  //fill files vectors for every list

  FILE *f_dumps;
  f_dumps = fopen(inputList_APDPN.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading inputFile_LC: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_APDPN.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_alpha.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading inputFile_alpha: " << dumps << std::endl;
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
    std::cout << "Reading inputFile_ADC: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_ADC.push_back(DUMPS);
  }  

  f_dumps = fopen(inputList_ChStat.c_str(),"r");
  while(fscanf(f_dumps,"%s \n", dumps) !=EOF ){
    std::cout << "Reading inputFile_ChStat: " << dumps << std::endl;
    std::string DUMPS = std::string(dumps);
    if(DUMPS.find("#") != std::string::npos) continue;
    inputFiles_ChStat.push_back(DUMPS);
  }  
  //fclose

  //make maps RawId --> ieta, iphi, iz  (ix, iy, iz)

  std::map<int,int> mapRawId_ieta;
  std::map<int,int> mapRawId_iphi;
  std::map<int,int> mapRawId_iz;

  FILE *fMap;                                                                 
  fMap = fopen("dump_EcalPedestals_hlt__since_00208206_till_00208206.dat","r");

  int ieta,iphi,iz;
  float f1,f2,f3,f4,f5,f6;
  long int RawId;
  fscanf(fMap, "%s \n", dumps); //skip the first line that begins with #
  //std::string DUMPS = std::string(dumps);
  //if(DUMPS.find("#") != std::string::npos) continue;
  int scan = 0;
  while(scan!=EOF){
    scan=fscanf(fMap,"%d %d %d %f %f %f %f %f %f %d \n", &ieta,&iphi,&iz,&f1,&f2,&f3;6f4,&f5,&f6,&RawId);
    mapRawId_ieta[RawId] = ieta;
    mapRawId_iphi[RawId] = iphi;
    mapRawId_iz[RawId] = iz; 
  }
  for (auto i : mapRawId_ieta)
    cout << "ieta = " <<  i.first << " ; RawId = " << i.second << endl;  //print one map to test
  
  //variables

  float ADCAmp_b = 8.;
  float ADCAmp_e = 12.;
  float eCut_b = 0.;
  float eCut_e = 0.;
  static const int  kEndcEtaRings  = 39;
  static const int  kBarlRings  = 85;
  int cutChStatus = 3;

  std::cout << "LOOP OVER THE APDPN FILES LIST";

  for(unsigned int ii = 0; ii < inputFiles_APDPN(size); ii++) { //loop over the list

  std::cout << "Loop number: " << ii << std::endl;    
  std::cout << "Reading inputFile_APDPN: " << inputFiles_APDPN.at(ii) << std::endl;
  std::cout << "Reading inputFile_alpha: " << inputFiles_alpha.at(ii) << std::endl;
  std::cout << "Reading inputFile_ADC: " << inputFiles_ADC.at(ii) << std::endl;
  std::cout << "Reading inputFile_IC: " << inputFiles_IC.at(ii) << std::endl;
  std::cout << "Reading inputFile_ChStat: " << inputFiles_ChStat.at(ii) << std::endl;

  //create the output file
  ostringstream t;
  t << "Espectra_" << ii+1;
  TFile *outputFile = new TFile(t + ".root").c_str(),"RECREATE");
  
  //create spectra
  std::vector<TH1F*> eCut_spectrum_b_histos;
  std::vector<TH1F*> eCut_spectrum_e_histos;

  eCut_spectrum_b_histos.resize(kBarlRings);
  eCut_spectrum_e_histos.resize(kEndcEtaRings);

  ostringstream t;
  for(Int_t i=0;i<kBarlRings;i++)
    {
      t << "eCut_spectrum_b_" << i+1;
      eCut_spectrum_b_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",50,0.,500.);
      t.str("");
    }
  for(Int_t i=0;i<kEndcEtaRings;i++)
    {
      t << "eCut_spectrum_e_" << i+1;
      eCut_spectrum_e_histos[i]=new TH1F(t.str().c_str(),";E_{CUT} [MeV]",75,0.,1500.);
      t.str("");
    }


  //read the Channel Status file and create ChStatus map: coord-->channel status
  std::map<int,std::map<int,std::map<int,int> > > ichStatus;
  std::ifstream infileChStatus(inputFiles_ChStat.at(ii).c_str());
  int ieta, iphi, iz, chStatus;
  while(infileChStatus >> ieta >> iphi >> iz >> chStatus)
    {
      ichStatus[ieta][iphi][iz]=chStatus; //for EC ieta =ix, iphi=iy 
    }


  //get ADCToGeV corrections (costants)
  std::string nameFile(inputFiles_ADC.at(ii).c_str());
  FILE *fADC;
  fADC = fopen (nameFile,"r");
  std::string s1,s2;
  float ADCToGeV_b, ADCToGeV_e;
  fscanf(fADC,"%s %f %s %f \n",&s1,&ADCToGeV_b,&s2,&ADCToGeV_e);
  cout << "ADCToGeV constant for EB = " <<  ADCToGeV_b << " ; ADCToGeV constant for EC = " << ADCToGeV_e << endl;
  
  //read the IC file and create IC map: coord-->IC
  std::map<int,std::map<int,std::map<int,int> > > ICmap;
  std::string nameFile(inputFiles_IC.at(ii).c_str()); 
  std::ifstream infileIC(nameFile);
  int ieta, iphi, iz;
  float IC;
  while(infileIC >> ieta >> iphi >> iz >> IC)
    {
      ICmap[ieta][iphi][iz]=IC; //for EC ieta =ix, iphi=iy 
    }

  //read the alpha file and create alpha map: coord -->alpha
  std::map<int,std::map<int,std::map<int,int> > > AlphaMap;
  std::string nameFile(inputFiles_alpha.at(ii).c_str()); 
  std::ifstream infileAlpha(nameFile);
  float alpha;
  while(infileAlpha >> ieta >> iphi >> iz >> alpha)
    {
      AlphaMap[ieta][iphi][iz]=alpha; //for EC ieta =ix, iphi=iy 
    }


  //read APDPN file and get apdpnratio
  std::string nameFile(inputFiles_APDPN.at(ii).c_str());
  FILE *fAPDPN;
  fLC = fopen (nameFile,"r");
  int BUFSIZ = 1000;
  char line [ BUFSIZ ];
  long int rawid;
  float apdpnratio;

  while ( fgets(line, sizeof line, file) != NULL ) { //reading the lines
   
    int n;
    int count = 0;
    char field [ BUFSIZ / 2 ], *ptr = std::strtok(line, "\n");
    printf("line  = \"%s\"\n", line);
    while ( sscanf(ptr, "%1999[^ ]%n", field, &n) == 1 ) { //splitting the line into fields (space separated)
      count++;
      ptr += n;
      if ( *ptr != ' ' ) break;
      if ( *field == 'T' ) break;
      if (count == 2) rawid = atoi(field); //get my rawid
      if (count == 4) apdpnratio = atof(field); //get the apdnpratio
      ++ptr;
    }
    putchar('\n');

    //get coordinates for my rawid from the map
    auto iter = mapRawId_ieta.find(rawid);
    if (iter != mapRawId_ieta.end())
      ieta = iter->second;
    auto iter = mapRawId_iphi.find(rawid);
    if (iter != mapRawId_iphi.end())
      iphi = iter->second;
    auto iter = mapRawId_iz.find(rawid);
    if (iter != mapRawId_iz.end())
      iz = iter->second;
  
    float a = AlphaMap[ieta][iphi][iz]; //get alpha from the map
    float LC = (apdpnratio)^a; //sintassi?? //get LC coefficient
    float IC = ICmap[ieta][iphi][iz]; //get IC coefficient from the map

    if(ichStatus[ieta][iphi][iz] < cutChStatus) continue;

    if (/*condizione per essere nel barrel*/) {
      eCut_b = ADCAmp_b * LC * IC * ADCToGeV_b;
      eCut_spectrum_b_histos[abs(ieta)-1]->Fill(eCut_b*1000.); //controllare come sono binnati gli histo 
    }
    else if (/*condizione per essere nell'endcap*/) {
      eCut_e = ADCAmp_e * LC * IC * ADCToGeV_e;
      eCut_spectrum_e_histos[abs(ieta)-1]->Fill(eCut_e*1000.); 
    }
    
  }

  //write ouput root file
  for(int i=0;i<kBarlRings;i++){
    eCut_spectrum_b_histos[i]->Write();
  }
  for(int i=0;i<kEndcEtaRings;i++){
    eCut_spectrum_e_histos[i]->Write();
  }
  fOutput.Close();
   
}//loop over the list


}//main 

std::vector<std::string> splitString(std::string& inputName, char split_char)
{
  std::vector<std::string> tokens;
  std::istringstream split(inputName);
  for(std::string each; getline(split, each, split_char); tokens.push_back(each));
  return tokens;
}

