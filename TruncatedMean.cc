#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TMath.h>

using std::string;

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  double recHitTime;
  double recHitTimeOld;
} RecHitInfo;


bool sortRecHitsInAscendingTime(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitTime < recHit2.recHitTime);
}

std::string itoa(int i) 
{
  char res[50];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

bool sameVal(double a, double b)
{
   return fabs(a - b) < 3.000e-03;
}

bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}


int EnergyWeighting(std::string infile, std::string outfile)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("ana/HGCTiming");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t            run;
  Int_t            event;
  Int_t            lumi;
  vector<double>*  Energydiff; 
  vector<double>*  Timediff;
  vector<double>*  TimeOld;
  vector<double>*  Time;

  Energydiff = 0;
  Timediff = 0;
  TimeOld = 0;
  Time = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("event", &(event));
  tree->SetBranchAddress("Energydiff", &(Energydiff));
  tree->SetBranchAddress("Timediff", &(Timediff));
  tree->SetBranchAddress("TimeOld", &(TimeOld));
  tree->SetBranchAddress("Time", &(Time));

  TH1F *h_TimeAverage = new TH1F("h_TimeAverage", "h_TimeAverage; Average time [ns]; Events", 200000,  -10.0, 10.0); h_TimeAverage->Sumw2();
  TH1F *h_TimeAverageOld = new TH1F("h_TimeAverageOld", "h_TimeAverageOld; Average time [ns]; Events", 200000,  -10.0, 10.0); h_TimeAverageOld->Sumw2();
  TH1F *h_smear = new TH1F("h_smear", "h_smear; smear; Entries", 200000,  -10.0, 10.0); h_smear->Sumw2();
  TH1F *h_nHits = new TH1F("h_nHits", "h_nHits; Number of hits; Events", 1000, -0.5, 999.5); h_nHits->Sumw2();
  TH1F *h_recHitTime = new TH1F("h_recHitTime", "h_recHitTime; recHitTime [ns]; Entries", 200000, -10.0, 10.0); h_recHitTime->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int ievent=0; ievent<nEvents; ++ievent)
  {
    tree->GetEvent(ievent);
    //filling the recHit info in a vector of structs
    std::vector<RecHitInfo> recHits;
    int hitNumber = 0.0;
 
    for (unsigned int k=0; k<Energydiff->size(); k++)
    {
      RecHitInfo recHit;
      recHit.recHitTimeOld = TimeOld->at(k)-1;
      recHit.recHitTime = Time->at(k);
      if(TimeOld->at(k) != -1 and Time->at(k) != -1) recHits.push_back(recHit);
    }

    std::sort(recHits.begin(), recHits.end(), sortRecHitsInAscendingTime);

    std::cout << "recHits.size() = " << recHits.size() << std::endl; 
   
    int removeElements = 0.80*recHits.size();

    for(int i=0; i<removeElements; i++)
    {
      recHits.pop_back();
    }
     
    std::cout << "reduced recHits.size() = " << recHits.size() << std::endl;
    
    double sumTimeAverage = 0.0;
    double totalSmearing = 0.0;
    double sumTimeAverageOld = 0.0;
    int nhits = 0;
    for (unsigned int j=0; j<recHits.size(); j++)
    {
      double smearing = recHits.at(j).recHitTime/recHits.at(j).recHitTimeOld - 1;
      if(recHits.at(j).recHitTimeOld > 0.0) h_recHitTime->Fill(recHits.at(j).recHitTime);
      if(recHits.at(j).recHitTimeOld > 0.0) sumTimeAverage += recHits.at(j).recHitTime*(1.0/(smearing*smearing));
      if(recHits.at(j).recHitTimeOld > 0.0) totalSmearing += 1.0/(smearing*smearing);
      if(recHits.at(j).recHitTimeOld > 0.0) nhits++;
      if(recHits.at(j).recHitTimeOld > 0.0) sumTimeAverageOld += recHits.at(j).recHitTimeOld + gRandom->Gaus(0, 0.05);
      if(recHits.at(j).recHitTimeOld > 0.0) h_smear->Fill(1.0/(smearing)*(smearing)); 
    }
    if(totalSmearing>0) h_smear->Fill(totalSmearing);
    //if(totalSmearing>0) h_nHits->Fill(nhits);
    h_nHits->Fill(recHits.size());
    //9GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-1.63557434456684803e-02);
    //8 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-2.07130106784667589e-02);
    //7 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-2.70781027686943787e-02);
    //6 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-3.69075531145384161e-02);
    //5 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-5.32697669522512740e-02);
    //4 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-8.35900898636587897e-02);
    //3 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-1.49993049482509022e-01);
    //2 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-3.46779644449904723e-01);
    //1 GeV
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing);
    //if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-1.63580068002320189e+00);
    if(totalSmearing>0) h_TimeAverage->Fill((sumTimeAverage/totalSmearing)-1.63580068002320189e+00);
    if(nhits>0) h_TimeAverageOld->Fill(sumTimeAverageOld/nhits);
  }//end of event loop

  
  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_TimeAverage->Write();
  h_TimeAverageOld->Write();
  h_smear->Write();
  h_nHits->Write();
  h_recHitTime->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0; 
}
