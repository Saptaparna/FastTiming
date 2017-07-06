#include <TF1.h>
#include <TH1D.h>
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
  double recHitE;
  double recHitX;
  double recHitY;
  double recHitZ;
  double recHitTime;
  double recHitSoverN;
  double recHitRmsTime;
  double recHitDeltaTime;
  double recHitSmearedTime; 
} RecHitInfo;

typedef struct
{
  double genParticleEta;
  double genParticlePhi;
} genParticleInfo;

bool sortRecHitsInDescendingE(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitE > recHit2.recHitE);
}

bool sortRecHitsInAscendingTime(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitTime < recHit2.recHitTime);
}

bool sortRecHitsInAscendingSmearedTime(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitSmearedTime < recHit2.recHitSmearedTime);
}

bool sortRecHitsInAscendingZ(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (fabs(recHit1.recHitZ) < fabs(recHit2.recHitZ));
}

bool sortRecHitsInAscendingE(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitE < recHit2.recHitE);
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


int TruncatedMean(std::string infile, std::string outfile, double fraction, double rhoCut)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("hgctiming/HGCTiming");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t           run;
  Int_t           event;
  Int_t           lumi;
  Int_t           reachedEE;
  Float_t         vertex_x;
  Float_t         vertex_y;
  Float_t         vertex_z;
  vector<double>   *recHit_energy;
  vector<double>   *recHit_x;
  vector<double>   *recHit_y;
  vector<double>   *recHit_z;
  vector<double>   *recHit_time;
  vector<double>   *recHit_smearedTime;
  vector<double>   *recHit_soverN;
  vector<double>   *recHit_rmsTime;
  vector<double>   *recHit_deltaTime;
  vector<double>   *cluster_x;
  vector<double>   *cluster_y;
  vector<double>   *cluster_z;
  vector<double>   *cluster_time;
  vector<double>   *cluster_eta;
  vector<double>   *simCluster_eta;
  vector<double>   *simCluster_phi;
  vector<double>   *genParticle_eta;
  vector<double>   *genParticle_phi;
  
  recHit_energy = 0;
  recHit_x = 0;
  recHit_y = 0;
  recHit_z = 0;
  recHit_time = 0;
  recHit_smearedTime = 0;
  recHit_soverN = 0;
  recHit_rmsTime = 0;
  recHit_deltaTime = 0;
  cluster_x = 0;
  cluster_y = 0;
  cluster_z = 0;
  cluster_time = 0;
  cluster_eta = 0;
  simCluster_eta = 0;
  simCluster_phi = 0;
  genParticle_eta = 0;
  genParticle_phi = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("event", &(event));
  tree->SetBranchAddress("reachedEE", &(reachedEE));
  tree->SetBranchAddress("vertex_x", &(vertex_x));
  tree->SetBranchAddress("vertex_y", &(vertex_y));
  tree->SetBranchAddress("vertex_z", &(vertex_z));
  tree->SetBranchAddress("recHit_energy", &(recHit_energy));
  tree->SetBranchAddress("recHit_x", &(recHit_x));
  tree->SetBranchAddress("recHit_y", &(recHit_y));
  tree->SetBranchAddress("recHit_z", &(recHit_z));
  tree->SetBranchAddress("recHit_time", &(recHit_time));
  tree->SetBranchAddress("recHit_smearedTime", &(recHit_smearedTime));
  tree->SetBranchAddress("recHit_soverN", &(recHit_soverN));
  tree->SetBranchAddress("recHit_rmsTime", &(recHit_rmsTime));
  tree->SetBranchAddress("recHit_deltaTime", &(recHit_deltaTime));
  tree->SetBranchAddress("cluster_x", &(cluster_x));
  tree->SetBranchAddress("cluster_y", &(cluster_y));
  tree->SetBranchAddress("cluster_z", &(cluster_z));
  tree->SetBranchAddress("cluster_time", &(cluster_time));
  tree->SetBranchAddress("cluster_eta", &(cluster_eta));
  tree->SetBranchAddress("simCluster_eta", &(simCluster_eta));
  tree->SetBranchAddress("simCluster_phi", &(simCluster_phi));
  tree->SetBranchAddress("genParticle_eta", &(genParticle_eta));
  tree->SetBranchAddress("genParticle_phi", &(genParticle_phi));

  TH1D *h_TimeAverage = new TH1D("h_TimeAverage", "h_TimeAverage; Average time [ns]; Events", 5000,  -5.0, 5.0); h_TimeAverage->Sumw2();
  TH1D *h_nHits_afterCut = new TH1D("h_nHits_afterCut", "h_nHits_afterCut; Number of hits; Events", 500, 0.0, 500.0); h_nHits_afterCut->Sumw2(); 
  TH1D *h_nHits_FullSize = new TH1D("h_nHits_FullSize", "h_nHits_FullSize; Number of hits before removal", 500, 0.0, 500.0); h_nHits_FullSize->Sumw2();
  TH1D *h_recHit_Time =  new TH1D("h_recHit_Time", "h_recHit_Time; recHit time [ns]; Events", 100, 0.0, 0.1);h_recHit_Time->Sumw2();
  TH1D *h_soverN = new TH1D("h_soverN", "h_soverN; SoverN; Entries", 40, -10.0, 1000.0); h_soverN->Sumw2(); 
  TH1D *h_rmsTime = new TH1D("h_rmsTime", "h_rmsTime; rmsTime; Entries", 100, 0.0, 0.1); h_rmsTime->Sumw2();
  TH1D *h_deltaTime = new TH1D("h_deltaTime", "h_deltaTime; rmsTime; Entries", 10000, -5.0, 5.0); h_deltaTime->Sumw2();
  TH1D *h_rho = new TH1D("h_rho", "h_rho; rho [cm]; Entries", 100, 0.0, 100.0); h_rho->Sumw2();
  TH1D *h_recHitX = new TH1D("h_recHitX", "h_recHitX; recHit X [cm]; Entries", 150, 0.0, 150.0);h_recHitX->Sumw2();
  TH1D *h_recHitY = new TH1D("h_recHitY", "h_recHitY; recHit Y [cm]; Entries", 100, 50.0, -50.0);h_recHitY->Sumw2();
  TH1D *h_recHitZ = new TH1D("h_recHitZ", "h_recHitZ; recHit Z [cm]; Entries", 100, 300.0, 400.0);h_recHitZ->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int ievent=0; ievent<nEvents; ++ievent)
  {
    tree->GetEvent(ievent);

    if(reachedEE==0) continue;

    std::vector<RecHitInfo> recHits;
    double sumOfEnergy = 0.0;

    for (unsigned int k=0; k<recHit_energy->size(); k++)
    {
      RecHitInfo recHit;
      recHit.recHitE = recHit_energy->at(k);
      recHit.recHitX = recHit_x->at(k);
      recHit.recHitY = recHit_y->at(k);
      recHit.recHitZ = recHit_z->at(k);
      recHit.recHitTime = recHit_time->at(k);
      recHit.recHitSoverN = recHit_soverN->at(k);
      recHit.recHitRmsTime = recHit_rmsTime->at(k);
      recHit.recHitSmearedTime = recHit_smearedTime->at(k);
      recHit.recHitDeltaTime = recHit_deltaTime->at(k);
      if(recHit_time->at(k) > -2.0) sumOfEnergy += recHit_energy->at(k);
      if(recHit_time->at(k) > -2.0) recHits.push_back(recHit);
    }

    std::sort(recHits.begin(), recHits.end(), sortRecHitsInAscendingSmearedTime);

    //std::cout << "recHits.size() = " << recHits.size() << std::endl; 
   
    h_nHits_FullSize->Fill(recHits.size());
  
    int removeElements = fraction*recHits.size(); //try fixed removal first

    for(int i=0; i<removeElements; i++)
    {
      recHits.pop_back();
    }

    //std::cout << "reduced recHits.size() = " << recHits.size() << std::endl;

    double sumTimeAverage = 0.0;
    double sumResolution = 0.0;
    int nhits = 0;
    for (unsigned int j=0; j<recHits.size(); j++)
    {
      if(recHits.size() >= 3) h_recHitX->Fill(recHits.at(j).recHitX);
      if(recHits.size() >= 3) h_recHitY->Fill(recHits.at(j).recHitY);
      if(recHits.size() >= 3) h_recHitZ->Fill(recHits.at(j).recHitZ);
      double rho = sqrt(pow((recHits.at(j).recHitX-h_recHitX->GetMean()), 2) + pow(recHits.at(j).recHitY, 2));
      if(recHits.size() >= 3 and rho<rhoCut) h_rho->Fill(rho);
      if(recHits.size() >= 3 and rho<rhoCut) sumTimeAverage += recHits.at(j).recHitSmearedTime*(1.0/(pow(recHits.at(j).recHitDeltaTime, 2)));
      if(recHits.size() >= 3 and rho<rhoCut) sumResolution += 1.0/(pow(recHits.at(j).recHitDeltaTime, 2));
      if(recHits.size() >= 3 and rho<rhoCut) nhits++;
      if(recHits.size() >= 3 and rho<rhoCut) h_rmsTime->Fill(recHits.at(j).recHitRmsTime);
      if(recHits.size() >= 3 and rho<rhoCut) h_soverN->Fill(recHits.at(j).recHitSoverN);
      if(recHits.size() >= 3 and rho<rhoCut) h_deltaTime->Fill(recHits.at(j).recHitSmearedTime-recHits.at(j).recHitTime);
      if(recHits.size() >= 3 and rho<rhoCut) h_recHit_Time->Fill(recHits.at(j).recHitSmearedTime);
    }
    if(nhits>0) h_TimeAverage->Fill(sumTimeAverage/sumResolution);
    if(nhits>0) h_nHits_afterCut->Fill(recHits.size());
 
  }//end of event loop

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_rho->Write();
  h_recHitX->Write();
  h_recHitY->Write();
  h_recHitZ->Write();
  h_TimeAverage->Write();
  h_recHit_Time->Write();
  h_nHits_afterCut->Write();
  h_nHits_FullSize->Write();
  h_rmsTime->Write();
  h_soverN->Write();
  h_deltaTime->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0; 
}
