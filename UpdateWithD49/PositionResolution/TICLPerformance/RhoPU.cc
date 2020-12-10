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
  double recHitEta;
  double recHitLayer;
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

bool sortByTimeOrder(double i, double j)
{
  return (i < j);
}

// go from HGC layer to z in cmf
float layerToZ( int layer, float eta ) 
{
    map<int, float> lToZ;
    lToZ[0]  = 320.75;
    lToZ[1]  = 321.50;
    lToZ[2]  = 322.73;
    lToZ[3]  = 323.48;
    lToZ[4]  = 324.71;
    lToZ[5]  = 325.46;
    lToZ[6]  = 326.69;
    lToZ[7]  = 327.44;
    lToZ[8]  = 328.67;
    lToZ[9]  = 329.42; //first set
    lToZ[10] = 330.73;
    lToZ[11] = 331.60;
    lToZ[12] = 332.91;
    lToZ[13] = 333.78;
    lToZ[14] = 335.09;
    lToZ[15] = 335.96;
    lToZ[16] = 337.27;
    lToZ[17] = 338.14;
    lToZ[18] = 339.45;
    lToZ[19] = 340.32; //second set
    lToZ[20] = 341.77;
    lToZ[21] = 342.84;
    lToZ[22] = 344.29;
    lToZ[23] = 345.36;
    lToZ[24] = 346.81;
    lToZ[25] = 347.88;
    lToZ[26] = 349.33;
    lToZ[27] = 350.40; //third set
    lToZ[28] = 356.33;
    lToZ[29] = 361.01;
    lToZ[30] = 365.69;
    lToZ[31] = 370.37;
    lToZ[32] = 375.05;
    lToZ[33] = 379.73;
    lToZ[34] = 384.41;
    lToZ[35] = 389.09;
    lToZ[36] = 393.77;
    lToZ[37] = 398.45;
    lToZ[38] = 403.13;
    lToZ[39] = 407.81; //fourth set
    
    float z = lToZ[ layer ];
    if( eta < 0 ) z *= -1.;
    return z;
}

void layerIntersection(double to[3], double from[3], double fromB[3])
{
  to[0]=(from[0]-fromB[0]) / (from[2] - fromB[2]) * (to[2] - from[2]) + from[0];
  to[1]=(from[1]-fromB[1]) / (from[2] - fromB[2]) * (to[2] - from[2]) + from[1];
}


int RhoPU(std::string infile, double rhoCut)
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
  vector<double>   *recHit_eta;
  vector<double>   *recHit_layer;
  vector<double>   *cluster_x;
  vector<double>   *cluster_y;
  vector<double>   *cluster_z;
  vector<double>   *cluster_time;
  vector<double>   *cluster_eta;
  vector<double>   *simCluster_eta;
  vector<double>   *simCluster_phi;
  vector<double>   *genParticle_eta;
  vector<double>   *genParticle_phi;
 
  std::cout<<"declarations "<<std::endl;
 
  recHit_energy = 0;
  recHit_x = 0;
  recHit_y = 0;
  recHit_z = 0;
  recHit_time = 0;
  recHit_eta = 0;
  recHit_layer = 0;
  cluster_x = 0;
  cluster_y = 0;
  cluster_z = 0;
  cluster_time = 0;
  cluster_eta = 0;
  simCluster_eta = 0;
  simCluster_phi = 0;
  genParticle_eta = 0;
  genParticle_phi = 0;

  std::cout<<"null pointer "<<std::endl;

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
  tree->SetBranchAddress("recHit_eta", &(recHit_eta));
  tree->SetBranchAddress("recHit_layer", &(recHit_layer));
  tree->SetBranchAddress("cluster_x", &(cluster_x));
  tree->SetBranchAddress("cluster_y", &(cluster_y));
  tree->SetBranchAddress("cluster_z", &(cluster_z));
  tree->SetBranchAddress("cluster_time", &(cluster_time));
  tree->SetBranchAddress("cluster_eta", &(cluster_eta));
  tree->SetBranchAddress("simCluster_eta", &(simCluster_eta));
  tree->SetBranchAddress("simCluster_phi", &(simCluster_phi));
  tree->SetBranchAddress("genParticle_eta", &(genParticle_eta));
  tree->SetBranchAddress("genParticle_phi", &(genParticle_phi));

  std::cout<<"branch addresses "<<std::endl;

  TH1D *h_TimeAverage = new TH1D("h_TimeAverage", "h_TimeAverage; Average time [ns]; Events", 5000,  -5.0, 5.0); h_TimeAverage->Sumw2();
  TH1D *h_nHits_afterCut = new TH1D("h_nHits_afterCut", "h_nHits_afterCut; Number of hits; Events", 80000, 0.0, 80000.0); h_nHits_afterCut->Sumw2(); 
  TH1D *h_nHits_FullSize = new TH1D("h_nHits_FullSize", "h_nHits_FullSize; Number of hits before removal", 200000, 0.0, 200000.0); h_nHits_FullSize->Sumw2();
  TH1D *h_recHit_Time =  new TH1D("h_recHit_Time", "h_recHit_Time; recHit time [ns]; Events", 2000, -10.0, 10.0);h_recHit_Time->Sumw2();
  TH1D *h_soverN = new TH1D("h_soverN", "h_soverN; SoverN; Entries", 40, -10.0, 1000.0); h_soverN->Sumw2(); 
  TH1D *h_rmsTime = new TH1D("h_rmsTime", "h_rmsTime; rmsTime; Entries", 100, 0.0, 0.1); h_rmsTime->Sumw2();
  TH1D *h_deltaTime = new TH1D("h_deltaTime", "h_deltaTime; rmsTime; Entries", 10000, -5.0, 5.0); h_deltaTime->Sumw2();
  TH1D *h_rho = new TH1D("h_rho", "h_rho; rho [cm]; Entries", 100, 0.0, 100.0); h_rho->Sumw2();
  TH1D *h_recHitX = new TH1D("h_recHitX", "h_recHitX; recHit X [cm]; Entries", 150, 0.0, 150.0);h_recHitX->Sumw2();
  TH1D *h_recHitY = new TH1D("h_recHitY", "h_recHitY; recHit Y [cm]; Entries", 100, 50.0, -50.0);h_recHitY->Sumw2();
  TH1D *h_recHitZ = new TH1D("h_recHitZ", "h_recHitZ; recHit Z [cm]; Entries", 100, 300.0, 400.0);h_recHitZ->Sumw2();

  TFile *outputFile;
  TTree *outputTree;
  int nHitsAfterRho;
  std::vector<double> v_recHitTime;

  std::cout<<"output filename "<<std::endl;

  std::string outputtreename=(infile+"_outputTree.root").c_str();
  outputFile = new TFile((outputtreename).c_str(),"RECREATE"); 
  outputTree=new TTree("PU", "PU");
  outputTree->Branch("v_recHitTime", &v_recHitTime);

  std::cout<<"output tree name "<<std::endl;

  std::vector<double> averageTime;
  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  int nEventsPassed = 0;
  int nEventsTimePassed = 0;
  int nEventsRhoPassed = 0;
  for (int ievent=0; ievent<nEvents; ++ievent)
  {
    tree->GetEvent(ievent);

    if(ievent % 100 == 0) std::cout << "Processing event " << ievent << std::endl;

    v_recHitTime.clear();

    if(reachedEE==0) continue;
    nEventsPassed++;
    double vtx[3] = {vertex_x, vertex_y, vertex_z};    

    double etaGen = genParticle_eta->at(0);
    double phiGen = genParticle_phi->at(0);
    double xGen = 5.0;
    double yGen = 0.0;
    double zGen = 5.0*TMath::SinH(etaGen);
    double fromAxis[3] = {xGen, yGen, zGen}; 
 
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
      recHit.recHitEta = recHit_eta->at(k);
      recHit.recHitLayer = recHit_layer->at(k);
      h_recHit_Time->Fill(recHit.recHitTime);
      double rhZ = layerToZ(recHit.recHitLayer, recHit.recHitEta);
      double toRH[3] = {0., 0., rhZ};
      layerIntersection(toRH, fromAxis, vtx);
      double Radius_RhGen = sqrt(pow(toRH[0]-recHit.recHitX, 2) + pow(toRH[1]-recHit.recHitY, 2));
      if(Radius_RhGen<rhoCut and recHit_time->at(k) > -6) h_rho->Fill(Radius_RhGen);
      if(Radius_RhGen<rhoCut and recHit_time->at(k) > -6) sumOfEnergy += recHit_energy->at(k);
      if(Radius_RhGen<rhoCut and recHit_time->at(k) > -6) recHits.push_back(recHit);
    }

    std::sort(recHits.begin(), recHits.end(), sortRecHitsInAscendingTime);

    h_nHits_FullSize->Fill(recHits.size());
  
    if(recHits.size()<1) continue;
    nEventsTimePassed++;

    int recHitsFullSize = recHits.size();
    double sumTimeAverage = 0.0;
    double sumResolution = 0.0;
    int nhits = 0;
    for (unsigned int j=0; j<recHits.size(); j++)
    {
      v_recHitTime.push_back(recHits.at(j).recHitTime);
      if(recHits.size() >= 3) h_recHitX->Fill(recHits.at(j).recHitX);
      if(recHits.size() >= 3) h_recHitY->Fill(recHits.at(j).recHitY);
      if(recHits.size() >= 3) h_recHitZ->Fill(recHits.at(j).recHitZ);
      sumTimeAverage += recHits.at(j).recHitTime;
      nhits++;
    }
    nHitsAfterRho = nhits;
    if(nhits>0) averageTime.push_back(sumTimeAverage/nhits);
    if(nhits>0) h_TimeAverage->Fill(sumTimeAverage/nhits);
    outputTree->Fill();
  }//end of event loop

  outputTree->Write();
  std::sort(averageTime.begin(), averageTime.end(), sortByTimeOrder);
  int removeElementsLow = 0.1586*averageTime.size();
  int removeElementsHigh = 0.1586*averageTime.size();

  std::cout << "averageTime.at(removeElementsLow) = " << averageTime.at(removeElementsLow) << std::endl;
  //std::cout << "averageTime.at(removeElementsHigh) = " << averageTime.at(removeElementsLow+0.68*averageTime.size()) << std::endl;
  //std::cout << "Width = " << (averageTime.at(removeElementsLow+0.68*averageTime.size()) - averageTime.at(removeElementsLow))/2.0 << std::endl;
  std::cout << "nEventsPassed = " << nEventsPassed << std::endl;
  std::cout << "nEventsTimePassed = " << nEventsTimePassed << std::endl;
  

  std::string histfilename=(infile+"_outputHist.root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_rho->Write();
  h_recHitX->Write();
  h_recHitY->Write();
  h_recHitZ->Write();
  h_TimeAverage->Write();
  h_recHit_Time->Write();
  h_nHits_afterCut->Write();
  h_nHits_FullSize->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0; 
}
