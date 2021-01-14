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

TLorentzVector fillTLorentzVector(float pT, float eta, float phi, float E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  float rechitE;
  float rechitX;
  float rechitXSub;
  float rechitY;
  float rechitZ;
  float rechitTime;
  float rechitpT;
  float rechitEta;
  int rechitLayer;
} RecHitInfo;

typedef struct
{
  float genParticleEta;
  float genParticlePhi;
} genParticleInfo;

typedef struct
{
  float multiClusX;
  float multiClusY;
  float multiClusZ;
  float multiClusE;
} multiClusInfo;


bool sortByTimeOrder(float i, float j)
{
  return (i < j);
}

bool sortByPositionOrder(float i, float j)
{
  return (i < j);
}

bool sortRecHitsInDescendingE(RecHitInfo rechit1, RecHitInfo rechit2)
{
  return (rechit1.rechitE > rechit2.rechitE);
}

bool sortRecHitsInAscendingTime(RecHitInfo rechit1, RecHitInfo rechit2)
{
  return (rechit1.rechitTime < rechit2.rechitTime);
}

bool sortRecHitsInAscendingPositionDifference(RecHitInfo rechit1, RecHitInfo rechit2)
{
  return (fabs(rechit1.rechitXSub) < fabs(rechit2.rechitXSub));
}

bool sortRecHitsInAscendingZ(RecHitInfo rechit1, RecHitInfo rechit2)
{
  return (fabs(rechit1.rechitZ) < fabs(rechit2.rechitZ));
}

bool sortRecHitsInAscendingE(RecHitInfo rechit1, RecHitInfo rechit2)
{
  return (rechit1.rechitE < rechit2.rechitE);
}

bool sortSimImpactPt(float i, float j)
{
  return (i > j);
}

std::string itoa(int i) 
{
  char res[50];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

bool sameVal(float a, float b)
{
   return fabs(a - b) < 3.000e-03;
}

bool sameValHighPrecision(float a, float b)
{
   return fabs(a - b) < 1.000e-08;
}

bool sortMultiClusInDescendingE(multiClusInfo multiClus1, multiClusInfo multiClus2)
{
  return (multiClus1.multiClusE > multiClus2.multiClusE);
}

void layerIntersection(float to[3], float from[3], float fromB[3])
{
  to[0]=(from[0]-fromB[0]) / (from[2] - fromB[2]) * (to[2] - from[2]) + from[0];
  to[1]=(from[1]-fromB[1]) / (from[2] - fromB[2]) * (to[2] - from[2]) + from[1];
}

int CopyTree(std::string infile, float pt)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("ana/hgc");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  ULong64_t       event;
  UInt_t          lumi;
  UInt_t          run;
  Float_t         vtx_x;
  Float_t         vtx_y;
  Float_t         vtx_z;  
  vector<float>   *rechit_eta;
  vector<float>   *rechit_phi;
  vector<float>   *rechit_pt;
  vector<float>   *rechit_energy;
  vector<float>   *rechit_x;
  vector<float>   *rechit_y;
  vector<float>   *rechit_z;
  vector<float>   *rechit_time;
  vector<float>   *rechit_thickness;
  vector<int>     *rechit_layer;
  vector<float>   *genpart_eta;
  vector<float>   *genpart_phi;
  vector<float>   *genpart_pt;
  vector<float>   *cluster2d_eta;
  vector<float>   *cluster2d_phi;
  vector<float>   *cluster2d_pt;
  vector<float>   *cluster2d_energy;
  vector<float>   *cluster2d_x;
  vector<float>   *cluster2d_y;
  vector<float>   *cluster2d_z;
  vector<int>     *cluster2d_layer;
  vector<int>     *cluster2d_nhitCore;
  vector<int>     *cluster2d_nhitAll;
  vector<int>     *cluster2d_multicluster;
  vector<vector<unsigned int> > *cluster2d_rechits;
  vector<int>     *genpart_reachedEE;
  vector<float>   *simcluster_eta;
  vector<float>   *simcluster_phi;
  vector<float>   *simcluster_pt;
  vector<float>   *simcluster_energy;
  vector<int>     *simcluster_pid;
  vector<float>   *simcluster_simEnergy;
  vector<vector<unsigned int> > *simcluster_hits;
  vector<vector<int> > *simcluster_hits_indices;
  vector<double>  *simHit_time;
  vector<float>   *simcluster_time_impactPoint;
  vector<float>   *simcluster_x_impactPoint;
  vector<float>   *simcluster_y_impactPoint;
  vector<float>   *simcluster_z_impactPoint;
  vector<float>   *simcluster_eta_impactPoint;
  vector<float>   *simcluster_phi_impactPoint;
  vector<float>   *simcluster_pt_impactMomentumPt;
  vector<float>   *simcluster_pt_impactMomentumEta;
  vector<float>   *simcluster_pt_impactMomentumPhi;
  vector<float>   *simcluster_pt_impactMomentumE;
  vector<float>   *calopart_eta;
  vector<float>   *calopart_phi;
  vector<float>   *calopart_pt;
  vector<float>   *calopart_energy;
  vector<float>   *calopart_simEnergy;
  vector<float>   *multiclus_eta;
  vector<float>   *multiclus_phi;
  vector<float>   *multiclus_pt;
  vector<float>   *multiclus_energy;
  vector<float>   *multiclus_z;
  vector<float>   *multiclus_slopeX;
  vector<float>   *multiclus_slopeY;
  vector<float>   *gen_eta;
  vector<float>   *gen_phi;
  vector<float>   *gen_pt;
  vector<float>   *gen_energy;
  vector<int>     *gen_charge;
  vector<int>     *gen_pdgid;
  vector<int>     *gen_status;
  vector<vector<int> > *gen_daughters;

  rechit_eta = 0;
  rechit_phi = 0;
  rechit_pt = 0;
  rechit_x = 0;
  rechit_y = 0;
  rechit_z = 0;
  rechit_time = 0;
  rechit_thickness = 0;
  rechit_layer = 0;
  rechit_energy = 0;
  genpart_eta = 0;
  genpart_phi = 0;
  genpart_pt = 0;
  cluster2d_eta = 0;
  cluster2d_phi= 0;
  cluster2d_pt = 0;
  cluster2d_energy = 0;
  cluster2d_x = 0;
  cluster2d_y = 0;
  cluster2d_z = 0;
  cluster2d_layer = 0;
  cluster2d_nhitCore = 0;
  cluster2d_nhitAll = 0;
  cluster2d_multicluster = 0;
  cluster2d_rechits = 0;
  genpart_reachedEE = 0;
  simcluster_eta = 0;
  simcluster_phi = 0;
  simcluster_pt = 0;
  simcluster_energy = 0;
  simcluster_pid = 0;
  simcluster_simEnergy = 0;
  simcluster_hits = 0;
  simcluster_hits_indices = 0;
  simHit_time = 0;
  simcluster_time_impactPoint = 0;
  simcluster_x_impactPoint = 0;
  simcluster_y_impactPoint = 0;
  simcluster_z_impactPoint = 0;
  simcluster_eta_impactPoint = 0;
  simcluster_phi_impactPoint = 0;
  simcluster_pt_impactMomentumPt = 0;
  simcluster_pt_impactMomentumPhi = 0; 
  simcluster_pt_impactMomentumEta = 0;
  simcluster_pt_impactMomentumE = 0; 
  calopart_eta = 0;
  calopart_phi = 0;
  calopart_pt = 0;
  calopart_energy = 0;
  calopart_simEnergy = 0;
  multiclus_eta = 0;
  multiclus_phi = 0;
  multiclus_pt = 0;
  multiclus_energy = 0;
  multiclus_z = 0;
  multiclus_slopeX = 0;
  multiclus_slopeY = 0;
  gen_eta = 0;
  gen_phi = 0;
  gen_pt = 0;
  gen_energy = 0;
  gen_charge = 0;
  gen_pdgid = 0;
  gen_status = 0;
  gen_daughters = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("event", &(event));
  tree->SetBranchAddress("vtx_x", &(vtx_x));
  tree->SetBranchAddress("vtx_y", &(vtx_y));
  tree->SetBranchAddress("vtx_z", &(vtx_z));
  tree->SetBranchAddress("rechit_energy", &(rechit_energy));
  tree->SetBranchAddress("rechit_x", &(rechit_x));
  tree->SetBranchAddress("rechit_y", &(rechit_y));
  tree->SetBranchAddress("rechit_z", &(rechit_z));
  tree->SetBranchAddress("rechit_time", &(rechit_time));
  tree->SetBranchAddress("rechit_pt", &(rechit_pt));
  tree->SetBranchAddress("genpart_eta", &(genpart_eta));
  tree->SetBranchAddress("genpart_phi", &(genpart_phi));
  tree->SetBranchAddress("rechit_eta", &(rechit_eta));
  tree->SetBranchAddress("rechit_layer", &(rechit_layer));
  tree->SetBranchAddress("cluster2d_eta", &(cluster2d_eta));
  tree->SetBranchAddress("cluster2d_phi", &(cluster2d_phi));
  tree->SetBranchAddress("cluster2d_energy", &(cluster2d_energy));
  tree->SetBranchAddress("cluster2d_x", &(cluster2d_x));
  tree->SetBranchAddress("cluster2d_y", &(cluster2d_y));
  tree->SetBranchAddress("cluster2d_z", &(cluster2d_z));
  tree->SetBranchAddress("cluster2d_layer", &(cluster2d_layer));
  tree->SetBranchAddress("cluster2d_nhitCore", &(cluster2d_nhitCore));
  tree->SetBranchAddress("cluster2d_nhitAll", &(cluster2d_nhitAll));
  tree->SetBranchAddress("cluster2d_multicluster", &(cluster2d_multicluster));
  tree->SetBranchAddress("cluster2d_rechits", &(cluster2d_rechits));
  tree->SetBranchAddress("genpart_reachedEE", &(genpart_reachedEE)); 
  tree->SetBranchAddress("simcluster_eta", &(simcluster_eta));
  tree->SetBranchAddress("simcluster_phi", &(simcluster_phi));
  tree->SetBranchAddress("simcluster_pt", &(simcluster_pt));
  tree->SetBranchAddress("simcluster_energy", &(simcluster_energy));
  tree->SetBranchAddress("simcluster_pid", &(simcluster_pid));
  tree->SetBranchAddress("simcluster_simEnergy", &(simcluster_simEnergy));
  tree->SetBranchAddress("simcluster_hits", &(simcluster_hits));
  tree->SetBranchAddress("simcluster_hits_indices", &(simcluster_hits_indices));
  tree->SetBranchAddress("simHit_time", &(simHit_time));
  tree->SetBranchAddress("simcluster_time_impactPoint", &(simcluster_time_impactPoint));
  tree->SetBranchAddress("simcluster_x_impactPoint", &(simcluster_x_impactPoint));
  tree->SetBranchAddress("simcluster_y_impactPoint", &(simcluster_y_impactPoint));
  tree->SetBranchAddress("simcluster_z_impactPoint", &(simcluster_z_impactPoint));
  tree->SetBranchAddress("simcluster_eta_impactPoint", &(simcluster_eta_impactPoint));
  tree->SetBranchAddress("simcluster_phi_impactPoint", &(simcluster_phi_impactPoint));
  tree->SetBranchAddress("simcluster_pt_impactMomentumPt", &(simcluster_pt_impactMomentumPt));
  tree->SetBranchAddress("simcluster_pt_impactMomentumPhi", &(simcluster_pt_impactMomentumPhi));
  tree->SetBranchAddress("simcluster_pt_impactMomentumEta", &(simcluster_pt_impactMomentumEta));
  tree->SetBranchAddress("simcluster_pt_impactMomentumE", &(simcluster_pt_impactMomentumE));
  tree->SetBranchAddress("calopart_eta", &(calopart_eta));
  tree->SetBranchAddress("calopart_phi", &(calopart_phi));
  tree->SetBranchAddress("calopart_pt", &(calopart_pt));
  tree->SetBranchAddress("calopart_energy", &(calopart_energy));
  tree->SetBranchAddress("calopart_simEnergy", &(calopart_simEnergy));
  tree->SetBranchAddress("multiclus_eta", &(multiclus_eta));
  tree->SetBranchAddress("multiclus_phi", &(multiclus_phi));
  tree->SetBranchAddress("multiclus_pt", &(multiclus_pt));
  tree->SetBranchAddress("multiclus_energy", &(multiclus_energy));
  tree->SetBranchAddress("multiclus_z", &(multiclus_z));
  tree->SetBranchAddress("multiclus_slopeX", &(multiclus_slopeX));
  tree->SetBranchAddress("multiclus_slopeY", &(multiclus_slopeY));
  tree->SetBranchAddress("gen_eta", &(gen_eta));
  tree->SetBranchAddress("gen_phi", &(gen_phi));
  tree->SetBranchAddress("gen_pt", &(gen_pt));
  tree->SetBranchAddress("gen_energy", &(gen_energy)); 
  tree->SetBranchAddress("gen_charge", &(gen_charge));
  tree->SetBranchAddress("gen_pdgid", &(gen_pdgid));
  tree->SetBranchAddress("gen_status", &(gen_status));
  tree->SetBranchAddress("gen_daughters", &(gen_daughters));

  TFile *outputFile;
  TTree *outputTree;
  int nHitsAfterRho;
  std::vector<float> v_recHitTime;
  std::vector<float> v_recHitX;
  std::vector<float> v_recHitY;
  std::vector<float> v_recHitZ;
  std::vector<float> v_recHitE;
  std::vector<float> v_recHitLayer;
  std::vector<float> v_recHitEta;
  std::vector<float> v_simPosX;
  std::vector<float> v_simPosY;
  std::vector<float> v_simPosZ;
  std::vector<float> v_simTime;
  std::vector<float> v_simID;
  std::vector<float> v_simEta;
  std::vector<float> v_simPhi;
  std::vector<float> v_simE;
  std::vector<std::vector<int> > v_simHitsIndices;
  std::cout<<"output filename "<<std::endl;

  std::string outputtreename=(infile+"_outputTree.root").c_str();
  outputFile = new TFile((outputtreename).c_str(),"RECREATE"); 
  outputTree=new TTree("ROI", "ROI");
  outputTree->Branch("v_recHitTime", &v_recHitTime);  
  outputTree->Branch("v_recHitX", &v_recHitX);
  outputTree->Branch("v_recHitY", &v_recHitY);
  outputTree->Branch("v_recHitZ", &v_recHitZ);
  outputTree->Branch("v_recHitE", &v_recHitE);
  outputTree->Branch("nHitsAfterRho", &nHitsAfterRho);
  outputTree->Branch("v_simPosX", &v_simPosX);
  outputTree->Branch("v_simPosY", &v_simPosY);
  outputTree->Branch("v_simPosZ", &v_simPosZ);
  outputTree->Branch("v_simTime", &v_simTime);
  outputTree->Branch("v_simID", &v_simID);
  outputTree->Branch("v_simEta", &v_simEta);
  outputTree->Branch("v_simPhi", &v_simPhi);
  outputTree->Branch("v_simE", &v_simE);
  outputTree->Branch("v_simHitsIndices", &v_simHitsIndices);  
  outputTree->Branch("v_recHitLayer", &v_recHitLayer);
  outputTree->Branch("v_recHitEta", &v_recHitEta);

  std::vector<float> averagePosition;
  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  int nEvents_passed=0;
  int rho_passed=0;
  int nEvents_passed_1hit=0;
  int n_eplus, n_eminus;
  n_eplus = n_eminus = 0;
  for (int ievent=0; ievent<nEvents; ++ievent)
  {
    tree->GetEvent(ievent);
    float vtx[3] = {0.0, 0.0, 0.0}; //consider replacing with original vertex
    float etaGen = gen_eta->at(0);
    float phiGen = gen_phi->at(0);
    float xGen = vtx_x;
    float yGen = vtx_y;
    float zGen = xGen*TMath::SinH(etaGen);
    float fromAxis[3] = {xGen, yGen, zGen};
    v_recHitTime.clear();
    v_recHitX.clear();
    v_recHitY.clear();
    v_recHitZ.clear();
    v_recHitE.clear();
    v_simPosX.clear();
    v_simPosY.clear();
    v_simPosZ.clear();
    v_simTime.clear();
    v_simID.clear();
    v_simEta.clear();
    v_simPhi.clear();
    v_simE.clear();
    v_simHitsIndices.clear();
    v_recHitLayer.clear();
    v_recHitEta.clear();
    std::vector<RecHitInfo> rechits;
    for (unsigned int k=0; k<rechit_energy->size(); k++)
    {
      RecHitInfo rechit;
      rechit.rechitE = rechit_energy->at(k);
      rechit.rechitX = rechit_x->at(k);
      rechit.rechitXSub = rechit_x->at(k) - xGen;
      rechit.rechitY = rechit_y->at(k);
      rechit.rechitZ = rechit_z->at(k);
      rechit.rechitTime = rechit_time->at(k) - 5;
      rechit.rechitpT = rechit_pt->at(k);
      rechit.rechitEta = rechit_eta->at(k);
      rechit.rechitLayer = rechit_layer->at(k);
      float toRH[3] = {0., 0., rechit.rechitZ};
      layerIntersection(toRH, fromAxis, vtx);
      float Radius_RhGen = sqrt(pow(toRH[0]-rechit.rechitX, 2) + pow(toRH[1]-rechit.rechitY, 2));
      if(Radius_RhGen < 5000.0) rechits.push_back(rechit);
    }   

    //implement a smart sorting algorithm
    std::sort(rechits.begin(), rechits.end(), sortRecHitsInAscendingPositionDifference); 

    for (unsigned int j=0; j<rechits.size(); j++)
    {
      if(rechits.size() >= 3)
      {
        v_recHitTime.push_back(rechits.at(j).rechitTime);
        v_recHitX.push_back(rechits.at(j).rechitX);
        v_recHitY.push_back(rechits.at(j).rechitY);
        v_recHitZ.push_back(rechits.at(j).rechitZ);
        v_recHitE.push_back(rechits.at(j).rechitE);
        v_recHitLayer.push_back(rechits.at(j).rechitLayer);
        v_recHitEta.push_back(rechits.at(j).rechitEta);
      }
    }
    for(unsigned int i=0; i<simcluster_x_impactPoint->size(); i++)
    {
      v_simPosX.push_back(simcluster_x_impactPoint->at(i));
      v_simPosY.push_back(simcluster_y_impactPoint->at(i));
      v_simPosZ.push_back(simcluster_z_impactPoint->at(i));
      v_simTime.push_back(simcluster_time_impactPoint->at(i));
    }
    for(unsigned int k=0; k<simcluster_eta->size(); k++)
    {
      v_simID.push_back(simcluster_pid->at(k));
      v_simEta.push_back(simcluster_eta->at(k));
      v_simPhi.push_back(simcluster_phi->at(k));
      v_simE.push_back(simcluster_energy->at(k));
      v_simHitsIndices.push_back(simcluster_hits_indices->at(k));
    }

    outputTree->Fill();
  }//end of event loop

  outputTree->Write();

  outputFile->Close();
  std::cout<<"Wrote output tree "<< outputtreename << std::endl;
  return 0; 
}
