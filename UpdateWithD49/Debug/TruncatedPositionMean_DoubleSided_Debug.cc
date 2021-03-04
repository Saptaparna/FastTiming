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
#include <TPrincipal.h>

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

void layerIntersection_v2(float to[3], float from[3], float fromB[3])
{
  to[0]=(from[0]-fromB[0]) / (fabs(from[2]) - fromB[2]) * (to[2] - fabs(from[2])) + from[0];
  to[1]=(from[1]-fromB[1]) / (fabs(from[2]) - fromB[2]) * (to[2] - fabs(from[2])) + from[1];
}

int TruncatedPositionMean_DoubleSided(std::string infile, std::string outfile, float pt, float rhoCut, float fractionLow, float fractionHigh)
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
  vector<vector<unsigned int> > *multiclus_cluster2d;

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
  multiclus_cluster2d = 0;

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
  tree->SetBranchAddress("multiclus_cluster2d", &(multiclus_cluster2d));
  tree->SetBranchAddress("gen_eta", &(gen_eta));
  tree->SetBranchAddress("gen_phi", &(gen_phi));
  tree->SetBranchAddress("gen_pt", &(gen_pt));
  tree->SetBranchAddress("gen_energy", &(gen_energy)); 
  tree->SetBranchAddress("gen_charge", &(gen_charge));
  tree->SetBranchAddress("gen_pdgid", &(gen_pdgid));
  tree->SetBranchAddress("gen_status", &(gen_status));
  tree->SetBranchAddress("gen_daughters", &(gen_daughters));

  TH1D *h_rho = new TH1D("h_rho", "h_rho; rho [cm]; Entries", 1000, 0.0, 100.0); h_rho->Sumw2();
  TH1D *h_meanXposition = new TH1D("h_meanXposition", "h_meanXposition; mean X position; Events", 300, 90, 120.0); h_meanXposition->Sumw2();
  TH1D *h_meanYposition = new TH1D("h_meanYposition", "h_meanYposition; mean Y position; Events", 100, -5.0, 5.0); h_meanYposition->Sumw2();
  TH1D *h_meanZposition = new TH1D("h_meanZposition", "h_meanZposition; mean Z position; Events", 2000, 300, 500.0); h_meanZposition->Sumw2();
  
  /*TH1D *h_meanXposition_multiclus = new TH1D("h_meanXposition_multiclus", "h_meanXposition_multiclus; mean X position; Events", 300, 90, 120.0); h_meanXposition_multiclus->Sumw2();
  TH1D *h_meanYposition_multiclus = new TH1D("h_meanYposition_multiclus", "h_meanYposition_multiclus; mean Y position; Events", 100, -5.0, 5.0); h_meanYposition_multiclus->Sumw2();
  TH1D *h_meanZposition_multiclus = new TH1D("h_meanZposition_multiclus", "h_meanZposition_multiclus; mean Z position; Events", 2000, 300, 500.0); h_meanZposition_multiclus->Sumw2();
  */
  TH1D *h_meanXposition_multiclus = new TH1D("h_meanXposition_multiclus", "h_meanXposition_multiclus; mean X position; Events", 6000, -300, 300.0); h_meanXposition_multiclus->Sumw2();
  TH1D *h_meanYposition_multiclus = new TH1D("h_meanYposition_multiclus", "h_meanYposition_multiclus; mean Y position; Events", 6000, -300, 300.0); h_meanYposition_multiclus->Sumw2();
  TH1D *h_meanZposition_multiclus = new TH1D("h_meanZposition_multiclus", "h_meanZposition_multiclus; mean Z position; Events", 10000, -500, 500.0); h_meanZposition_multiclus->Sumw2();

  TH1D *h_PositionAverage_X = new TH1D("h_PositionAverage_X", "h_PositionAverage_X;  mean X position; Events", 1200, -30, 30.0); h_PositionAverage_X->Sumw2();
  TH1D *h_PositionAverage_Y = new TH1D("h_PositionAverage_Y", "h_PositionAverage_Y;  mean Y position; Events", 1200, -30, 30.0); h_PositionAverage_Y->Sumw2();

  TH1D *h_nMult = new TH1D("h_nMult", "h_nMult; number of multiclusters; Events", 20, -0.5, 19.5); h_nMult->Sumw2();
  TH1D *h_nMult_AfterCuts = new TH1D("h_nMult_AfterCuts", "h_nMult_AfterCuts; number of multiclusters; Events", 10, -0.5, 9.5); h_nMult_AfterCuts->Sumw2();
  TH1D *h_res_Xposition = new TH1D("h_res_Xposition", "h_res_Xposition; X position [pred - truth [cm]]; Events", 100000, -10, 10); h_res_Xposition->Sumw2();
  TH1D *h_res_Yposition = new TH1D("h_res_Yposition", "h_res_Yposition; Y position [pred - truth [cm]]; Events", 100000, -10, 10); h_res_Yposition->Sumw2();
  //TH2D *h_reco_gen_X = new TH2D("h_reco_gen_X", "h_reco_gen_X; reco X cm; gen X cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_X->Sumw2();
  //TH2D *h_reco_gen_Y = new TH2D("h_reco_gen_Y", "h_reco_gen_Y; reco Y cm; gen Y cm", 1000, -0.5, 0.5, 1000, -0.5, 0.5); h_reco_gen_Y->Sumw2();
  TH2D *h_reco_gen_X = new TH2D("h_reco_gen_X", "h_reco_gen_X; gen X cm; reco X cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_X->Sumw2();
  TH2D *h_reco_gen_Y = new TH2D("h_reco_gen_Y", "h_reco_gen_Y; gen Y cm; reco Y cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_Y->Sumw2();

  TH2D *h_gen_X_Y = new TH2D("h_gen_X_Y", "h_gen_X_Y; gen X cm; gen Y cm", 6000, -300, 300.0, 6000, -300, 300.0); h_gen_X_Y->Sumw2();
  TH2D *h_reco_X_Y = new TH2D("h_reco_X_Y", "h_reco_X_Y; reco X cm; reco Y cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_X_Y->Sumw2();

  TH2D *h_reco_gen_X_perhit = new TH2D("h_reco_gen_X_perhit", "h_reco_gen_X_perhit; gen X cm; reco X cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_X_perhit->Sumw2();
  TH2D *h_reco_gen_Y_perhit = new TH2D("h_reco_gen_Y_perhit", "h_reco_gen_Y_perhit; gen Y cm; reco Y cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_Y_perhit->Sumw2();

  TH2D *h_reco_gen_X_perhit_beforeCuts = new TH2D("h_reco_gen_X_perhit_beforeCuts", "h_reco_gen_X_perhit_beforeCuts; gen X cm; reco X cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_X_perhit_beforeCuts->Sumw2();
  TH2D *h_reco_gen_Y_perhit_beforeCuts = new TH2D("h_reco_gen_Y_perhit_beforeCuts", "h_reco_gen_Y_perhit_beforeCuts; gen X cm; reco X cm", 6000, -300, 300.0, 6000, -300, 300.0); h_reco_gen_Y_perhit_beforeCuts->Sumw2();

  TH1D *h_nHits_afterCut = new TH1D("h_nHits_afterCut", "h_nHits_afterCut; Number of hits; Events", 500, 0.0, 500.0); h_nHits_afterCut->Sumw2(); 

  TH1D *h_rechitX = new TH1D("h_rechitX", "h_rechitX; rechit X [cm]; Entries", 150, 0.0, 150.0);h_rechitX->Sumw2();
  TH1D *h_rechitY = new TH1D("h_rechitY", "h_rechitY; rechit Y [cm]; Entries", 100, 50.0, -50.0);h_rechitY->Sumw2();
  TH1D *h_rechitZ = new TH1D("h_rechitZ", "h_rechitZ; rechit Z [cm]; Entries", 100, 300.0, 400.0);h_rechitZ->Sumw2();
  TH1D *h_rechitX_BeforeTruncation = new TH1D("h_rechitX_BeforeTruncation", "h_rechitX_BeforeTruncation; rechit X [cm]; Entries", 150, 0.0, 150.0);h_rechitX_BeforeTruncation->Sumw2();
  TH2D *h_rechit_rho_energy = new TH2D("h_rechit_rho_energy", "h_rechit_rho_energy; rechit #rho; rechit energy", 1000, 0.0, 100.0, 100, 0.0, 1.0);h_rechit_rho_energy->Sumw2();
  TH2D *h_multiclus_clus2D_energy = new TH2D("h_multiclus_clus2D_energy", "h_multiclus_clus2D_energy; Number of 2D clusters; Multicluster energy", 15, -0.5, 14.5, 100, 0.0, 50.0); h_multiclus_clus2D_energy->Sumw2(); 

  TH1D *h_nHits_beforeCut = new TH1D("h_nHits_beforeCut", "h_nHits_beforeCut; Number of hits; Events", 100000, -0.5, 99999.5); h_nHits_beforeCut->Sumw2();

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
  //int nEvents=1;//50;
  std::cout << "nEvents= " << nEvents << std::endl;
  int nEvents_passed=0;
  int rho_passed=0;
  int nEvents_passed_1hit=0;
  int n_eplus, n_eminus;
  n_eplus = n_eminus = 0;
  //TPrincipal* principal = new TPrincipal(3,"ND");
  for (int ievent=0; ievent<nEvents; ++ievent)
  {
    tree->GetEvent(ievent);
    float vtx[3] = {vtx_x, vtx_y, vtx_z};
    //float vtx[3] = {0.0, 0.0, 0.0}; //consider replacing with original vertex
    float etaGen = gen_eta->at(0);
    float phiGen = gen_phi->at(0);
    float xGen = gen_pt->at(0)*TMath::Cos(phiGen);//vtx_x;
    float yGen = gen_pt->at(0)*TMath::Sin(phiGen);//vtx_y;
    float zGen = gen_pt->at(0)*TMath::SinH(etaGen);
    //float zGen = xGen*TMath::SinH(etaGen);
    //double tanthetav = TMath::Sqrt(pow(vtx_x, 2) + pow(vtx_y, 2))/vtx_z;
    /*std::cout << "vtx_x = " << vtx_x << std::endl;
    std::cout << "vtx_y = " << vtx_y << std::endl;
    std::cout << "phiGen = " << phiGen << std::endl;
    std::cout << "etaGen = " << etaGen << std::endl;
    std::cout << "xGen = " << xGen << std::endl;
    std::cout << "yGen = " << yGen << std::endl;
    std::cout << "zGen = " << zGen << std::endl;
    */
    float rGen =  TMath::Sqrt(pow(xGen, 2) + pow(yGen, 2));
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
    h_nHits_beforeCut->Fill(rechit_energy->size());
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
      //float rh = (rechit.rechitZ - zGen)*tanthetav + rGen;
      //float xh = rh*TMath::Cos(phiGen);
      //float yh = rh*TMath::Sin(phiGen);
      layerIntersection(toRH, fromAxis, vtx);
      //if(rechit.rechitE > 0.02) std::cout << rechit.rechitX << "," << rechit.rechitY << "," << rechit.rechitZ << "," << rechit.rechitTime  << "," << rechit.rechitE << "," << grade << std::endl;
      //if(rechit.rechitE > 0.02) std::cout << rechit.rechitX << "," << rechit.rechitY << "," << rechit.rechitZ << std::endl;
      //if(rechit.rechitE > 0.05) std::cout << toRH[0]  << "," << toRH[1]  << "," << toRH[2] << std::endl;
      //principal->AddRow(rechit.rechitX);
      //principal->AddRow(rechit.rechitY);
      //principal->AddRow(rechit.rechitZ);
      if(rechit.rechitE > 0.05 and rechit.rechitTime > 0.) //0.05
      {
        float Radius_RhGen = sqrt(pow(toRH[0]-rechit.rechitX, 2) + pow(toRH[1]-rechit.rechitY, 2));
        h_rho->Fill(Radius_RhGen);
        h_rechit_rho_energy->Fill(Radius_RhGen, rechit.rechitE); 
        /*std::cout << "toRH[0] = " << toRH[0] << std::endl;
        std::cout << "toRH[1] = " << toRH[1] << std::endl;
        std::cout << "rechit.rechitX = " << rechit.rechitX << std::endl;
        std::cout << "rechit.rechitY = " << rechit.rechitY << std::endl;
        std::cout << "rechit.rechitZ = " << rechit.rechitZ << std::endl;
        */
        h_reco_X_Y->Fill(rechit.rechitX, rechit.rechitY);
        //h_gen_X_Y->Fill(xh, yh);
        h_gen_X_Y->Fill(toRH[0], toRH[1]);
        h_reco_gen_X_perhit_beforeCuts->Fill(toRH[0], rechit.rechitX);
        h_reco_gen_Y_perhit_beforeCuts->Fill(toRH[1], rechit.rechitY);
        //std::cout << toRH[0]  << "," << toRH[1]  << "," << toRH[2] << std::endl;
        if(Radius_RhGen < rhoCut) rechits.push_back(rechit);
        if(Radius_RhGen < rhoCut) h_rechitX_BeforeTruncation->Fill(rechit.rechitX);
        //if(Radius_RhGen < 2.0) rechits.push_back(rechit);
        //if(Radius_RhGen < 2.0) h_rechitX_BeforeTruncation->Fill(rechit.rechitX);
        //rechits.push_back(rechit);
        //h_rechitX_BeforeTruncation->Fill(rechit.rechitX);
      }
    }   
    //std::cout << "principal->Print() = " << principal->Print() << std::endl;
    //implement a smart sorting algorithm
    //std::sort(rechits.begin(), rechits.end(), sortRecHitsInAscendingPositionDifference); 

    double sumPositionAverageX = 0.0;
    double sumPositionAverageY = 0.0;
    double sumResolution = 0.0;
    double sumEnergy= 0.0;
    int nhits = 0;
    double toRH_average_x = 0;
    double toRH_average_y = 0;
    for (unsigned int j=0; j<rechits.size(); j++)
    {
      //if(rechits.at(j).rechitLayer < 28.0) continue;
      if(rechits.size() >= 3) h_rechitX->Fill(rechits.at(j).rechitX);
      if(rechits.size() >= 3) h_rechitY->Fill(rechits.at(j).rechitY);
      if(rechits.size() >= 3) h_rechitZ->Fill(rechits.at(j).rechitZ);
      double rho = sqrt(pow((rechits.at(j).rechitX-h_rechitX->GetMean()), 2) + pow(rechits.at(j).rechitY, 2));
      if(rechits.size() >= 3) sumPositionAverageX += rechits.at(j).rechitX*rechits.at(j).rechitE;
      if(rechits.size() >= 3) sumPositionAverageY += rechits.at(j).rechitY*rechits.at(j).rechitE;
      if(rechits.size() >= 3) nhits++;
      if(rechits.size() >= 3) sumEnergy += rechits.at(j).rechitE;
      float toRH[3] = {0., 0., rechits.at(j).rechitZ};
      layerIntersection(toRH, fromAxis, vtx);
      if(rechits.size() >= 3) toRH_average_x += toRH[0];
      h_reco_gen_X_perhit->Fill(toRH[0], rechits.at(j).rechitX);
      h_reco_gen_Y_perhit->Fill(toRH[1], rechits.at(j).rechitY);
      if(rechits.size() >= 3) toRH_average_y += toRH[1];
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
    if(nhits>0) h_PositionAverage_X->Fill(sumPositionAverageX/sumEnergy - toRH_average_x/nhits);
    if(nhits>0) h_PositionAverage_Y->Fill(sumPositionAverageY/sumEnergy - toRH_average_y/nhits);
    if(nhits>0) h_reco_gen_X->Fill(toRH_average_x/nhits, sumPositionAverageX/sumEnergy);
    if(nhits>0) h_reco_gen_Y->Fill(toRH_average_y/nhits, sumPositionAverageY/sumEnergy); 
    //h_reco_gen_X->Fill(toRH_average_x, sumPositionAverageX/sumEnergy);
    //h_reco_gen_Y->Fill(toRH_average_y, sumPositionAverageY/sumEnergy); 
    if(nhits>0) h_nHits_afterCut->Fill(rechits.size());
    if(nhits>0) averagePosition.push_back(sumPositionAverageX/sumEnergy);
    if(nhits>0) rho_passed++;  
    if(nhits>0)
    {
      nHitsAfterRho = nhits;
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
    }

    double sumXPosition = 0.0;
    double sumYPosition = 0.0;
    double sumZPosition = 0.0;
    for(unsigned int i=0; i<cluster2d_energy->size(); i++)
    {
      sumXPosition += cluster2d_x->at(i)*cluster2d_energy->at(i);
      sumYPosition += cluster2d_y->at(i)*cluster2d_energy->at(i);
      sumZPosition += cluster2d_z->at(i)*cluster2d_energy->at(i);
      sumEnergy += cluster2d_energy->at(i);
    }
    h_meanXposition->Fill(sumXPosition/sumEnergy);
    h_meanYposition->Fill(sumYPosition/sumEnergy);
    h_meanZposition->Fill(sumZPosition/sumEnergy); 
    double sumXPosition_multiClus = 0.0;
    double sumYPosition_multiClus = 0.0;
    double sumZPosition_multiClus = 0.0;
    double sumEnergy_multiClus = 0.0;
    int nMult = 0;
    h_nMult->Fill(multiclus_pt->size());
    double genXdir = 0.0;
    double genYdir = 0.0;
    //std::cout << "multiclus_cluster2d->size() = " << multiclus_cluster2d->size() << std::endl;
    for(unsigned int i=0; i<multiclus_pt->size(); i++)
    {
      
      if(multiclus_energy->at(i) > 40.0)
      { 
        nMult += 1;
        if(nMult!=1) continue;
       
        sumXPosition_multiClus += multiclus_slopeX->at(i)*multiclus_energy->at(i);
        sumYPosition_multiClus += multiclus_slopeY->at(i)*multiclus_energy->at(i);
        sumZPosition_multiClus += multiclus_z->at(i)*multiclus_energy->at(i);
        sumEnergy_multiClus += multiclus_energy->at(i);
        float toMC[3] = {0., 0., zGen};
        layerIntersection(toMC, fromAxis, vtx);
        
        genXdir = toMC[0];
        genYdir = toMC[1];
        float Radius_MCGen = sqrt(pow(toMC[0]-multiclus_slopeX->at(i), 2) + pow(toMC[1]-multiclus_slopeY->at(i), 2));//toMC[0];//sqrt(pow(toMC[0]-multiclus_slopeX->at(i), 2) + pow(toMC[1]-multiclus_slopeY->at(i), 2));
        /*if(Radius_MCGen < 0.5)
        {
          sumXPosition_multiClus += multiclus_slopeX->at(i)*multiclus_energy->at(i);
          sumYPosition_multiClus += multiclus_slopeY->at(i)*multiclus_energy->at(i);
          sumZPosition_multiClus += multiclus_z->at(i)*multiclus_energy->at(i);
          sumEnergy_multiClus += multiclus_energy->at(i);
        }*/
      }
    }
    h_nMult_AfterCuts->Fill(nMult); 
       
    //std::cout << "nMult = " << nMult << std::endl;  
    //h_meanXposition_multiclus->Fill(sumXPosition_multiClus/sumEnergy_multiClus);
    //h_meanYposition_multiclus->Fill(sumYPosition_multiClus/sumEnergy_multiClus);
    //h_meanZposition_multiclus->Fill(sumZPosition_multiClus/sumEnergy_multiClus);
    //h_res_Xposition->Fill(sumXPosition_multiClus/sumEnergy_multiClus - genXdir);
    //h_res_Yposition->Fill(sumYPosition_multiClus/sumEnergy_multiClus - genYdir);
    
    std::vector<multiClusInfo> v_multiClus;
    float genXdir_mc = 0.0;
    float genYdir_mc = 0.0;
    std::vector<int> flatten_multiclus_cluster2d;
    flatten_multiclus_cluster2d.clear();
    for (unsigned int k=0; k<multiclus_pt->size(); k++)
    {
      multiClusInfo multiClus;
      multiClus.multiClusX = multiclus_slopeX->at(k);
      multiClus.multiClusY = multiclus_slopeY->at(k);   
      multiClus.multiClusZ = multiclus_z->at(k);
      multiClus.multiClusE = multiclus_energy->at(k);
      v_multiClus.push_back(multiClus);
      //if(multiClus.multiClusE > 40) std::cout << "multiclus_cluster2d->at(k) = " << multiclus_cluster2d->at(k)[0] << std::endl;
      //flatten_multiclus_cluster2d.push_back(multiclus_cluster2d->at(k));
    }

    std::sort(v_multiClus.begin(), v_multiClus.end(), sortMultiClusInDescendingE);
    float toMC_MC[3];
    if(v_multiClus.size() > 0) 
    {
      toMC_MC[0] = 0.0;
      toMC_MC[1] = 0.0;
      toMC_MC[2] = v_multiClus.at(0).multiClusZ;
      //std::cout << "v_multiClus.at(0).multiClusZ = " << v_multiClus.at(0).multiClusZ << std::endl;
      layerIntersection(toMC_MC, fromAxis, vtx);
      genXdir_mc = toMC_MC[0];//xGen;//toMC_MC[0];
      genYdir_mc = toMC_MC[1];
      //std::cout << "v_multiClus.at(0).multiClusX = " << v_multiClus.at(0).multiClusX << std::endl;
      //std::cout << "toMC_MC[0] = " << toMC_MC[0] << std::endl;
      //std::cout << "v_multiClus.at(0).multiClusE = " << v_multiClus.at(0).multiClusE << std::endl;
      //std::cout << "genXdir_mc = " << genXdir_mc << std::endl;
      //std::cout << "genYdir_mc = " << genYdir_mc << std::endl;
      h_res_Xposition->Fill(v_multiClus.at(0).multiClusX - genXdir_mc);
      h_res_Yposition->Fill(v_multiClus.at(0).multiClusY - genYdir_mc);
      h_meanXposition_multiclus->Fill(v_multiClus.at(0).multiClusX);
      h_meanYposition_multiclus->Fill(v_multiClus.at(0).multiClusY);
      h_meanZposition_multiclus->Fill(v_multiClus.at(0).multiClusZ); 
      //h_multiclus_clus2D_energy->Fill(multiclus_cluster2d->size(), v_multiClus.at(0).multiClusE);
      h_multiclus_clus2D_energy->Fill(flatten_multiclus_cluster2d.size(), v_multiClus.at(0).multiClusE);
    }
    outputTree->Fill();
  }//end of event loop

  outputTree->Write();
  //std::cout << "nEvents_passed = " << nEvents_passed << std::endl;  
  //std::cout << "nEvents_passed_1hit = " << nEvents_passed_1hit << std::endl;
  //std::cout << "rho_passed = " << rho_passed << std::endl;

  //std::sort(averagePosition.begin(), averagePosition.end(), sortByPositionOrder);

  int removeElementsLow = 0.1586*averagePosition.size();
  /*
  std::cout << "averagePosition.at(removeElementsLow) = " << averagePosition.at(removeElementsLow) << std::endl;
  std::cout << "averagePosition.at(removeElementsHigh) = " << averagePosition.at(removeElementsLow+0.68*averagePosition.size()) << std::endl;
  std::cout << "68% CI = " << (averagePosition.at(removeElementsLow+0.68*averagePosition.size()) - averagePosition.at(removeElementsLow))/2.0 << std::endl;
  */
  int fracPercentLow = fractionLow*100;
  std::string strFracPercentLow = std::to_string(fracPercentLow);

  int fracPercentHigh = fractionHigh*100;
  std::string strFracPercentHigh = std::to_string(fracPercentHigh);

  std::string histfilename=(outfile+"_FractionLow_"+strFracPercentLow+"_FractionHigh_"+strFracPercentHigh+"_DoubleSided.root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_nHits_afterCut->Write();
  h_nHits_beforeCut->Write();
  h_rho->Write();
  h_meanXposition->Write();
  h_meanYposition->Write();
  h_meanZposition->Write();
  h_meanXposition_multiclus->Write();
  h_meanYposition_multiclus->Write();
  h_meanZposition_multiclus->Write();
  h_nMult_AfterCuts->Write();
  h_nMult->Write();
  h_res_Xposition->Write();
  h_res_Yposition->Write();
  h_PositionAverage_X->Write();
  h_PositionAverage_Y->Write();
  h_rechitX->Write();
  h_rechitY->Write();
  h_rechitZ->Write();
  h_rechitX_BeforeTruncation->Write();
  h_rechit_rho_energy->Write();
  h_multiclus_clus2D_energy->Write();
  h_reco_gen_X->Write();
  h_reco_gen_Y->Write();
  h_reco_gen_X_perhit->Write();
  h_reco_gen_Y_perhit->Write();
  h_reco_gen_X_perhit_beforeCuts->Write();
  h_reco_gen_Y_perhit_beforeCuts->Write();
  h_gen_X_Y->Write();
  h_reco_X_Y->Write();
  tFile->Close();
  outputFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0; 
}
