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

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double E)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiE(pT, eta, phi, E);
  return object_p4;
}

typedef struct
{
  float recHitE;
  float recHitX;
  float recHitY;
  float recHitZ;
} RecHitInfo;

bool sortRecHitsInDescendingE(RecHitInfo recHit1, RecHitInfo recHit2)
{
  return (recHit1.recHitE > recHit2.recHitE);
}

int ReadHGCTiming_Tree(std::string infile, std::string outfile)
{
  std::string inputfilename=(infile+".root").c_str();
  TChain *tree=new TChain("hgctiming/HGCTiming");
  tree->Add(inputfilename.c_str());
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  Int_t           run;
  Int_t           event;
  Int_t           lumi;
  Float_t         vertex_x;
  Float_t         vertex_y;
  Float_t         vertex_z;
  vector<float>   *recHit_energy;
  vector<float>   *recHit_x;
  vector<float>   *recHit_y;
  vector<float>   *recHit_z;
  TLorentzVector  *GenVertex;

  recHit_energy = 0;
  recHit_x = 0;
  recHit_y = 0;
  recHit_z = 0;
  GenVertex = 0;

  tree->SetBranchAddress("run", &(run));
  tree->SetBranchAddress("lumi", &(lumi));
  tree->SetBranchAddress("event", &(event));
  tree->SetBranchAddress("vertex_x", &(vertex_x));
  tree->SetBranchAddress("vertex_y", &(vertex_y));
  tree->SetBranchAddress("vertex_z", &(vertex_z));
  tree->SetBranchAddress("recHit_energy", &(recHit_energy));
  tree->SetBranchAddress("recHit_x", &(recHit_x));
  tree->SetBranchAddress("recHit_y", &(recHit_y));
  tree->SetBranchAddress("recHit_z", &(recHit_z));
  tree->SetBranchAddress("GenVertex", &(GenVertex));

  TH1F *h_vertex_x=new TH1F("h_vertex_x", "h_vertex_x; Vertex [mm]; Events", 100, -50.0, 50.0); h_vertex_x->Sumw2();
  TH1F *h_vertex_y=new TH1F("h_vertex_y", "h_vertex_y; Vertex [mm]; Events", 100, -50.0, 50.0); h_vertex_y->Sumw2();
  TH1F *h_vertex_z=new TH1F("h_vertex_z", "h_vertex_z; Vertex [mm]; Events", 100, -50.0, 50.0); h_vertex_z->Sumw2();
  TH1F *h_recHit_HighestEnergy = new TH1F("h_recHit_HighestEnergy", "h_recHit_HighestEnergy; Energy [keV]; Events/keV", 1000, 0, 10000); h_recHit_HighestEnergy->Sumw2();
  TH1F *h_recHit_X = new TH1F("h_recHit_X", "h_recHit_X; recHit x coordinate; Events", 200, -100, 100); h_recHit_X->Sumw2();
  TH1F *h_recHit_Y = new TH1F("h_recHit_Y", "h_recHit_Y; recHit y coordinate; Events", 200, -100, 100); h_recHit_Y->Sumw2();
  TH1F *h_recHit_Z = new TH1F("h_recHit_Z", "h_recHit_Z; recHit z coordinate; Events",1000, -500, 500); h_recHit_Z->Sumw2();
  TH1F *h_recHit_Phi = new TH1F("h_recHit_Phi", "h_recHit_Phi; recHit #phi; Events", 80,  -4,   4); h_recHit_Phi->Sumw2();
  TH1F *h_recHit_Eta = new TH1F("h_recHit_Eta", "h_recHit_Eta; recHit #eta; Events", 120,  -6,  6); h_recHit_Eta->Sumw2();

  int nEvents=tree->GetEntries();
  std::cout << "nEvents= " << nEvents << std::endl;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    //filling vertex info
    h_vertex_x->Fill(vertex_x);
    h_vertex_y->Fill(vertex_y);
    h_vertex_z->Fill(vertex_z);
    //filling the recHit info in a vector of structs
    std::vector<RecHitInfo> recHits;
    for (unsigned int j=0; j<recHit_energy->size(); j++)
    {
      RecHitInfo recHit;
      recHit.recHitE = recHit_energy->at(j)/(55.1);
      recHit.recHitX = recHit_x->at(j);
      recHit.recHitY = recHit_y->at(j);
      recHit.recHitZ = recHit_z->at(j);  
      recHits.push_back(recHit);
    }
    std::sort(recHits.begin(), recHits.end(), sortRecHitsInDescendingE);
    TVector3 recHitV;
    recHitV.SetXYZ(-9999.0, -9999.0, -9999.0);

    if(recHits.size() > 0) 
    {
      h_recHit_HighestEnergy->Fill(recHits.at(0).recHitE);
      h_recHit_X->Fill(recHits.at(0).recHitX);
      h_recHit_Y->Fill(recHits.at(0).recHitY);
      h_recHit_Z->Fill(recHits.at(0).recHitZ);
      recHitV.SetXYZ(recHits.at(0).recHitX, recHits.at(0).recHitY, recHits.at(0).recHitZ);
      h_recHit_Phi->Fill(recHitV.Phi());
      h_recHit_Eta->Fill(recHitV.Eta());
    }

  }//end of event loop

  std::string histfilename=(outfile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_vertex_x->Write();
  h_vertex_y->Write();
  h_vertex_z->Write();
  h_recHit_HighestEnergy->Write();
  h_recHit_X->Write();
  h_recHit_Y->Write();
  h_recHit_Z->Write();
  h_recHit_Phi->Write();
  h_recHit_Eta->Write();
  tFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;

  return 0; 
}
