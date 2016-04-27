#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include <TCanvas.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TLegend.h"
#include <TPad.h>
#include <sstream>
#include "TVectorD.h"
#include "TGraph.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

void plotVertexZ()
{
  gROOT->SetStyle("Plain");
  TCanvas c1("c1","Z coordinate of Vertex", 10, 10, 600, 400);
  
  TFile* file1 = TFile::Open("output_HGCTiming_SinglePhoton_Pt35GeV.root");
  TH1F *h1 = (TH1F*)file1->Get("h_vertex_z");
  h1->Rebin(5);
  h1->SetLineColor(kBlack);
  h1->SetLineWidth(2);
  h1->SetStats(kFALSE);
  h1->Draw("hist");

  TFile* file2 = TFile::Open("output_HGCTiming_SinglePhoton_Pt50GeV.root");
  TH1F *h2 = (TH1F*)file2->Get("h_vertex_z");
  h2->Rebin(5);
  h2->SetLineColor(kRed);
  h2->SetLineWidth(2);
  h2->SetStats(kFALSE);
  h2->Draw("hist same");

  TFile* file3 = TFile::Open("output_HGCTiming_SinglePhoton_Pt65GeV.root");
  TH1F *h3 = (TH1F*)file3->Get("h_vertex_z");
  h3->Rebin(5);
  h3->SetLineColor(kBlue);
  h3->SetLineWidth(2);
  h3->SetStats(kFALSE);
  h3->Draw("hist same");

  TFile* file4 = TFile::Open("output_HGCTiming_SinglePhoton_Pt80GeV.root");
  TH1F *h4 = (TH1F*)file4->Get("h_vertex_z");
  h4->Rebin(5);
  h4->SetLineColor(kGreen);
  h4->SetLineWidth(2);
  h4->SetStats(kFALSE);
  h4->Draw("hist same");
  
  TLegend *leg1 = new TLegend(0.5907383,0.6609651,0.8959732,0.8579088,NULL,"brNDC");
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.03);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(0);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->AddEntry(h1,"Single Photon p_{T} 35 GeV","l");
  leg1->AddEntry(h2,"Single Photon p_{T} 50 GeV","l");
  leg1->AddEntry(h3,"Single Photon p_{T} 65 GeV","l");
  leg1->AddEntry(h4,"Single Photon p_{T} 80 GeV","l");
  leg1->Draw();

  c1.SaveAs("Vertex_z_SinglePhoton.pdf");
  c1.SaveAs("Vertex_z_SinglePhoton.png");
}
