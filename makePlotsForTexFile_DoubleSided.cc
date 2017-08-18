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
#include "Riostream.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TLatex.h"

int makePlotsForTexFile_DoubleSided(std::string infile, float fractionLow, float fractionHigh)
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);

  std::string inputfilename=(infile+"_DoubleSided.root").c_str();
  TFile* file1 = TFile::Open(inputfilename.c_str());

  TH1F *h4 = (TH1F*)file1->Get("h_nHits_afterCut");
  h4->Rebin(2);
  h4->GetXaxis()->SetRangeUser(0, 300);
  h4->SetLineColor(kBlue);
  h4->SetLineWidth(2);
  h4->GetXaxis()->SetTitle("Number of hits");
  h4->GetYaxis()->SetTitle("Events");
  h4->GetYaxis()->SetTitleOffset(1.3);
  h4->SetTitle("");
  h4->Draw("HIST");

  c1.SaveAs((infile+"_h_nHits_afterCut.png").c_str());
  c1.SaveAs((infile+"_h_nHits_afterCut.pdf").c_str());

  gStyle->SetOptFit(1);
  TH1F *hfit = (TH1F*)file1->Get("h_TimeAverage");
  hfit->Rebin(5);
  hfit->GetXaxis()->SetRangeUser(-0.1, 0.4); 
  hfit->SetTitle("");
  TF1 *fa1 = new TF1("fa1", "gaus", -0.01, 0.01);
  fa1->SetLineColor(kRed);
  fa1->SetLineWidth(2);
  hfit->Fit("fa1", "R+");
  hfit->Fit("fa1");
  hfit->Draw();

  std::cout << "Constant = " << fa1->GetParameter(0) << std::endl;//this is the constant
  std::cout << "Bias = " << fa1->GetParameter(1) << std::endl;//this is the mean
  std::cout << "Variance = " << fa1->GetParameter(2) << std::endl;//this is the sigma

  int fracPercentLow = fractionLow*100;
  std::string strFracPercentLow = std::to_string(fracPercentLow);
  
  int fracPercentHigh = fractionHigh*100;
  std::string strFracPercentHigh = std::to_string(fracPercentHigh);

  TLatex *texf = new TLatex(0.10,1000.0,("Fraction Low_"+strFracPercentLow+" Fraction High_"+strFracPercentHigh).c_str());
  texf->SetTextSize(0.05);
  texf->Draw();
    
  c1.SaveAs((infile+"_h_TimeAverage.png").c_str());
  c1.SaveAs((infile+"_h_TimeAverage.pdf").c_str());
  c1.SaveAs("test.C");
  return 0;
} 

