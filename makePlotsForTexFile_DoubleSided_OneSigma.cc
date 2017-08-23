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

int makePlotsForTexFile_DoubleSided_OneSigma(std::string infile, float fractionLow, float fractionHigh)
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

  TH1F *h1 = (TH1F*)file1->Get("h_TimeAverage");
  h1->Rebin(1);
  h1->GetXaxis()->SetRangeUser(-0.1, 0.4); 
  double mean = h1->GetMean();
  std::cout << mean << std::endl;
  double meanBin = h1->FindBin(mean);
  int binLeft = -1;
  int binRight = -1;
  double OneSigma = 0.6827*100;
  std::cout << "h1->GetNbinsX() = " << h1->GetNbinsX() << std::endl;
  std::cout << (100-OneSigma)/2.0 << std::endl;
  std::cout << (100+OneSigma)/2.0 << std::endl;

  for(int ibin=1; ibin < h1->GetNbinsX()+1; ibin++)
  {
    double fractionalArea = (double)(h1->Integral(1, ibin))*100.0/(double)h1->Integral();
    if(fractionalArea <= 15.86)
    {
      binLeft = ibin;
    }
    if(fractionalArea <= 84.14)
    {
      binRight = ibin;
    }
  }
  std::cout << "binLeft = " << binLeft << std::endl;
  std::cout << "binRight = " << binRight << std::endl;
  std::cout << "Integral(binLeft, binRight) = " << h1->Integral(binLeft, binRight) << std::endl;
  std::cout << "Integral() = " << h1->Integral() << std::endl;
  double effectiveRMS = fabs(h1->GetBinCenter(binLeft)-h1->GetBinCenter(binRight))/2.0;
  std::cout << "One Sigma = " << effectiveRMS << std::endl;

  TH1F *hfit = (TH1F*)file1->Get("h_TimeAverage");
  hfit->Rebin(5);
  hfit->GetXaxis()->SetRangeUser(-0.1, 0.4); //for k-long
  //hfit->GetXaxis()->SetRangeUser(-0.1, 0.2);//for gamma
  hfit->SetTitle("");
  hfit->Draw();

  int fracPercentLow = fractionLow;
  std::string strFracPercentLow = std::to_string(fracPercentLow);

  int fracPercentHigh = fractionHigh;
  std::string strFracPercentHigh = std::to_string(fracPercentHigh);

  TLatex TL;
  TL.SetTextAlign(11);
  TL.SetTextSize(0.05);
  TL.SetTextFont(22);
  TL.DrawLatexNDC(0.1,0.94,("Fraction Low: "+strFracPercentLow+"% Fraction High: "+strFracPercentHigh+"%").c_str());

  c1.SaveAs((infile+"_h_TimeAverage.png").c_str());
  c1.SaveAs((infile+"_h_TimeAverage.pdf").c_str());

  return 0;
} 

