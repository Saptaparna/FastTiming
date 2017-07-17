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

bool sameVal(double a, double b)
{
   return fabs(a - b) < 2.000e-00;
}

bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}

typedef struct
{
  int leftEdge;
  int rightEdge;
  double binWidth;
} binInfo;

bool sortBinsInBinWidth (binInfo bin1, binInfo bin2) 
{ 
  return (bin1.binWidth < bin2.binWidth); 
}

void computeEffectiveRMS(std::string infile)
{
  std::string inputfilename=(infile+".root").c_str();
  TFile *rootfile = TFile::Open((infile+".root").c_str());
  TH1F *h1 = (TH1F*) rootfile->Get("h_recHit_Time");
  h1->GetXaxis()->SetRangeUser(-2.0, 5.0);
  h1->Rebin(1);
  double mean = h1->GetMean();
  double meanBin = h1->FindBin(mean);
  int binLeft = -1;
  int binRight = -1;  
  double OneSigma = 0.6827*100;
  std::vector<binInfo> bins; 
 
  for(int ibin=1; ibin < h1->GetNbinsX()+1-1; ibin++)
  {
    for(int jbin=ibin+1; jbin < h1->GetNbinsX()+1; jbin++)
    {
      double fractionalArea = (double)(h1->Integral(ibin, jbin))*100.0/(double)h1->Integral();  
      if(sameVal(fractionalArea, OneSigma))
      {
        binInfo bin;
        bin.leftEdge=ibin;
        bin.rightEdge=jbin;
        bin.binWidth=jbin-ibin;
        bins.push_back(bin);
      }
    }
  }
  std::sort(bins.begin(), bins.end(), sortBinsInBinWidth);
  std::cout << "bins.at(0).leftEdge = " << bins.at(0).leftEdge << std::endl;
  std::cout << "bins.at(0).rightEdge = " << bins.at(0).rightEdge << std::endl;
  std::cout << "bins.at(0).binWidth = " << bins.at(0).binWidth << std::endl;
}
