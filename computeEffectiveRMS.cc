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
   return fabs(a - b) < 1.000e-01;
}

bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}


double computeEffectiveRMS(std::string infile)
{
  std::string inputfilename=(infile+".root").c_str();
  TFile *rootfile = TFile::Open((infile+".root").c_str());
  TH1F *h1 = (TH1F*) rootfile->Get("h_TimeAverage");
  h1->Rebin(1);
  double mean = h1->GetMean();
  double meanBin = h1->FindBin(mean);
  int binLeft = -1;
  int binRight = -1;  
  double OneSigma = 0.6827*100;

  for(int ibin=1; ibin < h1->GetNbinsX()+1; ibin++)
  {
    int binLeftTemp = meanBin - ibin;
    int binRightTemp = meanBin + ibin;
    double fractionalArea = (double)(h1->Integral(binLeftTemp, binRightTemp))*100.0/(double)h1->Integral();  
    if(sameVal(fractionalArea, OneSigma))
    {
      binLeft = binLeftTemp;
      binRight = binRightTemp;
    }
  }

  double effectiveRMS = fabs(h1->GetBinCenter(binLeft)-h1->GetBinCenter(binRight));
  return effectiveRMS;

}
