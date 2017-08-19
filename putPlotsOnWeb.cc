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
#include <fstream>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>

using std::string;

void readTextFiles()
{
  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 1200, 800);

  int perCentTruncationLow[9] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
  int perCentTruncationHigh[9] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
  std::string line;
  TH2D *h_Variance = new TH2D("h_Variance", "h_Variance", 10, 0.0, 100.0, 10, 0.0, 100.0);h_Variance->Sumw2();

  for(unsigned int i=0; i<sizeof(perCentTruncationLow)/sizeof(int); i++)
  {
    for(unsigned int j=0; j<sizeof(perCentTruncationHigh)/sizeof(int); j++)
    {
      if(perCentTruncationLow[i]+perCentTruncationHigh[j] > 90) continue;
      std::string strFracPercentLow = std::to_string(perCentTruncationLow[i]);
      std::string strFracPercentHigh = std::to_string(perCentTruncationHigh[j]);
      ifstream inFile(("OutputFiles/SoverN5000_Floor30ps_DoubleSided/output_Fit_HGCTiming_kLong_Pt10_SoverN5000ps_Floor30ps_FractionLow_"+strFracPercentLow+"_FractionHigh_"+strFracPercentHigh+".txt").c_str());
      if (inFile.is_open())
      {
        while (!inFile.eof())
        {
          getline(inFile, line);
          if(line.find("Variance = ")!= string::npos) 
          {
            //std::string::size_type sz;
            //std::cout << sz << std::endl;
            //double variance = 0.0;
            //double variance = std::stod (line.substr(sz));
            //std::cout << variance << std::endl;
            std::string number = line.substr(11);
            std::cout << i << " " << j  << " " << perCentTruncationLow[i] << " " << perCentTruncationHigh[j] << " " << number << std::endl;
            double variance = std::stod(number);
            h_Variance->SetBinContent(i+1, j+1, variance);
          }
        }
      }
    }
  }
  h_Variance->GetXaxis()->SetTitle("Fraction Low (%)");
  h_Variance->GetYaxis()->SetTitle("Fraction High (%)");
  h_Variance->SetStats(kFALSE);
  h_Variance->SetTitle("");
  
  h_Variance->Draw("COLZ TEXT");
  gStyle->SetPaintTextFormat("4.4f");

  TLatex TL;
  TL.SetTextAlign(11);
  TL.SetTextSize(0.05);
  TL.SetTextFont(22);
  TL.DrawLatexNDC(0.1,0.94,"Variance: K^{0}_{L} with p_{T} = 10 GeV and eta = 1.8");

  c1.SaveAs("OutputFiles/SoverN5000_Floor30ps_DoubleSided/h_Variance.png");
  c1.SaveAs("OutputFiles/SoverN5000_Floor30ps_DoubleSided/h_Variance.pdf");
}


void putPlotsOnWeb()
{
  int PT = 10;
  int perCentTruncationLow[9] = {10, 20, 30, 40, 50, 60, 70, 80, 90};
  int perCentTruncationHigh[9] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

  std::ofstream outfile("OutputFiles/SoverN5000_Floor30ps_DoubleSided/index.html");
  outfile<<"<html>"<<std::endl;
  outfile<<"<head>"<<std::endl;
  outfile<<"</head>"<<std::endl;
  outfile<<"<body>"<<std::endl;
  outfile<<"<h1 align='center'> Bias Variance Study kLong pT 10 GeV</h1>"<<std::endl;
  outfile<<"<table border='1'>"<<std::endl;

  // Write to HTML file
  outfile<<"<tr>"<<std::endl;
  outfile<<" <td>"<<std::endl;
  for(unsigned int i=0; i<sizeof(perCentTruncationLow)/sizeof(int); i++)
  {
    for(unsigned int j=0; j<sizeof(perCentTruncationHigh)/sizeof(int); j++)
    {
      if(perCentTruncationLow[i]+perCentTruncationHigh[j] > 90) continue;
      //std::cout << "perCentTruncationLow[" << i << "] = " << perCentTruncationLow[i] << std::endl;
      //std::cout << "perCentTruncationHigh[" << j << "] = " << perCentTruncationHigh[j] << std::endl;
      outfile<<"  <img src=\"output_HGCTiming_kLong_Pt"<<PT<<"_SoverN5000ps_Floor30ps_FractionLow_"<<perCentTruncationLow[i]<<"_FractionHigh_"<<perCentTruncationHigh[j]<<"_h_TimeAverage"<<".png\"/>" << std::endl;
    }
  }
  outfile<<"<br/>"<<std::endl;
  outfile<<"  <img src=\"h_Variance.png\"/>" << std::endl; 
  outfile<<" </td>"<<std::endl;

  outfile<<"</body>"<<std::endl;
  outfile<<"</html>"<<std::endl;


}
 
