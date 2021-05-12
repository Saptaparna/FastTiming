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

int makeComparison()
{

  gROOT->SetStyle("Plain");
  TCanvas c1("c1","PlotsForTexFile", 10, 10, 600, 400);
  gStyle->SetOptFit(1);
  //c1.SetLogx();

  TFile* file_new = TFile::Open("output_hgcalNtuple_April25_20k_May5.root");
  TH1F *hEnergyTotalRadiusTwo = (TH1F*)file_new->Get("h_EnergyTotalRadiusTwo");
  hEnergyTotalRadiusTwo->Rebin(2);
  hEnergyTotalRadiusTwo->GetXaxis()->SetRangeUser(0.0, 70.0);
  hEnergyTotalRadiusTwo->SetTitle("");
  double max_new = hEnergyTotalRadiusTwo->GetMaximum();
  hEnergyTotalRadiusTwo->SetLineColor(kBlue);
  hEnergyTotalRadiusTwo->SetLineWidth(2);
  hEnergyTotalRadiusTwo->SetFillColor(kBlue);  
  hEnergyTotalRadiusTwo->SetFillStyle(3004);
  hEnergyTotalRadiusTwo->SetMaximum(1.2*max_new);
  hEnergyTotalRadiusTwo->GetXaxis()->SetTitle("Total energy in a cynlinder with timing requirement [GeV]");
  hEnergyTotalRadiusTwo->Draw("HIST SAMES");
  c1.Modified();
  c1.Update();
  TPaveStats *stats2 = (TPaveStats*)hEnergyTotalRadiusTwo->GetListOfFunctions()->FindObject("stats");
  stats2->SetName("Nominal");
  stats2->SetY1NDC(0.7);
  stats2->SetY2NDC(0.9);
  stats2->SetTextColor(kBlue);
  hEnergyTotalRadiusTwo->Draw("HIST SAMES");
  c1.Update();  
  
  TH1F *hEnergyTotalRadiusFive = (TH1F*)file_new->Get("h_EnergyTotalRadiusFive");
  hEnergyTotalRadiusFive->Rebin(2);
  hEnergyTotalRadiusFive->GetXaxis()->SetRangeUser(0.0, 100.0);
  hEnergyTotalRadiusFive->SetTitle("");
  double max_old = hEnergyTotalRadiusFive->GetMaximum();
  hEnergyTotalRadiusFive->SetLineColor(kRed);
  hEnergyTotalRadiusFive->SetLineWidth(2);
  hEnergyTotalRadiusFive->SetFillColor(kRed);
  hEnergyTotalRadiusFive->SetFillStyle(3004);
  hEnergyTotalRadiusFive->Draw("HIST SAMES");
  c1.Draw();
  c1.Modified();
  c1.Update();
  TPaveStats *stats = (TPaveStats*)hEnergyTotalRadiusFive->GetListOfFunctions()->FindObject("stats");
  stats->SetName("Displaced");
  stats->SetY1NDC(0.7);
  stats->SetY2NDC(0.5);
  stats->SetTextColor(kRed);
  hEnergyTotalRadiusFive->Draw("HIST SAMES");
  c1.Update();
  
  TH1F *hEnergyTotalRadiusTen = (TH1F*)file_new->Get("h_EnergyTotalRadiusTen");
  hEnergyTotalRadiusTen->Rebin(2);
  hEnergyTotalRadiusTen->GetXaxis()->SetRangeUser(0.0, 100.0);
  hEnergyTotalRadiusTen->SetTitle("");
  hEnergyTotalRadiusTen->SetLineColor(kOrange-7);
  hEnergyTotalRadiusTen->SetLineWidth(2);
  hEnergyTotalRadiusTen->SetFillColor(kOrange-7);
  hEnergyTotalRadiusTen->SetFillStyle(3004);
  hEnergyTotalRadiusTen->Draw("HIST SAMES");
  c1.Draw();
  c1.Modified();
  c1.Update();
  TPaveStats *stats3 = (TPaveStats*)hEnergyTotalRadiusTen->GetListOfFunctions()->FindObject("stats");
  stats3->SetName("Displaced");
  stats3->SetY1NDC(0.5);
  stats3->SetY2NDC(0.3);
  stats3->SetTextColor(kOrange-7);
  hEnergyTotalRadiusTen->Draw("HIST SAMES");
  c1.Update();

  TH1F *hEnergyTotalRadiusFifty = (TH1F*)file_new->Get("h_EnergyTotalRadiusFifty");
  hEnergyTotalRadiusFifty->Rebin(2);
  hEnergyTotalRadiusFifty->GetXaxis()->SetRangeUser(0.0, 100.0);
  hEnergyTotalRadiusFifty->SetTitle("");
  hEnergyTotalRadiusFifty->SetLineColor(kGreen+3);
  hEnergyTotalRadiusFifty->SetLineWidth(2);
  hEnergyTotalRadiusFifty->SetFillColor(kGreen+3);
  hEnergyTotalRadiusFifty->SetFillStyle(3004);
  hEnergyTotalRadiusFifty->Draw("HIST SAMES");
  c1.Draw();
  c1.Modified();
  c1.Update();
  TPaveStats *stats4 = (TPaveStats*)hEnergyTotalRadiusFifty->GetListOfFunctions()->FindObject("stats");
  stats4->SetName("Displaced");
  stats4->SetY1NDC(0.3);
  stats4->SetY2NDC(0.1);
  stats4->SetTextColor(kGreen+3);
  hEnergyTotalRadiusFifty->Draw("HIST SAMES");
  c1.Update();
 
  TLatex TL;
  TL.SetTextAlign(11);
  TL.SetTextSize(0.05);
  TL.SetTextFont(22);
  //TL.DrawLatexNDC(0.1,0.94,("Fraction Low: 10% Fraction High: 10 %"));
  
  TLatex TL_new;
  TL_new.SetTextAlign(11);
  TL_new.SetTextSize(0.03);
  TL_new.SetTextFont(22);
  TL_new.SetTextColor(2);
  int hfitNewIntegral = hEnergyTotalRadiusTwo->Integral();
  std::string strhfitNewIntegral = std::to_string(hfitNewIntegral);
  //TL_new.DrawLatexNDC(0.2,0.30,("Integral: "+strhfitNewIntegral).c_str());
  
  TLatex TL_old;
  TL_old.SetTextAlign(11);
  TL_old.SetTextSize(0.03);
  TL_old.SetTextFont(22);
  TL_old.SetTextColor(4);
  int hfitOldIntegral = hEnergyTotalRadiusFive->Integral();
  std::string strhfitOldIntegral = std::to_string(hfitOldIntegral);
  //TL_old.DrawLatexNDC(0.2,0.20,("Integral: "+strhfitOldIntegral).c_str());
 
  TLegend *leg = new TLegend(0.10,0.65,0.30,0.80,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetLineColor(1);
  leg->SetLineStyle(0);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(0);
  leg->AddEntry(hEnergyTotalRadiusTwo,"Radius < 2 cm","l");
  leg->AddEntry(hEnergyTotalRadiusFive,"Radius < 5 cm","l");
  leg->AddEntry(hEnergyTotalRadiusTen,"Radius < 10 cm","l");
  leg->AddEntry(hEnergyTotalRadiusFifty,"Radius < 50 cm","l");
  leg->Draw();

  c1.SaveAs("h_ShowerAxis_Energy_TimeCut.png");
  c1.SaveAs("h_ShowerAxis_Energy_TimeCut.pdf");
  return 0;

}
