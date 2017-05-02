// user include files
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include "TRandom3.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "RecoNtuples/HGCalAnalysis/plugins/RecHiTimeEstimator.h"

RecHiTimeEstimator::RecHiTimeEstimator(const edm::ParameterSet& ps){
  //  std::cout << " >>> constructor " << std::endl;
  //section type
  keV2fC[0] =  ps.getParameter<double>("HGCEE_keV2fC");
  keV2fC[1] =  ps.getParameter<double>("HGCHEF_keV2fC");
  keV2MIP = ps.getParameter<double>("HGCHB_keV2MIP");

  noisefC[0] = (ps.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
  noisefC[1] = (ps.getParameter<std::vector<double> >("HGCEF_noisefC")).at(1);
  noiseMIP = ps.getParameter<double>("HGCBH_noiseMIP");

  //cell type
  fCPerMIP[0] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(0);
  fCPerMIP[1] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(1);
  fCPerMIP[2] =  (ps.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(2);

  const auto& rcorr = ps.getParameter<std::vector<double> >("thicknessCorrection");
  scaleCorrection.clear();
  scaleCorrection.push_back(1.f);
  for( auto corr : rcorr ) {
    scaleCorrection.push_back(1.0/corr);
  }

  const auto& dweights = ps.getParameter<std::vector<double> >("dEdXweights");
  for( auto weight : dweights ) {
    weights.push_back(weight);
  }

  keV2GeV = 1e-6;
  keV2MeV = 1e-3;

  //100um
  paramA[0] = 1.00;
  parErrA[0] = 0.01;
  paramC[0] = 0.009;
  parErrC[0] = 0.001;

  //200um
  paramA[1] = 1.06;
  parErrA[1] = 0.02;
  paramC[1] = 0.008;
  parErrC[1] = 0.001;

  //300um
  paramA[2] = 1.11;
  parErrA[2] = 0.02;
  paramC[2] = 0.010;
  parErrC[2] = 0.001;

  timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);
}


void RecHiTimeEstimator::setEventSetup(const edm::EventSetup& es){
  es.get<CaloGeometryRecord>().get(pG);
}


const HGCalDDDConstants* RecHiTimeEstimator::get_ddd(const CaloSubdetectorGeometry* geom){
  const HGCalGeometry* hg = static_cast<const HGCalGeometry*>(geom);
  const HGCalDDDConstants* ddd = &(hg->topology().dddConstants());
  if( nullptr == ddd ) {
    throw cms::Exception("hgcal::RecHitTools")
      << "DDDConstants not accessibl to hgcal::RecHitTools!";
  }
  return ddd;
}


double RecHiTimeEstimator::getTimeHit(int thick, double SoverN){
  //  std::cout << " >>> getTimeHit  " << std::endl;
  timeResolution->SetParameters(paramA[thick], paramC[thick]);

  double sigma = 0.2;
  if(SoverN > 1000) return paramC[thick];
  else sigma = timeResolution->Eval(SoverN);

  TRandom3* rand = new TRandom3(); 
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}


/*
double RecHiTimeEstimator::getTimeHitFixThr(){
  //  std::cout << " >>> getTimeHit  " << std::endl;

  double sigma = 0.05;

  TRandom3* rand = new TRandom3(); 
  double smearing = rand->Gaus(0., sigma);
  return smearing;
}
*/

void RecHiTimeEstimator::correctTime(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){
  //  std::cout << " >>> correctTime  " << std::endl;
  const CaloGeometry* caloGeom = pG.product();

  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const HGCalDetId detid = it_hit->detid();
    auto ddd = get_ddd(caloGeom->getSubdetectorGeometry(detid));
    int thick = ddd->waferTypeL(detid.wafer()) - 1;
    unsigned int layer = detid.layer();

    int sectionType = 2;
    if(layer < 29) sectionType = 0;
    else if(layer < 41) sectionType = 1;

    double sigmaNoiseMIP = noiseMIP;
    if(sectionType != 2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];

    // double keVtoMip = keV2fC[sectionType]/fCPerMIP[thick];
    // if(sectionType == 2) keVtoMip = keV2MIP;

    // double sigmaNoise = sigmaNoiseMIP / keV2MIP * (weights.at(layer-1)/keV2MeV*keV2MIP) * keV2GeV * scaleCorrection[thick];
    // if(sectionType != 2) sigmaNoise = sigmaNoiseMIP * (weights.at(layer-1)/keV2MeV*keVtoMip) / keVtoMip * keV2GeV * scaleCorrection[thick];


    HGCRecHit myrechit(*it_hit);
    float energy = myrechit.energy();

    int energyMIP = 0.;
    if(sectionType == 2) energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer-1)/keV2MeV );
    else energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer-1)/keV2MeV);


    if(energyMIP > 3. && myrechit.time() >= 0.){
      float SoverN = energyMIP / sigmaNoiseMIP;
      double smearedTime = getTimeHit(thick, SoverN);
      //      std::cout << " smeared time  = " << smearedTime << " SoverN = " << SoverN << " thick = " << thick << std::endl;
      //cross-check with Sapta for selection on original time
      myrechit.setTime((myrechit.time() - 1.) * (1+smearedTime));
    }
    else myrechit.setTime(-1.); 
    Newrechits->push_back(myrechit);
  }
  return;
}

/*
void RecHiTimeEstimator::correctTimeFixThr(const HGCRecHitCollection& rechits, HGCRecHitCollection* Newrechits){
  //  std::cout << " >>> correctTime  " << std::endl;
  for(HGCRecHitCollection::const_iterator it_hit = rechits.begin(); it_hit < rechits.end(); ++it_hit) {
    const HGCalDetId detid = it_hit->detid();
    unsigned int layer = detid.layer();

    int sectionType = 2;
    if(layer < 29) sectionType = 0;
    else if(layer < 41) sectionType = 1;

    HGCRecHit myrechit(*it_hit);
    float energy = myrechit.energy();

    float energyCharge = 0.;
    // from SimCalorimetry/HGCalSimProducers/src/HGCHEbackDigitizer.cc L 58
    if(sectionType == 2) energyCharge = energy * 1.e6 * 1.; 
    else energyCharge = energy * 1.e6 * keV2fC[sectionType];

    if(energyCharge > 60 && myrechit.time() >= 0.){
      double smearedTime = getTimeHitFixThr();
      myrechit.setTime((myrechit.time() - 1.) * (1+smearedTime));
    }
    else myrechit.setTime(-1.); 
    Newrechits->push_back(myrechit);
  }
  return;
}*/
