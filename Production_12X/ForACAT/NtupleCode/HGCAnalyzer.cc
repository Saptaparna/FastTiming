// -*- C++ -*-
//
// Package:    HGCAnalyzer/HGCAnalyzer
// Class:      HGCAnalyzer
//
/**\class HGCAnalyzer HGCAnalyzer.cc HGCAnalyzer/HGCAnalyzer/plugins/HGCAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saptaparna Bhattacharya
//         Created:  Tue, 16 Nov 2021 15:27:10 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TTree.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
//#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
//#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
//#include "Geometry/HGCalCommonData/interface/HGCalGeometryMode.h"
//#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
//#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

//#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
//#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class HGCAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit HGCAnalyzer(const edm::ParameterSet&);
  ~HGCAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillRecHit(const DetId &detid, const float &fraction, const int &cluster_index_ = -1);
  // ----parameters------
  bool rawRecHits_;

  // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> recHitsEE_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsFH_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitsBH_;

  std::vector<float> rechit_energy_;
  std::vector<float> rechit_x_;
  std::vector<float> rechit_y_;
  std::vector<float> rechit_z_;
  std::vector<float> rechit_time_;
  std::vector<unsigned int> rechit_detid_;

  TTree *t_;
  std::set<DetId> storedRecHits_;
  std::map<DetId, const HGCRecHit *> hitmap_;
  unsigned int rechit_index_;
  std::map<DetId, unsigned int> detIdToRecHitIndexMap_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HGCAnalyzer::HGCAnalyzer(const edm::ParameterSet& iConfig)
    : rawRecHits_(iConfig.getParameter<bool>("rawRecHits"))
{
    recHitsEE_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCEERecHits"));
    recHitsFH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEFRecHits"));
    recHitsBH_ = consumes<HGCRecHitCollection>(edm::InputTag("HGCalRecHit", "HGCHEBRecHits")); 
}

HGCAnalyzer::~HGCAnalyzer() {
}


//
// member functions
//

// ------------ method called for each event  ------------
void HGCAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  using namespace edm;

  rechit_energy_.clear();
  rechit_x_.clear();
  rechit_y_.clear();
  rechit_z_.clear();
  rechit_time_.clear();
  rechit_detid_.clear();
  storedRecHits_.clear();
  hitmap_.clear();

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  iEvent.getByToken(recHitsEE_, recHitHandleEE);
  iEvent.getByToken(recHitsFH_, recHitHandleFH);
  iEvent.getByToken(recHitsBH_, recHitHandleBH);
  const auto &rechitsEE = *recHitHandleEE;
  const auto &rechitsFH = *recHitHandleFH;
  const auto &rechitsBH = *recHitHandleBH;
  for (unsigned int i = 0; i < rechitsEE.size(); ++i) 
  {
    hitmap_[rechitsEE[i].detid()] = &rechitsEE[i];
  }
  for (unsigned int i = 0; i < rechitsFH.size(); ++i) 
  {
    hitmap_[rechitsFH[i].detid()] = &rechitsFH[i];
  }
  for (unsigned int i = 0; i < rechitsBH.size(); ++i) 
  {
    hitmap_[rechitsBH[i].detid()] = &rechitsBH[i];
  }

  if (rawRecHits_)
  {
    const HGCRecHitCollection &rechitsEE = *recHitHandleEE;
    const HGCRecHitCollection &rechitsFH = *recHitHandleFH;
    const HGCRecHitCollection &rechitsBH = *recHitHandleBH; 

    for (HGCRecHitCollection::const_iterator it_hit = rechitsEE.begin(); it_hit < rechitsEE.end(); ++it_hit) 
    {
      const HGCalDetId detid = it_hit->detid();
      if (storedRecHits_.find(detid) == storedRecHits_.end()) 
      {
        fillRecHit(detid, -1);
      }
    }
    for (HGCRecHitCollection::const_iterator it_hit = rechitsFH.begin(); it_hit < rechitsFH.end();  ++it_hit) 
    {
      const HGCalDetId detid = it_hit->detid();
      if (storedRecHits_.find(detid) == storedRecHits_.end()) 
      {
        fillRecHit(detid, -1);
      }
    }
    for (HGCRecHitCollection::const_iterator it_hit = rechitsBH.begin(); it_hit < rechitsBH.end();  ++it_hit) 
    {
      const HGCalDetId detid = it_hit->detid();
      if (storedRecHits_.find(detid) == storedRecHits_.end()) 
      {
        fillRecHit(detid, -1);
      }
    }
  }//rawRecHits_
  t_->Fill();
}

void HGCAnalyzer::fillRecHit(const DetId &detid, const float &fraction, const int &cluster_index_)
{
  const HGCRecHit *hit = hitmap_[detid];
  //edm::ESHandle<HGCalGeometry> geoHandle_ee;
  //iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandle_ee);
  //edm::ESHandle<HGCalGeometry> geoHandle_hef;
  //iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",geoHandle_hef);
  //const GlobalPoint position = geoHandle_ee->getPosition(detid);
  rechit_energy_.push_back(hit->energy());  
  rechit_detid_.push_back(detid);
  //rechit_x_.push_back(position.x());
  //rechit_y_.push_back(position.y());
  //rechit_z_.push_back(position.z());
  rechit_time_.push_back(hit->time());
  storedRecHits_.insert(detid);
  detIdToRecHitIndexMap_[detid] = rechit_index_;
  ++rechit_index_;
}

// ------------ method called once each job just before starting event loop  ------------
void HGCAnalyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void HGCAnalyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HGCAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCAnalyzer);
