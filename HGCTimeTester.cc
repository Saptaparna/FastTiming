// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"



#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2F.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//#include "RecoLocalCalo/HGCalRecAlgos/interface/ClusterTools.h"
//#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalDepthPreClusterer.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"



#include "RecoNtuples/HGCalAnalysis/plugins/RecHiTimeEstimator.h"

#include <string>
#include <map>

class HGCTimeTester : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
//
// constructors and destructor
//
  explicit HGCTimeTester(const edm::ParameterSet&);
  ~HGCTimeTester();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

 // ----------member data ---------------------------

  edm::EDGetTokenT<HGCRecHitCollection> _recHitsEE;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsFH;
  edm::EDGetTokenT<HGCRecHitCollection> _recHitsBH;

  RecHiTimeEstimator* myEstimator;

  Int_t run_,event_,lumi_;
  std::vector<double> Energydiff;
  std::vector<double> Timediff;
  std::vector<double> TimeOld;
  std::vector<double> Time;
  TBranch *branch_;
  TTree *tree_;
  //TH1F* hEnergydiff;
  //TH1F* hTimediff;
  //TH1F* hTimeOld;
  //TH1F* hTime;

};  



HGCTimeTester::HGCTimeTester(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  _recHitsEE = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCEEInput"));
  _recHitsFH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCFHInput"));
  _recHitsBH = consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("HGCBHInput"));

  myEstimator = new RecHiTimeEstimator(iConfig);

  edm::Service<TFileService> fs;


  //hEnergydiff = fs->make<TH1F>("hEnergydiff", "", 50000, -5, 5);
  //hTimediff = fs->make<TH1F>("hTimediff", "", 50000, -10, 10);
  //hTimeOld = fs->make<TH1F>("hTimeOld", "", 50000, -10, 10);
  //hTime = fs->make<TH1F>("hTime", "", 50000, -10, 10);



}

HGCTimeTester::~HGCTimeTester()
{
  //delete myEstimator;
  //std::cout<<"myEstimator = "<< myEstimator<<std::endl;
}


void
HGCTimeTester::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout << " >>> analyzer " << std::endl;
  using namespace edm;
  Energydiff.clear();
  Timediff.clear();
  TimeOld.clear();
  Time.clear();

  run_   = iEvent.id().run();
  event_ = iEvent.id().event();
  lumi_  = iEvent.luminosityBlock();
  myEstimator->setEventSetup(iSetup);

  Handle<HGCRecHitCollection> recHitHandleEE;
  Handle<HGCRecHitCollection> recHitHandleFH;
  Handle<HGCRecHitCollection> recHitHandleBH;

  iEvent.getByToken(_recHitsEE,recHitHandleEE);
  iEvent.getByToken(_recHitsFH,recHitHandleFH);
  iEvent.getByToken(_recHitsBH,recHitHandleBH);
  
  const HGCRecHitCollection& rechitsEEOld = *recHitHandleEE;
  HGCRecHitCollection NewrechitsEE; 
  myEstimator->correctTime(rechitsEEOld, &NewrechitsEE);
  std::cout << "EE old size = " << rechitsEEOld.size() << " new size = " << NewrechitsEE.size() << std::endl; 

  /*const HGCRecHitCollection& rechitsFHOld = *recHitHandleFH;
  HGCRecHitCollection NewrechitsFH; 
  myEstimator->correctTime(rechitsFHOld, &NewrechitsFH);*/
  //std::cout << "FH old size = " << rechitsFHOld.size() << " new size = " << NewrechitsFH.size() << std::endl; 


  // loop over EE new RecHits                                                                                                                                                 
  HGCRecHitCollection::const_iterator it_hitOld = rechitsEEOld.begin();
  for(HGCRecHitCollection::const_iterator it_hit = NewrechitsEE.begin(); it_hit < NewrechitsEE.end(); ++it_hit) {
    float energyN = it_hit->energy();
    float energyO = it_hitOld->energy();
    Energydiff.push_back(energyN - energyO);
    Timediff.push_back(it_hit->time() - it_hitOld->time());
    TimeOld.push_back(it_hitOld->time());
    Time.push_back(it_hit->time());
    //if(it_hit->time() > 0.0) std::cout << "it_hit->time() = " << it_hit->time() << std::endl;
    /*hEnergydiff->Fill(energyN - energyO);
    hTimediff->Fill(it_hit->time() - it_hitOld->time());
    hTimeOld->Fill(it_hitOld->time());
    hTime->Fill(it_hit->time());*/
    ++it_hitOld;
  }
  tree_->Fill();
}

void
HGCTimeTester::beginJob()
{
  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("HGCTiming","HGCTiming");
  branch_=tree_->Branch("run",   &run_,   "run/I");
  branch_=tree_->Branch("event", &event_, "event/I");
  branch_=tree_->Branch("lumi",  &lumi_,  "lumi/I");
  branch_=tree_->Branch("Energydiff", &Energydiff);
  branch_=tree_->Branch("Timediff", &Timediff);
  branch_=tree_->Branch("TimeOld", &TimeOld);
  branch_=tree_->Branch("Time", &Time);
}

// ------------ method called once each job just after ending the event loop  ------------
void
HGCTimeTester::endJob()
{
  //delete myEstimator;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCTimeTester::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCTimeTester);
