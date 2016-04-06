#include <memory>
#include "HGCTimingAnalyzer/HGCTimingAnalyzer/interface/HGCTimingAnalyzer.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimG4CMS/Calo/interface/CaloHitID.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DetectorDescription/OfflineDBLoader/interface/GeometryInfoDump.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

HGCTimingAnalyzer::HGCTimingAnalyzer(const edm::ParameterSet& iConfig)
{
  //genSource_        = iConfig.getUntrackedParameter< std::string >("genSource");
  //genCandsFromSimTracksSource_ = iConfig.getUntrackedParameter< std::string >("genCandsFromSimTracksSource");

}


HGCTimingAnalyzer::~HGCTimingAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HGCTimingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run_   = iEvent.id().run();
   event_ = iEvent.id().event();
   lumi_  = iEvent.luminosityBlock();

   vertex_x = 0.0;
   vertex_y = 0.0;
   vertex_z = 0.0;
   vertex_nTracks = 0;
   vertex_pt = 0.0;
   vertex_nChi2 = 0.0;
   recHit_x.clear();
   recHit_y.clear();
   recHit_z.clear();
   recHit_energy.clear();
   recHit_time.clear();
   //gen particles
   //edm::Handle<std::vector<SimTrack> > simTk;
   //iEvent.getByLabel("g4SimHits", simTk);
   //edm::Handle<std::vector<SimVertex> > simVtx;
   //iEvent.getByLabel("g4SimHits", simVtx);
   edm::Handle<std::vector<int> > genParticlePDGId;
   iEvent.getByLabel("genParticles", genParticlePDGId);
   
   std::cout << "genParticlePDGId->size() = " << genParticlePDGId->size() << std::endl;


   /*for (unsigned int j=0; j<simTk->size(); j++) 
   {
     double pt = simTk->at(j).momentum().pt();
     if (pt<2.5) continue;
     double eta = simTk->at(j).momentum().eta();
     if (abs(eta)>2.5) continue;
     if (simTk->at(j).noVertex()) continue;
     int vertIndex = simTk->at(j).vertIndex();
     vertex_x = simVtx->at(vertIndex).position().x();
     vertex_y = simVtx->at(vertIndex).position().y();
     vertex_z = simVtx->at(vertIndex).position().z();
     //vertex_nTracks = simVtx->at(vertIndex)->tracksSize();
     //vertex_pt = simVtx->at(vertIndex).pt();
     //vertex_nChi2 = vtx.normalizedChi2();
   }*/
   //if(simVtx.isValid()) std::cout << "simVtx->size() = " << simVtx->size() << std::endl;
   
   // get the gen particles
   //edm::Handle<edm::View<reco::GenParticle> > genParticles;
   //iEvent.getByLabel("genParticles", genParticles);
   //edm::Handle<reco::GenParticleCollection> genCandsFromSimTracks;
   //iEvent.getByLabel("genCandsFromSimTracks", genCandsFromSimTracks);
 
   // get the vertex info
   /*edm::Handle<edm::HepMCProduct>  hepmcevent;
   iEvent.getByLabel("generator", hepmcevent);
   const HepMC::GenEvent& genevt = hepmcevent->getHepMCData();
   genVertex_->SetXYZT(0.,0.,0.,0.);
   if( genevt.vertices_size() > 0.0)
   {
    HepMC::FourVector temp = (*genevt.vertices_begin())->position() ;
    genVertex_->SetXYZT(0.1*temp.x(),0.1*temp.y(),0.1*temp.z(),temp.t()/299.792458); // convert positions to cm and time to ns (it's in mm to start)
   }
   */
   //for(size_t i = 0; i < genParticles->size(); ++ i)
   //{
     //const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>((*genParticles)[i]);
     //if(p.mother()) std::cout << " p.mother()->pdgId() = " << p.mother()->pdgId() << std::endl;
   //}
   tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HGCTimingAnalyzer::beginJob()
{
  outputFile_ = new TFile("HGCTiming.root","RECREATE");
  tree_=new TTree("HGCTiming","HGCTiming");
  branch_=tree_->Branch("run",   &run_,   "run/I");
  branch_=tree_->Branch("event", &event_, "event/I");
  branch_=tree_->Branch("lumi",  &lumi_,  "lumi/I");
  branch_=tree_->Branch("vertex_x", &vertex_x, "vertex_x/F");
  branch_=tree_->Branch("vertex_y", &vertex_y, "vertex_y/F");
  branch_=tree_->Branch("vertex_z", &vertex_z, "vertex_z/F");
  /*branch_=tree_->Branch("vertex_nTracks", &vertex_nTracks, "vertex_nTracks/I");
  branch_=tree_->Branch("vertex_pt", &vertex_pt, "vertex_pt/F");
  branch_=tree_->Branch("vertex_nChi2", &vertex_nChi2, "vertex_nChi2/F");*/
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HGCTimingAnalyzer::endJob() 
{ 
  outputFile_->Write();
  outputFile_->Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCTimingAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCTimingAnalyzer);
