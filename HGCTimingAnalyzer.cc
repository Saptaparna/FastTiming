#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HGCTimingAnalyzer/HGCTimingAnalyzer/interface/HGCTimingAnalyzer.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/Math/interface/deltaR.h"
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
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <iostream> 

HGCTimingAnalyzer::HGCTimingAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

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

   edm::Handle<std::vector<SimTrack> > simTk;
   iEvent.getByLabel("g4TracksSource", simTk);
   edm::Handle<std::vector<SimVertex> > simVtx;
   iEvent.getByLabel("g4VerticesSource", simVtx);
   edm::Handle<std::vector<int> > genParticlePDGId;
   iEvent.getByLabel("genParticles", genParticlePDGId);
   edm::Handle<reco::VertexCollection> recoVtx;
   iEvent.getByLabel("offlinePrimaryVertices", recoVtx);

   for(size_t iv=0; iv<recoVtx->size(); iv++)
   { 
     const reco::Vertex &vtx=recoVtx->at(iv);
     if(!vtx.isValid()) continue;
     if(vtx.isFake()) continue;
     vertex_x = vtx.x();
     vertex_y = vtx.y();
     vertex_z = vtx.z();
     vertex_nTracks = vtx.nTracks();
     vertex_pt = vtx.p4().pt();
     vertex_nChi2 = vtx.normalizedChi2();
   } 

   // get the gen particles
   edm::Handle<edm::View<reco::Candidate> > genParticles;
   iEvent.getByLabel("genParticles", genParticles);
   edm::Handle<reco::GenParticleCollection> genCandsFromSimTracks;
   iEvent.getByLabel("genCandsFromSimTracks", genCandsFromSimTracks);
  
   // get the vertex info
   edm::Handle<edm::HepMCProduct>  hepmcevent;
   iEvent.getByLabel("generator", hepmcevent);
   const HepMC::GenEvent& genevt = hepmcevent->getHepMCData();
   genVertex_->SetXYZT(0.,0.,0.,0.);
   if( genevt.vertices_size() > 0.0) 
   {
    HepMC::FourVector temp = (*genevt.vertices_begin())->position() ;
    genVertex_->SetXYZT(0.1*temp.x(),0.1*temp.y(),0.1*temp.z(),temp.t()/299.792458); // convert positions to cm and time to ns (it's in mm to start)
   }

   for(size_t i = 0; i < genParticles->size(); ++ i)
   {
     const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>((*genParticles)[i]); 
     if(p.mother()) std::cout << " p.mother()->pdgId() = " << p.mother()->pdgId() << std::endl;
   }  

   //EE hits
   edm::Handle<edm::PCaloHitContainer> eeSimHits;
   iEvent.getByLabel(edm::InputTag("g4SimHits", "eeSimHitsSource"), eeSimHits);
   edm::Handle<HGCRecHitCollection> eeRecHits;
   iEvent.getByLabel(edm::InputTag("HGCalRecHit","HGCEERecHits"),eeRecHits);
   edm::ESHandle<HGCalGeometry> eeGeom;
   iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",eeGeom);
   const HGCalTopology &topo=eeGeom->topology();
   const HGCalDDDConstants &dddConst=topo.dddConstants();

   if(eeRecHits.isValid())
   {

     for(HGCRecHitCollection::const_iterator hit_it=eeRecHits->begin(); hit_it!=eeRecHits->end(); hit_it++)
     {
        uint32_t recoDetId(hit_it->id());
        const GlobalPoint pos( std::move( eeGeom->getPosition(recoDetId) ) );
        recHit_x.push_back(pos.x());
        recHit_y.push_back(pos.y());
        recHit_z.push_back(pos.z());
        recHit_energy.push_back(hit_it->energy());
        recHit_time.push_back(hit_it->time());
        std::cout << "dddConst.getFirstModule(true)->cellSize = " << dddConst.getFirstModule(true)->cellSize << std::endl;
     } 
   }

   if(eeSimHits.isValid())
   {
     for(edm::PCaloHitContainer::const_iterator hit_it=eeSimHits->begin(); hit_it!=eeSimHits->end(); hit_it++)
     {
        HGCalDetId simId(hit_it->id());
        int layer(simId.layer()),cell(simId.cell());
        std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
        cell  = recoLayerCell.first;
        layer = recoLayerCell.second;
 
     }
   }


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
  branch_=tree_->Branch("vertex_nTracks", &vertex_nTracks, "vertex_nTracks/I");
  branch_=tree_->Branch("vertex_pt", &vertex_pt, "vertex_pt/F");
  branch_=tree_->Branch("vertex_nChi2", &vertex_nChi2, "vertex_nChi2/F");
  branch_=tree_->Branch("recHit_x", &recHit_x);
  branch_=tree_->Branch("recHit_y", &recHit_y);
  branch_=tree_->Branch("recHit_z", &recHit_z);
  branch_=tree_->Branch("GenVertex", "TLorentzVector", &genVertex_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HGCTimingAnalyzer::endJob() 
{
  outputFile_->Write();
  outputFile_->Close();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
HGCTimingAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HGCTimingAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HGCTimingAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HGCTimingAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

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
