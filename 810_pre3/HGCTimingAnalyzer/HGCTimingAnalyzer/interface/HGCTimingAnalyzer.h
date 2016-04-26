// -*- C++ -*-
//
// Package:    HGCTimingAnalyzer/HGCTimingAnalyzer
// Class:      HGCTimingAnalyzer
// 
/**\class HGCTimingAnalyzer HGCTimingAnalyzer.cc HGCTimingAnalyzer/HGCTimingAnalyzer/plugins/HGCTimingAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Saptaparna Bhattacharya
//         Created:  Sun, 24 Apr 2016 01:18:28 GMT
//
//


// system include files
#ifndef HGCTIMINGANALYZER_H
#define HGCTIMINGANALYZER_H
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
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
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class HGCTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HGCTimingAnalyzer(const edm::ParameterSet&);
      ~HGCTimingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      Int_t run_,event_,lumi_,vertex_nTracks;
      Float_t vertex_x, vertex_y, vertex_z, vertex_pt, vertex_nChi2;
      std::vector<float> recHit_x;
      std::vector<float> recHit_y; 
      std::vector<float> recHit_z;
      std::vector<float> recHit_energy; 
      std::vector<float> recHit_time;
      std::vector<uint32_t> recHit_recoDetId; 
      TLorentzVector *genVertex_;
      //TFile *outputFile_;
      TBranch *branch_;      
      TTree *tree_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > srcGenParticles_;      
      edm::EDGetTokenT<std::vector<SimTrack> > srcSimTracks_;
      edm::EDGetTokenT<std::vector<SimVertex> > srcSimVertices_; 
      edm::EDGetTokenT<edm::HepMCProduct> srcHepmcevent_;
      edm::EDGetTokenT<std::vector<reco::PFRecHit> > srcPFRecHit_;
      };
      //
      // //
      // // constants, enums and typedefs
      // //
      //
      // //
      // // static data member definitions
      // //
      //
      // //
      // // constructors and destructor
      // //
      //
#endif
