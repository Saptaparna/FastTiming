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
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalRecHitWorkerBaseClass.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/HGCalRecHitSimpleAlgo.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/HGCalCommonData/interface/HGCalParameters.h"
#include "Geometry/HGCalCommonData/interface/HGCalDDDConstants.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include <iterator>
#include <iostream>
#include <algorithm>
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerFECodecBase.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerBackendProcessor.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

namespace 
{  
  template<typename T>
  struct array_deleter
  {
    void operator () (T* arr) { delete [] arr; }
  };
}

class HGCTimingAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HGCTimingAnalyzer(const edm::ParameterSet&);
      ~HGCTimingAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void setTreeCellSize(const size_t);
      // ----------member data ---------------------------
      Int_t run_,event_,lumi_,vertex_nTracks;
      Float_t vertex_x, vertex_y, vertex_z, vertex_pt, vertex_nChi2;
      std::vector<float> recHit_x;
      std::vector<float> recHit_y; 
      std::vector<float> recHit_z;
      std::vector<float> recHit_energy; 
      std::vector<float> recHit_time;
      std::vector<float> recHit_pt;
      std::vector<float> recHit_eta;
      std::vector<float> recHit_phi;
      std::vector<float> genParticle_eta;
      std::vector<float> genParticle_phi;
      std::vector<float> genParticle_pdgId;
      std::vector<float> cluster_x;
      std::vector<float> cluster_y;
      std::vector<float> cluster_z;
      std::vector<float> cluster_energy;
      std::vector<float> cluster_time;
      std::vector<float> cluster_layer;
      std::vector<float> cluster_eta;
      std::vector<float> simCluster_eta;
      std::vector<float> simCluster_phi;
      std::vector<float> simCluster_energy;
      std::vector<float> simCluster_pt;
      std::vector<float> simCluster_simEnergy;
      std::vector<float> recHit_layer;
      std::vector<float> v_isIn3x3;
      std::vector<float> v_isIn5x5;
      std::vector<float> v_isIn7x7;
      std::vector<float> pcaShowerPos_x;
      std::vector<float> pcaShowerPos_y;
      std::vector<float> pcaShowerPos_z;
      std::vector<float> pcaShowerDir_x;
      std::vector<float> pcaShowerDir_y;
      std::vector<float> pcaShowerDir_z;
      TLorentzVector *genVertex_;
      TBranch *branch_;      
      TTree *tree_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > srcGenParticles_;      
      edm::EDGetTokenT<std::vector<SimTrack> > srcSimTracks_;
      edm::EDGetTokenT<std::vector<SimVertex> > srcSimVertices_; 
      edm::EDGetTokenT<std::vector<reco::PFRecHit> > srcPFRecHit_;
      edm::EDGetTokenT<std::vector<PCaloHit> > srcCaloHit_;
      HGCalTriggerGeometryBase::es_info es_info_;
      std::unique_ptr<HGCalTriggerGeometryBase> triggerGeometry_;
      edm::EDGetTokenT<std::vector<reco::PFCluster> > srcPFCluster_; 
      edm::EDGetTokenT<std::vector<CaloParticle> > srcCaloParticle_;
      edm::EDGetToken srcDigiee_, srcDigifh_, srcDigibh_;
      edm::EDGetTokenT<edm::SortedCollection<HGCRecHit> > srcRecHitEE_;
      edm::EDGetTokenT<edm::SortedCollection<HGCRecHit> > srcRecHitHEB_;
      edm::EDGetTokenT<edm::SortedCollection<HGCRecHit> > srcRecHitHEF_;
      hgcal::RecHitTools recHitTools;
      int nhits_ee;
      int nhits_hf;
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
