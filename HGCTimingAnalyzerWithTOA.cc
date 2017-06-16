#include "HGCTimingAnalyzer/HGCTimingAnalyzer/interface/HGCTimingAnalyzerWithTOA.h"

// Returns the resolution of the recHit time.
double HGCTimingAnalyzerWithTOA::getTimeHit(int thick, double SoverN){
  timeResolution->SetParameters(paramA[thick], paramC[thick]);
  double sigma = 0.2;
  if (SoverN > 1000) sigma = paramC[thick];
  else if (SoverN <= 1) sigma = timeResolution->Eval(1.);
  else sigma = timeResolution->Eval(SoverN);
  return sigma;
}

HGCTimingAnalyzerWithTOA::HGCTimingAnalyzerWithTOA(const edm::ParameterSet& iConfig)
:   srcGenParticles_ (consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag> ("srcGenParticles"))),
    srcSimTracks_ (consumes<std::vector<SimTrack> >(iConfig.getParameter<edm::InputTag> ("srcSimTracks"))), 
    srcSimVertices_ (consumes<std::vector<SimVertex> >(iConfig.getParameter<edm::InputTag> ("srcSimVertices"))),
    srcPFRecHit_ (consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag> ("srcPFRecHit"))),
    srcPFCluster_ (consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag> ("srcPFCluster"))),
    srcCaloParticle_ (consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag> ("srcCaloParticle"))),
    srcRecHitEE_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitEE"))),
    srcRecHitHEF_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitHEF"))),
    srcRecHitBH_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitBH")))
{
   usesResource("TFileService");
   const edm::ParameterSet& geometryConfig = iConfig.getParameterSet("TriggerGeometry");
   const std::string& trigGeomName = geometryConfig.getParameter<std::string>("TriggerGeometryName");
   nhits_ee=0;
   nhits_hf=0;
   nhits_bh=0;
   keV2fC[0] =  iConfig.getParameter<double>("HGCEE_keV2fC");
   keV2fC[1] =  iConfig.getParameter<double>("HGCHEF_keV2fC");
   keV2MIP   =  iConfig.getParameter<double>("HGCHB_keV2MIP");

   noisefC[0] = (iConfig.getParameter<std::vector<double> >("HGCEE_noisefC")).at(0);
   noisefC[1] = (iConfig.getParameter<std::vector<double> >("HGCEF_noisefC")).at(1);
   noiseMIP = iConfig.getParameter<double>("HGCBH_noiseMIP");

   //cell type
   fCPerMIP[0] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(0);
   fCPerMIP[1] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(1);
   fCPerMIP[2] =  (iConfig.getParameter<std::vector<double> >("HGCEE_fCPerMIP")).at(2);

   const auto& rcorr = iConfig.getParameter<std::vector<double> >("thicknessCorrection");
   scaleCorrection.clear();
   scaleCorrection.push_back(1.f);
   for( auto corr : rcorr ) {
     scaleCorrection.push_back(1.0/corr);
   }
   const auto& dweights = iConfig.getParameter<std::vector<double> >("dEdXweights");
   for( auto weight : dweights ) {
     weights.push_back(weight);
     //std::cout << "weight = " << weight << std::endl;
   }

   keV2GeV = 1e-6;
   keV2MeV = 1e-3;
/*
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
*/

   //100um
   paramA[0] = 1.50;//Jim's corrected number
   //paramA[0] = 3.00;//Jim's point 2//change as suggested by Michael
   parErrA[0] = 0.01;
   paramC[0] = 0.020;//Jim's point 1 case 1
   //paramC[0] = 0.030;//Jim's point 1 case 2
   parErrC[0] = 0.001;

   //200um
   paramA[1] = 1.8;//Jim's point 2
   parErrA[1] = 0.02;
   paramC[1] = 0.020;//Jim's point 1, case 1
   //paramC[1] = 0.030;//Jim's point 1 case 2
   parErrC[1] = 0.001;

   //300um
   paramA[2] = 1.2;//Jim's point 2
   parErrA[2] = 0.02;
   paramC[2] = 0.020;//Jim's point 1, case 1
   //paramC[2] = 0.030;//Jim's point 1 case 2
   parErrC[2] = 0.001;
/*
   //100um
   eneryDegradation[0] = 0.50;
   //200um
   eneryDegradation[1] = 0.50;
   //300um
   eneryDegradation[2] = 0.70;
*/
 
   eneryDegradation[0] = 1.0;
   eneryDegradation[1] = 1.0;
   eneryDegradation[2] = 1.0;

   timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);                     
}


HGCTimingAnalyzerWithTOA::~HGCTimingAnalyzerWithTOA()
{
}

// ------------ method called for each event  ------------
void
HGCTimingAnalyzerWithTOA::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   iSetup.get<CaloGeometryRecord>().get(pG);
   const CaloGeometry* caloGeom = pG.product();
   recHitTools.getEventSetup(iSetup);
   run_   = iEvent.id().run();
   event_ = iEvent.id().event();
   lumi_  = iEvent.luminosityBlock();
   vertex_x = -999.0;
   vertex_y = -999.0;
   vertex_z = -999.0;
   vertex_nTracks = 0;
   vertex_pt = 0.0;
   vertex_nChi2 = 0.0;
   recHit_x.clear();
   recHit_y.clear();
   recHit_z.clear();
   cluster_x.clear();
   cluster_y.clear();
   cluster_z.clear();
   cluster_eta.clear();
   cluster_energy.clear();
   cluster_time.clear();
   recHit_energy.clear();
   recHit_time.clear();
   recHit_smearedTime.clear();
   recHit_layer.clear();
   recHit_soverN.clear();
   recHit_rmsTime.clear();
   recHit_deltaTime.clear();
   cluster_layer.clear();
   genParticle_eta.clear();
   genParticle_phi.clear();
   genParticle_pdgId.clear();
   simCluster_eta.clear();
   simCluster_phi.clear();
   simCluster_energy.clear();
   simCluster_pt.clear();
   simCluster_simEnergy.clear();
   recHit_pt.clear();
   recHit_eta.clear();
   recHit_phi.clear();

   edm::Handle<std::vector<SimTrack> > simTk;
   iEvent.getByToken(srcSimTracks_, simTk);
   edm::Handle<std::vector<SimVertex> > simVtx;
   iEvent.getByToken(srcSimVertices_, simVtx);
   edm::Handle<std::vector<reco::GenParticle> > genParticles;
   iEvent.getByToken(srcGenParticles_, genParticles);
 
   for (unsigned int j=0; j<simTk->size(); j++) 
   {
     //double pt = simTk->at(j).momentum().pt();
     //if(pt < 2.0) continue;
     //double eta = simTk->at(j).momentum().eta();
     //if (simTk->at(j).noVertex()) continue;
     int vertIndex = simTk->at(j).vertIndex();
     vertex_x = simVtx->at(vertIndex).position().x();
     vertex_y = simVtx->at(vertIndex).position().y();
     vertex_z = simVtx->at(vertIndex).position().z();
   }
   
   for(size_t i = 0; i < genParticles->size(); ++ i)
   {
     const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>((*genParticles)[i]);
     genParticle_eta.push_back(p.eta());
     genParticle_phi.push_back(p.phi());
     genParticle_pdgId.push_back(p.pdgId());
   }
 
   edm::Handle<edm::SortedCollection<HGCRecHit> > srcRecHitEE;
   iEvent.getByToken(srcRecHitEE_, srcRecHitEE);
   edm::Handle<edm::SortedCollection<HGCRecHit> > srcRecHitHEF;
   iEvent.getByToken(srcRecHitHEF_, srcRecHitHEF);
   edm::Handle<edm::SortedCollection<HGCRecHit> > srcRecHitBH;
   iEvent.getByToken(srcRecHitBH_, srcRecHitBH);
   edm::Handle<std::vector<reco::PFCluster> > pfCluster;
   iEvent.getByToken(srcPFCluster_, pfCluster); 
   edm::Handle<std::vector<CaloParticle> > caloParticle;
   iEvent.getByToken(srcCaloParticle_, caloParticle);
   edm::ESHandle<HGCalGeometry> geoHandle_ee;
   iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",geoHandle_ee);
   edm::ESHandle<HGCalGeometry> geoHandle_hef;
   iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",geoHandle_hef);

   //EE 
   for(unsigned int l=0; l<pfCluster->size(); l++) // Iterating over sim clusters
   {
       cluster_x.push_back(pfCluster->at(l).position().X());
       cluster_y.push_back(pfCluster->at(l).position().Y());
       cluster_z.push_back(pfCluster->at(l).position().Z());
       cluster_energy.push_back(pfCluster->at(l).energy());
       cluster_time.push_back(pfCluster->at(l).time());
       cluster_layer.push_back(pfCluster->at(l).layer());
       cluster_eta.push_back(pfCluster->at(l).eta());   
       for(unsigned ifrac = 0; ifrac < pfCluster->at(l).hitsAndFractions().size(); ifrac++)
       {
         uint32_t id = pfCluster->at(l).hitsAndFractions().at(ifrac).first;
         edm::SortedCollection<HGCRecHit>::const_iterator hgcHitEE = srcRecHitEE->find(id);
         if(hgcHitEE != srcRecHitEE->end())
         {
           const HGCalDetId detId_recHit_ee = hgcHitEE->detid();
           //DetId detId_recHit_ee = hgcHitEE->detid();
           if(id==detId_recHit_ee.rawId())// and hgcHitEE->time() > 0)
           {
             nhits_ee++;
             //unsigned int layerSecond = recHitTools.getLayerWithOffset(id);
             //recHit_layer.push_back(layerSecond);
             const GlobalPoint& recHitPos_ee = geoHandle_ee->getPosition(id); 
             recHit_phi.push_back(recHitTools.getPhi(recHitPos_ee));
             recHit_eta.push_back(recHitTools.getEta(recHitPos_ee, vertex_z));
             recHit_pt.push_back(recHitTools.getPt(recHitPos_ee, hgcHitEE->energy(), vertex_z));
             recHit_x.push_back(recHitPos_ee.x());
             recHit_y.push_back(recHitPos_ee.y());
             recHit_z.push_back(recHitPos_ee.z());
             recHit_energy.push_back(hgcHitEE->energy());
             auto ddd = get_ddd(caloGeom->getSubdetectorGeometry(detId_recHit_ee));
             int thick = ddd->waferTypeL(detId_recHit_ee.wafer()) - 1;
             unsigned int layer = detId_recHit_ee.layer();
             recHit_layer.push_back(layer);
             int sectionType;
             sectionType = 2;
             if(layer < 29) sectionType = 0;
             else if(layer < 41) sectionType = 1;
             double sigmaNoiseMIP = 0.0;
             sigmaNoiseMIP = noiseMIP;
             if(sectionType != 2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];
             float energyMIP;
             float energy = hgcHitEE->energy()*eneryDegradation[thick];
             if(sectionType == 2) energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV );
             else energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV);
             if(energyMIP > 3.) 
             {
               float SoverN = energyMIP / sigmaNoiseMIP;
               double rmsTime = getTimeHit(thick, SoverN);               
	       TRandom3* rand = new TRandom3();
	       rand->SetSeed(0);
	       double deltaTime = rand->Gaus(0., rmsTime);
	       double smearedTime = hgcHitEE->time()-1. + deltaTime;
	       recHit_time.push_back(hgcHitEE->time()-1);
               recHit_smearedTime.push_back(smearedTime);
               recHit_soverN.push_back(SoverN);
               recHit_rmsTime.push_back(rmsTime);
               recHit_deltaTime.push_back(deltaTime); 
             }
             else 
             {
               recHit_time.push_back(-99.0);
               recHit_smearedTime.push_back(-99.0);
               recHit_soverN.push_back(-99.0);
               recHit_rmsTime.push_back(-99.0);
               recHit_deltaTime.push_back(-99.0);
             }
           }
         }
       }
       
       for(unsigned ifrac = 0; ifrac < pfCluster->at(l).hitsAndFractions().size(); ifrac++)
       {
         uint32_t id = pfCluster->at(l).hitsAndFractions().at(ifrac).first;
         edm::SortedCollection<HGCRecHit>::const_iterator hgcHitHEF = srcRecHitHEF->find(id);
         if(hgcHitHEF != srcRecHitHEF->end())
         {
           //DetId detId_recHit_hef = hgcHitHEF->detid();
           const HGCalDetId detId_recHit_hef = hgcHitHEF->detid();
           if(id==detId_recHit_hef.rawId())// and hgcHitHEF->time() > 0)
           { 
             nhits_hf++;
             const GlobalPoint& recHitPos_hef = geoHandle_hef->getPosition(id);
             //unsigned int layer = recHitTools.getLayerWithOffset(id);//alternate way to access the layer
             recHit_x.push_back(recHitPos_hef.x());
             recHit_y.push_back(recHitPos_hef.y());
             recHit_z.push_back(recHitPos_hef.z());
             recHit_energy.push_back(hgcHitHEF->energy());
             auto ddd = get_ddd(caloGeom->getSubdetectorGeometry(detId_recHit_hef));
             int thick = ddd->waferTypeL(detId_recHit_hef.wafer()) - 1;
             unsigned int layer = detId_recHit_hef.layer();
             recHit_layer.push_back(layer);
             int sectionType;
             sectionType = 2;
             if(layer < 29) sectionType = 0;
             else if(layer < 41) sectionType = 1;
             double sigmaNoiseMIP = 0.0;
             sigmaNoiseMIP = noiseMIP;
             if(sectionType !=2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];
             float energyMIP;
             //float energy = hgcHitHEF->energy();
             float energy = hgcHitHEF->energy()*eneryDegradation[thick];
             if(sectionType == 2) energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV );
             else energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV);
             if(energyMIP > 3.)
             {
               float SoverN = energyMIP / sigmaNoiseMIP;
               double rmsTime = getTimeHit(thick, SoverN);
               TRandom3* rand = new TRandom3();
               rand->SetSeed(0);
               double deltaTime = rand->Gaus(0., rmsTime);
               double smearedTime = hgcHitHEF->time()-1. + deltaTime;
               recHit_time.push_back(hgcHitHEF->time()-1);
               recHit_smearedTime.push_back(smearedTime);
               recHit_soverN.push_back(SoverN);
               recHit_rmsTime.push_back(rmsTime);
               recHit_deltaTime.push_back(deltaTime);
             }
             else
             {
               recHit_time.push_back(-99.0);
               recHit_smearedTime.push_back(-99.0);
               recHit_soverN.push_back(-99.0);
               recHit_rmsTime.push_back(-99.0);
               recHit_deltaTime.push_back(-99.0);
             }
           }
         }
       }
    
     for(unsigned ifrac = 0; ifrac < pfCluster->at(l).hitsAndFractions().size(); ifrac++)
     {
         uint32_t id = pfCluster->at(l).hitsAndFractions().at(ifrac).first;
         edm::SortedCollection<HGCRecHit>::const_iterator hgcHitBH = srcRecHitBH->find(id);
         if(hgcHitBH != srcRecHitBH->end())
         {
           const HGCalDetId detId_recHit_bh = hgcHitBH->detid();
           if(id==detId_recHit_bh.rawId())// and hgcHitHEF->time() > 0)
           {
             nhits_bh++;
             //recHit_time.push_back(hgcHitBH->time());
             std::cout << "hgcHitBH->time() = " << hgcHitBH->time() << std::endl;
           }
         }
     }  
 
   }
   // loop over caloParticles based on https://github.com/CMS-HGCAL/reco-ntuples/blob/master/HGCalAnalysis/plugins/HGCalAnalysis.cc
   for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticle->begin(); it_caloPart != caloParticle->end(); ++it_caloPart) 
   {
     const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
     std::vector<uint32_t> simClusterIndex;
     simCluster_eta.push_back(it_caloPart->eta());
     simCluster_phi.push_back(it_caloPart->phi());
     simCluster_energy.push_back(it_caloPart->energy());
     simCluster_pt.push_back(it_caloPart->pt());
     simCluster_simEnergy.push_back(it_caloPart->simEnergy());
   }
tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void HGCTimingAnalyzerWithTOA::beginJob()
{
  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("HGCTiming","HGCTiming");
  branch_=tree_->Branch("run",   &run_,   "run/I");
  branch_=tree_->Branch("event", &event_, "event/I");
  branch_=tree_->Branch("lumi",  &lumi_,  "lumi/I");
  branch_=tree_->Branch("vertex_x", &vertex_x, "vertex_x/F");
  branch_=tree_->Branch("vertex_y", &vertex_y, "vertex_y/F");
  branch_=tree_->Branch("vertex_z", &vertex_z, "vertex_z/F");
  branch_=tree_->Branch("recHit_energy", &recHit_energy);
  branch_=tree_->Branch("recHit_x", &recHit_x);
  branch_=tree_->Branch("recHit_y", &recHit_y);
  branch_=tree_->Branch("recHit_z", &recHit_z);
  branch_=tree_->Branch("recHit_pt", &recHit_pt);
  branch_=tree_->Branch("recHit_eta", &recHit_eta);
  branch_=tree_->Branch("recHit_phi", &recHit_phi);
  branch_=tree_->Branch("cluster_x", &cluster_x);
  branch_=tree_->Branch("cluster_y", &cluster_y);
  branch_=tree_->Branch("cluster_z", &cluster_z);
  branch_=tree_->Branch("cluster_energy", &cluster_energy);
  branch_=tree_->Branch("cluster_time", &cluster_time);
  branch_=tree_->Branch("recHit_time", &recHit_time);
  branch_=tree_->Branch("recHit_smearedTime", &recHit_smearedTime);
  branch_=tree_->Branch("recHit_soverN", &recHit_soverN);
  branch_=tree_->Branch("recHit_rmsTime", &recHit_rmsTime);
  branch_=tree_->Branch("recHit_deltaTime", &recHit_deltaTime);
  branch_=tree_->Branch("recHit_soverN", &recHit_soverN);
  branch_=tree_->Branch("recHit_layer", &recHit_layer);
  branch_=tree_->Branch("cluster_layer", &cluster_layer);
  branch_=tree_->Branch("cluster_eta", &cluster_eta);
  branch_=tree_->Branch("genParticle_eta", &genParticle_eta);
  branch_=tree_->Branch("genParticle_phi", &genParticle_phi);
  branch_=tree_->Branch("genParticle_pdgId", &genParticle_pdgId);
  branch_=tree_->Branch("simCluster_eta", &simCluster_eta);
  branch_=tree_->Branch("simCluster_phi", &simCluster_phi);
  branch_=tree_->Branch("simCluster_energy", &simCluster_energy);
  branch_=tree_->Branch("simCluster_pt", &simCluster_pt);
  branch_=tree_->Branch("simCluster_simEnergy", &simCluster_simEnergy);
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
HGCTimingAnalyzerWithTOA::endJob() 
{
  std::cout << "nhits_ee = " << nhits_ee << std::endl;
  std::cout << "nhits_hf = " << nhits_hf << std::endl;
  std::cout << "nhits_bh = " << nhits_bh << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCTimingAnalyzerWithTOA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

const HGCalDDDConstants* HGCTimingAnalyzerWithTOA::get_ddd(const CaloSubdetectorGeometry* geom){
  const HGCalGeometry* hg = static_cast<const HGCalGeometry*>(geom);
  const HGCalDDDConstants* ddd = &(hg->topology().dddConstants());
  if( nullptr == ddd ) {
    throw cms::Exception("hgcal::RecHitTools")
      << "DDDConstants not accessibl to hgcal::RecHitTools!";
  }
  return ddd;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCTimingAnalyzerWithTOA);
