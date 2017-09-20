#include "HGCTimingAnalyzer/HGCTimingAnalyzer/interface/HGCTimingAnalyzerWithTOAPU.h"

bool sameValCheck(double a, double b)
{
	return fabs(a - b) < 1.000e-02;
}

// Returns the resolution of the recHit time.
double HGCTimingAnalyzerWithTOAPU::getTimeHit(int thick, double SoverN){
	timeResolution->SetParameters(paramA[thick], paramC[thick]);
	double sigma = 0.2;
	if (SoverN > 1000) sigma = paramC[thick];
	else if (SoverN <= 1) sigma = timeResolution->Eval(1.);
	else sigma = timeResolution->Eval(SoverN);
	return sigma;
}

HGCTimingAnalyzerWithTOAPU::HGCTimingAnalyzerWithTOAPU(const edm::ParameterSet& iConfig)
	:   srcGenParticles_ (consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag> ("srcGenParticles"))),
	srcSimTracks_ (consumes<std::vector<SimTrack> >(iConfig.getParameter<edm::InputTag> ("srcSimTracks"))), 
	srcSimVertices_ (consumes<std::vector<SimVertex> >(iConfig.getParameter<edm::InputTag> ("srcSimVertices"))),
	srcPFRecHit_ (consumes<std::vector<reco::PFRecHit> >(iConfig.getParameter<edm::InputTag> ("srcPFRecHit"))),
	srcPFCluster_ (consumes<std::vector<reco::PFCluster> >(iConfig.getParameter<edm::InputTag> ("srcPFCluster"))),
	srcCaloParticle_ (consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag> ("srcCaloParticle"))),
	srcRecHitEE_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitEE"))),
	srcRecHitHEF_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitHEF"))),
	srcRecHitBH_ (consumes<edm::SortedCollection<HGCRecHit> >(iConfig.getParameter<edm::InputTag> ("srcRecHitBH"))),
	_part (consumes<std::vector<TrackingParticle> >(iConfig.getParameter<edm::InputTag> ("srcPartHandle")))  
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
	}

	keV2GeV = 1e-6;
	keV2MeV = 1e-3;
	/*
	//100um
	energyDegradation[0] = 0.50;
	//200um
	energyDegradation[1] = 0.50;
	//300um
	energyDegradation[2] = 0.70;
	*/
	//100um
	paramA[0]  =  iConfig.getParameter<double>("SlopeParameterA");
	parErrA[0] =  0.01;
        paramC[0]  =  iConfig.getParameter<double>("FloorParameterC");
	parErrC[0] =  0.001;
	//200um
	paramA[1]  =  iConfig.getParameter<double>("SlopeParameterA");
	parErrA[1] =  0.01;
        paramC[1]  =  iConfig.getParameter<double>("FloorParameterC");
	parErrC[1] =  0.001;
	//300um
	paramA[2]  =  iConfig.getParameter<double>("SlopeParameterA");
	parErrA[2] =  0.01;
        paramC[2]  =  iConfig.getParameter<double>("FloorParameterC");
	parErrC[2] =  0.001;
        
        /*paramA[0]  =  4.0;
        parErrA[0] =  0.01;
        paramC[0]  =  0.020;
        parErrC[0] =  0.001;

        paramA[1]  =  4.0;
        parErrA[1] =  0.01;
        paramC[1]  =  0.020;
        parErrC[1] =  0.001;        

        paramA[2]  =  4.0;
        parErrA[2] =  0.01;
        paramC[2]  =  0.020;
        parErrC[2] =  0.001;     
        */
	energyDegradation[0] = 1.0;
	energyDegradation[1] = 1.0;
	energyDegradation[2] = 1.0;
	/*
	energyDegradation[0] = 0.50;
	energyDegradation[1] = 0.50;
	energyDegradation[2] = 0.70;
        */	   
	//timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x/sqrt(2.), 2) + pow([1], 2) )", 1., 1000.);                    
	timeResolution = new TF1("timeSi100", "sqrt(pow([0]/x, 2) + pow([1], 2) )", 1., 1000.);  
}


HGCTimingAnalyzerWithTOAPU::~HGCTimingAnalyzerWithTOAPU()
{
}

// ------------ method called for each event  ------------
	void
HGCTimingAnalyzerWithTOAPU::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
	recHit_energyMIP.clear();
	recHit_id.clear();

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
	//for looking into non-interacting and non-converting particles
	edm::Handle<std::vector<TrackingParticle> > partHandle;
	const std::vector<TrackingParticle> *part;
	iEvent.getByToken(_part, partHandle);
	part = &(*partHandle);

	unsigned int npart = part->size();
	//std::cout << "npart = " << npart << std::endl;
	int reachedEE=-999;
	for(unsigned int i=0; i<npart; ++i) 
	{
		reachedEE=-500;
		// event=0 is the hard scattering (all PU have event()>=1)
		// bunchCrossing == 0 intime, buncCrossing!=0 offtime, standard generation has [-12,+3]
		if((*part)[i].eventId().event() ==0 and (*part)[i].eventId().bunchCrossing()==0) 
		{
			// default values for decay position is outside detector, i.e. ~stable
			reachedEE=1;
			if((*part)[i].decayVertices().size()>=1) 
			{ //they can be delta rays, in this case you have multiple decay verices
				if ((*part)[i].decayVertices()[0]->inVolume()) reachedEE=0; //if it decays inside the tracker volume
				break;
			}
		} 
	}
	reachedEE_=reachedEE;
        double axisX = 0;
        double axisY = 0;
        double axisZ = 0;
        double sumEnergyToNorm = 0;
        GlobalPoint showerAxis;
        for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticle->begin(); it_caloPart != caloParticle->end(); ++it_caloPart)
        {
          const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
          int bc = it_caloPart->eventId().bunchCrossing();
          int ev = it_caloPart->eventId().event(); 
          //std::cout << "bc = " << bc << std::endl;
          //std::cout << "ev = " << ev << std::endl;
          //std::cout << "it_caloPart->eta() = " << it_caloPart->eta() << std::endl;
          //std::cout << " simClusterRefVector.size() = " << simClusterRefVector.size() << std::endl; 
          for (CaloParticle::sc_iterator it_sc = simClusterRefVector.begin(); it_sc != simClusterRefVector.end(); ++it_sc) 
          {
            const std::vector<std::pair<uint32_t,float> > hits_and_fractions = (*it_sc)->hits_and_fractions();      
            for (std::vector<std::pair<uint32_t,float> >::const_iterator it_haf = hits_and_fractions.begin(); it_haf != hits_and_fractions.end(); ++it_haf) 
            {
              DetId id = (it_haf->first);
              if(bc==0 and ev==0) std::cout << "it_haf->second signal = " << it_haf->second << std::endl;
              else if(bc!=0)  std::cout << "it_haf->second PU = " << it_haf->second << std::endl;
              //std::cout << "id = " << id.rawId() << std::endl; 

              edm::SortedCollection<HGCRecHit>::const_iterator hgcHitEE = srcRecHitEE->find(id);
              if(hgcHitEE != srcRecHitEE->end())
              {
                DetId detId_recHit_ee = hgcHitEE->detid();
                if(id==detId_recHit_ee.rawId())
                {
                  sumEnergyToNorm += hgcHitEE->energy()*it_haf->second;
	          axisX += hgcHitEE->energy()*it_haf->second * recHitTools.getPosition(id).x();
                  axisY += hgcHitEE->energy()*it_haf->second * recHitTools.getPosition(id).y();
                  axisZ += hgcHitEE->energy()*it_haf->second * recHitTools.getPosition(id).z();
                  //std::cout << "energy = " << energyFraction << " x = " << x << " y = " << y << " z = " << z << std::endl;
                }
              }
            }
          }
        } 
        showerAxis = GlobalPoint(axisX/sumEnergyToNorm, axisY/sumEnergyToNorm, (axisZ - vertex_z)/sumEnergyToNorm);
        /*std::cout << "showerAxis.x() = " << showerAxis.x() << std::endl;
        std::cout << "showerAxis.y() = " << showerAxis.y() << std::endl;
        std::cout << "showerAxis.z() = " << showerAxis.z() << std::endl;
        std::cout << "showerAxis.eta() = " << showerAxis.eta() << std::endl;
        std::cout << "showerAxis.phi() = " << showerAxis.phi() << std::endl;
	*/
        for(unsigned int l=0; l<pfCluster->size(); l++) // Iterating over sim clusters
        {
          for(unsigned ifrac = 0; ifrac < pfCluster->at(l).hitsAndFractions().size(); ifrac++)
          {
            uint32_t id = pfCluster->at(l).hitsAndFractions().at(ifrac).first;
            edm::SortedCollection<HGCRecHit>::const_iterator hgcHitEE = srcRecHitEE->find(id);
            if(hgcHitEE != srcRecHitEE->end())
            {
              const HGCalDetId detId_recHit_ee = hgcHitEE->detid();
              //std::cout << "id recHit = " << id << std::endl; 
              //DetId detId_recHit_ee = hgcHitEE->detid(); 
              if(id==detId_recHit_ee.rawId())
              {
                //std::cout << "id recHit = " << id << std::endl;
              }
            }
          }
        }
        
        for(unsigned int l=0; l<pfCluster->size(); l++) // Iterating over sim clusters
        {
          for(unsigned ifrac = 0; ifrac < pfCluster->at(l).hitsAndFractions().size(); ifrac++)
          {
            uint32_t id = pfCluster->at(l).hitsAndFractions().at(ifrac).first;
            edm::SortedCollection<HGCRecHit>::const_iterator hgcHitHEF = srcRecHitHEF->find(id);
            if(hgcHitHEF != srcRecHitHEF->end())
            {
              const HGCalDetId detId_recHit_hef = hgcHitHEF->detid();
              //std::cout << "id recHit = " << id << std::endl;
              //DetId detId_recHit_hef = hgcHitHEF->detid(); 
              if(id==detId_recHit_hef.rawId())
              {
                //std::cout << "id recHit = " << id << std::endl;
              }
            }
          }
        } 
        //EE
	/*
	for(unsigned int l=0; l<pfCluster->size(); l++) // Iterating over sim clusters
	{
		//for (std::vector<CaloParticle>::const_iterator it_caloPart = caloParticle->begin(); it_caloPart != caloParticle->end(); ++it_caloPart)
		//{
                        //const SimClusterRefVector simClusterRefVector = it_caloPart->simClusters();
                        //if(simClusterRefVector.size() > 1) continue;
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
						const GlobalPoint& recHitPos_ee = geoHandle_ee->getPosition(id); 
						recHit_phi.push_back(recHitTools.getPhi(recHitPos_ee));
						recHit_eta.push_back(recHitTools.getEta(recHitPos_ee, vertex_z));
						recHit_pt.push_back(recHitTools.getPt(recHitPos_ee, hgcHitEE->energy(), vertex_z));
						recHit_x.push_back(recHitPos_ee.x());
						recHit_y.push_back(recHitPos_ee.y());
						recHit_z.push_back(recHitPos_ee.z());
						recHit_energy.push_back(hgcHitEE->energy());
                                                unsigned int layer = recHitTools.getLayerWithOffset(id);
                                                recHit_id.push_back(id);
						auto ddd = get_ddd(caloGeom->getSubdetectorGeometry(detId_recHit_ee));
						int thick = ddd->waferTypeL(detId_recHit_ee.wafer()) - 1;
						//unsigned int layer = detId_recHit_ee.layer();
						recHit_layer.push_back(layer);
						int sectionType;
						sectionType = 2;
						if(layer < 29) sectionType = 0;
						else if(layer < 41) sectionType = 1;
						double sigmaNoiseMIP = 0.0;
						sigmaNoiseMIP = noiseMIP;
						if(sectionType != 2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];
						double energyMIP;
						double energy = hgcHitEE->energy()*energyDegradation[thick];
						if(sectionType == 2) energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV );
						else energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV);
						if(energyMIP > 3.) 
						{
							double SoverN = energyMIP / sigmaNoiseMIP;
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
							recHit_energyMIP.push_back(energyMIP); 
						}
						else 
						{
							recHit_time.push_back(-99.0);
							recHit_smearedTime.push_back(-99.0);
							recHit_soverN.push_back(-99.0);
							recHit_rmsTime.push_back(-99.0);
							recHit_deltaTime.push_back(-99.0);
							recHit_energyMIP.push_back(-99.0);
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
						unsigned int layer = recHitTools.getLayerWithOffset(id);//alternate way to access the layer
						recHit_phi.push_back(recHitTools.getPhi(recHitPos_hef));
						recHit_eta.push_back(recHitTools.getEta(recHitPos_hef, vertex_z));
						recHit_pt.push_back(recHitTools.getPt(recHitPos_hef, hgcHitHEF->energy(), vertex_z));
						recHit_x.push_back(recHitPos_hef.x());
						recHit_y.push_back(recHitPos_hef.y());
						recHit_z.push_back(recHitPos_hef.z());
                                                recHit_energy.push_back(hgcHitHEF->energy());
						recHit_id.push_back(id);
						auto ddd = get_ddd(caloGeom->getSubdetectorGeometry(detId_recHit_hef));
						int thick = ddd->waferTypeL(detId_recHit_hef.wafer()) - 1;
						int thick2 = recHitTools.getSiThickness(detId_recHit_hef) / 100. - 1.;
						if(thick!=thick2) std::cout << "thick = " << thick << ", thick2 = " << thick2 << std::endl;
						recHit_layer.push_back(layer);
						int sectionType;
						sectionType = 2;
						if(layer < 29) sectionType = 0;
						else if(layer < 41) sectionType = 1;
						double sigmaNoiseMIP = 0.0;
						sigmaNoiseMIP = noiseMIP;
						if(sectionType !=2) sigmaNoiseMIP = noisefC[sectionType]/fCPerMIP[thick];
						double energyMIP;
						double energy = hgcHitHEF->energy()*energyDegradation[thick];
						if(sectionType == 2) energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV );
						else energyMIP = energy/scaleCorrection[thick]/keV2GeV / (weights.at(layer)/keV2MeV);
						if(energyMIP > 3.)
						{
							double SoverN = energyMIP / sigmaNoiseMIP;
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
							recHit_energyMIP.push_back(energyMIP);
						}
						else
						{
							recHit_time.push_back(-99.0);
							recHit_smearedTime.push_back(-99.0);
							recHit_soverN.push_back(-99.0);
							recHit_rmsTime.push_back(-99.0);
							recHit_deltaTime.push_back(-99.0);
							recHit_energyMIP.push_back(-99.0);
						}
					}
				}
			}
	        //}//eta matching
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
					//std::cout << "hgcHitBH->time() = " << hgcHitBH->time() << std::endl;
				}
			}
		}  
	}*/

	tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void HGCTimingAnalyzerWithTOAPU::beginJob()
{
	edm::Service<TFileService> fs;
	tree_=fs->make<TTree>("HGCTiming","HGCTiming");
	branch_=tree_->Branch("run",   &run_,   "run/I");
	branch_=tree_->Branch("event", &event_, "event/I");
	branch_=tree_->Branch("lumi",  &lumi_,  "lumi/I");
	branch_=tree_->Branch("reachedEE", &reachedEE_, "reachedEE/I");
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
	branch_=tree_->Branch("recHit_energyMIP", &recHit_energyMIP);
	branch_=tree_->Branch("recHit_id", &recHit_id);
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
HGCTimingAnalyzerWithTOAPU::endJob() 
{
	std::cout << "nhits_ee = " << nhits_ee << std::endl;
	std::cout << "nhits_hf = " << nhits_hf << std::endl;
	std::cout << "nhits_bh = " << nhits_bh << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HGCTimingAnalyzerWithTOAPU::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

const HGCalDDDConstants* HGCTimingAnalyzerWithTOAPU::get_ddd(const CaloSubdetectorGeometry* geom){
	const HGCalGeometry* hg = static_cast<const HGCalGeometry*>(geom);
	const HGCalDDDConstants* ddd = &(hg->topology().dddConstants());
	if( nullptr == ddd ) {
		throw cms::Exception("hgcal::RecHitTools")
			<< "DDDConstants not accessibl to hgcal::RecHitTools!";
	}
	return ddd;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HGCTimingAnalyzerWithTOAPU);
