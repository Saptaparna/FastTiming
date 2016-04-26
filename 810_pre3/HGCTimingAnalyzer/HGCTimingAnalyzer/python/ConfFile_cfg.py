import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCTiming")

process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part1.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part2.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part3.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part4.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part5.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part6.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part7.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part8.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part9.root',
        'file:/uscms_data/d2/sapta/work/HighGranularityCalorimeter/TimingStudies_8X/CMSSW_8_1_0_pre3/src/RecoLocalCalo/HGCalRecProducers/test/SinglePhotonGunPt50_n1000_part10.root'
    )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCTiming_SinglePhoton_Pt50GeV_1.root'))

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.hgctiming = cms.EDAnalyzer('HGCTimingAnalyzer',
                                   srcGenParticles = cms.InputTag('genParticles'),
                                   srcSimTracks = cms.InputTag('g4SimHits'),
                                   srcSimVertices = cms.InputTag('g4SimHits'),
                                   srcHepmcevent  = cms.InputTag('generator', 'unsmeared'),
                                   srcPFRecHit = cms.InputTag('particleFlowRecHitHGC', '')
                                  )
                                   
process.p = cms.Path(process.hgctiming)
