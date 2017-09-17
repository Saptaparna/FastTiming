import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as cms


from RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi import dEdX_weights, HGCalRecHit
from RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi import HGCalUncalibRecHit
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer, hgchefrontDigitizer, hgchebackDigitizer


process = cms.Process("HGCTimingWithTOA")

process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D13_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
# get uncalibrechits with weights method
process.load("RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi")

# get rechits e.g. from the weights
process.load("RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi")


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'root://cmsxrootd.fnal.gov//store/user/sapta/TimingStudies_9X/kLongTime_step3_PU/step3_kLong_Pt10_n40_part1_directional_PU.root',
    'root://cmsxrootd.fnal.gov//store/user/sapta/TimingStudies_9X/kLongTime_step3_PU/step3_kLong_Pt10_n40_part2_directional_PU.root',
   )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCTiming_kLong_Pt10_SoverN5000ps_Floor30ps_PU_Part1.root'))

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.hgctiming = cms.EDAnalyzer('HGCTimingAnalyzerWithTOA',
                                   #SlopeParameterA = cms.double(1.0),
                                   #SlopeParameterA = cms.double(4.0),
                                   SlopeParameterA = cms.double(5.0),
                                   #FloorParameterC = cms.double(0.020),
                                   FloorParameterC = cms.double(0.030),
                                   HGCEE_keV2fC  = hgceeDigitizer.digiCfg.keV2fC,
                                   HGCHEF_keV2fC = hgchefrontDigitizer.digiCfg.keV2fC,
                                   HGCHB_keV2MIP = hgchebackDigitizer.digiCfg.keV2MIP,
                                   dEdXweights = cms.vdouble(dEdX_weights),
                                   thicknessCorrection = cms.vdouble(HGCalRecHit.thicknessCorrection),
                                   HGCEE_fCPerMIP = cms.vdouble(HGCalUncalibRecHit.HGCEEConfig.fCPerMIP),
                                   HGCEE_noisefC = cms.vdouble(hgceeDigitizer.digiCfg.noise_fC),
                                   HGCEF_noisefC = cms.vdouble(hgchefrontDigitizer.digiCfg.noise_fC),
                                   HGCBH_noiseMIP = hgchebackDigitizer.digiCfg.noise_MIP, 
                                   srcGenParticles = cms.InputTag('genParticles'),
                                   srcSimTracks = cms.InputTag('g4SimHits'),
                                   srcSimVertices = cms.InputTag('g4SimHits'),
                                   srcPFRecHit = cms.InputTag('particleFlowRecHitHGC', 'Cleaned'),
                                   srcPFCluster = cms.InputTag('particleFlowClusterHGCal'),
                                   srcRecHitEE = cms.InputTag('HGCalRecHit', 'HGCEERecHits'),
                                   srcRecHitHEF = cms.InputTag('HGCalRecHit', 'HGCHEFRecHits'),
                                   srcRecHitBH = cms.InputTag('HGCalRecHit', 'HGCHEBRecHits'),
                                   srcCaloParticle = cms.InputTag('mix', 'MergedCaloTruth'),
                                   srcPartHandle = cms.InputTag('mix','MergedTrackTruth'),
                                   TriggerGeometry = cms.PSet(
                                       TriggerGeometryName = cms.string('HGCalTriggerGeometryImp1'),
                                       L1TCellsMapping = cms.FileInPath("L1Trigger/L1THGCal/data/cellsToTriggerCellsMap.txt"),
                                       ),
                                  )

process.p = cms.Path(process.hgctiming)
