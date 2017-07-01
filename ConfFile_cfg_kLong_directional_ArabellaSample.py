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
'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/GammaTime_step3/step3_Gamma_Pt10_n1000_part1_directional.root',
'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/GammaTime_step3/step3_Gamma_Pt10_n1000_part2_directional.root',
'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/GammaTime_step3/step3_Gamma_Pt10_n1000_part3_directional.root',
'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/GammaTime_step3/step3_Gamma_Pt10_n1000_part4_directional.root'
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part1_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part2_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part3_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part4_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part5_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part6_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part7_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part8_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part9_directional.root',
#'file:/uscms_data/d1/sapta/work/HighGranularityCalorimeter/TimingStudies_9X/CMSSW_9_1_0_pre3/src/kLongTime_step3/step3_kLong_Pt10_n1000_part10_directional.root'
   )
)

process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCTiming_Gamma_Pt10_SoverN1000ps_Floor20ps_EE_FH_Test.root'))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('HGCTiming_kLong_Pt10_SoverN1000ps_Floor20ps_EE_FH_Test.root'))

process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.hgctiming = cms.EDAnalyzer('HGCTimingAnalyzerWithTOA',
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
