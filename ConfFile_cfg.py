import FWCore.ParameterSet.Config as cms

process = cms.Process("HGCTiming")

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuonReco_cff')
process.load('Configuration.Geometry.GeometryExtended2023HGCalMuon_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/s/sapta/HighGranularityCalorimeter/TimingStudies/CMSSW_6_2_0_SLHC26_patch3/src/UserCode/HGCanalysis/Events_22_8_57.root'
    )
)

process.hgctiming = cms.EDAnalyzer('HGCTimingAnalyzer')

process.p = cms.Path(process.hgctiming)
