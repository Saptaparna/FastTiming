# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCUETP8M1_cfi --conditions auto:phase2_realistic -n 10 --era Phase2C4 --eventcontent FEVTDEBUG --relval 9000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

#process = cms.Process('SIM',eras.Phase2C4)
#process = cms.Process('SIM',eras.Phase2C4_timing)
process = cms.Process('SIM',eras.Phase2_timing)

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgchefrontDigitizer
hgchefrontDigitizer.tofDelay = cms.double(3.)
#hgchefrontDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(5., 5., 5.)
#hgchefrontDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02)
#hgchefrontDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(3.75, 7.71, 11.64)
hgchefrontDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(25., 25., 25.)
hgchefrontDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.0004, 0.0004, 0.0004)
hgchefrontDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(12.0, 12.0, 12.0)
hgchefrontDigitizer.digiCfg.feCfg.toaLSB_ns = cms.double(0.0244)

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer
hgceeDigitizer.tofDelay = cms.double(3.)
#hgceeDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(5., 5., 5.)
#hgceeDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02)
#hgceeDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(3.75, 7.71, 11.64)
hgceeDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(25., 25., 25.)
hgceeDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.0004, 0.0004, 0.0004)
hgceeDigitizer.digiCfg.feCfg.tdcForToAOnset_fC  = cms.vdouble(12.0, 12.0, 12.0)
hgceeDigitizer.digiCfg.feCfg.toaLSB_ns = cms.double(0.0244)
#cms.vdouble(1.25,2.57,3.88)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D28_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
#process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC14TeV_cfi')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1)
    input = cms.untracked.int32(1000)
)

# random seeds
process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(1+1000)
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = cms.untracked.uint32(1+1000)

# Input source
process.source = cms.Source("EmptySource")
process.source.firstLuminosityBlock = cms.untracked.uint32(1)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('MultipleParticleGun_cfi'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1_Gamma_Pt10_n1000_part1_directional_D17File.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    AddAntiParticle = cms.bool(False),
    PGunParameters = cms.PSet(
        PartID = cms.vint32(22),
        MinPhi = cms.double(0.0),
        MaxPhi = cms.double(0.0),
        MinPt = cms.double(10),
        MaxPt = cms.double(10),
        MinEta = cms.double(1.79),
        MaxEta = cms.double(1.81)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('multiple particles predefined pT/E eta 1p479 to 3')
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 


# Customisation from command line
#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(4)
process.options.numberOfStreams=cms.untracked.uint32(0)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
