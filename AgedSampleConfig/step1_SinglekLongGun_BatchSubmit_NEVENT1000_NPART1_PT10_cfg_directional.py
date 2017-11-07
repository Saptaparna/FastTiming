# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: SingleGammaPt35_cfi --conditions auto:run2_mc -n 10 --era Phase2C2 --eventcontent FEVTDEBUGHLT --relval 10000,100 -s GEN,SIM,DIGI:pdigi_valid,L1,DIGI2RAW,HLT:@fake --datatier GEN-SIM-DIGI-RAW --beamspot Realistic50ns13TeVCollision --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023tilted --geometry Extended2023D3 --no_exec --fileout file:step1_DIGIRAW.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Phase2)
from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgchefrontDigitizer
hgchefrontDigitizer.tofDelay = cms.double(3.)
#hgchefrontDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(5., 5., 5.)
#hgchefrontDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02)
#hgchefrontDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(3.75, 7.71, 11.64)
hgchefrontDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(25., 25., 25.)
hgchefrontDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.0004, 0.0004, 0.0004)
hgchefrontDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(12.0, 12.0, 12.0)
hgchefrontDigitizer.digiCfg.feCfg.toaLSB_ns = cms.double(0.0244)
hgchefrontDigitizer.digiCfg.chargeCollectionEfficiencies = cms.vdouble(0.5, 0.5, 0.7)
hgchefrontDigitizer.digiCfg.noise_fC = cms.vdouble(0.384, 0.360, 0.280)

from SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi import hgceeDigitizer
hgceeDigitizer.tofDelay = cms.double(3.)
#hgceeDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(5., 5., 5.)
#hgceeDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.02, 0.02, 0.02)
#hgceeDigitizer.digiCfg.feCfg.tdcForToAOnset_fC = cms.vdouble(3.75, 7.71, 11.64)
hgceeDigitizer.digiCfg.feCfg.jitterNoise_ns = cms.vdouble(25., 25., 25.)
hgceeDigitizer.digiCfg.feCfg.jitterConstant_ns = cms.vdouble(0.0004, 0.0004, 0.0004)
hgceeDigitizer.digiCfg.feCfg.tdcForToAOnset_fC  = cms.vdouble(12.0, 12.0, 12.0)
hgceeDigitizer.digiCfg.feCfg.toaLSB_ns = cms.double(0.0244)
hgceeDigitizer.digiCfg.chargeCollectionEfficiencies = cms.vdouble(0.5, 0.5, 0.7)
hgceeDigitizer.digiCfg.noise_fC = cms.vdouble(0.384, 0.360, 0.280)
#cms.vdouble(1.25,2.57,3.88)


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.VtxSmearedNoSmear_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.L1TrackTrigger_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_Fake2_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
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

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
                                              SelectEvents = cms.untracked.PSet(
                                                                                SelectEvents = cms.vstring('generation_step')
                                                                                ),
                                              dataset = cms.untracked.PSet(
                                                                           dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
                                                                           filterName = cms.untracked.string('')
                                                                           ),
                                              eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
                                              fileName = cms.untracked.string('file:step1_kLong_Pt10_n1000_part1_directional.root'),
                                              outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
                                              splitLevel = cms.untracked.int32(0)
                                              )

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")

process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

gunmode = 'default'

if gunmode == 'default':
    process.generator = cms.EDProducer("FlatRandomPtGunProducer",
                                       AddAntiParticle = cms.bool(False),
                                       PGunParameters = cms.PSet(
                                                                 MaxEta = cms.double(1.81),
                                                                 MaxPhi = cms.double(0.0),
                                                                 MaxPt = cms.double(10),
                                                                 MinEta = cms.double(1.79),
                                                                 MinPhi = cms.double(0.0),
                                                                 MinPt = cms.double(10),
                                                                 PartID = cms.vint32(130)
                                                                 ),
                                       Verbosity = cms.untracked.int32(0),
                                       psethack = cms.string('multiple particles predefined pT/E eta 1p479 to 3')
                                       )
elif gunmode == 'pythia8':
    process.generator = cms.EDFilter("FlatRandomPtGunProducer",
                                     maxEventsToPrint = cms.untracked.int32(1),
                                     pythiaPylistVerbosity = cms.untracked.int32(1),
                                     pythiaHepMCVerbosity = cms.untracked.bool(True),
                                     PGunParameters = cms.PSet(
                                                               ParticleID = cms.vint32(130),
                                                               AddAntiParticle = cms.bool(False),
                                                               MinPhi = cms.double(0.0),
                                                               MaxPhi = cms.double(0.0),
                                                               MinPt = cms.double(10),
                                                               MaxPt = cms.double(10),
                                                               MinEta = cms.double(1.79),
                                                               MaxEta = cms.double(1.81)
                                                               ),
                                     PythiaParameters = cms.PSet(parameterSets = cms.vstring())
                                     )

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
# filter all path with the production filter sequence
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq

#Setup FWK for multithreaded
#process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
