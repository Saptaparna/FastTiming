# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TTbar_14TeV_TuneCP5_cfi --conditions auto:phase2_realistic_T15 -n 10 --era Phase2C9 --eventcontent FEVTDEBUG -s GEN,SIM,DIGI:pdigi_valid,L1,L1TrackTrigger,DIGI2RAW,HLT:@fake2 --datatier GEN-SIM --beamspot NoSmear --geometry Extended2026D49 --pileup AVE_200_BX_25ns --procModifier fineCalo --no_exec --python_filename=GSDfineCalo_fragment.py
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.ProcessModifiers.fineCalo_cff import fineCalo

process = cms.Process('HLT',Phase2C9, fineCalo)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
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
    input = cms.untracked.int32(NEVENT),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

process.RandomNumberGeneratorService.generator.initialSeed = cms.untracked.uint32(NPART+NEVENT)
process.RandomNumberGeneratorService.VtxSmeared.initialSeed = cms.untracked.uint32(NPART+NEVENT)
process.RandomNumberGeneratorService.mix.initialSeed = cms.untracked.uint32(NPART+NEVENT)

# Input source
process.source = cms.Source("EmptySource")
process.source.firstLuminosityBlock = cms.untracked.uint32(NPART)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('singleGamma_TuneCP5_cfi nevts:10'),
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
    fileName = cms.untracked.string('file:step1_Gamma_nNEVENT_partNPART.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.FEVTDEBUGoutput.outputCommands.append("keep *_*G4*_*_*")
process.FEVTDEBUGoutput.outputCommands.append("keep SimClustersedmAssociation_mix_*_*")
process.FEVTDEBUGoutput.outputCommands.append("keep CaloParticlesedmAssociation_mix_*_*")


# Additional output definition

# Other statements
#process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
#process.mix.bunchspace = cms.int32(25)
#process.mix.minBunch = cms.int32(-3)
#process.mix.maxBunch = cms.int32(3)
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.mix.digitizers = cms.PSet(process.theDigitizersValid)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

#process.generator = cms.EDProducer("CloseByParticleGunProducer",
#    PGunParameters = cms.PSet(PartID = cms.vint32(22), # Since I'm shooting antiparticles, I'm putting gamma and K0 twice in the vector
#        EnMin = cms.double(20.0),
#        EnMax = cms.double(200.0),
#        RMin = cms.double(60),
#        RMax = cms.double(120),
#        ZMin = cms.double(320),#320#-410
#        ZMax = cms.double(321),
#        Delta = cms.double(0),
#        Pointing = cms.bool(True),
#        Overlapping = cms.bool(False),
#        RandomShoot = cms.bool(False),
#        NParticles = cms.int32(1),
#        MaxEta = cms.double(2.7),
#        MinEta = cms.double(1.7),#-2.7
#        MaxPhi = cms.double(3.14159265359),
#        MinPhi = cms.double(-3.14159265359)
#    ),
#    Verbosity = cms.untracked.int32(0),
#    firstRun = cms.untracked.uint32(1),
#    AddAntiParticle = cms.bool(False),
#    psethack = cms.string('Random particles in front of HGCAL')
#)
#

process.generator = cms.EDProducer("FlatEtaRangeGunProducer",
    # particle ids
    particleIDs=cms.vint32(22), # 22, 11,-11,211,-211,13,-13, 310, 130, 111, 311, 321, -321),
    # max number of particles to shoot at a time
    nParticles=cms.int32(1),
    # shoot exactly the particles defined in particleIDs in that order
    exactShoot=cms.bool(False),
    # randomly shoot [1, nParticles] particles, each time randomly drawn from particleIDs
    randomShoot=cms.bool(False),
    # energy range
    eMin=cms.double(20),
    eMax=cms.double(200.0),
    # phi range
    phiMin=cms.double(-3.14159265359),
    phiMax=cms.double(3.14159265359),
    # eta range
    etaMin=cms.double(1.52),
    etaMax=cms.double(3.00),

    debug=cms.untracked.bool(True),
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi_valid)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.L1TrackTrigger_step = cms.Path(process.L1TrackTrigger)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.L1TrackTrigger_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.FEVTDEBUGoutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.ProductionFilterSequence)

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
