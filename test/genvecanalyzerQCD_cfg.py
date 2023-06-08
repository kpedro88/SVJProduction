import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing("analysis")
options.register("sample", "", VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.parseArguments()

process = cms.Process("demo")

process.load("Configuration.StandardSequences.Services_cff")
process.load("SVJ.Production.GenVecAnalyzer_cfi")

process.GenVecAnalyzer.signal = False

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedGenParticles"),
    cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
)

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets as _ak4GenJets
process.ak8GenJetsNoNu = _ak4GenJets.clone(
    rParam = 0.8,
    src = "packedGenParticlesForJetsNoNu",
)

process.GenVecAnalyzer.JetTag = cms.InputTag("ak8GenJetsNoNu")

process.genParticlesForMETAllVisible = cms.EDProducer("InputGenJetsParticleSelector",
    src = cms.InputTag("packedGenParticles"),
    partonicFinalState = cms.bool(False),
    excludeResonances = cms.bool(False),
    excludeFromResonancePids = cms.vuint32(),
    tausAsJets = cms.bool(False),

    ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 12, 14, 16
    )
)

process.genMetTrue = cms.EDProducer("GenMETProducer",
    src = cms.InputTag("genParticlesForMETAllVisible"), ## Input  product label       
    alias = cms.string('GenMETAllVisible'), ## Alias  for FWLite          
    onlyFiducialParticles = cms.bool(False), ## Use only fiducial GenParticles
    globalThreshold = cms.double(0.0), ## Global Threshold for input objects
    usePt   = cms.bool(True), ## using Pt instead Et
    applyFiducialThresholdForFractions   = cms.bool(False),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)
with open('CMS_MonteCarlo2016_RunIISummer16MiniAODv2_{}_MINIAODSIM_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_file_index.txt'.format(options.sample),'r') as infile:
	process.source.fileNames = [line.rstrip() for line in infile]

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("genvec_{}_n-all.root".format(options.sample))
)

process.p1 = cms.Path(process.packedGenParticlesForJetsNoNu*process.ak8GenJetsNoNu*process.genParticlesForMETAllVisible*process.genMetTrue*process.GenVecAnalyzer)
