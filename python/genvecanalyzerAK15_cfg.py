import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.load("Configuration.StandardSequences.Services_cff")
process.load("SVJ.Production.GenVecAnalyzer_cfi")

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets as _ak4GenJets
process.ak15GenJetsNoNu = _ak4GenJets.clone(
    rParam = 1.5,
    src = "genParticlesForJetsNoNu",
)

process.GenVecAnalyzer.JetTag = cms.InputTag("ak15GenJetsNoNu")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("genvecanalysis.root")
)

process.p1 = cms.Path(process.ak15GenJetsNoNu*process.GenVecAnalyzer)
