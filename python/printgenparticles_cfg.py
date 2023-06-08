import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")

process.load("Configuration.StandardSequences.Services_cff")
process.load("SVJ.Production.PrintGenParticles_cfi")
# fileName is set directly in process.PrintGenParticles

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.p1 = cms.Path(process.PrintGenParticles)
