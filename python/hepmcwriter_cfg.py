import FWCore.ParameterSet.Config as cms

process = cms.Process("write")

process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)

process.writer = cms.EDAnalyzer("HepMCEventWriterNew",
    hepMCProduct = cms.InputTag('generator','unsmeared'),
    fileName = cms.string("GenEvent_ASCII.dat"),
)
process.writer_path = cms.EndPath(process.writer)

# Path and EndPath definitions
process.write_step = cms.Path(process.writer)
