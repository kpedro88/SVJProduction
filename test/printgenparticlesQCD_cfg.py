import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing("analysis")
options.register("sample", "", VarParsing.multiplicity.singleton, VarParsing.varType.string)
options.parseArguments()

process = cms.Process("demo")

process.load("Configuration.StandardSequences.Services_cff")
process.load("SVJ.Production.PrintGenParticles_cfi")

process.PrintGenParticles.PartTag = cms.InputTag("packedGenParticlesForJetsNoNu")
process.PrintGenParticles.fileName = cms.string("vecs_{}_n-all".format(options.sample))

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedGenParticles"),
    cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
)
with open('CMS_MonteCarlo2016_RunIISummer16MiniAODv2_{}_MINIAODSIM_PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_file_index.txt'.format(options.sample),'r') as infile:
	process.source.fileNames = [line.rstrip() for line in infile]

process.p1 = cms.Path(process.packedGenParticlesForJetsNoNu*process.PrintGenParticles)
