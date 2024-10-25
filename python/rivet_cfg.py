# based on https://twiki.cern.ch/twiki/bin/view/CMS/RivetontoAODSIM
import FWCore.ParameterSet.Config as cms
process = cms.Process("runRivetAnalysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.source = cms.Source("PoolSource",  fileNames = cms.untracked.vstring('file:AODSIM.root') )
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.HepMCCollection = cms.InputTag("generator","unsmeared")
process.rivetAnalyzer.AnalysisNames = cms.vstring('ATLAS_2023_I2663256')
process.p = cms.Path(process.rivetAnalyzer)
