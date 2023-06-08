import FWCore.ParameterSet.Config as cms

from genvecanalyzerQCD_cfg import process, options

process.ak8GenJetsNoNu.rParam = 1.5

process.TFileService.fileName = cms.string("genvecAK15_{}_n-all.root".format(options.sample))

