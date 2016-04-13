#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("LHE")

#@#
names_of_files = ''
for i in range(200):
	names_of_files += 'file:unweighted_events_%s.lhe'%i
#@#

#@#process.source = cms.Source("LHESource",
#@#	fileNames = cms.untracked.vstring('file:unweighted_events.lhe'
#@#))
process.sorce = cms.Source("LHESource",
	fileNames = cms.untracked.vstring(names_of_files
))

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string('ttbar')
)

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cout = cms.untracked.PSet( threshold = cms.untracked.string('INFO') )

process.LHE = cms.OutputModule("PoolOutputModule",
	dataset = cms.untracked.PSet(dataTier = cms.untracked.string('LHE')),
	fileName = cms.untracked.string('lhe.root')
)

#process.lhedump = cms.EDAnalyzer("DummyLHEAnalyzer",
#                                 src = cms.InputTag("source")
#                                 )

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.outpath = cms.EndPath(process.LHE)
