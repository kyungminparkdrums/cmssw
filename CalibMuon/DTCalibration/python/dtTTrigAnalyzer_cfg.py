import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("DTTTrigAnalyzer",eras.Run3)

process.load("CondCore.CondDB.CondDB_cfi")

process.source = cms.Source("EmptySource",
    numberEventsInRun = cms.untracked.uint32(1),
    firstRun = cms.untracked.uint32(1)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.dtTTrigAnalyzer = cms.EDAnalyzer("DTTTrigAnalyzer",
    dbLabel = cms.untracked.string(""),
    rootFileName = cms.untracked.string("") 
)

process.p = cms.Path(process.dtTTrigAnalyzer)
