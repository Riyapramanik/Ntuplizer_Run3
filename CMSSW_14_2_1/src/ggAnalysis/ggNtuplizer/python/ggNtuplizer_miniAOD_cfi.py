import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import *
# add by me                                                                                                                                                                    
from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import cleanPatJets 

ggNtuplizer = cms.EDAnalyzer("ggNtuplizer",
                             doGenParticles       = cms.bool(True),
                             runOnParticleGun     = cms.bool(False),
                             runOnSherpa          = cms.bool(False),
                             runL1ECALPrefire     = cms.bool(False),
                             dumpPFPhotons        = cms.bool(True), 
                             dumpPDFSystWeight    = cms.bool(False),
                             development          = cms.bool(False),
                             addFilterInfoMINIAOD = cms.bool(False),
                             store_electrons = cms.untracked.bool(True),
                             store_muons = cms.untracked.bool(True),
                             store_photons = cms.untracked.bool(True),
                             store_ak4jets = cms.untracked.bool(True),
                             store_CHS_met = cms.untracked.bool(True),
                             store_PUPPI_met = cms.untracked.bool(True),
                             year                 = cms.int32(2022),
                             
                             patTriggerResults    = cms.InputTag("TriggerResults", "", "RECO"),
                             genParticleSrc       = cms.InputTag("prunedGenParticles"),
                             generatorLabel       = cms.InputTag("generator"),
                             LHEEventLabel        = cms.InputTag("externalLHEProducer"),
                             pileupCollection     = cms.InputTag("slimmedAddPileupInfo"),
                             VtxLabel             = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             rhoLabel             = cms.InputTag("fixedGridRhoFastjetAll"),
                             electronSrc          = cms.InputTag("slimmedElectrons"),
                             photonSrc            = cms.InputTag("slimmedPhotons"),
                             Muons              = cms.InputTag("slimmedMuons"),
                             gsfTrackSrc          = cms.InputTag("reducedEgamma", "reducedGsfTracks"),
                             reducedEcalRecHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                             reducedEcalRecHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                             reducedEcalRecHitsES = cms.InputTag("reducedEgamma", "reducedESRecHits"),
                             
                             PFJetsAK4                 = cms.InputTag("slimmedJets"),
                             GENJetAK4            = cms.InputTag("slimmedGenJets"),
                             PFRho                = cms.InputTag("fixedGridRhoFastjetAll"),
                             packedPFCands             = cms.InputTag("packedPFCandidates"),
                             ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),
                             phoWP80MapToken = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIIWinter22-v1-wp80"),
                             phoWP90MapToken = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIIWinter22-v1-wp90"),
                             minJetPt = cms.untracked.double(25.),
                             maxEta = cms.untracked.double(3.),
                             
                             PFJetAK4Collection = cms.InputTag("slimmedJetsPuppi"),
                             #MET
                             PFMet = cms.InputTag("slimmedMETsUpdated"),          #"Updated" comes from postfix in runMetCorAndUncFromMiniAOD
                             PuppiMet = cms.InputTag("slimmedMETsPuppi"),  #"Updated" comes from postfix in runMetCorAndUncFromMiniAOD
                             MET_Filters = cms.InputTag("TriggerResults","","RECO"),
                             
                             #trigger parameters
                             bits = cms.InputTag("TriggerResults","","HLT"),
                             prescales = cms.InputTag("patTrigger","","RECO"),
                             TriggerObjects = cms.InputTag("slimmedPatTrigger"),
                            

                             # JEC/JER file parameters:
                             jecL1FastFileAK4 = cms.string(""),
                             jecL2RelativeFileAK4 = cms.string(""),
                             jecL3AbsoluteFileAK4 = cms.string(""),
                             jecL2L3ResidualFileAK4 = cms.string(""),
                             JECUncFileAK4 = cms.string(""),
                             PtResoFileAK4 = cms.string(""),
                             PtSFFileAK4 = cms.string(""),
                             JetVetoMap = cms.string(""),
                             
                            # EGM correction parameters:
                             dataYear = cms.int32(2022),
                             dataPeriod = cms.string("C"),
                             isData = cms.bool(True),
                             applyEGMCorrections = cms.bool(False),

                             YEAR = cms.untracked.string("2023"),
                             Data = cms.untracked.bool(True),
                             MonteCarlo = cms.untracked.bool(False),
                             UltraLegacy = cms.untracked.bool(False),
                             isRun3 = cms.untracked.bool(True),
                             
                             
)
