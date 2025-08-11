import FWCore.ParameterSet.Config as cms
process = cms.Process('ggKit')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
# For Jets
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('RecoJets.JetProducers.PileupJetIDParams_cfi')
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

#===== MET Filters ==

process.load('RecoMET.METFilters.primaryVertexFilter_cfi')
process.primaryVertexFilter.vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilterResultProducerNoMinZ = process.HBHENoiseFilterResultProducer.clone(minZeros = cms.int32(99999))

process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices") 
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonDzFilter_cfi')
process.BadPFMuonDzFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonDzFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonDzFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices")
process.BadPFMuonDzFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
#process.eeBadScFilter.EERecHitSource = cms.InputTag('reducedEgamma','reducedEERecHits')

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

process.load('RecoMET.METFilters.hfNoisyHitsFilter_cfi')

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('eventoutput.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


from RecoJets.Configuration.GenJetParticles_cff import *

from Configuration.AlCa.GlobalTag import GlobalTag
from RecoJets.Configuration.GenJetParticles_cff import *
IsDATA = True
YEAR = "2022"  
ERA = "A"

if IsDATA:
    if YEAR=="2022":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2')
        JEC_tag = "Summer22_22Sep2023_RunCD_V2_DATA"
        JER_tag = 'Summer22_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22_23Sep2023_RunCD_v1'
    elif YEAR=="2022EE":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1')
        JEC_tag = "Summer22EE_22Sep2023_Run"+ERA+"_V2_DATA" 
        JER_tag = 'Summer22EE_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22EE_23Sep2023_RunEFG_v1'
    elif YEAR=="2023":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1')
        if ERA == "Cv4":
            JEC_tag = "Summer23Prompt23_RunCv4_V1_DATA"
        else:
            JEC_tag = "Summer23Prompt23_RunCv123_V1_DATA"
        JER_tag = 'Summer23Prompt23_RunCv123_JRV1_MC'       
        JetVeto_tag = "Summer23Prompt23_RunC_v1"
    elif YEAR=="2023BPiX":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1')
        JEC_tag = "Summer23BPixPrompt23_RunD_V1_DATA" 
        JER_tag = 'Summer23BPixPrompt23_RunD_JRV1_MC'       
        JetVeto_tag = "Summer23BPixPrompt23_RunD_v1"
else:
    if YEAR=="2022":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v5')
        JEC_tag = 'Summer22_22Sep2023_V2_MC'
        JER_tag = 'Summer22_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22_23Sep2023_RunCD_v1'
    elif YEAR=="2022EE":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_postEE_v6')
        JEC_tag = 'Summer22EE_22Sep2023_V2_MC'
        JER_tag = 'Summer22EE_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22EE_23Sep2023_RunEFG_v1'
    elif YEAR=="2023":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14')
        JEC_tag = "Summer23Prompt23_V1_MC"
        if ERA=="Cv4":
            JER_tag = 'Summer23Prompt23_RunCv4_JRV1_MC'
        else:
            JER_tag = 'Summer23Prompt23_RunCv123_JRV1_MC'   
        JetVeto_tag = "Summer23Prompt23_RunC_v1"
    elif YEAR=="2023BPiX":
        process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_postBPix_v2')
        JEC_tag = "Summer23BPixPrompt23_V1_MC"
        JER_tag = 'Summer23BPixPrompt23_RunD_JRV1_MC'     
        JetVeto_tag = "Summer23BPixPrompt23_RunD_v1"


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                               '/store/data/Run2022A/EGamma/MINIAOD/22Sep2023-v1/30000/19531535-6c96-4871-aa56-8cf0b9581e6d.root'                               
        ))

process.cleanJets = cms.EDProducer("PATJetCleaner",
                                   src = cms.InputTag("selectedUpdatedPatJetsWithPNetInfo"),
                                   preselection = cms.string(''),
                                   checkOverlaps = cms.PSet(
                                    muons = cms.PSet(
                                       src          = cms.InputTag("slimmedMuons"),
                                       algorithm    = cms.string("byDeltaR"),
                                       preselection = cms.string(""),
                                       deltaR       = cms.double(0.4),
                                       checkRecoComponents = cms.bool(False),
                                       pairCut             = cms.string(""),
                                       requireNoOverlaps   = cms.bool(True)
                                    ),
                                      electrons = cms.PSet(
                                         src          = cms.InputTag("slimmedElectrons"),
                                         algorithm    = cms.string("byDeltaR"),
                                         preselection = cms.string(""),
                                         deltaR       = cms.double(0.4),
                                         checkRecoComponents = cms.bool(False),
                                         pairCut             = cms.string(""),
                                         requireNoOverlaps   = cms.bool(True)
                                      )
                                   ),
                                   finalCut = cms.string('')
                                   )

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
input_tags = dict()
switchOnVIDPhotonIdProducer(process, dataFormat)
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = [
'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_RunIIIWinter22_122X_V1_cff',
'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Winter22_122X_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_iso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff'
]

for idmod in my_id_modules:
   setupAllVIDIdsInModule(process, idmod, setupVIDElectronSelection)
   setupAllVIDIdsInModule(process, idmod, setupVIDPhotonSelection)

process.TFileService = cms.Service("TFileService", fileName = cms.string('ggtree_data.root'))
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")

# Basic configuration
process.ggNtuplizer.year = cms.int32(2022)
process.ggNtuplizer.dataYear = cms.int32(2022)
process.ggNtuplizer.dataPeriod = cms.string("A")
process.ggNtuplizer.useETDependentCorrections = cms.bool(True)
process.ggNtuplizer.applyEGMCorrections = cms.bool(True)
process.ggNtuplizer.isData = cms.bool(True)

# Data/MC flags
process.ggNtuplizer.doGenParticles = cms.bool(False)  # No gen particles for data
process.ggNtuplizer.runOnParticleGun = cms.bool(False)
process.ggNtuplizer.runOnSherpa = cms.bool(False)
process.ggNtuplizer.runOnMC = cms.bool(False)  # Important for data

# Control what to dump
process.ggNtuplizer.dumpPFPhotons = cms.bool(True)
process.ggNtuplizer.dumpJets = cms.bool(True)
process.ggNtuplizer.dumpAK8Jets = cms.bool(False)
process.ggNtuplizer.dumpTaus = cms.bool(False)
process.ggNtuplizer.dumpHFElectrons = cms.bool(False)
process.ggNtuplizer.dumpPDFSystWeight = cms.bool(False)
process.ggNtuplizer.dumpSoftDrop = cms.bool(True)

# Flags
process.ggNtuplizer.runL1ECALPrefire = cms.bool(False)
process.ggNtuplizer.development = cms.bool(False)
process.ggNtuplizer.addFilterInfoMINIAOD = cms.bool(False)
process.ggNtuplizer.doNoHFMET = cms.bool(False)

# Trigger cuts
process.ggNtuplizer.trgFilterDeltaPtCut = cms.double(0.5)
process.ggNtuplizer.trgFilterDeltaRCut = cms.double(0.3)

# Input collections - updated for 2022 MINIAOD
process.ggNtuplizer.triggerEvent = cms.InputTag("slimmedPatTrigger", "", "")
process.ggNtuplizer.triggerResults = cms.InputTag("TriggerResults", "", "HLT")
process.ggNtuplizer.patTriggerResults = cms.InputTag("TriggerResults", "", "RECO")
process.ggNtuplizer.genParticleSrc = cms.InputTag("prunedGenParticles")
process.ggNtuplizer.generatorLabel = cms.InputTag("generator")
process.ggNtuplizer.LHEEventLabel = cms.InputTag("externalLHEProducer")
process.ggNtuplizer.pileupCollection = cms.InputTag("slimmedAddPileupInfo")
process.ggNtuplizer.VtxLabel = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ggNtuplizer.VtxBSLabel = cms.InputTag("offlinePrimaryVerticesWithBS")

# Rho and MET
process.ggNtuplizer.rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
process.ggNtuplizer.rhoCentralLabel = cms.InputTag("fixedGridRhoFastjetCentralNeutral")
process.ggNtuplizer.PFRho = cms.InputTag("fixedGridRhoFastjetAll")
process.ggNtuplizer.pfMETLabel = cms.InputTag("slimmedMETs")  # Updated for MINIAOD
process.ggNtuplizer.puppiMETLabel = cms.InputTag("slimmedMETsPuppi")

# Lepton collections
process.ggNtuplizer.electronSrc = cms.InputTag("slimmedElectrons")
process.ggNtuplizer.calibelectronSrc = cms.InputTag("Electrons")
process.ggNtuplizer.photonSrc = cms.InputTag("slimmedPhotons") 
process.ggNtuplizer.calibphotonSrc = cms.InputTag("Photons")
process.ggNtuplizer.muonSrc = cms.InputTag("cleanedMu") 

# Track and other collections
process.ggNtuplizer.gsfTrackSrc = cms.InputTag("reducedEgamma", "reducedGsfTracks")
process.ggNtuplizer.TrackLabel = cms.InputTag("generalTracks")
process.ggNtuplizer.gsfElectronLabel = cms.InputTag("gsfElectrons")
process.ggNtuplizer.PFAllCandidates = cms.InputTag("particleFlow")
process.ggNtuplizer.packedPFCands = cms.InputTag("packedPFCandidates")

# ECAL RecHits
process.ggNtuplizer.ebReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits")
process.ggNtuplizer.eeReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits")
process.ggNtuplizer.esReducedRecHitCollection = cms.InputTag("reducedEgamma", "reducedESRecHits")
process.ggNtuplizer.reducedEcalRecHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits")
process.ggNtuplizer.reducedEcalRecHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits")
process.ggNtuplizer.reducedEcalRecHitsES = cms.InputTag("reducedEgamma", "reducedESRecHits")

# Other collections
process.ggNtuplizer.recoPhotonSrc = cms.InputTag("reducedEgamma", "reducedGedPhotonCores")

# Jet collections
process.ggNtuplizer.PFJetAK4Collection = cms.InputTag("slimmedJets")  # From cfi
process.ggNtuplizer.ak8JetSrc = cms.InputTag("slimmedJetsAK8")
process.ggNtuplizer.GENJetAK4 = cms.InputTag("slimmedGenJets")
process.ggNtuplizer.tauSrc = cms.InputTag("slimmedTaus")

# Isolation producers (if available)
process.ggNtuplizer.elePFClusEcalIsoProducer = cms.InputTag("electronEcalPFClusterIsolationProducer")
process.ggNtuplizer.elePFClusHcalIsoProducer = cms.InputTag("electronHcalPFClusterIsolationProducer")

# Filters
process.ggNtuplizer.ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter")

# Photon ID maps
process.ggNtuplizer.phoWP80MapToken = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIIWinter22-v1-wp80")
process.ggNtuplizer.phoWP90MapToken = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIIWinter22-v1-wp90")

# Jet and object selection cuts
process.ggNtuplizer.minJetPt = cms.untracked.double(25.)
process.ggNtuplizer.maxEta = cms.untracked.double(3.)

# Jet veto map
process.ggNtuplizer.JetVetoMap = cms.string('JetVetoMaps/'+JetVeto_tag+'.root')

# New particles for BSM searches (from cfi)
process.ggNtuplizer.newParticles = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35)

# Storage flags (from cfi file)
process.ggNtuplizer.store_electrons = cms.untracked.bool(True)
process.ggNtuplizer.store_muons = cms.untracked.bool(True)
process.ggNtuplizer.store_photons = cms.untracked.bool(True)
process.ggNtuplizer.store_ak4jets = cms.untracked.bool(True)
process.ggNtuplizer.store_CHS_met = cms.untracked.bool(True)
process.ggNtuplizer.store_PUPPI_met = cms.untracked.bool(True)
process.ggNtuplizer.store_electron_scalnsmear = cms.untracked.bool(False) 
process.ggNtuplizer.store_photon_scalnsmear = cms.untracked.bool(False)   
process.ggNtuplizer.store_electron_idSF = cms.untracked.bool(True)
process.ggNtuplizer.store_photon_idSF = cms.untracked.bool(True)

# Define the process path
process.p = cms.Path(
    process.egmPhotonIDSequence * 
    process.ggNtuplizer
    )


