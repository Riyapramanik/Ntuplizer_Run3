import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.patSequences_cff import *
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import *

process = cms.Process('ggKit')


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(allowUnscheduled=cms.untracked.bool(True))
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.EventContent.EventContent_cff")
#process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");


#!!!!!!!!!!!!!!!!!!!!!!
#Tags for diff year
#!!!!!!!!!!!!!!!!!!!!!!
from RecoJets.Configuration.GenJetParticles_cff import *
#take this from argument
IsDATA = True
YEAR = "2022"
ERA = "C"

if IsDATA:
    if YEAR=="2022":
        process.GlobalTag.globaltag = "130X_dataRun3_v2"
        JEC_tag = "Summer22_22Sep2023_RunCD_V2_DATA"
        JER_tag = 'Summer22_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22_23Sep2023_RunCD_v1'
    elif YEAR=="2022EE":
        process.GlobalTag.globaltag = "130X_dataRun3_PromptAnalysis_v1"
        JEC_tag = "Summer22EE_22Sep2023_Run"+ERA+"_V2_DATA" 
        JER_tag = 'Summer22EE_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22EE_23Sep2023_RunEFG_v1'
    elif YEAR=="2023":
        process.GlobalTag.globaltag = "130X_dataRun3_PromptAnalysis_v1"
        if ERA == "Cv4":
            JEC_tag = "Summer23Prompt23_RunCv4_V1_DATA"
        else:
            JEC_tag = "Summer23Prompt23_RunCv123_V1_DATA"
        JER_tag = 'Summer23Prompt23_RunCv123_JRV1_MC'       
        JetVeto_tag = "Summer23Prompt23_RunC_v1"
    elif YEAR=="2023BPiX":
        process.GlobalTag.globaltag = "130X_dataRun3_PromptAnalysis_v1"
        JEC_tag = "Summer23BPixPrompt23_RunD_V1_DATA" 
        JER_tag = 'Summer23BPixPrompt23_RunD_JRV1_MC'       
        JetVeto_tag = "Summer23BPixPrompt23_RunD_v1"
else:
    if YEAR=="2022":
        process.GlobalTag.globaltag = "130X_mcRun3_2022_realistic_v5"
        JEC_tag = 'Summer22_22Sep2023_V2_MC'
        JER_tag = 'Summer22_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22_23Sep2023_RunCD_v1'
    elif YEAR=="2022EE":
        process.GlobalTag.globaltag = "130X_mcRun3_2022_realistic_postEE_v6"
        JEC_tag = 'Summer22EE_22Sep2023_V2_MC'
        JER_tag = 'Summer22EE_22Sep2023_JRV1_MC'
        JetVeto_tag = 'Summer22EE_23Sep2023_RunEFG_v1'
    elif YEAR=="2023":
        process.GlobalTag.globaltag = "130X_mcRun3_2023_realistic_v14"
        JEC_tag = "Summer23Prompt23_V1_MC"
        if ERA=="Cv4":
            JER_tag = 'Summer23Prompt23_RunCv4_JRV1_MC'
        else:
            JER_tag = 'Summer23Prompt23_RunCv123_JRV1_MC'   
        JetVeto_tag = "Summer23Prompt23_RunC_v1"
    elif YEAR=="2023BPiX":
        process.GlobalTag.globaltag = "130X_mcRun3_2023_realistic_postBPix_v2"
        JEC_tag = "Summer23BPixPrompt23_V1_MC"
        JER_tag = 'Summer23BPixPrompt23_RunD_JRV1_MC'     
        JetVeto_tag = "Summer23BPixPrompt23_RunD_v1"
    else:
        process.GlobalTag.globaltag = "130X_mcRun3_2022_realistic_v5"
        JEC_tag = "Summer22_22Sep2023_RunCD_V2_DATA" 
        JER_tag = 'Summer22_22Sep2023_JRV1_MC'       
        


#!!!!!!!!!!!!!!
#! Can make it as argument
#!!!!!!!!!!!!!!
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_v2')
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/data/Run2022C/EGamma/MINIAOD/22Sep2023-v1/50000/6015928f-0763-4e7f-9d48-ff0f2ddaf12e.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/mc/Run3Summer22EEMiniAODv4/ZGto2NuG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_postEE_v6-v3/2560000/0adadef6-3fb4-453d-ae24-6bc79b4338cb.root'))
print(process.source)

#!!!!!!!!!!!!!!
#PAT Conditions
#!!!!!!!!!!!!!!
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

#output root file
process.TFileService = cms.Service("TFileService", fileName=cms.string("data.root"))

#!!!!!!!!!!!!!!!!!!!!!!!!
#PAT Sequences Loading
#!!!!!!!!!!!!!!!!!!!!!!!
from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#!!!!!!!!!!!!!!!!!!!!!!!!
# VID For Run3
#!!!!!!!!!!!!!!!!!!!!!!!!
# define which Electron IDs want to produce
#needed for electron and photon ID mapping 
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *  
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
el_id_modules = [
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_iso_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff",
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    ]
#Add them to the VID producer
for iModule in el_id_modules:
    setupAllVIDIdsInModule(process, iModule, setupVIDElectronSelection)

#Photon ID variables
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
pho_id_modules = [
    'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_RunIIIWinter22_122X_V1_cff',
    'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Winter22_122X_V1_cff',
    ]

for iModule in pho_id_modules:
    setupAllVIDIdsInModule(process, iModule, setupVIDPhotonSelection)

#!!!!!!!!!!!!!!!!!!!
#Jet Tagger
#For QG likelihood 
from RecoJets.JetProducers.QGTagger_cfi import  QGTagger
process.qgtagger = QGTagger.clone(
    srcJets = cms.InputTag("slimmedJetsPuppi"),
    srcVertexCollection="offlineSlimmedPrimaryVertices"
)

#For pileup jet ID 
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId
process.pileupJetID= pileupJetId.clone(
    jets = cms.InputTag('slimmedJetsPuppi'),
    inputIsCorrected = False,
    applyJec = False,
    vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
)

#Creating a new AK4 jet collection adding QG likelihood and pileup jet ID values
process.slimmedJetsPuppiWithInfo = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("slimmedJetsPuppi"),
    userFloats = cms.PSet(
        qgLikelihood = cms.InputTag('qgtagger:qgLikelihood'),
        pileupJetId_fullDiscriminant = cms.InputTag('pileupJetID:fullDiscriminant'),
    )
)

#!!!!!!!!!!!!!!!!!
# MET 
# For CHS MET
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
    process,
    isData  = IsDATA,
    postfix = "Updated"
)
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True )
# For Puppi
runMetCorAndUncFromMiniAOD(process,
                           isData=IsDATA,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )
process.puppi.useExistingWeights = True

#Muon information
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import patMuons #.miniIsoParams #as _miniIsoParams

process.slimmedMuonsUpdated = cms.EDProducer("PATMuonUpdater",
    src = cms.InputTag("slimmedMuons"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    computeMiniIso = cms.bool(False),
    fixDxySign = cms.bool(True),
    pfCandsForMiniIso = cms.InputTag("packedPFCandidates"),
    miniIsoParams = patMuons.miniIsoParams, #PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.miniIsoParams, # so they're in sync
    recomputeMuonBasicSelectors = cms.bool(False),
)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if YEAR == "2022":
    year_int = 2022
elif YEAR == "2022EE":
    year_int = 2022
elif YEAR == "2023":
    year_int = 2023
elif YEAR == "2023BPiX":
    year_int = 2023
else:
    year_int = 2022
period_str = ERA

#===== MET Filters ==
process.load('RecoMET.METFilters.primaryVertexFilter_cfi')
process.primaryVertexFilter.vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
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
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
process.load('RecoMET.METFilters.hfNoisyHitsFilter_cfi')
process.allMetFilterPaths = cms.Sequence(
    process.primaryVertexFilter *
    process.globalSuperTightHalo2016Filter *
    process.BadPFMuonFilter *
    process.BadPFMuonDzFilter #*
#    process.eeBadScFilter *
#    process.ecalBadCalibFilter *
#    process.hfNoisyHitsFilter
)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#!!!!!!!!!!!!!!!!!!!!!!!!!
#Analyser!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!
#Load the Ntuplizer
process.load("ggAnalysis.ggNtuplizer.ggNtuplizer_miniAOD_cfi")

#Required parameters
process.ggNtuplizer.Data = cms.untracked.bool(IsDATA)
process.ggNtuplizer.MonteCarlo = cms.untracked.bool(not IsDATA)
process.ggNtuplizer.YEAR = cms.untracked.string(YEAR)

process.ggNtuplizer.dataYear = cms.int32(year_int)
process.ggNtuplizer.dataPeriod = cms.string(period_str) 
process.ggNtuplizer.isData = cms.bool(IsDATA)

process.ggNtuplizer.PFJetsAK4 = cms.InputTag("slimmedJetsPuppiWithInfo")  # updated collection
process.ggNtuplizer.GENJetAK4 = cms.InputTag("slimmedGenJets")
process.ggNtuplizer.minJetPt = cms.untracked.double(25.)
process.ggNtuplizer.maxEta = cms.untracked.double(3.)

process.ggNtuplizer.Muons = cms.InputTag("slimmedMuonsUpdated") #updated muon collection

process.ggNtuplizer.electronSrc = cms.InputTag("slimmedElectrons")
process.ggNtuplizer.photonSrc = cms.InputTag("slimmedPhotons")
process.ggNtuplizer.PFMet = cms.InputTag("slimmedMETs")
process.ggNtuplizer.PuppiMet = cms.InputTag("slimmedMETsPuppi")

process.ggNtuplizer.genParticleSrc = cms.InputTag("prunedGenParticles")
process.ggNtuplizer.packedPFCands = cms.InputTag("packedPFCandidates")  
process.ggNtuplizer.VtxLabel = cms.InputTag("offlineSlimmedPrimaryVertices")
process.ggNtuplizer.pileupCollection = cms.InputTag("slimmedAddPileupInfo")
process.ggNtuplizer.rhoLabel = cms.InputTag("fixedGridRhoFastjetAll")
process.ggNtuplizer.PFRho = cms.InputTag("fixedGridRhoFastjetAll")

process.ggNtuplizer.bits = cms.InputTag("TriggerResults","","HLT")
process.ggNtuplizer.TriggerObjects = cms.InputTag("slimmedPatTrigger")
process.ggNtuplizer.prescales = cms.InputTag("patTrigger","","RECO") 
process.ggNtuplizer.MET_Filters = cms.InputTag("TriggerResults","","RECO") # For 2024 data

#JER and JEC Files                                                                                                         
process.ggNtuplizer.jecL1FastFileAK4          = cms.string('JECfiles/'+JEC_tag+'/'+JEC_tag+'_L1FastJet_AK4PFPuppi.txt')
process.ggNtuplizer.jecL2RelativeFileAK4      = cms.string('JECfiles/'+JEC_tag+'/'+JEC_tag+'_L2Relative_AK4PFPuppi.txt')
process.ggNtuplizer.jecL3AbsoluteFileAK4      = cms.string('JECfiles/'+JEC_tag+'/'+JEC_tag+'_L3Absolute_AK4PFPuppi.txt')
process.ggNtuplizer.jecL2L3ResidualFileAK4    = cms.string('JECfiles/'+JEC_tag+'/'+JEC_tag+'_L2L3Residual_AK4PFPuppi.txt')
process.ggNtuplizer.JECUncFileAK4 = cms.string('JECfiles/'+JEC_tag+'/'+JEC_tag+'_UncertaintySources_AK4PFPuppi.txt')

process.ggNtuplizer.PtResoFileAK4  = cms.string('JERfiles/'+JER_tag+'/'+JER_tag+'_PtResolution_AK4PFPuppi.txt')
process.ggNtuplizer.PtSFFileAK4 = cms.string('JERfiles/'+JER_tag+'/'+JER_tag+'_SF_AK4PFPuppi.txt')

process.ggNtuplizer.JetVetoMap = cms.string('JetVetoMaps/'+JetVeto_tag+'.root')

process.ggNtuplizer.UltraLegacy = cms.untracked.bool(False)
process.ggNtuplizer.isRun3 = cms.untracked.bool(True)
process.ggNtuplizer.doGenParticles = cms.bool(False)
process.ggNtuplizer.isData = cms.bool(True)
process.ggNtuplizer.store_electrons = cms.untracked.bool(True)
process.ggNtuplizer.store_muons = cms.untracked.bool(True)  
process.ggNtuplizer.store_photons = cms.untracked.bool(True)  
process.ggNtuplizer.store_ak4jets = cms.untracked.bool(True)
process.ggNtuplizer.store_CHS_met = cms.untracked.bool(True)    
process.ggNtuplizer.store_PUPPI_met = cms.untracked.bool(True)
process.ggNtuplizer.applyEGMCorrections = cms.bool(True)

process.edTask = cms.Task()
for key in process.__dict__.keys():
    if(type(getattr(process,key)).__name__=='EDProducer' or type(getattr(process,key)).__name__=='EDFilter'):
        process.edTask.add(getattr(process,key))

process.p = cms.Path(
    process.egmPhotonIDSequence 
    *process.allMetFilterPaths              # NEEDED for MET filter flags
    *process.egmGsfElectronIDSequence       # NEEDED for electron ID
    *process.electronMVAValueMapProducer
    *process.slimmedMuonsUpdated            # NEEDED for updated muons
    *process.qgtagger                       # ADD THIS
    *process.pileupJetID                    # ADD THIS
    *process.slimmedJetsPuppiWithInfo       # ADD THIS
    *process.ggNtuplizer
)

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
process.patAlgosToolsTask = getPatAlgosToolsTask(process)
process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)

process.schedule = cms.Schedule(
#    process.pathRunPatAlgos,
    process.p
)

# Options
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.options.allowUnscheduled = cms.untracked.bool(True)

# Add message logger settings:
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.ggNtuplizer = cms.untracked.PSet(limit=cms.untracked.int32(-1))
