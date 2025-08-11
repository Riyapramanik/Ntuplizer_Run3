from CRABClient.UserUtilities import config
import sys

config = config()


#**************************submit function***********************
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#****************************************************************


config.General.workArea     = 'crab_projects_ntuples'
config.JobType.pluginName   = 'Analysis' # PrivateMC
config.JobType.psetName     = 'run_mc2022_133X.py'
config.Data.inputDBS        = 'global'    
config.Data.splitting       = 'FileBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.totalUnits      = -1
config.Data.publication     = False
#config.Site.storageSite     = 'T2_CH_CERN'
#config.Site.storageSite = 'T2_IN_TIFR'
config.Site.storageSite ='T3_CH_CERNBOX'

# dataset dependent configuration
#############G+jet samples
config.Data.unitsPerJob    = 10
#config.Data.outLFNDirBase  = '/store/user/shilpi/highPtPhotonID/GJetPtBin/v2'
config.Data.outLFNDirBase ='/store/user/rpramani'
#config.Data.outLFNDirBase ='/eos/user/r/rpramani'

###
config.General.requestName = 'DY_MC'
config.Data.inputDataset   = '/DYTo2L_MLL-50_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v1/MINIAODSIM'
