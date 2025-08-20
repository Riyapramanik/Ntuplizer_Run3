from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_ZGto2NuG-1Jets_PTG-200to400'
config.General.workArea = 'crab_2023_MC'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run.py'
config.JobType.inputFiles = ['JECfiles','JERfiles', 'JetVetoMaps', 'ElectronJson', 'PhotonJson']
config.JobType.pyCfgParams = ['--YEAR','2023','--ERA','C','--IsRun3']
config.JobType.outputFiles = ['ntupls.root']
config.JobType.maxJobRuntimeMin = 2700
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/ZGto2NuG-1Jets_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 1 #for test job

config.Data.outLFNDirBase = '/store/user/rpramani/run3/2023/MC'
config.Data.publication = False
config.Data.outputDatasetTag = 'ZGto2NuG-1Jets_PTG-200to400'

config.Site.storageSite = 'T3_CH_CERNBOX'
