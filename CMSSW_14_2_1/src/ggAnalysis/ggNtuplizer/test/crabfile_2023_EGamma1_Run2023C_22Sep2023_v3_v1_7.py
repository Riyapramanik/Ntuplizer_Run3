from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_EGamma1_Run2023C_22Sep2023_v3_v1_7'
config.General.workArea = 'crab_2023_DATA_C'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run.py'
config.JobType.inputFiles = ['JECfiles','JERfiles', 'JetVetoMaps', 'ElectronJson', 'PhotonJson']
config.JobType.pyCfgParams = ['--YEAR','2023','--ERA','C','--IsDATA','--IsRun3']
config.JobType.outputFiles = ['ntupls.root']
config.JobType.maxJobRuntimeMin = 2700
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = '/EGamma1/Run2023C-22Sep2023_v3-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 1 #for test job

config.Data.outLFNDirBase = '/store/user/rpramani/run3/2023/data_era_C'
config.Data.publication = False
config.Data.outputDatasetTag = 'EGamma1_Run2023C_22Sep2023_v3_v1_7'

config.Site.storageSite = 'T3_CH_CERNBOX'
