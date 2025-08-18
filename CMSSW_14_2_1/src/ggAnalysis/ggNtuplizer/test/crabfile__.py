from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'crab_'
config.General.workArea = 'crab_'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = ''
config.JobType.inputFiles = ['JECfiles','JERfiles', 'JetVetoMaps', 'ElectronJson', 'PhotonJson']
config.JobType.pyCfgParams = ['--YEAR','','--ERA','','--IsRun3']
config.JobType.outputFiles = ['ntupls.root']
config.JobType.maxJobRuntimeMin = 2700
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = ''
config.Data.inputDBS = ''
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 1 #for test job

config.Data.outLFNDirBase = '/eos/user/r/rpramani/run3/2023/data_era_C'
config.Data.publication = False
config.Data.outputDatasetTag = ''

config.Site.storageSite = ''
