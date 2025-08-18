#!/bin/sh
production_tag=$1
workarea=$2
config=$3
Dataset=$4
publication=$5
site=$6
DBS=$7
YEAR=$8
IsDATA=$9
ERA=${10}

temp=crabfile_${YEAR}_${production_tag}.py

echo "IsDATA" ${IsDATA}

golden_json=""
if [ $YEAR == "2022" ] || [ $YEAR == "2022EE" ]; then
	golden_json='/eos/user/c/cmsdqm/www/CAF/certification/Collisions22/Cert_Collisions2022_355100_362760_Golden.json'
elif [ $YEAR == "2023" ] || [ $YEAR == "2023BPiX" ]; then
	golden_json='/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'
else
	golden_json='/eos/user/c/cmsdqm/www/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
fi

# Clear the file using truncate (original method)
truncate -s 0 $temp

# Write the basic config using echo with proper redirection
echo "from CRABClient.UserUtilities import config" >> $temp
echo "config = config()" >> $temp
echo "" >> $temp
echo "config.General.requestName = 'crab_${production_tag}'" >> $temp
echo "config.General.workArea = 'crab_${workarea}'" >> $temp
echo "config.General.transferOutputs = True" >> $temp
echo "config.General.transferLogs = True" >> $temp
echo "" >> $temp
echo "config.JobType.pluginName = 'Analysis'" >> $temp
echo "config.JobType.psetName = '${config}'" >> $temp
echo "config.JobType.inputFiles = ['JECfiles','JERfiles', 'JetVetoMaps', 'ElectronJson', 'PhotonJson']" >> $temp

# Add parameters based on data/MC
if [[ "$IsDATA" == "1" ]]; then
	echo "config.JobType.pyCfgParams = ['--YEAR','${YEAR}','--ERA','${ERA}','--IsDATA','--IsRun3']" >> $temp
else
	echo "config.JobType.pyCfgParams = ['--YEAR','${YEAR}','--ERA','${ERA}','--IsRun3']" >> $temp
fi

# Continue with common config
echo "config.JobType.outputFiles = ['ntupls.root']" >> $temp
echo "config.JobType.maxJobRuntimeMin = 2700" >> $temp
echo "config.JobType.maxMemoryMB = 4000" >> $temp
echo "config.JobType.allowUndistributedCMSSW = True" >> $temp
echo "" >> $temp
echo "config.Data.inputDataset = '$Dataset'" >> $temp
echo "config.Data.inputDBS = '$DBS'" >> $temp

# Data-specific or MC-specific splitting
if [[ "$IsDATA" == "1" ]]; then
	echo "config.Data.splitting = 'LumiBased'" >> $temp
	echo "config.Data.lumiMask = '${golden_json}'" >> $temp
	echo "config.Data.unitsPerJob = 10" >> $temp
else
	echo "config.Data.splitting = 'FileBased'" >> $temp
	echo "config.Data.unitsPerJob = 1" >> $temp
fi

# Final config
echo "#config.Data.totalUnits = 1 #for test job" >> $temp
echo "" >> $temp
echo "config.Data.outLFNDirBase = '/store/user/rpramani/run3/2023/data_era_C'" >> $temp
echo "config.Data.publication = False" >> $temp
echo "config.Data.outputDatasetTag = '${production_tag}'" >> $temp
echo "" >> $temp
echo "config.Site.storageSite = '$site'" >> $temp
