#!/bin/bash

source /cvmfs/cms.cern.ch/cmsset_default.sh;
presentDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )
cd $presentDir
eval 'scramv1 project CMSSW CMSSW_10_6_24'
source /cvmfs/cms.cern.ch/cmsset_default.sh;
cd ${presentDir}/CMSSW_10_6_24/src
eval `scramv1 runtime -sh`;
git cms-init
git cms-addpkg CondFormats/JetMETObjects
git cms-addpkg JetMETCorrections/Modules;
git cms-addpkg RecoBTag/Combined;
git cms-addpkg JetMETCorrections/Type1MET;
git cms-addpkg RecoMET/METFilters;
git cms-addpkg RecoEgamma/PhotonIdentification/;
git cms-addpkg RecoEgamma/ElectronIdentification/;
git cms-addpkg CommonTools/PileupAlgos;
git cms-merge-topic -u cms-tau-pog:CMSSW_10_6_X_tau-pog_MVAdm;
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_102X_v3;
git clone git@github.com:cms-egamma/EgammaPostRecoTools.git EgammaUser/EgammaPostRecoTools;
cd EgammaUser/EgammaPostRecoTools; git checkout master; cd -;
git clone -b CMSSW_10_6_24 git@github.com:wadud001/ggAnalysis.git;
scram b -j$(nproc);