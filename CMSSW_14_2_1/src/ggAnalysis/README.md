cmsrel CMSSW_10_6_24 <br>
cd CMSSW_10_6_24/src <br>
cmsenv <br>
git cms-init<br>
git cms-addpkg CondFormats/JetMETObjects<br>
git cms-addpkg JetMETCorrections/Modules;<br>
git cms-addpkg RecoBTag/Combined;  <br> 
git cms-addpkg JetMETCorrections/Type1MET;  <br> 
git cms-addpkg RecoMET/METFilters;  <br> 
git cms-addpkg RecoEgamma/PhotonIdentification/;  <br> 
git cms-addpkg RecoEgamma/ElectronIdentification/;  <br> 
git cms-addpkg CommonTools/PileupAlgos;  <br>
git cms-merge-topic -u cms-tau-pog:CMSSW_10_6_X_tau-pog_MVAdm;   <br> 
git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_102X_v3;  <br>
git cms-addpkg RecoEgamma/EgammaTools;  <br>
git clone https://github.com/cms-egamma/EgammaPostRecoTools.git;  <br>
mv EgammaPostRecoTools/python/EgammaPostRecoTools.py RecoEgamma/EgammaTools/python/;  <br>
rm -rf EgammaPostRecoTools;  <br>
git clone https://github.com/cms-data/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/ ;  <br>
git cms-addpkg EgammaAnalysis/ElectronTools;  <br>
git clone -b CMSSW_10_6_24 git@github.com:wadud001/ggAnalysis.git;<br>
scram b -j$(nproc); <br>
