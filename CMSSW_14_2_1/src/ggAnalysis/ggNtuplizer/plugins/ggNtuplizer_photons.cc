#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/EGMIDSFManager.h"

using namespace std;
Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoESEnP1_;
vector<float>  phoESEnP2_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<float>  phoSCBrem_;
vector<int>    phohasPixelSeed_;
vector<int>    phoEleVeto_;
vector<float>  phoR9_;
vector<float>  phoHoverE_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE2x2Full5x5_;
vector<float>  phoE5x5Full5x5_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoPFChWorstIso_;
vector<float>  phoSeedBCE_;
vector<float>  phoSeedBCEta_;
vector<float>  phoIDMVA_;
vector<float>  phoSeedTime_;
vector<float>  phoSeedEnergy_;
vector<float>  phoMIPTotEnergy_;
//vector<float>  phoSeedTimeFull5x5_;
vector<float>  phoMIPChi2_;
vector<float>  phoMIPSlope_;
vector<float>  phoMIPIntercept_;
vector<float>  phoMIPNhitCone_;
vector<float>  phoMIPIsHalo_;
vector<UShort_t> phoIDbit_;

vector<Float_t>      phoPFClusEcalIso_;
vector<Float_t>      phoPFClusHcalIso_;
vector<UChar_t>      nPhoTrkHollowConeDR04_;
vector<UChar_t>      nPhoTrkSolidConeDR04_;
vector<Float_t>      phoTrkSumPtSolidConeDR04_;
vector<Float_t>      phoTrkSumPtHollowConeDR04_;
vector<UChar_t>      nPhoTrkHollowConeDR03_;
vector<UChar_t>      nPhoTrkSolidConeDR03_;
vector<Float_t>      phoTrkSumPtSolidConeDR03_;
vector<Float_t>      phoTrkSumPtHollowConeDR03_;
vector<Float_t>      phoECALIso_;
vector<Float_t>      phoHCALIso_;

vector<Float_t>      phoE3x3Full5x5_;
vector<Float_t>      phoMaxEnergyXtal;
vector<Float_t>      phoE2ndFull5x5_;
vector<Float_t>      phoE1x3Full5x5_;
vector<Float_t>      phoE1x5Full5x5_;
vector<Float_t>      phoE2x5Full5x5_;
vector<Float_t>      phoE2x5BottomFull5x5_;
vector<Float_t>      phoE2x5LeftFull5x5_;
vector<Float_t>      phoE2x5MaxFull5x5_;
vector<Float_t>      phoE2x5RightFull5x5_;
vector<Float_t>      phoE2x5TopFull5x5_;
vector<Float_t>      phoEBottomFull5x5_;
vector<Float_t>      phoELeftFull5x5_;
vector<Float_t>      phoERightFull5x5_;
vector<Float_t>      phoETopFull5x5_;

//Scale and smearing
vector<float>  phoEnergyPostCorr_;
vector<float>  phoenergyScaleValue_;
vector<float>  phoenergySigmaValue_;
vector<float>  phoScale_unc_up_;
vector<float>  phoScale_unc_dn_;
vector<float>  phoSigma_unc_up_;
vector<float>  phoSigma_unc_dn_;

//ID scale factor
vector<float>  phoIDSF_Tight_;
vector<float>  phoIDSFUp_Tight_;
vector<float>  phoIDSFDown_Tight_;

vector<float> phoHaloTaggerMVAVal_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesPhotons(TTree* tree) {
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
  tree->Branch("phoSCE",                  &phoSCE_);
  tree->Branch("phoSCRawE",               &phoSCRawE_);
  tree->Branch("phoESEnP1",               &phoESEnP1_);
  tree->Branch("phoESEnP2",               &phoESEnP2_);
  tree->Branch("phoSCEta",                &phoSCEta_);
  tree->Branch("phoSCPhi",                &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",               &phoSCBrem_);
  tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
  tree->Branch("phoEleVeto",              &phoEleVeto_);
  tree->Branch("phoR9",                   &phoR9_);
  tree->Branch("phoHoverE",               &phoHoverE_);
  tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
  tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  tree->Branch("phoSeedBCE",              &phoSeedBCE_);
  tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
  tree->Branch("phoIDMVA",                &phoIDMVA_);
  tree->Branch("phoSeedTime",             &phoSeedTime_);
  tree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
  tree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);

  tree->Branch("phoPFClusEcalIso",          &phoPFClusEcalIso_);
  tree->Branch("phoPFClusHcalIso",          &phoPFClusHcalIso_);
  tree->Branch("nPhoTrkSolidConeDR03",      &nPhoTrkSolidConeDR03_);
  tree->Branch("nPhoTrkHollowConeDR03",     &nPhoTrkHollowConeDR03_);
  tree->Branch("phoTrkSumPtSolidConeDR03",  &phoTrkSumPtSolidConeDR03_);
  tree->Branch("phoTrkSumPtHollowConeDR03", &phoTrkSumPtHollowConeDR03_);
  tree->Branch("nPhoTrkSolidConeDR04",      &nPhoTrkSolidConeDR04_);
  tree->Branch("nPhoTrkHollowConeDR04",     &nPhoTrkHollowConeDR04_);
  tree->Branch("phoTrkSumPtSolidConeDR04",   &phoTrkSumPtSolidConeDR04_);
  tree->Branch("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04_);
  tree->Branch("phoECALIso",                &phoECALIso_);
  tree->Branch("phoHCALIso",                &phoHCALIso_);
  //tree->Branch("phoSeedTimeFull5x5",              &phoSeedTimeFull5x5_);
  tree->Branch("phoMIPChi2",                      &phoMIPChi2_);
  tree->Branch("phoMIPSlope",                     &phoMIPSlope_);
  tree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
  tree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
  tree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);
  tree->Branch("phoIDbit",         &phoIDbit_);

  tree->Branch("phoE2x5BottomFull5x5",      &phoE2x5BottomFull5x5_);
  tree->Branch("phoE2x5LeftFull5x5",        &phoE2x5LeftFull5x5_);
  tree->Branch("phoE2x5MaxFull5x5",         &phoE2x5MaxFull5x5_);
  tree->Branch("phoE2x5RightFull5x5",       &phoE2x5RightFull5x5_);
  tree->Branch("phoE2x5TopFull5x5",         &phoE2x5TopFull5x5_);
  tree->Branch("phoEBottomFull5x5",         &phoEBottomFull5x5_);
  tree->Branch("phoELeftFull5x5",           &phoELeftFull5x5_);
  tree->Branch("phoERightFull5x5",          &phoERightFull5x5_);
  tree->Branch("phoETopFull5x5",            &phoETopFull5x5_);

  tree->Branch("phoE3x3Full5x5",            &phoE3x3Full5x5_);
  tree->Branch("phoMaxEnergyXtal",          &phoMaxEnergyXtal);
  tree->Branch("phoE2ndFull5x5",            &phoE2ndFull5x5_);
  tree->Branch("phoE1x3Full5x5",            &phoE1x3Full5x5_);
  tree->Branch("phoE1x5Full5x5",            &phoE1x5Full5x5_);
  tree->Branch("phoE2x5Full5x5",            &phoE2x5Full5x5_);

  tree->Branch("phoecalTrkEnergyPostCorr",      &phoEnergyPostCorr_);
  tree->Branch("phoenergyScaleValue",            &phoenergyScaleValue_);
  tree->Branch("phoenergySigmaValue_",        &phoenergySigmaValue_);
  tree->Branch("phoScale_unc_up_",            &phoScale_unc_up_);
  tree->Branch("phoScale_unc_dn_",            &phoScale_unc_dn_);
  tree->Branch("phoSigma_unc_up_",            &phoSigma_unc_up_);
  tree->Branch("phoSigma_unc_dn_",            &phoSigma_unc_dn_);  
 
  tree->Branch("phoIDSF_Tight",      &phoIDSF_Tight_);
  tree->Branch("phoIDSFUp_Tight",    &phoIDSFUp_Tight_);
  tree->Branch("phoIDSFDown_Tight",  &phoIDSFDown_Tight_);
 
  
  tree->Branch("phoHaloTaggerMVAVal",       &phoHaloTaggerMVAVal_);
  
}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {

  // cleanup from previous execution
  phoE_                   .clear();
  phoEt_                  .clear();
  phoEta_                 .clear();
  phoPhi_                 .clear();
  phoSCE_                 .clear();
  phoSCRawE_              .clear();
  phoESEnP1_              .clear();
  phoESEnP2_              .clear();
  phoSCEta_               .clear();
  phoSCPhi_               .clear();
  phoSCEtaWidth_          .clear();
  phoSCPhiWidth_          .clear();
  phoSCBrem_              .clear();
  phohasPixelSeed_        .clear();
  phoEleVeto_             .clear();
  phoR9_                  .clear();
  phoHoverE_              .clear();
  phoESEffSigmaRR_        .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE2x2Full5x5_         .clear();
  phoE5x5Full5x5_         .clear();
  phoR9Full5x5_           .clear();
  phoPFChIso_             .clear();
  phoPFPhoIso_            .clear();
  phoPFNeuIso_            .clear();
  phoPFChWorstIso_        .clear();
  phoSeedBCE_           .clear();
  phoSeedBCEta_         .clear();
  phoIDMVA_               .clear();
  phoSeedTime_            .clear();
  phoSeedEnergy_          .clear();
  phoMIPTotEnergy_        .clear();


  phoPFClusEcalIso_               . clear();
  phoPFClusHcalIso_               . clear();
  nPhoTrkHollowConeDR03_          . clear();
  nPhoTrkSolidConeDR03_           . clear();
  phoTrkSumPtSolidConeDR03_       . clear();
  phoTrkSumPtHollowConeDR03_      . clear();
  nPhoTrkHollowConeDR04_          . clear();
  nPhoTrkSolidConeDR04_           . clear();
  phoTrkSumPtSolidConeDR04_       . clear();
  phoTrkSumPtHollowConeDR04_      . clear();
  phoECALIso_                     . clear();
  phoHCALIso_                     . clear();
	
  //phoSeedTimeFull5x5_   .clear();
  phoMIPChi2_           .clear();
  phoMIPSlope_          .clear();
  phoMIPIntercept_      .clear();
  phoMIPNhitCone_       .clear();
  phoMIPIsHalo_         .clear();
  
  phoIDbit_.clear();

  phoE2x5BottomFull5x5_.clear();
  phoE2x5LeftFull5x5_.clear();
  phoE2x5MaxFull5x5_.clear();
  phoE2x5RightFull5x5_.clear();
  phoE2x5TopFull5x5_.clear();
  phoEBottomFull5x5_.clear();
  phoELeftFull5x5_.clear();
  phoERightFull5x5_.clear();
  phoETopFull5x5_.clear();
  phoE3x3Full5x5_.clear();
  phoMaxEnergyXtal.clear();
  phoE2ndFull5x5_.clear();
  phoE1x3Full5x5_.clear();
  phoE1x5Full5x5_.clear();
  phoE2x5Full5x5_.clear();

  phoEnergyPostCorr_.clear();
  phoenergyScaleValue_.clear();
  phoenergySigmaValue_.clear();
  phoScale_unc_up_.clear();
  phoScale_unc_dn_.clear();
  phoSigma_unc_up_.clear();
  phoSigma_unc_dn_.clear();

  phoIDSF_Tight_       .clear();
  phoIDSFUp_Tight_     .clear();
  phoIDSFDown_Tight_   .clear();
  
  phoHaloTaggerMVAVal_.clear();

  nPho_ = 0;
  
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  if (photonHandle.isValid()) {
    
  //edm::Handle<reco::PhotonCollection> recoPhotonHandle;
  //e.getByToken(recophotonCollection_, recoPhotonHandle);

  edm::Handle<EcalRecHitCollection> RecHitsEB;
  e.getByToken(ebReducedRecHitCollection_, RecHitsEB);

  edm::Handle<EcalRecHitCollection> RecHitsEE;
  e.getByToken(eeReducedRecHitCollection_, RecHitsEE);
  CaloTopology const *p_topology = &es.getData(caloTop);  // get calo topology
  const EcalRecHitCollection *ebRecHits = RecHitsEB.product();
  const EcalRecHitCollection *eeRecHits = RecHitsEE.product();
  EcalClusterLazyTools       lazyTool    (e, ecalClusterToolsESGetTokens_.get(es), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, ecalClusterToolsESGetTokens_.get(es), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  if(store_photons){
  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

    phoE_             .push_back(iPho->energy());
    phoEt_            .push_back(iPho->et());
    phoEta_           .push_back(iPho->eta());
    phoPhi_           .push_back(iPho->phi());
    phoSCE_           .push_back((*iPho).superCluster()->energy());
    phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
    phoESEnP1_        .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
    phoESEnP2_        .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
    phoSCEta_         .push_back((*iPho).superCluster()->eta());
    phoSCPhi_         .push_back((*iPho).superCluster()->phi());
    phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
    phoSCBrem_        .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
    phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
    phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
    phoR9_            .push_back(iPho->r9());
    phoHoverE_        .push_back(iPho->hadTowOverEm());
    phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));

    phoPFChIso_       .push_back(iPho->chargedHadronIso());
    phoPFPhoIso_      .push_back(iPho->photonIso());
    phoPFNeuIso_      .push_back(iPho->neutralHadronIso());
    phoPFChWorstIso_  .push_back(iPho->chargedHadronWorstVtxIso());
    
    phoPFClusEcalIso_          . push_back(iPho->ecalPFClusterIso());
    phoPFClusHcalIso_          . push_back(iPho->hcalPFClusterIso());
    phoTrkSumPtSolidConeDR04_  . push_back(iPho->trkSumPtSolidConeDR04() );
    phoTrkSumPtSolidConeDR03_  . push_back(iPho->trkSumPtSolidConeDR03());
    phoTrkSumPtHollowConeDR04_ . push_back(iPho->trkSumPtHollowConeDR04());
    phoTrkSumPtHollowConeDR03_ . push_back(iPho->trkSumPtHollowConeDR03());
    
    nPhoTrkHollowConeDR04_     . push_back(iPho->nTrkHollowConeDR04());
    nPhoTrkSolidConeDR04_      . push_back(iPho->nTrkSolidConeDR04());
    nPhoTrkHollowConeDR03_     . push_back(iPho->nTrkHollowConeDR03());
    nPhoTrkSolidConeDR03_      . push_back(iPho->nTrkSolidConeDR03());
    phoECALIso_                . push_back(iPho->ecalIso());
    phoHCALIso_                . push_back(iPho->hcalIso());


    phoE2x5BottomFull5x5_      . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5Bottom);
    phoE2x5LeftFull5x5_        . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5Left);
    phoE2x5MaxFull5x5_         . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5Max);
    phoE2x5RightFull5x5_       . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5Right);
    phoE2x5TopFull5x5_         . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5Top);
    phoEBottomFull5x5_         . push_back(iPho->full5x5_showerShapeVariables()                       . eBottom);
    phoELeftFull5x5_           . push_back(iPho->full5x5_showerShapeVariables()                       . eLeft);
    phoERightFull5x5_          . push_back(iPho->full5x5_showerShapeVariables()                       . eRight);
    phoETopFull5x5_            . push_back(iPho->full5x5_showerShapeVariables()                       . eTop);
    phoMaxEnergyXtal           . push_back(iPho->full5x5_maxEnergyXtal());
    phoE2ndFull5x5_            . push_back(iPho->full5x5_showerShapeVariables()                       . e2nd);
    phoE1x3Full5x5_            . push_back(iPho->full5x5_showerShapeVariables()                       . e1x3);
    phoE1x5Full5x5_            . push_back(iPho->full5x5_showerShapeVariables()                       . e1x5);
    phoE2x5Full5x5_            . push_back(iPho->full5x5_showerShapeVariables()                       . e2x5);
    phoE3x3Full5x5_          . push_back(iPho->full5x5_showerShapeVariables()                                    . e3x3);
    
    phoHaloTaggerMVAVal_       .push_back(iPho->haloTaggerMVAVal());

    //For scale and smearing
    float originalPt = iPho->pt();
    float scEta = iPho->superCluster()->eta();
    float r9 = iPho->r9();
    DetId seedDetId = iPho->superCluster()->seed()->seed();
    int seedGain = EGMCorrectionManager::GetSeedGain(seedDetId, e, ebReducedRecHitCollection_, eeReducedRecHitCollection_);
    int run = isData ? e.id().run() : 1;
    double randomNum = isData ? 0.0 : normalDistribution_(randomGenerator_);
    double correctedPt = egmCorrectionManager_->applyCorrectedPhotonPt(originalPt, run, scEta, r9, seedGain, isData, randomNum);
    phoEnergyPostCorr_.push_back(correctedPt); 
    
    if (isData) {
      double scale = egmCorrectionManager_->getPhotonScale(run, scEta, r9, originalPt, seedGain);
      double scaleUnc = egmCorrectionManager_->getPhotonScaleUnc(originalPt, r9, std::abs(scEta));
      
      phoenergyScaleValue_.push_back(scale);
      phoScale_unc_up_.push_back(scale + scaleUnc);
      phoScale_unc_dn_.push_back(scale - scaleUnc);
      
      phoenergySigmaValue_.push_back(0.0);
      phoSigma_unc_up_.push_back(0.0);
      phoSigma_unc_dn_.push_back(0.0);
      
    } else {
      double smear = egmCorrectionManager_->getPhotonSmear(originalPt, r9, std::abs(scEta));
      double smearUnc = egmCorrectionManager_->getPhotonSmearUnc(originalPt, r9, std::abs(scEta));
      phoenergyScaleValue_.push_back(1.0);
      phoScale_unc_up_.push_back(0.0);
      phoScale_unc_dn_.push_back(0.0);
      
      phoenergySigmaValue_.push_back(smear);
      phoSigma_unc_up_.push_back(smear + smearUnc);
      phoSigma_unc_dn_.push_back(smear - smearUnc);
      }
    double pt = iPho->pt();
    double eta = iPho->eta();
    double phi = iPho->phi();
    phoIDSF_Tight_.push_back(egmIDSFManager_->getPhotonIDSF("Tight", pt, eta, phi));
    phoIDSFUp_Tight_.push_back(egmIDSFManager_->getPhotonIDSFUncUp("Tight", pt, eta, phi));
    phoIDSFDown_Tight_.push_back(egmIDSFManager_->getPhotonIDSFUncDown("Tight", pt, eta, phi));
    
    ///Run 3
    ///cut values here: https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/PhotonIdentification/python/Identification/mvaPhotonID_Winter22_122X_V1_cff.py#L25-L36

    phoIDMVA_         .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIIWinter22v1Values"));
    
    // VID decisions     
    UShort_t tmpphoIDbit = 0;
    ///Run 3
    bool isPassLoose  = iPho->photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-loose");
    if (isPassLoose)  setbit(tmpphoIDbit, 0);   
    bool isPassMedium = iPho->photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-medium");
    if (isPassMedium) setbit(tmpphoIDbit, 1);    
    bool isPassTight  = iPho->photonID("cutBasedPhotonID-RunIIIWinter22-122X-V1-tight");
    if (isPassTight)  setbit(tmpphoIDbit, 2);

    phoIDbit_.push_back(tmpphoIDbit);
    DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    
    //For Run3
    const EcalRecHitCollection* rechits = isBarrel ? ebRecHits : eeRecHits;
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {

      phoSeedTime_  .push_back((*theSeedHit).time());     
      phoSeedEnergy_.push_back((*theSeedHit).energy());
    } else{
      phoSeedTime_  .push_back(-99.);
      phoSeedEnergy_.push_back(-99.);
    }
    const auto &vCov  = (isBarrel?EcalClusterTools::localCovariances(*((*iPho).superCluster()->seed()), ebRecHits,  p_topology) : EcalClusterTools::localCovariances(*((*iPho).superCluster()->seed()), eeRecHits,  p_topology));
    
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
    phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());
    phoMIPTotEnergy_         .push_back(iPho->mipTotEnergy());
    phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
    phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
    //phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
    //std::cout<<"seed time from supercluster "<<lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster()))<<std::endl;
    phoMIPChi2_        .push_back(iPho->mipChi2());
    phoMIPSlope_       .push_back(iPho->mipSlope());
    phoMIPIntercept_   .push_back(iPho->mipIntercept());
    phoMIPNhitCone_    .push_back(iPho->mipNhitCone());
    phoMIPIsHalo_      .push_back(iPho->mipIsHalo());
    nPho_++;
  }
  }
  }

}

