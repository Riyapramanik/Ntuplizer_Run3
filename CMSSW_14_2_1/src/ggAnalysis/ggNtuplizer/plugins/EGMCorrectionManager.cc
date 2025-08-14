#include "ggAnalysis/ggNtuplizer/interface/EGMCorrectionManager.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "FWCore/Framework/interface/Event.h"
// Constructor - just store the year and period, then load the corrections
EGMCorrectionManager::EGMCorrectionManager(int year, const std::string& period)
    : year_(year), period_(period) {
    initializeCorrections();
}

void EGMCorrectionManager::initializeCorrections() {
    std::string electronFile = getElectronJSONFile();
    std::string photonFile = getPhotonJSONFile();
    
    electronCorrectionSet_ = correction::CorrectionSet::from_file(electronFile);
    photonCorrectionSet_ = correction::CorrectionSet::from_file(photonFile);
    
    setupElectronEvaluators();
    setupPhotonEvaluators();
    
}

std::string EGMCorrectionManager::getElectronJSONFile() {

  std::string basePath = "/eos/user/r/rpramani/run3_ntuplizer/CMSSW_14_2_1/src/ggAnalysis/ggNtuplizer/test/ElectronJson/";
  if (year_ == 2022){
    if (period_ == "B" || period_ == "C" || period_ == "D"){
      return basePath + "electronSS_EtDependent_22_BCD.json";
        }
    else if (period_ == "E" || period_ == "F" || period_ == "G"){
      return basePath + "electronSS_EtDependent_22_EFG.json";
        }
  }

  else if (year_ == 2023) {
        if (period_ == "C") {
            return basePath + "electronSS_EtDependent_23_C.json";
        }
        else if (period_ == "D") {
            return basePath + "electronSS_EtDependent_23_D.json";
        }
    }
  return "";
}


std::string EGMCorrectionManager::getPhotonJSONFile() {
    std::string basePath = "/eos/user/r/rpramani/run3_ntuplizer/CMSSW_14_2_1/src/ggAnalysis/ggNtuplizer/test/PhotonJson/";
    
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D" || 
            period_ == "preEE" || period_.empty()) {
            return basePath + "photonSS_EtDependent_22_BCD.json";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G" || period_ == "postEE") {
            return basePath + "photonSS_EtDependent_22_EFG.json";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C") {
            return basePath + "photonSS_EtDependent_23_C.json";
        }
        else if (period_ == "D") {
            return basePath + "photonSS_EtDependent_23_D.json";
        }
    }
    return "";
}

std::string EGMCorrectionManager::getYearSuffix() {
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D") {
            return "_2022preEE";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G") {
            return "_2022postEE";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C" || period_ == "preBPIX") {
            return "_2023preBPIX";
        }
        else if (period_ == "D") {
            return "_2023postBPIX";
        }
    }

    return ""; 
    
}

void EGMCorrectionManager::setupPhotonEvaluators() {
    std::string yearSuffix = getYearSuffix();
    std::string scaleName = "EGMScale_Compound_Pho" + yearSuffix;
    //std::string smearName = "EGMSmearAndSyst_PhoEtaR9" + yearSuffix;
    std::string smearName = "EGMSmearAndSyst_PhoPTsplit" + yearSuffix;
    photonScaleEvaluator_ = photonCorrectionSet_->compound().at(scaleName);    
    photonSmearEvaluator_ = photonCorrectionSet_->at(smearName);
}

void EGMCorrectionManager::setupElectronEvaluators() {
    std::string yearSuffix = getYearSuffix();
    
    std::string scaleName = "EGMScale_Compound_Ele" + yearSuffix;
    std::string smearName;
    //if (yearSuffix == "_2022preEE") {
    smearName = "EGMSmearAndSyst_ElePTsplit"+yearSuffix; 
	/*} else {
        smearName = "EGMSmearAndSyst_EleEtaR9" + yearSuffix;  
	}*/
   
    electronScaleEvaluator_ = electronCorrectionSet_->compound().at(scaleName);
    electronSmearEvaluator_ = electronCorrectionSet_->at(smearName);
}

double EGMCorrectionManager::getElectronScale(int run, double scEta, double r9, double pt, int seedGain) {
  return electronScaleEvaluator_->evaluate({"scale", static_cast<double>(run), static_cast<double>(scEta), static_cast<double>(r9), static_cast<double>(std::abs(scEta)), static_cast<double>(pt), static_cast<double>(seedGain)});
}

double EGMCorrectionManager::getElectronSmear(double pt, double r9, double scEta) {
  return electronSmearEvaluator_->evaluate({"smear",static_cast<double>(pt),static_cast<double>(r9),static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::getElectronScaleUnc(double pt, double r9, double scEta) {
  return electronSmearEvaluator_->evaluate({"escale",static_cast<double>(pt), static_cast<double>(r9), static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::getElectronSmearUnc(double pt, double r9, double scEta) {
  return electronSmearEvaluator_->evaluate({"esmear", static_cast<double>(pt), static_cast<double>(r9),static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::getPhotonScale(int run, double scEta, double r9, double pt, int seedGain) {
  return photonScaleEvaluator_->evaluate({"scale", static_cast<double>(run), static_cast<double>(scEta),static_cast<double>(r9), static_cast<double>(std::abs(scEta)),static_cast<double>(pt), static_cast<double>(seedGain)});
}

double EGMCorrectionManager::getPhotonSmear(double pt, double r9, double scEta) {
  return photonSmearEvaluator_->evaluate({"smear", static_cast<double>(pt),static_cast<double>(r9), static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::getPhotonScaleUnc(double pt, double r9, double scEta) {
  return photonSmearEvaluator_->evaluate({"escale", static_cast<double>(pt), static_cast<double>(r9), static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::getPhotonSmearUnc(double pt, double r9, double scEta) {
  return photonSmearEvaluator_->evaluate({"esmear", static_cast<double>(pt), static_cast<double>(r9), static_cast<double>(std::abs(scEta))});
}

double EGMCorrectionManager::applyCorrectedElectronPt(double originalPt, int run, double scEta, double r9, int seedGain, bool isData, double randomNum) {
    if (isData) {
        double scale = getElectronScale(run, scEta, r9, originalPt, seedGain);
        return scale * originalPt;
    } else {
        double smear = getElectronSmear(originalPt, r9, std::abs(scEta));
        return originalPt * (1.0 + smear * randomNum);
    }
}

double EGMCorrectionManager::applyCorrectedPhotonPt(double originalPt, int run, double scEta, double r9, int seedGain, bool isData, double randomNum) {

  if (isData) {
        double scale = getPhotonScale(run, scEta, r9, originalPt, seedGain);
        return scale * originalPt;
    } else {
        double smear = getPhotonSmear(originalPt, r9, std::abs(scEta));
        return originalPt * (1.0 + smear * randomNum);
    }
}

//Getting the seedGain
int EGMCorrectionManager::GetSeedGain(const DetId& seedDetId, const edm::Event& e,
                                     const edm::EDGetTokenT<EcalRecHitCollection>& ebToken,
                                     const edm::EDGetTokenT<EcalRecHitCollection>& eeToken) {
    int gain = 12; // Default ECAL gain
    if (seedDetId.subdetId() == EcalBarrel) {
            edm::Handle<EcalRecHitCollection> ebRecHits;
            e.getByToken(ebToken, ebRecHits);
            
            if (ebRecHits.isValid()) {
                EcalRecHitCollection::const_iterator rechit = ebRecHits->find(seedDetId);
                if (rechit != ebRecHits->end()) {
                    if (rechit->checkFlag(EcalRecHit::kHasSwitchToGain1)) {
                        gain = 1;
                    } else if (rechit->checkFlag(EcalRecHit::kHasSwitchToGain6)) {
                        gain = 6;
                    } else {
                        gain = 12;
                    }
                }
            }
        }
        else if (seedDetId.subdetId() == EcalEndcap) {
            edm::Handle<EcalRecHitCollection> eeRecHits;
            e.getByToken(eeToken, eeRecHits);
            
            if (eeRecHits.isValid()) {
                EcalRecHitCollection::const_iterator rechit = eeRecHits->find(seedDetId);
                if (rechit != eeRecHits->end()) {
                    if (rechit->checkFlag(EcalRecHit::kHasSwitchToGain1)) {
                        gain = 1;
                    } else if (rechit->checkFlag(EcalRecHit::kHasSwitchToGain6)) {
                        gain = 6;
                    } else {
                        gain = 12;
                    }
                }
            }
        }
    
    return gain;
}
