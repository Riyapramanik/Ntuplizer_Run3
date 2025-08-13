#include "ggAnalysis/ggNtuplizer/interface/EGMIDSFManager.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include "ggAnalysis/ggNtuplizer/interface/EGMIDSFManager.h"

EGMIDSFManager::EGMIDSFManager(int year, const std::string& period)
    : year_(year), period_(period) {
    initializeCorrections();
}


void EGMIDSFManager::initializeCorrections() {

  std::cout<<"======initializeCorrections========="<<std::endl;
  
    std::string electronIDFile = getElectronIDJSONFile();
    std::string photonIDFile = getPhotonIDJSONFile();
    
    electronIDCorrectionSet_ = correction::CorrectionSet::from_file(electronIDFile);
    photonIDCorrectionSet_ = correction::CorrectionSet::from_file(photonIDFile);
    
    setupElectronIDEvaluators();
    setupPhotonIDEvaluators();
}


std::string EGMIDSFManager::getElectronIDJSONFile() {

  std::cout<<"======EGMIDSFManager.cc========="<<std::endl;
  std::string basePath = "/eos/user/r/rpramani/run3_ntuplizer/CMSSW_14_2_1/src/ggAnalysis/ggNtuplizer/test/ElectronJson/";

  if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D") {
            return basePath + "electronID_highPt_22_BCD.json";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G") {
            return basePath + "electronID_highPt_22_EFG.json";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C") {
            return basePath + "electronID_highPt_23_C.json";
        }
        else if (period_ == "D") {
            return basePath + "electronID_highPt_23_D.json";
        }
    }
    
    return ""; 
}


std::string EGMIDSFManager::getPhotonIDJSONFile() {
    std::string basePath = "/eos/user/r/rpramani/run3_ntuplizer/CMSSW_14_2_1/src/ggAnalysis/ggNtuplizer/test/PhotonJson/";
    
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D") {
            return basePath + "photonID_highPt_22_BCD.json";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G") {
            return basePath + "photonID_highPt_22_EFG.json";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C") {
            return basePath + "photonID_highPt_23_C.json";
        }
        else if (period_ == "D") {
            return basePath + "photonID_highPt_23_D.json";
        }
    }
    
    return "";  
}


void EGMIDSFManager::setupElectronIDEvaluators() {
    std::vector<std::string> idTypes = {
        "Loose", "Medium", "Tight", "wp80iso", "wp90iso"
    };

    std::cout<<"****************setupElectronIDEvaluators*******"<<std::endl;
    
    std::string correctionName = "Electron-ID-SF";
    auto correctionRef = electronIDCorrectionSet_->at(correctionName);
        for (const auto& idType : idTypes) {
            electronIDEvaluators_[idType] = correctionRef;
	}
}

void EGMIDSFManager::setupPhotonIDEvaluators() {
    std::vector<std::string> idTypes = {
        "Loose", "Medium", "Tight", "wp80", "wp90"
    };

    std::cout<<"****************setupPhotonIDEvaluators*******"<<std::endl;
    std::cout << "Available corrections in photon ID file:" << std::endl;
    for (const auto& [name, corr] : *photonIDCorrectionSet_) {
        std::cout << "  Available correction: " << name << std::endl;
    }
    
    std::string correctionName = "Photon-ID-SF";
    auto correctionRef = photonIDCorrectionSet_->at(correctionName);
        for (const auto& idType : idTypes) {
            photonIDEvaluators_[idType] = correctionRef;
        }
}

std::string EGMIDSFManager::getYearValue() {
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D") {
            return "2022Re-recoBCD";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G") {
            return "2022Re-recoE+PromptFG";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C") {
            return "2023PromptC";
        }
        else if (period_ == "D" ) {
            return "2023PromptD";
        }
    }
    
    return ""; // Default
}

double EGMIDSFManager::getElectronIDSF(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sf", idType, eta, pt});
}

double EGMIDSFManager::getElectronIDSFUncUp(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sfup", idType, eta, pt});
}

double EGMIDSFManager::getElectronIDSFUncDown(const std::string& idType, double pt, double eta, double phi) {
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sfdown", idType, eta, pt});
}

double EGMIDSFManager::getPhotonIDSF(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    std::cout << "***getPhotonIDSF***"<<std::endl;
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sf", idType, eta, pt});
}

double EGMIDSFManager::getPhotonIDSFUncUp(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sfup", idType, eta, pt});
}

double EGMIDSFManager::getPhotonIDSFUncDown(const std::string& idType, double pt, double eta, double phi) {   
    std::string yearValue = getYearValue();
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sfdown", idType, eta, pt});
}
