#include "ggAnalysis/ggNtuplizer/interface/EGMIDSFManager.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cstdlib>
#include "ggAnalysis/ggNtuplizer/interface/EGMIDSFManager.h"

EGMIDSFManager::EGMIDSFManager(int year, const std::string& period, bool useHTT)
    : year_(year), period_(period), useHTT_(useHTT) {
    initializeCorrections();
}


void EGMIDSFManager::initializeCorrections() {
    std::string electronIDFile = getElectronIDJSONFile();
    std::string photonIDFile = getPhotonIDJSONFile();
    
    unzipFileIfNeeded(electronIDFile);
    unzipFileIfNeeded(photonIDFile);
    
    electronIDCorrectionSet_ = correction::CorrectionSet::from_file(electronIDFile);
    photonIDCorrectionSet_ = correction::CorrectionSet::from_file(photonIDFile);
    
    setupElectronIDEvaluators();
    setupPhotonIDEvaluators();
}


std::string EGMIDSFManager::getElectronIDJSONFile() {

  std::cout<<"======EGMIDSFManager.cc========="<<std::endl;
  
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D" || 
            period_ == "preEE" || period_.empty()) {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoBCD/electron.json.gz";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G" || period_ == "postEE") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoE+PromptFG/electron.json.gz";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C" || period_ == "preBPIX") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2023_PromptCD/2023PromptC/electron.json.gz";
        }
        else if (period_ == "D" || period_ == "postBPIX") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2023_PromptCD/2023PromptD/electron.json.gz";
        }
    }
    
    throw std::runtime_error("Year " + std::to_string(year_) + " period " + period_ + " not supported for electron ID");
}


std::string EGMIDSFManager::getPhotonIDJSONFile() {
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D" || 
            period_ == "preEE" || period_.empty()) {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoBCD/photon.json.gz";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G" || period_ == "postEE") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2022_recoBCDE_PromptFG/2022Re-recoE+PromptFG/photon.json.gz";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C" || period_ == "preBPIX") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2023_PromptCD/2023PromptC/photon.json.gz";
        }
        else if (period_ == "D" || period_ == "postBPIX") {
            return "/eos/cms/store/group/phys_egamma/correctionlibJSONs/Run3_2023_PromptCD/2023PromptD/photon.json.gz";
        }
    }
    
    throw std::runtime_error("Year " + std::to_string(year_) + " period " + period_ + " not supported for photon ID");
}

void EGMIDSFManager::unzipFileIfNeeded(const std::string& filePath) {
    std::ifstream file(filePath);
    if (file.good()) {
        return;
    }
    
    std::string gzFilePath = filePath + ".gz";
    std::string unzipCmd = "gunzip -k " + gzFilePath;
    int result = system(unzipCmd.c_str());
    
    if (result != 0) {
        throw std::runtime_error("Failed to unzip file: " + gzFilePath);
    }
}


void EGMIDSFManager::setupElectronIDEvaluators() {
  std::vector<std::string> idTypes = {
        "Loose", "Medium", "Tight", "wp80iso", "wp90iso"
    };
    
    std::string correctionName = "Electron-ID-SF";
    
    try {
        auto correctionRef = electronIDCorrectionSet_->at(correctionName);
        for (const auto& idType : idTypes) {
            electronIDEvaluators_[idType] = correctionRef;
        }
    } catch (const std::exception& e) {
        std::cout << "Warning: Could not find " << correctionName << " in electron ID file: " << e.what() << std::endl;
    }
}


void EGMIDSFManager::setupPhotonIDEvaluators() {
    // Setup evaluators for photon ID working points
    std::vector<std::string> idTypes = {
        "Loose", "Medium", "Tight", "wp80", "wp90"
    };
    
    std::string correctionName = "Photon-ID-SF";
    
    try {
        auto correctionRef = photonIDCorrectionSet_->at(correctionName);
        for (const auto& idType : idTypes) {
            photonIDEvaluators_[idType] = correctionRef;
        }
    } catch (const std::exception& e) {
        std::cout << "Warning: Could not find " << correctionName << " in photon ID file: " << e.what() << std::endl;
    }
}

std::string EGMIDSFManager::getYearValue() {
    if (year_ == 2022) {
        if (period_ == "B" || period_ == "C" || period_ == "D" || 
            period_ == "preEE" || period_.empty()) {
            return "2022Re-recoBCD";
        }
        else if (period_ == "E" || period_ == "F" || period_ == "G" || period_ == "postEE") {
            return "2022Re-recoE+PromptFG";
        }
    } 
    else if (year_ == 2023) {
        if (period_ == "C" || period_ == "preBPIX") {
            return "2023PromptC";
        }
        else if (period_ == "D" || period_ == "postBPIX") {
            return "2023PromptD";
        }
    }
    
    throw std::runtime_error("Year " + std::to_string(year_) + " period " + period_ + " not supported");
}

double EGMIDSFManager::getElectronIDSF(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sf", idType, eta, pt, phi});
}

double EGMIDSFManager::getElectronIDSFUncUp(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sfup", idType, eta, pt, phi});
}

double EGMIDSFManager::getElectronIDSFUncDown(const std::string& idType, double pt, double eta, double phi) {
    std::string yearValue = getYearValue();
    return electronIDEvaluators_[idType]->evaluate({yearValue, "sfdown", idType, eta, pt, phi});
}

// Photon ID Scale Factors
double EGMIDSFManager::getPhotonIDSF(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sf", idType, eta, pt, phi});
}

double EGMIDSFManager::getPhotonIDSFUncUp(const std::string& idType, double pt, double eta, double phi) {    
    std::string yearValue = getYearValue();
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sfup", idType, eta, pt, phi});
}

double EGMIDSFManager::getPhotonIDSFUncDown(const std::string& idType, double pt, double eta, double phi) {   
    std::string yearValue = getYearValue();
    return photonIDEvaluators_[idType]->evaluate({yearValue, "sfdown", idType, eta, pt, phi});
}
