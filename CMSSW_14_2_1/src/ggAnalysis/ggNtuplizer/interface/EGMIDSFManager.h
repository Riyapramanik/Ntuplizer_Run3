#ifndef EGMIDSFManager_h
#define EGMIDSFManager_h

#include <string>
#include <memory>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "correction.h"

class EGMIDSFManager {
public:
    // Constructor
    EGMIDSFManager(int year, const std::string& period);
    
    // Electron ID Scale Factors (always require phi for accuracy)
    double getElectronIDSF(const std::string& idType, double pt, double eta,double phi);
    double getElectronIDSFUncUp(const std::string& idType, double pt, double eta, double phi);
    double getElectronIDSFUncDown(const std::string& idType, double pt, double eta, double phi);
    
    // Photon ID Scale Factors (always require phi for accuracy)
    double getPhotonIDSF(const std::string& idType, double pt, double eta, double phi);
    double getPhotonIDSFUncUp(const std::string& idType, double pt, double eta, double phi);
    double getPhotonIDSFUncDown(const std::string& idType, double pt, double eta, double phi);

private:
    int year_;
    std::string period_;
    std::unique_ptr<correction::CorrectionSet> electronIDCorrectionSet_;
    std::unique_ptr<correction::CorrectionSet> photonIDCorrectionSet_;
    std::map<std::string, correction::Correction::Ref> electronIDEvaluators_;
    std::map<std::string, correction::Correction::Ref> photonIDEvaluators_;

    void initializeCorrections();
    std::string getElectronIDJSONFile();
    std::string getPhotonIDJSONFile();
    void unzipFileIfNeeded(const std::string& filePath);
    void setupElectronIDEvaluators();
    void setupPhotonIDEvaluators();
    std::string getYearValue();
};

#endif
