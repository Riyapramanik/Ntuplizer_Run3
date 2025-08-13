#ifndef EGMCorrectionManager_h
#define EGMCorrectionManager_h
#include "correction.h"
#include <string>
#include <memory>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h" 
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


class EGMCorrectionManager {
public:
    EGMCorrectionManager(int year, const std::string& period = "");
    ~EGMCorrectionManager() = default;
    
    // Apply corrections and get corrected pT
    double applyCorrectedElectronPt(double originalPt, int run, double scEta, double r9, 
                                   int seedGain, bool isData, double randomNum = 0.0);
    double applyCorrectedPhotonPt(double originalPt, int run, double scEta, double r9, 
                                 int seedGain, bool isData, double randomNum = 0.0);
    
    // Electron corrections
    double getElectronScale(int run, double scEta, double r9, double pt, int seedGain);
    double getElectronSmear(double pt, double r9, double absScEta);
    double getElectronScaleUnc(double pt, double r9, double absScEta);
    double getElectronSmearUnc(double pt, double r9, double absScEta);
    
    // Photon corrections  
    double getPhotonScale(int run, double scEta, double r9, double pt, int seedGain);
    double getPhotonSmear(double pt, double r9, double absScEta);
    double getPhotonScaleUnc(double pt, double r9, double absScEta);
    double getPhotonSmearUnc(double pt, double r9, double absScEta);

    //Getting seed Gain
     static int GetSeedGain(const DetId& seedDetId, const edm::Event& e, 
                          const edm::EDGetTokenT<EcalRecHitCollection>& ebToken,
                          const edm::EDGetTokenT<EcalRecHitCollection>& eeToken);

private:
    
    // Setup functions
    void initializeCorrections();
    void setupElectronEvaluators();
    void setupPhotonEvaluators();
    
    // File path builders (return uncompressed .json paths)
    std::string getElectronJSONFile();
    std::string getPhotonJSONFile();
    std::string getYearSuffix();
    
    // Utility (checks if uncompressed file exists, unzips .gz if needed)
    void unzipFileIfNeeded(const std::string& filePath);
        
    // Configuration
    int year_;
    std::string period_;
    bool useETDependent_;
    
    // Correction data loaded from JSON files
    std::unique_ptr<correction::CorrectionSet> electronCorrectionSet_;
    std::unique_ptr<correction::CorrectionSet> photonCorrectionSet_;
    
    // The actual correction functions
    correction::CompoundCorrection::Ref electronScaleEvaluator_;  // For data scale corrections
    correction::Correction::Ref electronSmearEvaluator_; // For MC smearing corrections  
    correction::CompoundCorrection::Ref photonScaleEvaluator_;   // For data scale correction
    correction::Correction::Ref photonSmearEvaluator_;   // For MC smearing corrections
};

#endif
