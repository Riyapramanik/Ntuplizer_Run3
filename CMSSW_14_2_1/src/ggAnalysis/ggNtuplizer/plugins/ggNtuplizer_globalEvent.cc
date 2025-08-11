#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
// Pileup reweight variables
Float_t     rho_;
Bool_t      isData_;
// Photon triggers
Bool_t      HLT_Photon165_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon175EB_TightID_TightIso_;
Bool_t      HLT_Photon175_;
Bool_t      HLT_Photon200EB_TightID_TightIso_;
Bool_t      HLT_Photon200_;

//MET triggers  
Bool_t      HLT_PFMET200_BeamHaloCleaned_;
Bool_t      HLT_PFMET200_NotCleaned_;
Bool_t      HLT_PFMET250_NotCleaned_;
Bool_t      HLT_PFMET300_NotCleaned_;

// Prescale information
Int_t       HLT_Photon165_R9Id90_HE10_IsoM_Prescale_;
Int_t       HLT_Photon175EB_TightID_TightIso_Prescale_;
Int_t       HLT_Photon175_Prescale_;
Int_t       HLT_Photon200EB_TightID_TightIso_Prescale_;
Int_t       HLT_Photon200_Prescale_;
Int_t       HLT_PFMET200_BeamHaloCleaned_Prescale_;
Int_t       HLT_PFMET200_NotCleaned_Prescale_;
Int_t       HLT_PFMET250_NotCleaned_Prescale_;
Int_t       HLT_PFMET300_NotCleaned_Prescale_;

const int nHLTmx = 9;
Bool_t booltrg[nHLTmx];
Int_t prescaletrg[nHLTmx];

const char *hlt_name[nHLTmx]  = {
  "HLT_Photon165_R9Id90_HE10_IsoM",        // index 0
  "HLT_Photon175EB_TightID_TightIso",      // index 1  
  "HLT_Photon175",                         // index 2
  "HLT_Photon200EB_TightID_TightIso",      // index 3
  "HLT_Photon200",                         // index 4
  "HLT_PFMET200_BeamHaloCleaned",          // index 5
  "HLT_PFMET200_NotCleaned",               // index 6
  "HLT_PFMET250_NotCleaned",               // index 7
  "HLT_PFMET300_NotCleaned"                // index 8
};

void ggNtuplizer::branchesGlobalEvent(TTree* tree) {

   // Pileup reweight branches
  tree->Branch("rho",                               &rho_);

   // Photon trigger branches
  tree->Branch("HLT_Photon165_R9Id90_HE10_IsoM",        &HLT_Photon165_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon175EB_TightID_TightIso",      &HLT_Photon175EB_TightID_TightIso_);
  tree->Branch("HLT_Photon175",                         &HLT_Photon175_);
  tree->Branch("HLT_Photon200EB_TightID_TightIso",      &HLT_Photon200EB_TightID_TightIso_);
  tree->Branch("HLT_Photon200",                         &HLT_Photon200_);

  // MET trigger branches
  tree->Branch("HLT_PFMET200_BeamHaloCleaned",          &HLT_PFMET200_BeamHaloCleaned_);
  tree->Branch("HLT_PFMET200_NotCleaned",               &HLT_PFMET200_NotCleaned_);
  tree->Branch("HLT_PFMET250_NotCleaned",               &HLT_PFMET250_NotCleaned_);
  tree->Branch("HLT_PFMET300_NotCleaned",               &HLT_PFMET300_NotCleaned_);

  // Prescale branches
  tree->Branch("HLT_Photon165_R9Id90_HE10_IsoM_Prescale",        &HLT_Photon165_R9Id90_HE10_IsoM_Prescale_);
  tree->Branch("HLT_Photon175EB_TightID_TightIso_Prescale",      &HLT_Photon175EB_TightID_TightIso_Prescale_);
  tree->Branch("HLT_Photon175_Prescale",                         &HLT_Photon175_Prescale_);
  tree->Branch("HLT_Photon200EB_TightID_TightIso_Prescale",      &HLT_Photon200EB_TightID_TightIso_Prescale_);
  tree->Branch("HLT_Photon200_Prescale",                         &HLT_Photon200_Prescale_);
  tree->Branch("HLT_PFMET200_BeamHaloCleaned_Prescale",          &HLT_PFMET200_BeamHaloCleaned_Prescale_);
  tree->Branch("HLT_PFMET200_NotCleaned_Prescale",               &HLT_PFMET200_NotCleaned_Prescale_);
  tree->Branch("HLT_PFMET250_NotCleaned_Prescale",               &HLT_PFMET250_NotCleaned_Prescale_);
  tree->Branch("HLT_PFMET300_NotCleaned_Prescale",               &HLT_PFMET300_NotCleaned_Prescale_);
}

void ggNtuplizer::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  // Pileup reweight information
  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  rho_ = *(rhoHandle.product());
  
  edm::Handle<edm::TriggerResults> trigRes;	
  e.getByToken(triggerBits_, trigRes);
  const edm::TriggerNames &names_ = e.triggerNames(*trigRes);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  e.getByToken(triggerObjects_, triggerObjects);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  e.getByToken(triggerPrescales_, triggerPrescales);
  
  const char* variab_trig;
  // Initialize trigger arrays                                                                                          
  for(int jk=0; jk<nHLTmx; jk++) {
    booltrg[jk] = false;
    prescaletrg[jk] = 1;
  }

  for (int jk=0; jk<nHLTmx; jk++) {
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names_.triggerName(ij);
      variab_trig = name.c_str(); 
      if (strstr(variab_trig, hlt_name[jk]) && ((strlen(variab_trig)-strlen(hlt_name[jk]))<5))
      {
        if ((trigRes->accept(ij))){   //||(isMC)) {
          booltrg[jk] = true; booltrg[nHLTmx] = true;
          // Get prescale
          if (triggerPrescales.isValid()) {
            prescaletrg[jk] = triggerPrescales->getPrescaleForIndex<double>(ij);
          }
          break;
        }
      }
    }//ij     
  }
  
  for(int jk=0; jk<nHLTmx; jk++) {
    // Photon triggers
    if(jk==0)      { HLT_Photon165_R9Id90_HE10_IsoM_ = booltrg[jk]; 
                     HLT_Photon165_R9Id90_HE10_IsoM_Prescale_ = prescaletrg[jk]; }
    else if(jk==1) { HLT_Photon175EB_TightID_TightIso_ = booltrg[jk]; 
                     HLT_Photon175EB_TightID_TightIso_Prescale_ = prescaletrg[jk]; }
    else if(jk==2) { HLT_Photon175_ = booltrg[jk]; 
                     HLT_Photon175_Prescale_ = prescaletrg[jk]; }
    else if(jk==3) { HLT_Photon200EB_TightID_TightIso_ = booltrg[jk]; 
                     HLT_Photon200EB_TightID_TightIso_Prescale_ = prescaletrg[jk]; }
    else if(jk==4) { HLT_Photon200_ = booltrg[jk]; 
                     HLT_Photon200_Prescale_ = prescaletrg[jk]; }
    
    // MET triggers
    else if(jk==5) { HLT_PFMET200_BeamHaloCleaned_ = booltrg[jk]; 
                     HLT_PFMET200_BeamHaloCleaned_Prescale_ = prescaletrg[jk]; }
    else if(jk==6) { HLT_PFMET200_NotCleaned_ = booltrg[jk]; 
                     HLT_PFMET200_NotCleaned_Prescale_ = prescaletrg[jk]; }
    else if(jk==7) { HLT_PFMET250_NotCleaned_ = booltrg[jk]; 
                     HLT_PFMET250_NotCleaned_Prescale_ = prescaletrg[jk]; }
    else if(jk==8) { HLT_PFMET300_NotCleaned_ = booltrg[jk]; 
                     HLT_PFMET300_NotCleaned_Prescale_ = prescaletrg[jk]; }
  }
}
