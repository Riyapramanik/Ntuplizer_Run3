#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
// Pileup reweight variables
Float_t     rho_;
Int_t       run_;
Long64_t    event_;
UShort_t    lumis_;
// Primary vertex variables
UChar_t     nVtx_;
UChar_t     nGoodVtx_;
Bool_t      isPVGood_;
Float_t     vtx_;
Float_t     vty_;
Float_t     vtz_;

// HLT decisions stored in 6 UShort_t variables (16 bits each)
UShort_t    HLTriggerWord0_;  // bits 0-15   (triggers 0-15)
UShort_t    HLTriggerWord1_;  // bits 16-31  (triggers 16-31)
UShort_t    HLTriggerWord2_;  // bits 32-47  (triggers 32-47)
UShort_t    HLTriggerWord3_;  // bits 48-63  (triggers 48-63)
UShort_t    HLTriggerWord4_;  // bits 64-79  (triggers 64-79)
UShort_t    HLTriggerWord5_;  // bits 80-93  (triggers 80-93, uses 14 bits)

const int nHLTmx = 94;
Bool_t booltrg[nHLTmx];

const char *hlt_name[nHLTmx] = {
  // photon triggers
  "HLT_Photon100EBHE10", "HLT_Photon110EB_TightID_TightIso", "HLT_Photon120_R9Id90_HE10_IsoM", "HLT_Photon120",
  "HLT_Photon130EB_TightID_TightIso", "HLT_Photon150EB_TightID_TightIso", "HLT_Photon150", "HLT_Photon165_R9Id90_HE10_IsoM",
  "HLT_Photon175EB_TightID_TightIso", "HLT_Photon175", "HLT_Photon200EB_TightID_TightIso", "HLT_Photon200",
  "HLT_Photon20_HoverELoose", "HLT_Photon300_NoHE", "HLT_Photon30EB_TightID_TightIso", "HLT_Photon30_HoverELoose",
  "HLT_Photon32_OneProng32_M50To105", "HLT_Photon33", "HLT_Photon35_TwoProngs35", "HLT_Photon50EB_TightID_TightIso",
  "HLT_Photon50_R9Id90_HE10_IsoM", "HLT_Photon50_TimeGt2p5ns", "HLT_Photon50_TimeLtNeg2p5ns", "HLT_Photon50",
  "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350", "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380",
  "HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400", "HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
  "HLT_Photon75EB_TightID_TightIso", "HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3",
  "HLT_Photon75_R9Id90_HE10_IsoM", "HLT_Photon75", "HLT_Photon90EB_TightID_TightIso", "HLT_Photon90_R9Id90_HE10_IsoM",
  "HLT_Photon90",

  // dielectron triggers
  "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", "HLT_DoubleEle33_CaloIdL_MW",

  // single electron triggers
  "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele135_CaloIdVT_GsfTrkIdT", "HLT_Ele30_WPTight_Gsf",
  "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
  "HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_Ele38_WPTight_Gsf", "HLT_Ele40_WPTight_Gsf",
  "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06", "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40",
  "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10", "HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40",
  "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", "HLT_Ele50_IsoVVVL_PFHT450",

  // PFJet triggers
  "HLT_PFJet110", "HLT_PFJet140", "HLT_PFJet200_TimeGt2p5ns", "HLT_PFJet200_TimeLtNeg2p5ns", "HLT_PFJet200",
  "HLT_PFJet260", "HLT_PFJet320", "HLT_PFJet400", "HLT_PFJet40_GPUvsCPU", "HLT_PFJet40",
  "HLT_PFJet450", "HLT_PFJet500", "HLT_PFJet550", "HLT_PFJet60", "HLT_PFJet80",

  // PFJet forward triggers
  "HLT_PFJetFwd140", "HLT_PFJetFwd200", "HLT_PFJetFwd260", "HLT_PFJetFwd320", "HLT_PFJetFwd400",
  "HLT_PFJetFwd440", "HLT_PFJetFwd450", "HLT_PFJetFwd500",

  // PFMET triggers
  "HLT_PFMET105_IsoTrk50", "HLT_PFMET120_PFMHT120_IDTight_PFHT60", "HLT_PFMET120_PFMHT120_IDTight",
  "HLT_PFMET130_PFMHT130_IDTight", "HLT_PFMET140_PFMHT140_IDTight", "HLT_PFMET200_BeamHaloCleaned",
  "HLT_PFMET200_NotCleaned", "HLT_PFMET250_NotCleaned", "HLT_PFMET300_NotCleaned",
  "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF",
  "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
  "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight",
  "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
  "HLT_PFMETTypeOne140_PFMHT140_IDTight", "HLT_PFMETTypeOne200_BeamHaloCleaned"
};

void ggNtuplizer::branchesGlobalEvent(TTree* tree) {
  
  tree->Branch("run",                    & run_);
  tree->Branch("event",                  & event_);
  tree->Branch("lumis",                  & lumis_);
  tree->Branch("nVtx",                   & nVtx_);
  tree->Branch("nGoodVtx",               & nGoodVtx_);
  tree->Branch("isPVGood",               & isPVGood_);
  tree->Branch("vtx",                    & vtx_);
  tree->Branch("vty",                    & vty_);
  tree->Branch("vtz",                    & vtz_);
  tree->Branch("rho",                    &rho_);

// HLT trigger decisions stored in 6 UShort_t variables (16 bits each)
  tree->Branch("HLTriggerWord0",         & HLTriggerWord0_);   // triggers 0-15
  tree->Branch("HLTriggerWord1",         & HLTriggerWord1_);   // triggers 16-31
  tree->Branch("HLTriggerWord2",         & HLTriggerWord2_);   // triggers 32-47
  tree->Branch("HLTriggerWord3",         & HLTriggerWord3_);   // triggers 48-63
  tree->Branch("HLTriggerWord4",         & HLTriggerWord4_);   // triggers 64-79
  tree->Branch("HLTriggerWord5",         & HLTriggerWord5_);   // triggers 80-93

}

void ggNtuplizer::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  nVtx_     = -1;
  nGoodVtx_ = -1;
  if (vtxHandle.isValid()) {
    nVtx_     = 0;
    nGoodVtx_ = 0;
    
    for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
      
      if (nVtx_ == 0) {
        vtx_     = v->x();
        vty_     = v->y();
        vtz_     = v->z();

        isPVGood_ = false;
        if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) isPVGood_ = true;
      }

      if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
      nVtx_++;

    }
  } else edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  rho_ = *(rhoHandle.product());
  
  edm::Handle<edm::TriggerResults> trigRes;	
  e.getByToken(triggerBits_, trigRes);
  const edm::TriggerNames &names_ = e.triggerNames(*trigRes);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  e.getByToken(triggerObjects_, triggerObjects);
  
  const char* variab_trig;

  // Initialize all HLT words to 0
  HLTriggerWord0_ = 0;
  HLTriggerWord1_ = 0;
  HLTriggerWord2_ = 0;
  HLTriggerWord3_ = 0;
  HLTriggerWord4_ = 0;
  HLTriggerWord5_ = 0;
  
  // Initialize trigger arrays                                                                                          
  for(int jk=0; jk<nHLTmx; jk++) {
    booltrg[jk] = false;
  }
  
  for (int jk=0; jk<nHLTmx; jk++) {
    
    for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
      std::string name = names_.triggerName(ij);
      variab_trig = name.c_str();
      if (strstr(variab_trig, hlt_name[jk]) && ((strlen(variab_trig)-strlen(hlt_name[jk]))<5))
      {	
        if ((trigRes->accept(ij))){   //||(isMC)) {
          booltrg[jk] = true;
	  if (jk < 16) {
            setbit(HLTriggerWord0_, jk);
          } else if (jk < 32) {
            setbit(HLTriggerWord1_, jk - 16);
          } else if (jk < 48) {
            setbit(HLTriggerWord2_, jk - 32);
          } else if (jk < 64) {
            setbit(HLTriggerWord3_, jk - 48);
          } else if (jk < 80) {
            setbit(HLTriggerWord4_, jk - 64);
          } else {
            setbit(HLTriggerWord5_, jk - 80);
          }
          break;
	}
      }
    }//ij   
  }
    
}

