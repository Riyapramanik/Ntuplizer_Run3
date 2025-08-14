#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
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

// Photon triggers
Bool_t      HLT_Photon100EBHE10_;
Bool_t      HLT_Photon110EB_TightID_TightIso_;
Bool_t      HLT_Photon120_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon120_;
Bool_t      HLT_Photon130EB_TightID_TightIso_;
Bool_t      HLT_Photon150EB_TightID_TightIso_;
Bool_t      HLT_Photon150_;
Bool_t      HLT_Photon165_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon175EB_TightID_TightIso_;
Bool_t      HLT_Photon175_;
Bool_t      HLT_Photon200EB_TightID_TightIso_;
Bool_t      HLT_Photon200_;
Bool_t      HLT_Photon20_HoverELoose_;
Bool_t      HLT_Photon300_NoHE_;
Bool_t      HLT_Photon30EB_TightID_TightIso_;
Bool_t      HLT_Photon30_HoverELoose_;
Bool_t      HLT_Photon32_OneProng32_M50To105_;
Bool_t      HLT_Photon33_;
Bool_t      HLT_Photon35_TwoProngs35_;
Bool_t      HLT_Photon50EB_TightID_TightIso_;
Bool_t      HLT_Photon50_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon50_TimeGt2p5ns_;
Bool_t      HLT_Photon50_TimeLtNeg2p5ns_;
Bool_t      HLT_Photon50_;
Bool_t      HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350_;
Bool_t      HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380_;
Bool_t      HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400_;
Bool_t      HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_;
Bool_t      HLT_Photon75EB_TightID_TightIso_;
Bool_t      HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_;
Bool_t      HLT_Photon75_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon75_;
Bool_t      HLT_Photon90EB_TightID_TightIso_;
Bool_t      HLT_Photon90_R9Id90_HE10_IsoM_;
Bool_t      HLT_Photon90_;

// Dielectron triggers
Bool_t      HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_;
Bool_t      HLT_DoubleEle33_CaloIdL_MW_;

// Single electron triggers
Bool_t      HLT_Ele115_CaloIdVT_GsfTrkIdT_;
Bool_t      HLT_Ele135_CaloIdVT_GsfTrkIdT_;
Bool_t      HLT_Ele30_WPTight_Gsf_;
Bool_t      HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_;
Bool_t      HLT_Ele32_WPTight_Gsf_L1DoubleEG_;
Bool_t      HLT_Ele32_WPTight_Gsf_;
Bool_t      HLT_Ele35_WPTight_Gsf_;
Bool_t      HLT_Ele38_WPTight_Gsf_;
Bool_t      HLT_Ele40_WPTight_Gsf_;
Bool_t      HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06_;
Bool_t      HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_;
Bool_t      HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10_;
Bool_t      HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_;
Bool_t      HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_;
Bool_t      HLT_Ele50_IsoVVVL_PFHT450_;

// PFJet triggers
Bool_t      HLT_PFJet110_;
Bool_t      HLT_PFJet140_;
Bool_t      HLT_PFJet200_TimeGt2p5ns_;
Bool_t      HLT_PFJet200_TimeLtNeg2p5ns_;
Bool_t      HLT_PFJet200_;
Bool_t      HLT_PFJet260_;
Bool_t      HLT_PFJet320_;
Bool_t      HLT_PFJet400_;
Bool_t      HLT_PFJet40_GPUvsCPU_;
Bool_t      HLT_PFJet40_;
Bool_t      HLT_PFJet450_;
Bool_t      HLT_PFJet500_;
Bool_t      HLT_PFJet550_;
Bool_t      HLT_PFJet60_;
Bool_t      HLT_PFJet80_;

// PFJet forward triggers
Bool_t      HLT_PFJetFwd140_;
Bool_t      HLT_PFJetFwd200_;
Bool_t      HLT_PFJetFwd260_;
Bool_t      HLT_PFJetFwd320_;
Bool_t      HLT_PFJetFwd400_;
Bool_t      HLT_PFJetFwd440_;
Bool_t      HLT_PFJetFwd450_;
Bool_t      HLT_PFJetFwd500_;

// PFMET triggers
Bool_t      HLT_PFMET105_IsoTrk50_;
Bool_t      HLT_PFMET120_PFMHT120_IDTight_PFHT60_;
Bool_t      HLT_PFMET120_PFMHT120_IDTight_;
Bool_t      HLT_PFMET130_PFMHT130_IDTight_;
Bool_t      HLT_PFMET140_PFMHT140_IDTight_;
Bool_t      HLT_PFMET200_BeamHaloCleaned_;
Bool_t      HLT_PFMET200_NotCleaned_;
Bool_t      HLT_PFMET250_NotCleaned_;
Bool_t      HLT_PFMET300_NotCleaned_;
Bool_t      HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_;
Bool_t      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_;
Bool_t      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_;
Bool_t      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_;
Bool_t      HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_;
Bool_t      HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_;
Bool_t      HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_;
Bool_t      HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_;
Bool_t      HLT_PFMETTypeOne140_PFMHT140_IDTight_;
Bool_t      HLT_PFMETTypeOne200_BeamHaloCleaned_;

const int nHLTmx = 94;
Bool_t booltrg[nHLTmx];
Int_t prescaletrg[nHLTmx];



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
    
// Photon trigger branches
  tree->Branch("HLT_Photon100EBHE10", &HLT_Photon100EBHE10_);
  tree->Branch("HLT_Photon110EB_TightID_TightIso", &HLT_Photon110EB_TightID_TightIso_);
  tree->Branch("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon120", &HLT_Photon120_);
  tree->Branch("HLT_Photon130EB_TightID_TightIso", &HLT_Photon130EB_TightID_TightIso_);
  tree->Branch("HLT_Photon150EB_TightID_TightIso", &HLT_Photon150EB_TightID_TightIso_);
  tree->Branch("HLT_Photon150", &HLT_Photon150_);
  tree->Branch("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon175EB_TightID_TightIso", &HLT_Photon175EB_TightID_TightIso_);
  tree->Branch("HLT_Photon175", &HLT_Photon175_);
  tree->Branch("HLT_Photon200EB_TightID_TightIso", &HLT_Photon200EB_TightID_TightIso_);
  tree->Branch("HLT_Photon200", &HLT_Photon200_);
  tree->Branch("HLT_Photon20_HoverELoose", &HLT_Photon20_HoverELoose_);
  tree->Branch("HLT_Photon300_NoHE", &HLT_Photon300_NoHE_);
  tree->Branch("HLT_Photon30EB_TightID_TightIso", &HLT_Photon30EB_TightID_TightIso_);
  tree->Branch("HLT_Photon30_HoverELoose", &HLT_Photon30_HoverELoose_);
  tree->Branch("HLT_Photon32_OneProng32_M50To105", &HLT_Photon32_OneProng32_M50To105_);
  tree->Branch("HLT_Photon33", &HLT_Photon33_);
  tree->Branch("HLT_Photon35_TwoProngs35", &HLT_Photon35_TwoProngs35_);
  tree->Branch("HLT_Photon50EB_TightID_TightIso", &HLT_Photon50EB_TightID_TightIso_);
  tree->Branch("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon50_TimeGt2p5ns", &HLT_Photon50_TimeGt2p5ns_);
  tree->Branch("HLT_Photon50_TimeLtNeg2p5ns", &HLT_Photon50_TimeLtNeg2p5ns_);
  tree->Branch("HLT_Photon50", &HLT_Photon50_);
  tree->Branch("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350_);
  tree->Branch("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380_);
  tree->Branch("HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400", &HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400_);
  tree->Branch("HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_);
  tree->Branch("HLT_Photon75EB_TightID_TightIso", &HLT_Photon75EB_TightID_TightIso_);
  tree->Branch("HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3", &HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_);
  tree->Branch("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon75", &HLT_Photon75_);
  tree->Branch("HLT_Photon90EB_TightID_TightIso", &HLT_Photon90EB_TightID_TightIso_);
  tree->Branch("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM_);
  tree->Branch("HLT_Photon90", &HLT_Photon90_);

  // Dielectron trigger branches
  tree->Branch("HLT_DiEle27_WPTightCaloOnly_L1DoubleEG", &HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_);
  tree->Branch("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW_);

  // Single electron trigger branches
  tree->Branch("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT_);
  tree->Branch("HLT_Ele135_CaloIdVT_GsfTrkIdT", &HLT_Ele135_CaloIdVT_GsfTrkIdT_);
  tree->Branch("HLT_Ele30_WPTight_Gsf", &HLT_Ele30_WPTight_Gsf_);
  tree->Branch("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned", &HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_);
  tree->Branch("HLT_Ele32_WPTight_Gsf_L1DoubleEG", &HLT_Ele32_WPTight_Gsf_L1DoubleEG_);
  tree->Branch("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf_);
  tree->Branch("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf_);
  tree->Branch("HLT_Ele38_WPTight_Gsf", &HLT_Ele38_WPTight_Gsf_);
  tree->Branch("HLT_Ele40_WPTight_Gsf", &HLT_Ele40_WPTight_Gsf_);
  tree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06_);
  tree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_);
  tree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10_);
  tree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40", &HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_);
  tree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_);
  tree->Branch("HLT_Ele50_IsoVVVL_PFHT450", &HLT_Ele50_IsoVVVL_PFHT450_);

  // PFJet trigger branches
  tree->Branch("HLT_PFJet110", &HLT_PFJet110_);
  tree->Branch("HLT_PFJet140", &HLT_PFJet140_);
  tree->Branch("HLT_PFJet200_TimeGt2p5ns", &HLT_PFJet200_TimeGt2p5ns_);
  tree->Branch("HLT_PFJet200_TimeLtNeg2p5ns", &HLT_PFJet200_TimeLtNeg2p5ns_);
  tree->Branch("HLT_PFJet200", &HLT_PFJet200_);
  tree->Branch("HLT_PFJet260", &HLT_PFJet260_);
  tree->Branch("HLT_PFJet320", &HLT_PFJet320_);
  tree->Branch("HLT_PFJet400", &HLT_PFJet400_);
  tree->Branch("HLT_PFJet40_GPUvsCPU", &HLT_PFJet40_GPUvsCPU_);
  tree->Branch("HLT_PFJet40", &HLT_PFJet40_);
  tree->Branch("HLT_PFJet450", &HLT_PFJet450_);
  tree->Branch("HLT_PFJet500", &HLT_PFJet500_);
  tree->Branch("HLT_PFJet550", &HLT_PFJet550_);
  tree->Branch("HLT_PFJet60", &HLT_PFJet60_);
  tree->Branch("HLT_PFJet80", &HLT_PFJet80_);

  // PFJet forward trigger branches
  tree->Branch("HLT_PFJetFwd140", &HLT_PFJetFwd140_);
  tree->Branch("HLT_PFJetFwd200", &HLT_PFJetFwd200_);
  tree->Branch("HLT_PFJetFwd260", &HLT_PFJetFwd260_);
  tree->Branch("HLT_PFJetFwd320", &HLT_PFJetFwd320_);
  tree->Branch("HLT_PFJetFwd400", &HLT_PFJetFwd400_);
  tree->Branch("HLT_PFJetFwd440", &HLT_PFJetFwd440_);
  tree->Branch("HLT_PFJetFwd450", &HLT_PFJetFwd450_);
  tree->Branch("HLT_PFJetFwd500", &HLT_PFJetFwd500_);

  // PFMET trigger branches
  tree->Branch("HLT_PFMET105_IsoTrk50", &HLT_PFMET105_IsoTrk50_);
  tree->Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_);
  tree->Branch("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight_);
  tree->Branch("HLT_PFMET130_PFMHT130_IDTight", &HLT_PFMET130_PFMHT130_IDTight_);
  tree->Branch("HLT_PFMET140_PFMHT140_IDTight", &HLT_PFMET140_PFMHT140_IDTight_);
  tree->Branch("HLT_PFMET200_BeamHaloCleaned", &HLT_PFMET200_BeamHaloCleaned_);
  tree->Branch("HLT_PFMET200_NotCleaned", &HLT_PFMET200_NotCleaned_);
  tree->Branch("HLT_PFMET250_NotCleaned", &HLT_PFMET250_NotCleaned_);
  tree->Branch("HLT_PFMET300_NotCleaned", &HLT_PFMET300_NotCleaned_);
  tree->Branch("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_);
  tree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_);
  tree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_);
  tree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_);
  tree->Branch("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_);
  tree->Branch("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight", &HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_);
  tree->Branch("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_);
  tree->Branch("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight", &HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_);
  tree->Branch("HLT_PFMETTypeOne140_PFMHT140_IDTight", &HLT_PFMETTypeOne140_PFMHT140_IDTight_);
  tree->Branch("HLT_PFMETTypeOne200_BeamHaloCleaned", &HLT_PFMETTypeOne200_BeamHaloCleaned_);
}

void ggNtuplizer::fillGlobalEvent(const edm::Event& e, const edm::EventSetup& es) {
  std::cout<<"***********************globalEvent*****************"<<std::endl;

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
  
// Assign trigger results to individual variables
  for(int jk=0; jk<nHLTmx; jk++) {
    // Photon triggers
    if(jk==0)      { HLT_Photon100EBHE10_ = booltrg[jk];}
    else if(jk==1) { HLT_Photon110EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==2) { HLT_Photon120_R9Id90_HE10_IsoM_ = booltrg[jk];}
    else if(jk==3) { HLT_Photon120_ = booltrg[jk];}
    else if(jk==4) { HLT_Photon130EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==5) { HLT_Photon150EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==6) { HLT_Photon150_ = booltrg[jk];}
    else if(jk==7) { HLT_Photon165_R9Id90_HE10_IsoM_ = booltrg[jk];}
    else if(jk==8) { HLT_Photon175EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==9) { HLT_Photon175_ = booltrg[jk];}
    else if(jk==10) { HLT_Photon200EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==11) { HLT_Photon200_ = booltrg[jk];}
    else if(jk==12) { HLT_Photon20_HoverELoose_ = booltrg[jk];}
    else if(jk==13) { HLT_Photon300_NoHE_ = booltrg[jk];}
    else if(jk==14) { HLT_Photon30EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==15) { HLT_Photon30_HoverELoose_ = booltrg[jk];}
    else if(jk==16) { HLT_Photon32_OneProng32_M50To105_ = booltrg[jk];}
    else if(jk==17) { HLT_Photon33_ = booltrg[jk];}
    else if(jk==18) { HLT_Photon35_TwoProngs35_ = booltrg[jk];}
    else if(jk==19) { HLT_Photon50EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==20) { HLT_Photon50_R9Id90_HE10_IsoM_ = booltrg[jk];}
    else if(jk==21) { HLT_Photon50_TimeGt2p5ns_ = booltrg[jk];}
    else if(jk==22) { HLT_Photon50_TimeLtNeg2p5ns_ = booltrg[jk];}
    else if(jk==23) { HLT_Photon50_ = booltrg[jk];}
    else if(jk==24) { HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT350_ = booltrg[jk];}
    else if(jk==25) { HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT380_ = booltrg[jk];}
    else if(jk==26) { HLT_Photon60_R9Id90_CaloIdL_IsoL_DisplacedIdL_PFHT400_ = booltrg[jk];}
    else if(jk==27) { HLT_Photon60_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_ = booltrg[jk];}
    else if(jk==28) { HLT_Photon75EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==29) { HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_ = booltrg[jk];}
    else if(jk==30) { HLT_Photon75_R9Id90_HE10_IsoM_ = booltrg[jk];}
    else if(jk==31) { HLT_Photon75_ = booltrg[jk];}
    else if(jk==32) { HLT_Photon90EB_TightID_TightIso_ = booltrg[jk];}
    else if(jk==33) { HLT_Photon90_R9Id90_HE10_IsoM_ = booltrg[jk];}
    else if(jk==34) { HLT_Photon90_ = booltrg[jk];}
    
    // Dielectron triggers
    else if(jk==35) { HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_ = booltrg[jk];}
    else if(jk==36) { HLT_DoubleEle33_CaloIdL_MW_ = booltrg[jk];}
    
    // Single electron triggers
    else if(jk==37) { HLT_Ele115_CaloIdVT_GsfTrkIdT_ = booltrg[jk]; }
    else if(jk==38) { HLT_Ele135_CaloIdVT_GsfTrkIdT_ = booltrg[jk]; }
    else if(jk==39) { HLT_Ele30_WPTight_Gsf_ = booltrg[jk]; }
    else if(jk==40) { HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_ = booltrg[jk]; }
    else if(jk==41) { HLT_Ele32_WPTight_Gsf_L1DoubleEG_ = booltrg[jk]; }
    else if(jk==42) { HLT_Ele32_WPTight_Gsf_ = booltrg[jk]; }
    else if(jk==43) { HLT_Ele35_WPTight_Gsf_ = booltrg[jk]; }
    else if(jk==44) { HLT_Ele38_WPTight_Gsf_ = booltrg[jk]; }
    else if(jk==45) { HLT_Ele40_WPTight_Gsf_ = booltrg[jk]; }
    else if(jk==46) { HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_PNetBB0p06_ = booltrg[jk]; }
    else if(jk==47) { HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet220_SoftDropMass40_ = booltrg[jk]; }
    else if(jk==48) { HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_PNetBB0p10_ = booltrg[jk]; }
    else if(jk==49) { HLT_Ele50_CaloIdVT_GsfTrkIdT_AK8PFJet230_SoftDropMass40_ = booltrg[jk]; }
    else if(jk==50) { HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_ = booltrg[jk]; }
    else if(jk==51) { HLT_Ele50_IsoVVVL_PFHT450_ = booltrg[jk]; }
    
    // PFJet triggers
    else if(jk==52) { HLT_PFJet110_ = booltrg[jk]; }
    else if(jk==53) { HLT_PFJet140_ = booltrg[jk]; }
    else if(jk==54) { HLT_PFJet200_TimeGt2p5ns_ = booltrg[jk]; }
    else if(jk==55) { HLT_PFJet200_TimeLtNeg2p5ns_ = booltrg[jk]; }
    else if(jk==56) { HLT_PFJet200_ = booltrg[jk]; }
    else if(jk==57) { HLT_PFJet260_ = booltrg[jk]; }
    else if(jk==58) { HLT_PFJet320_ = booltrg[jk]; }
    else if(jk==59) { HLT_PFJet400_ = booltrg[jk]; }
    else if(jk==60) { HLT_PFJet40_GPUvsCPU_ = booltrg[jk]; }
    else if(jk==61) { HLT_PFJet40_ = booltrg[jk]; }
    else if(jk==62) { HLT_PFJet450_ = booltrg[jk]; }
    else if(jk==63) { HLT_PFJet500_ = booltrg[jk]; }
    else if(jk==64) { HLT_PFJet550_ = booltrg[jk]; }
    else if(jk==65) { HLT_PFJet60_ = booltrg[jk]; }
    else if(jk==66) { HLT_PFJet80_ = booltrg[jk]; }
    
    // PFJet forward triggers
    else if(jk==67) { HLT_PFJetFwd140_ = booltrg[jk]; }
    else if(jk==68) { HLT_PFJetFwd200_ = booltrg[jk]; }
    else if(jk==69) { HLT_PFJetFwd260_ = booltrg[jk]; }
    else if(jk==70) { HLT_PFJetFwd320_ = booltrg[jk]; }
    else if(jk==71) { HLT_PFJetFwd400_ = booltrg[jk]; }
    else if(jk==72) { HLT_PFJetFwd440_ = booltrg[jk]; }
    else if(jk==73) { HLT_PFJetFwd450_ = booltrg[jk]; }
    else if(jk==74) { HLT_PFJetFwd500_ = booltrg[jk]; }
    
    // PFMET triggers
    else if(jk==75) { HLT_PFMET105_IsoTrk50_ = booltrg[jk]; }
    else if(jk==76) { HLT_PFMET120_PFMHT120_IDTight_PFHT60_ = booltrg[jk]; }
    else if(jk==77) { HLT_PFMET120_PFMHT120_IDTight_ = booltrg[jk]; }
    else if(jk==78) { HLT_PFMET130_PFMHT130_IDTight_ = booltrg[jk]; }
    else if(jk==79) { HLT_PFMET140_PFMHT140_IDTight_ = booltrg[jk]; }
    else if(jk==80) { HLT_PFMET200_BeamHaloCleaned_ = booltrg[jk]; }
    else if(jk==81) { HLT_PFMET200_NotCleaned_ = booltrg[jk]; }
    else if(jk==82) { HLT_PFMET250_NotCleaned_ = booltrg[jk]; }
    else if(jk==83) { HLT_PFMET300_NotCleaned_ = booltrg[jk]; }
    else if(jk==84) { HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_FilterHF_ = booltrg[jk]; }
    else if(jk==85) { HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_FilterHF_ = booltrg[jk]; }
    else if(jk==86) { HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_ = booltrg[jk]; }
    else if(jk==87) { HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_ = booltrg[jk]; }
    else if(jk==88) { HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_FilterHF_ = booltrg[jk]; }
    else if(jk==89) { HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_ = booltrg[jk]; }
    else if(jk==90) { HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_FilterHF_ = booltrg[jk]; }
    else if(jk==91) { HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_ = booltrg[jk]; }
    else if(jk==92) { HLT_PFMETTypeOne140_PFMHT140_IDTight_ = booltrg[jk]; }
    else if(jk==99) { HLT_PFMETTypeOne200_BeamHaloCleaned_ = booltrg[jk]; }

    
  }
}
