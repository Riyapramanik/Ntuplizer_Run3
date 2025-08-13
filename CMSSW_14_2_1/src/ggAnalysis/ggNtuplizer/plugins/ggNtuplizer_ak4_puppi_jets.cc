#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TROOT.h"
#include "TRandom.h"
#include "ggAnalysis/ggNtuplizer/interface/JetIDTools.h"

using namespace std;

UShort_t          nAK4PUPPIJet_;
vector<Char_t>    AK4PUPPIJet_Charge_;
vector<float>     AK4PUPPIJet_Pt_;
vector<float>     AK4PUPPIJet_En_;
vector<float>     AK4PUPPIJet_Eta_;
vector<float>     AK4PUPPIJet_Phi_;
vector<float>     AK4PUPPIJet_RawPt_;
vector<float>     AK4PUPPIJet_RawEn_;
vector<float>     AK4PUPPIJet_Mt_;
vector<float>     AK4PUPPIJet_Area_;
vector<float>     AK4PUPPIJet_CHF_;
vector<float>     AK4PUPPIJet_NHF_;
vector<float>     AK4PUPPIJet_CEF_;
vector<float>     AK4PUPPIJet_NEF_;
vector<UShort_t>  AK4PUPPIJet_NNP_;
vector<float>     AK4PUPPIJet_MUF_;

vector<float>     AK4PUPPIJet_JECUnc_;
vector<Float_t>   AK4PUPPIJet_JER_SF_;
vector<Float_t>   AK4PUPPIJet_JER_SFup_;
vector<Float_t>   AK4PUPPIJet_JER_SFdown_;
vector<ULong64_t> AK4PUPPIJet_FiredTrgs_;
//By me
std::vector<float> PFJetAK4_JEC_;
double min_pt_AK4jet;
double max_eta;

std::vector<float> jesup_AbsoluteStat_;
std::vector<float> jesup_AbsoluteScale_;
std::vector<float> jesup_AbsoluteMPFBias_;
std::vector<float> jesup_FlavorQCD_;
std::vector<float> jesup_Fragmentation_;
std::vector<float> jesup_PileUpDataMC_;
std::vector<float> jesup_PileUpPtBB_;
std::vector<float> jesup_PileUpPtEC1_;
std::vector<float> jesup_PileUpPtEC2_;
std::vector<float> jesup_PileUpPtRef_;
std::vector<float> jesup_RelativeFSR_;
std::vector<float> jesup_RelativeJEREC1_;
std::vector<float> jesup_RelativeJEREC2_;
std::vector<float> jesup_RelativePtBB_;
std::vector<float> jesup_RelativePtEC1_;
std::vector<float> jesup_RelativePtEC2_;
std::vector<float> jesup_RelativeBal_;
std::vector<float> jesup_RelativeSample_;
std::vector<float> jesup_RelativeStatEC_;
std::vector<float> jesup_RelativeStatFSR_;
std::vector<float> jesup_SinglePionECAL_;
std::vector<float> jesup_SinglePionHCAL_;
std::vector<float> jesup_TimePtEta_;
std::vector<float> jesup_Total_;

std::vector<float> jesdn_AbsoluteStat_;
std::vector<float> jesdn_AbsoluteScale_;
std::vector<float> jesdn_AbsoluteMPFBias_;
std::vector<float> jesdn_FlavorQCD_;
std::vector<float> jesdn_Fragmentation_;
std::vector<float> jesdn_PileUpDataMC_;
std::vector<float> jesdn_PileUpPtBB_;
std::vector<float> jesdn_PileUpPtEC1_;
std::vector<float> jesdn_PileUpPtEC2_;
std::vector<float> jesdn_PileUpPtRef_;
std::vector<float> jesdn_RelativeFSR_;
std::vector<float> jesdn_RelativeJEREC1_;
std::vector<float> jesdn_RelativeJEREC2_;
std::vector<float> jesdn_RelativePtBB_;
std::vector<float> jesdn_RelativePtEC1_;
std::vector<float> jesdn_RelativePtEC2_;
std::vector<float> jesdn_RelativeBal_;
std::vector<float> jesdn_RelativeSample_;
std::vector<float> jesdn_RelativeStatEC_;
std::vector<float> jesdn_RelativeStatFSR_;
std::vector<float> jesdn_SinglePionECAL_;
std::vector<float> jesdn_SinglePionHCAL_;
std::vector<float> jesdn_TimePtEta_;
std::vector<float> jesdn_Total_;
std::vector<float> jetID_;

std::vector<int>    jetveto_Flag_;
std::vector<int>    jetveto_eep_Flag_;
std::vector<int>    hadronflav_;
std::vector<int>    partonflav_;
std::vector<float>  qgl_;
std::vector<float>  PUID_;
std::vector<float>  charge_kappa_0p3_;
std::vector<float>  charge_kappa_0p6_;
std::vector<float>  charge_kappa_1p0_;
std::vector<float>  charged_ptsum_;
std::vector<bool> jetID_tightlepveto_;

TFile *file_jetvetomap;
TH2D *h_jetvetomap;
TH2D *h_jetvetomap_eep;

bool store_ak4jets, subtractLepton_fromAK4;
void Read_JEC(double &total_JEC,  double &tmprecpt, 
			  TLorentzVector jetp4, 
			  double Rho, 
			  bool isData,
			  pat::Jet jet,
			  FactorizedJetCorrector *jecL1Fast, FactorizedJetCorrector *jecL2Relative, FactorizedJetCorrector *jecL3Absolute, FactorizedJetCorrector*jecL2L3Residual)
{
	
    double jeteta = jetp4.Eta();
    double jetphi = jetp4.Phi();
	
    double total_cor =1;
      
    jecL1Fast->setJetPt(tmprecpt); jecL1Fast->setJetA(jet.jetArea()); jecL1Fast->setRho(Rho);jecL1Fast->setJetEta(jeteta); 
    double corFactorL1Fast = jecL1Fast->getCorrection();
    total_cor *= corFactorL1Fast;
    tmprecpt = tmprecpt * corFactorL1Fast;
      
    jecL2Relative->setJetPt(tmprecpt); jecL2Relative->setJetEta(jeteta); jecL2Relative->setJetPhi(jetphi);
    double corFactorL2Relative = jecL2Relative->getCorrection();
    total_cor *= corFactorL2Relative ;
    tmprecpt = tmprecpt * corFactorL2Relative;
      
    jecL3Absolute->setJetPt(tmprecpt); jecL3Absolute->setJetEta(jeteta);
    double corFactorL3Absolute = jecL3Absolute->getCorrection();
    total_cor *= corFactorL3Absolute ;
    tmprecpt = tmprecpt * corFactorL3Absolute;
      
    double corFactorL2L3Residual=1.;
      
    if(isData){
		jecL2L3Residual->setJetPt(tmprecpt); jecL2L3Residual->setJetEta(jeteta);
		corFactorL2L3Residual = jecL2L3Residual->getCorrection();
		total_cor*= corFactorL2L3Residual;
		tmprecpt *=corFactorL2L3Residual;
	}
	
	total_JEC = total_cor;
	
	return;     
}

void Read_JER(std::string mPtResoFile, std::string mPtSFFile, double tmprecpt, TLorentzVector pfjet4v, double Rho, edm::Handle<reco::GenJetCollection>  genjets, double dRcut, vector<double> &SFs)
{
 
	JME::JetResolution resolution;
	resolution = JME::JetResolution(mPtResoFile.c_str());
	JME::JetResolutionScaleFactor res_sf;
	res_sf = JME::JetResolutionScaleFactor(mPtSFFile.c_str());
	
	JME::JetParameters parameters_5 = {{JME::Binning::JetPt, tmprecpt}, {JME::Binning::JetEta, pfjet4v.Eta()}, {JME::Binning::Rho, Rho}};
	double rp = resolution.getResolution(parameters_5);
	double gaus_rp = gRandom->Gaus(0.,rp);
	double sf = res_sf.getScaleFactor(parameters_5, Variation::NOMINAL);
	double sf_up = res_sf.getScaleFactor(parameters_5, Variation::UP);
	double sf_dn = res_sf.getScaleFactor(parameters_5, Variation::DOWN);
	
	bool match = false;
	int match_gen = -1;
		
	for (unsigned get = 0; get<(genjets->size()); get++) {
		TLorentzVector genjet4v((*genjets)[get].px(),(*genjets)[get].py(),(*genjets)[get].pz(), (*genjets)[get].energy());
		if((deltaR2(pfjet4v.Rapidity(),pfjet4v.Phi(),genjet4v.Rapidity(),genjet4v.Phi()) < (dRcut)) &&(fabs(tmprecpt-genjet4v.Pt())<(3*fabs(rp)*tmprecpt))){
			match = true;
			match_gen = get;
			break;
		}
	}
		
	if(match && (match_gen>=0)){
	  
		SFs.push_back((sf-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_up-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_dn-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
	  
	}else{
	  
		SFs.push_back(sqrt(max(0.,(sf*sf-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_up*sf_up-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_dn*sf_dn-1))) * gaus_rp);
	}
      	
}


bool Assign_JetVeto(TLorentzVector p4, bool jetID, JetIDVars IDVars,  edm::Handle<edm::View<pat::Muon>>  Muon_collection, TH2D *h_map, float pt_cut=15, float EMF_cut=0.9, float dRmu_cut=0.2)
{
	
	//See: https://cms-talk.web.cern.ch/t/jet-veto-maps-for-run3-data/18444
	
	bool is_veto = false;
	
	TLorentzVector jet_p4;
	jet_p4.SetPtEtaPhiM(p4.Pt(),p4.Eta(), p4.Phi(), p4.M());
		
	if(jet_p4.Pt()>pt_cut && jetID && ((IDVars.NEMF+IDVars.CEMF)<EMF_cut)){
				  
			for(unsigned int iMuon = 0; iMuon < Muon_collection->size(); iMuon++ ) {  
			
				const auto &muon_cand = (*Muon_collection)[iMuon];
				TLorentzVector muon_cand_p4;
				muon_cand_p4.SetPtEtaPhiM(muon_cand.pt(),muon_cand.eta(),muon_cand.phi(),muon_cand.mass());
			
				if(muon_cand.isPFMuon() && jet_p4.DeltaR(muon_cand_p4)>dRmu_cut){
				
					int eta_bin_index = h_map -> GetXaxis() -> FindBin(jet_p4.Eta());
					int phi_bin_index = h_map -> GetYaxis() -> FindBin(jet_p4.Phi());
					if(h_map ->  GetBinContent(eta_bin_index, phi_bin_index) > 0 )  { is_veto = true; }
				
			}  
		}
	}
		
	return is_veto;
}

void ggNtuplizer::branchesAK4PUPPIJets(TTree* tree) {
	tree->Branch("nAK4PUPPIJet",                          & nAK4PUPPIJet_);
	tree->Branch("AK4PUPPIJet_Charge",                    & AK4PUPPIJet_Charge_);
	tree->Branch("AK4PUPPIJet_Pt",                        & AK4PUPPIJet_Pt_);
	tree->Branch("AK4PUPPIJet_En",                        & AK4PUPPIJet_En_);
	tree->Branch("AK4PUPPIJet_Eta",                       & AK4PUPPIJet_Eta_);
	tree->Branch("AK4PUPPIJet_Phi",                       & AK4PUPPIJet_Phi_);
	tree->Branch("AK4PUPPIJet_RawPt",                     & AK4PUPPIJet_RawPt_);
	tree->Branch("AK4PUPPIJet_RawEn",                     & AK4PUPPIJet_RawEn_);
	tree->Branch("AK4PUPPIJet_Mt",                        & AK4PUPPIJet_Mt_);
	tree->Branch("AK4PUPPIJet_Area",                      & AK4PUPPIJet_Area_);
       
	tree->Branch("AK4PUPPIJet_JECUnc",   & AK4PUPPIJet_JECUnc_);
	tree->Branch("AK4PUPPIJet_JER_SF",     & AK4PUPPIJet_JER_SF_);
	tree->Branch("AK4PUPPIJet_JER_SFup",   & AK4PUPPIJet_JER_SFup_);
	tree->Branch("AK4PUPPIJet_JER_SFdown", & AK4PUPPIJet_JER_SFdown_);
	tree->Branch("AK4PUPPIJet_FiredTrgs",  & AK4PUPPIJet_FiredTrgs_);
	tree->Branch("AK4PUPPIJet_CHF",        & AK4PUPPIJet_CHF_);
	tree->Branch("AK4PUPPIJet_NHF",        & AK4PUPPIJet_NHF_);
	tree->Branch("AK4PUPPIJet_CEF",        & AK4PUPPIJet_CEF_);
	tree->Branch("AK4PUPPIJet_NEF",        & AK4PUPPIJet_NEF_);
	tree->Branch("AK4PUPPIJet_NNP",        & AK4PUPPIJet_NNP_);
	tree->Branch("AK4PUPPIJet_MUF",        & AK4PUPPIJet_MUF_);

	tree_->Branch("AK4PUPPIJet_jetID",       &jetID_);
	tree_->Branch("AK4PUPPIJet_jesup_AbsoluteStat",       &jesup_AbsoluteStat_);
	tree_->Branch("AK4PUPPIJet_jesup_AbsoluteScale",      &jesup_AbsoluteScale_);
	tree_->Branch("AK4PUPPIJet_jesup_AbsoluteMPFBias",    &jesup_AbsoluteMPFBias_);
	tree_->Branch("AK4PUPPIJet_jesup_FlavorQCD",          &jesup_FlavorQCD_);
	tree_->Branch("AK4PUPPIJet_jesup_Fragmentation",      &jesup_Fragmentation_);
	tree_->Branch("AK4PUPPIJet_jesup_PileUpDataMC",       &jesup_PileUpDataMC_);
	tree_->Branch("AK4PUPPIJet_jesup_PileUpPtBB",         &jesup_PileUpPtBB_);
	tree_->Branch("AK4PUPPIJet_jesup_PileUpPtEC1",        &jesup_PileUpPtEC1_);
	tree_->Branch("AK4PUPPIJet_jesup_PileUpPtEC2",        &jesup_PileUpPtEC2_);
	tree_->Branch("AK4PUPPIJet_jesup_PileUpPtRef",        &jesup_PileUpPtRef_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeFSR",        &jesup_RelativeFSR_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeJEREC1",     &jesup_RelativeJEREC1_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeJEREC2",     &jesup_RelativeJEREC2_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativePtBB",       &jesup_RelativePtBB_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativePtEC1",      &jesup_RelativePtEC1_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativePtEC2",      &jesup_RelativePtEC2_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeBal",        &jesup_RelativeBal_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeSample",     &jesup_RelativeSample_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeStatEC",     &jesup_RelativeStatEC_);
	tree_->Branch("AK4PUPPIJet_jesup_RelativeStatFSR",    &jesup_RelativeStatFSR_);
	tree_->Branch("AK4PUPPIJet_jesup_SinglePionECAL",     &jesup_SinglePionECAL_);
	tree_->Branch("AK4PUPPIJet_jesup_SinglePionHCAL",     &jesup_SinglePionHCAL_);
	tree_->Branch("AK4PUPPIJet_jesup_TimePtEta",          &jesup_TimePtEta_);
	tree_->Branch("AK4PUPPIJet_jesup_Total",              &jesup_Total_);

	tree->Branch("AK4PUPPIJet_jesdn_AbsoluteStat",       &jesdn_AbsoluteStat_);
	tree->Branch("AK4PUPPIJet_jesdn_AbsoluteScale",      &jesdn_AbsoluteScale_);
	tree->Branch("AK4PUPPIJet_jesdn_AbsoluteMPFBias",    &jesdn_AbsoluteMPFBias_);
	tree->Branch("AK4PUPPIJet_jesdn_FlavorQCD",          &jesdn_FlavorQCD_);
	tree->Branch("AK4PUPPIJet_jesdn_Fragmentation",      &jesdn_Fragmentation_);
	tree->Branch("AK4PUPPIJet_jesdn_PileUpDataMC",       &jesdn_PileUpDataMC_);
	tree->Branch("AK4PUPPIJet_jesdn_PileUpPtBB",         &jesdn_PileUpPtBB_);
	tree->Branch("AK4PUPPIJet_jesdn_PileUpPtEC1",        &jesdn_PileUpPtEC1_);
	tree->Branch("AK4PUPPIJet_jesdn_PileUpPtEC2",        &jesdn_PileUpPtEC2_);
	tree->Branch("AK4PUPPIJet_jesdn_PileUpPtRef",        &jesdn_PileUpPtRef_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeFSR",        &jesdn_RelativeFSR_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeJEREC1",     &jesdn_RelativeJEREC1_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeJEREC2",     &jesdn_RelativeJEREC2_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativePtBB",       &jesdn_RelativePtBB_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativePtEC1",      &jesdn_RelativePtEC1_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativePtEC2",      &jesdn_RelativePtEC2_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeBal",        &jesdn_RelativeBal_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeSample",     &jesdn_RelativeSample_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeStatEC",     &jesdn_RelativeStatEC_);
	tree->Branch("AK4PUPPIJet_jesdn_RelativeStatFSR",    &jesdn_RelativeStatFSR_);
	tree->Branch("AK4PUPPIJet_jesdn_SinglePionECAL",     &jesdn_SinglePionECAL_);
	tree->Branch("AK4PUPPIJet_jesdn_SinglePionHCAL",     &jesdn_SinglePionHCAL_);
	tree->Branch("AK4PUPPIJet_jesdn_TimePtEta",          &jesdn_TimePtEta_);
	tree->Branch("AK4PUPPIJet_jesdn_Total",              &jesdn_Total_);

	tree_->Branch("AK4PUPPIJet_jetveto_Flag",         &jetveto_Flag_);
	tree_->Branch("AK4PUPPIJet_jetveto_eep_Flag",     &jetveto_eep_Flag_);
	tree_->Branch("AK4PUPPIJet_hadronflav",           &hadronflav_);
	tree_->Branch("AK4PUPPIJet_partonflav",           &partonflav_);
	tree_->Branch("AK4PUPPIJet_qgl",                  &qgl_);
	tree_->Branch("AK4PUPPIJet_PUID",                 &PUID_);
	tree_->Branch("AK4PUPPIJet_charge_kappa_0p3",     &charge_kappa_0p3_);
	tree_->Branch("AK4PUPPIJet_charge_kappa_0p6",     &charge_kappa_0p6_);
	tree_->Branch("AK4PUPPIJet_charge_kappa_1p0",     &charge_kappa_1p0_);
	tree_->Branch("AK4PUPPIJet_charged_ptsum",        &charged_ptsum_);
	tree_->Branch("AK4PUPPIJet_etID_tightlepveto", &jetID_tightlepveto_);

}


void ggNtuplizer::fillAK4PUPPIJets(const edm::Event& e, const edm::EventSetup& es) {
	AK4PUPPIJet_Charge_                       . clear();
	AK4PUPPIJet_Pt_                           . clear();
	AK4PUPPIJet_En_                           . clear();
	AK4PUPPIJet_Eta_                          . clear();
	AK4PUPPIJet_Phi_                          . clear();
	AK4PUPPIJet_RawPt_                        . clear();
	AK4PUPPIJet_RawEn_                        . clear();
	AK4PUPPIJet_Mt_                           . clear();
	AK4PUPPIJet_Area_                         . clear();	
	AK4PUPPIJet_JECUnc_                       . clear();
	AK4PUPPIJet_JER_SF_                       . clear();
	AK4PUPPIJet_JER_SFup_                     . clear();
	AK4PUPPIJet_JER_SFdown_                   . clear();
	AK4PUPPIJet_FiredTrgs_                    . clear();
	AK4PUPPIJet_CHF_                          . clear();
	AK4PUPPIJet_NHF_                          . clear();
	AK4PUPPIJet_CEF_                          . clear();
	AK4PUPPIJet_NEF_                          . clear();
	AK4PUPPIJet_NNP_                          . clear();
	AK4PUPPIJet_MUF_                          . clear();
 	nAK4PUPPIJet_ = 0;	
	PFJetAK4_JEC_                             . clear();

	jesup_AbsoluteStat_        .clear();
	jesup_AbsoluteScale_       .clear();
	jesup_AbsoluteMPFBias_     .clear();
	jesup_FlavorQCD_           .clear();
	jesup_Fragmentation_       .clear();
	jesup_PileUpDataMC_        .clear();
	jesup_PileUpPtBB_          .clear();
	jesup_PileUpPtEC1_         .clear();
	jesup_PileUpPtEC2_         .clear();
	jesup_PileUpPtRef_         .clear();
	jesup_RelativeFSR_         .clear();
	jesup_RelativeJEREC1_      .clear();
	jesup_RelativeJEREC2_      .clear();
	jesup_RelativePtBB_        .clear();
	jesup_RelativePtEC1_       .clear();
	jesup_RelativePtEC2_       .clear();
	jesup_RelativeBal_         .clear();
	jesup_RelativeSample_      .clear();
	jesup_RelativeStatEC_      .clear();
	jesup_RelativeStatFSR_     .clear();
	jesup_SinglePionECAL_      .clear();
	jesup_SinglePionHCAL_      .clear();
	jesup_TimePtEta_           .clear();
	jesup_Total_               .clear();

	jesdn_AbsoluteStat_       .clear();
	jesdn_AbsoluteScale_      .clear();
	jesdn_AbsoluteMPFBias_    .clear();
	jesdn_FlavorQCD_          .clear();
	jesdn_Fragmentation_      .clear();
	jesdn_PileUpDataMC_       .clear();
	jesdn_PileUpPtBB_         .clear();
	jesdn_PileUpPtEC1_        .clear();
	jesdn_PileUpPtEC2_        .clear();
	jesdn_PileUpPtRef_        .clear();
	jesdn_RelativeFSR_        .clear();
	jesdn_RelativeJEREC1_     .clear();
	jesdn_RelativeJEREC2_     .clear();
	jesdn_RelativePtBB_       .clear();
	jesdn_RelativePtEC1_      .clear();
	jesdn_RelativePtEC2_      .clear();
	jesdn_RelativeBal_        .clear();
	jesdn_RelativeSample_     .clear();
	jesdn_RelativeStatEC_     .clear();
	jesdn_RelativeStatFSR_    .clear();
	jesdn_SinglePionECAL_     .clear();
	jesdn_SinglePionHCAL_     .clear();
	jesdn_TimePtEta_          .clear();
	jesdn_Total_              .clear();
	jetID_                    .clear();

	jetveto_Flag_.clear();
	jetveto_eep_Flag_.clear();
	hadronflav_.clear();
	partonflav_.clear();
	qgl_.clear();
	PUID_.clear();
	charge_kappa_0p3_.clear();
	charge_kappa_0p6_.clear();
	charge_kappa_1p0_.clear();
	charged_ptsum_.clear();
	jetID_tightlepveto_.clear();

	std::cout<<"***************************Started Jets************************"<<std::endl;
	
	edm::Handle<edm::View<pat::Jet>> pfjetAK4s;
	e.getByToken(tok_pfjetAK4s_, pfjetAK4s);
	if(pfjetAK4s.isValid() && store_ak4jets){
	  for (edm::View<pat::Jet>::const_iterator iJet = pfjetAK4s->begin(); iJet != pfjetAK4s->end(); ++iJet) {
	    const pat::Jet &ak4jet = *iJet;

	    //Rho                                                                                                                                            
        edm::Handle<double> Rho_PF;
        e.getByToken(tok_Rho_,Rho_PF);

        edm::Handle<edm::View<pat::Muon>> muons;
        e.getByToken(tok_muons_, muons);

        double Rho = *Rho_PF;
	edm::Handle<reco::GenJetCollection> genjetAK4s;
	e.getByToken(tok_genjetAK4s_, genjetAK4s);
	    
	TLorentzVector pfjetAK4_4v(
				   ak4jet.correctedP4("Uncorrected").px(),
				   ak4jet.correctedP4("Uncorrected").py(),
				   ak4jet.correctedP4("Uncorrected").pz(),
				   ak4jet.correctedP4("Uncorrected").energy()
				   );
        double tmprecpt = pfjetAK4_4v.Pt();
        double total_cor = 1.0;
	//JEC orrection on Jet
	Read_JEC(total_cor,tmprecpt,pfjetAK4_4v,Rho,isData,ak4jet,jecL1FastAK4,jecL2RelativeAK4,jecL3AbsoluteAK4,jecL2L3ResidualAK4);
	 // basic selection conditions on JEC-corrected jets //                                                                                         
	if(tmprecpt<min_pt_AK4jet) continue;
	if(abs(ak4jet.eta())>max_eta) continue;
	
	//JEC corrected 4-momentum //                                                                                                                   
	TLorentzVector pfjetAK4_4v_jecor;
	pfjetAK4_4v_jecor = pfjetAK4_4v * total_cor;
	PFJetAK4_JEC_.push_back(total_cor);
	TLorentzVector rawP4 = pfjetAK4_4v;
	TLorentzVector corrP4 = pfjetAK4_4v_jecor;
	
	//JEC corrected Jet variable
	AK4PUPPIJet_Pt_.push_back(corrP4.Pt());
	AK4PUPPIJet_En_.push_back(corrP4.Energy());
	AK4PUPPIJet_Eta_.push_back(corrP4.Eta());
	AK4PUPPIJet_Phi_.push_back(corrP4.Phi());
	
	AK4PUPPIJet_Charge_.push_back(ak4jet.charge());
	AK4PUPPIJet_RawPt_.push_back(rawP4.Pt());
	AK4PUPPIJet_RawEn_.push_back(rawP4.Energy());
	AK4PUPPIJet_Mt_.push_back(ak4jet.mt());
	AK4PUPPIJet_Area_.push_back(ak4jet.jetArea());
	AK4PUPPIJet_CEF_.push_back(ak4jet.chargedEmEnergyFraction());
	AK4PUPPIJet_NEF_.push_back(ak4jet.neutralEmEnergyFraction());
	AK4PUPPIJet_CHF_.push_back(ak4jet.chargedHadronEnergyFraction());
	AK4PUPPIJet_NHF_.push_back(ak4jet.neutralHadronEnergyFraction());
	AK4PUPPIJet_NNP_.push_back(ak4jet.neutralMultiplicity());
	AK4PUPPIJet_MUF_.push_back(ak4jet.muonEnergyFraction());

	//JER for MC                                                                                                                                     
        if (isMC) {
          TLorentzVector tmp4v = pfjetAK4_4v_jecor;
          std::vector<double> SFs;
          Read_JER(mPtResoFileAK4, mPtSFFileAK4, tmprecpt, tmp4v, Rho, genjetAK4s, 0.5 * 0.4, SFs);

          AK4PUPPIJet_JER_SF_.push_back(SFs[0]);
          AK4PUPPIJet_JER_SFup_.push_back(SFs[1]);
          AK4PUPPIJet_JER_SFdown_.push_back(SFs[2]);
        }

	//JES uncertainty //
	static const int njetmx = 10;
	static const int nsrc = 24;
	const int njecmcmx = 2*nsrc + 1 ;

	for (int isrc = 0; isrc < njecmcmx; ++isrc) {
	  double sup = 1.0;
	  JetCorrectionUncertainty* jecUnc = nullptr;
	  if (isrc > 0 && isrc <= nsrc) {
	    jecUnc = vsrc[isrc - 1];
	    jecUnc->setJetEta(pfjetAK4_4v_jecor.Eta());
	    jecUnc->setJetPt(pfjetAK4_4v_jecor.Pt());
	    sup += jecUnc->getUncertainty(true);
	    
	    switch (isrc) {
            case 1:  jesup_AbsoluteStat_.push_back(sup); break;
            case 2:  jesup_AbsoluteScale_.push_back(sup); break;
            case 3:  jesup_AbsoluteMPFBias_.push_back(sup); break;
            case 4:  jesup_FlavorQCD_.push_back(sup); break;
            case 5:  jesup_Fragmentation_.push_back(sup); break;
            case 6:  jesup_PileUpDataMC_.push_back(sup); break;
            case 7:  jesup_PileUpPtBB_.push_back(sup); break;
            case 8:  jesup_PileUpPtEC1_.push_back(sup); break;
            case 9:  jesup_PileUpPtEC2_.push_back(sup); break;
            case 10: jesup_PileUpPtRef_.push_back(sup); break;
            case 11: jesup_RelativeFSR_.push_back(sup); break;
            case 12: jesup_RelativeJEREC1_.push_back(sup); break;
            case 13: jesup_RelativeJEREC2_.push_back(sup); break;
            case 14: jesup_RelativePtBB_.push_back(sup); break;
            case 15: jesup_RelativePtEC1_.push_back(sup); break;
            case 16: jesup_RelativePtEC2_.push_back(sup); break;
            case 17: jesup_RelativeBal_.push_back(sup); break;
            case 18: jesup_RelativeSample_.push_back(sup); break;
            case 19: jesup_RelativeStatEC_.push_back(sup); break;
            case 20: jesup_RelativeStatFSR_.push_back(sup); break;
            case 21: jesup_SinglePionECAL_.push_back(sup); break;
            case 22: jesup_SinglePionHCAL_.push_back(sup); break;
            case 23: jesup_TimePtEta_.push_back(sup); break;
            case 24: jesup_Total_.push_back(sup); break;
	    }
	  }
	  else if (isrc > nsrc) {
	    jecUnc = vsrc[isrc - 1 - nsrc];
	    jecUnc->setJetEta(pfjetAK4_4v_jecor.Eta());
	    jecUnc->setJetPt(pfjetAK4_4v_jecor.Pt());
	    sup -= jecUnc->getUncertainty(false);
	    
	    switch (isrc - nsrc) {
            case 1:  jesdn_AbsoluteStat_.push_back(sup); break;
            case 2:  jesdn_AbsoluteScale_.push_back(sup); break;
            case 3:  jesdn_AbsoluteMPFBias_.push_back(sup); break;
            case 4:  jesdn_FlavorQCD_.push_back(sup); break;
            case 5:  jesdn_Fragmentation_.push_back(sup); break;
            case 6:  jesdn_PileUpDataMC_.push_back(sup); break;
            case 7:  jesdn_PileUpPtBB_.push_back(sup); break;
            case 8:  jesdn_PileUpPtEC1_.push_back(sup); break;
            case 9:  jesdn_PileUpPtEC2_.push_back(sup); break;
            case 10: jesdn_PileUpPtRef_.push_back(sup); break;
            case 11: jesdn_RelativeFSR_.push_back(sup); break;
            case 12: jesdn_RelativeJEREC1_.push_back(sup); break;
            case 13: jesdn_RelativeJEREC2_.push_back(sup); break;
            case 14: jesdn_RelativePtBB_.push_back(sup); break;
            case 15: jesdn_RelativePtEC1_.push_back(sup); break;
            case 16: jesdn_RelativePtEC2_.push_back(sup); break;
            case 17: jesdn_RelativeBal_.push_back(sup); break;
            case 18: jesdn_RelativeSample_.push_back(sup); break;
            case 19: jesdn_RelativeStatEC_.push_back(sup); break;
            case 20: jesdn_RelativeStatFSR_.push_back(sup); break;
            case 21: jesdn_SinglePionECAL_.push_back(sup); break;
            case 22: jesdn_SinglePionHCAL_.push_back(sup); break;
            case 23: jesdn_TimePtEta_.push_back(sup); break;
            case 24: jesdn_Total_.push_back(sup); break;
	    }
	  }
	}
              
    // Jet ID variables
    JetIDVars AK4idvars;

    AK4idvars.NHF  = (ak4jet.isPFJet()) ? ak4jet.neutralHadronEnergyFraction() : 1;
    AK4idvars.NEMF = (ak4jet.isPFJet()) ? ak4jet.neutralEmEnergyFraction()     : 1;
    AK4idvars.MUF  = (ak4jet.isPFJet()) ? ak4jet.muonEnergyFraction()          : 1;
    AK4idvars.CHF  = (ak4jet.isPFJet()) ? ak4jet.chargedHadronEnergyFraction() : 1;
    AK4idvars.CEMF = (ak4jet.isPFJet()) ? ak4jet.chargedEmEnergyFraction()     : 1;
    AK4idvars.NumConst = (ak4jet.isPFJet()) ? (ak4jet.chargedMultiplicity() + ak4jet.neutralMultiplicity()) : 1;
    AK4idvars.NumNeutralParticle = (ak4jet.isPFJet()) ? ak4jet.neutralMultiplicity() : 1;
    AK4idvars.CHM  = (ak4jet.isPFJet()) ? ak4jet.chargedHadronMultiplicity() : 1;
    jetID_.push_back(getJetID(AK4idvars, "CHS", year, pfjetAK4_4v.Eta(), false, isUltraLegacy, isRun3));
    jetID_tightlepveto_.push_back(getJetID(AK4idvars, "CHS", year, pfjetAK4_4v.Eta(), true, isUltraLegacy, isRun3));
    
    // Jet veto flag
    bool vetoFlag = Assign_JetVeto(pfjetAK4_4v, jetID_.back(), AK4idvars, muons, h_jetvetomap);
    bool vetoEEP  = (year == "2022EE") ? Assign_JetVeto(pfjetAK4_4v, jetID_.back(), AK4idvars, muons, h_jetvetomap_eep, 30.) : false;
    jetveto_Flag_.push_back(vetoFlag);
    jetveto_eep_Flag_.push_back(vetoEEP);
    
    // Flavor and QGL/PUID
    hadronflav_.push_back(ak4jet.hadronFlavour());
    partonflav_.push_back(ak4jet.partonFlavour());
    qgl_.push_back(ak4jet.userFloat("qgLikelihood"));
    PUID_.push_back(ak4jet.userFloat("pileupJetId_fullDiscriminant"));
    
    // Charge-weighted pt sums
    std::vector<reco::CandidatePtr> daughters = ak4jet.daughterPtrVector();
    std::sort(daughters.begin(), daughters.end(), [](const reco::CandidatePtr &a, const reco::CandidatePtr &b) {
      return a->pt() > b->pt();
 });
    
    float sumptchg_kp3 = 0, sumptchg_kp6 = 0, sumptchg_k1 = 0, chg_ptsum = 0;
    for (const auto &cand : daughters) {
      sumptchg_kp3 += pow(cand->pt(), 0.3) * cand->charge();
      sumptchg_kp6 += pow(cand->pt(), 0.6) * cand->charge();
      sumptchg_k1  += pow(cand->pt(), 1.0) * cand->charge();
      if (std::abs(cand->charge()) > 0)
        chg_ptsum += cand->pt();
    }
    
    float jetPt = pfjetAK4_4v.Pt();
    charge_kappa_0p3_.push_back((daughters.size() > 0 && jetPt > 0) ? sumptchg_kp3 / pow(jetPt, 0.3) : 0);
    charge_kappa_0p6_.push_back((daughters.size() > 0 && jetPt > 0) ? sumptchg_kp6 / pow(jetPt, 0.6) : 0);
    charge_kappa_1p0_.push_back((daughters.size() > 0 && jetPt > 0) ? sumptchg_k1 / jetPt : 0);
    charged_ptsum_.push_back(chg_ptsum);

    if(nAK4PUPPIJet_ >= njetmx) { break;}
    nAK4PUPPIJet_++;
    
  }
  
	}

	std::cout<<"*****************END Jets***********************"<<std::endl;
}



