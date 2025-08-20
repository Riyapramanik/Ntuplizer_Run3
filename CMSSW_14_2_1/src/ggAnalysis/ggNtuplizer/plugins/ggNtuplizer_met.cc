#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"


using namespace std;

// MET filter booleans //
  
  bool Flag_goodVertices_;
  bool Flag_globalSuperTightHalo2016Filter_;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter_;
  bool Flag_BadPFMuonFilter_;
  bool Flag_BadPFMuonDzFilter_;
  bool Flag_hfNoisyHitsFilter_;
  bool Flag_eeBadScFilter_;
//bool Flag_ecalBadCalibFilter_;

float miset , misphi , sumEt, misetsig;
float miset_covXX, miset_covXY, miset_covYY;
float miset_UnclusEup, miset_UnclusEdn;
float misphi_UnclusEup, misphi_UnclusEdn;

float miset_PUPPI , misphi_PUPPI , sumEt_PUPPI, misetsig_PUPPI;
float miset_PUPPI_covXX, miset_PUPPI_covXY, miset_PUPPI_covYY;
float miset_PUPPI_JESup, miset_PUPPI_JESdn, miset_PUPPI_JERup, miset_PUPPI_JERdn, miset_PUPPI_UnclusEup, miset_PUPPI_UnclusEdn;
float misphi_PUPPI_JESup, misphi_PUPPI_JESdn, misphi_PUPPI_JERup, misphi_PUPPI_JERdn, misphi_PUPPI_UnclusEup, misphi_PUPPI_UnclusEdn;

void ggNtuplizer::branchesMET(TTree* tree) {
  //MET Filters
   tree->Branch("Flag_goodVertices",&Flag_goodVertices_,"Flag_goodVertices_/O");
  tree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter_,"Flag_globalSuperTightHalo2016Filter_/O");
  tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter_,"Flag_EcalDeadCellTriggerPrimitiveFilter_/O");
  tree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter_,"Flag_BadPFMuonFilter_/O");
  tree->Branch("Flag_BadPFMuonDzFilter",&Flag_BadPFMuonDzFilter_,"Flag_BadPFMuonDzFilter_/O");
  tree->Branch("Flag_hfNoisyHitsFilter",&Flag_hfNoisyHitsFilter_,"Flag_hfNoisyHitsFilter_/O");
  tree->Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter_,"Flag_eeBadScFilter_/O");
  //tree->Branch("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter_,"Flag_ecalBadCalibFilter_/O");
  
if(store_CHS_met){
  
  tree->Branch("CHSMET_pt",&miset,"miset/F") ;
  tree->Branch("CHSMET_phi",&misphi,"misphi/F") ;
  tree->Branch("CHSMET_sig",&misetsig,"misetsig/F");
  tree->Branch("CHSMET_sumEt",&sumEt,"sumEt/F");
  
  tree->Branch("CHSMET_covXX",&miset_covXX,"miset_covXX/F") ;
  tree->Branch("CHSMET_covXY",&miset_covXY,"miset_covXY/F") ;
  tree->Branch("CHSMET_covYY",&miset_covYY,"miset_covYY/F") ;
  
  tree->Branch("CHSMET_pt_UnclusEup",&miset_UnclusEup,"miset_CHS_UnclusEup/F") ;
  tree->Branch("CHSMET_pt_UnclusEdn",&miset_UnclusEdn,"miset_CHS_UnclusEdn/F") ;
  tree->Branch("CHSMET_phi_UnclusEup",&misphi_UnclusEup,"CHSMET_phi_UnclusEup/F") ;
  tree->Branch("CHSMET_phi_UnclusEdn",&misphi_UnclusEdn,"CHSMET_phi_UnclusEdn/F") ;
  
  }

 if(store_PUPPI_met){
  
  tree->Branch("PuppiMET_pt",&miset_PUPPI,"miset_PUPPI/F") ;
  tree->Branch("PuppiMET_phi",&misphi_PUPPI,"misphi_PUPPI/F") ;
  tree->Branch("PuppiMET_sig",&misetsig_PUPPI,"misetsig_PUPPI/F");
  tree->Branch("PuppiMET_sumEt",&sumEt_PUPPI,"sumEt_PUPPI/F");
  
  tree->Branch("PuppiMET_covXX",&miset_PUPPI_covXX,"miset_PUPPI_covXX/F") ;
  tree->Branch("PuppiMET_covXY",&miset_PUPPI_covXY,"miset_PUPPI_covXY/F") ;
  tree->Branch("PuppiMET_covYY",&miset_PUPPI_covYY,"miset_PUPPI_covYY/F") ;
  
  tree->Branch("PuppiMET_pt_JESup",&miset_PUPPI_JESup,"miset_PUPPI_JESup/F") ;
  tree->Branch("PuppiMET_pt_JESdn",&miset_PUPPI_JESdn,"miset_PUPPI_JESdn/F") ;
  tree->Branch("PuppiMET_pt_JERup",&miset_PUPPI_JERup,"miset_PUPPI_JERup/F") ;
  tree->Branch("PuppiMET_pt_JERdn",&miset_PUPPI_JERdn,"miset_PUPPI_JERdn/F") ;
  tree->Branch("PuppiMET_pt_UnclusEup",&miset_PUPPI_UnclusEup,"miset_PUPPI_UnclusEup/F") ;
  tree->Branch("PuppiMET_pt_UnclusEdn",&miset_PUPPI_UnclusEdn,"miset_PUPPI_UnclusEdn/F") ;
  tree->Branch("PuppiMET_phi_JESup",&misphi_PUPPI_JESup,"misphi_PUPPI_JESup/F") ;
  tree->Branch("PuppiMET_phi_JESdn",&misphi_PUPPI_JESdn,"misphi_PUPPI_JESdn/F") ;
  tree->Branch("PuppiMET_phi_JERup",&misphi_PUPPI_JERup,"misphi_PUPPI_JERup/F") ;
  tree->Branch("PuppiMET_phi_JERdn",&misphi_PUPPI_JERdn,"misphi_PUPPI_JERdn/F") ;
  tree->Branch("PuppiMET_phi_UnclusEup",&misphi_PUPPI_UnclusEup,"misphi_PUPPI_UnclusEup/F") ;
  tree->Branch("PuppiMET_phi_UnclusEdn",&misphi_PUPPI_UnclusEdn,"misphi_PUPPI_UnclusEdn/F") ;
  
  }
}

void ggNtuplizer::fillMET(const edm::Event& e, const edm::EventSetup& es) {
  //Initialize MET filter flags
  Flag_goodVertices_ = false;
  Flag_globalSuperTightHalo2016Filter_ = false;
  Flag_EcalDeadCellTriggerPrimitiveFilter_ = false;
  Flag_BadPFMuonFilter_ = false;
  Flag_BadPFMuonDzFilter_ = false;
  Flag_hfNoisyHitsFilter_ = false;
  Flag_eeBadScFilter_ = false;
  //Flag_ecalBadCalibFilter_ = false;
// MET uncertainty ids are taken from: https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/MET.h#L152-L158 //
  
  // CHS MET //
  if(store_CHS_met) {
  edm::Handle<pat::METCollection> pfmet_ ;
  e.getByToken(tok_mets_,pfmet_) ;
  
  if(pfmet_.isValid()){
    
    const pat::MET &met = pfmet_->front();
    
    miset = met.corPt(); //met.pt();
    misphi = met.corPhi();//met.phi();
    misetsig = met.metSignificance();
    sumEt = met.corSumEt();//sumEt();

    if(sumEt<30){return;}
            
    miset_covXX = met.getSignificanceMatrix().At(0,0);
    miset_covXY = met.getSignificanceMatrix().At(0,1);
    miset_covYY = met.getSignificanceMatrix().At(1,1);
    
    //MET uncertainty numbering scheme: https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/MET.h 
	    
    miset_UnclusEup = met.shiftedPt(pat::MET::UnclusteredEnUp);  // met.shiftedPt(pat::MET::METUncertainty(10));
    miset_UnclusEdn = met.shiftedPt(pat::MET::UnclusteredEnDown);// met.shiftedPt(pat::MET::METUncertainty(11));
	
    misphi_UnclusEup = met.shiftedPhi(pat::MET::UnclusteredEnUp);  //(pat::MET::METUncertainty(10));
    misphi_UnclusEdn = met.shiftedPhi(pat::MET::UnclusteredEnDown);//(pat::MET::METUncertainty(11));
  }
	    
  }//store_CHS_met
  
  // PUPPI MET //
  
  if(store_PUPPI_met){
  
  edm::Handle<pat::METCollection> pfmet_PUPPI_ ;
  e.getByToken(tok_mets_PUPPI_,pfmet_PUPPI_) ;
  
  if(pfmet_PUPPI_.isValid()){
    
    const pat::MET &met = pfmet_PUPPI_->front();
    
    miset_PUPPI = met.corPt(); 
    misphi_PUPPI = met.corPhi();
    misetsig_PUPPI = met.metSignificance();
    sumEt_PUPPI = met.corSumEt();

    if(sumEt_PUPPI<30){return;}
    
    miset_PUPPI_covXX = met.getSignificanceMatrix().At(0,0);
    miset_PUPPI_covXY = met.getSignificanceMatrix().At(0,1);
    miset_PUPPI_covYY = met.getSignificanceMatrix().At(1,1);
    
    //MET uncertainty numbering scheme: https://cmssdt.cern.ch/lxr/source/DataFormats/PatCandidates/interface/MET.h
    
    miset_PUPPI_JESup = met.shiftedPt(pat::MET::JetEnUp); //(pat::MET::METUncertainty(2));
    miset_PUPPI_JESdn = met.shiftedPt(pat::MET::JetEnDown); //(pat::MET::METUncertainty(3));
    miset_PUPPI_JERup = met.shiftedPt(pat::MET::JetResUp); //(pat::MET::METUncertainty(0));
    miset_PUPPI_JERdn = met.shiftedPt(pat::MET::JetResDown); //(pat::MET::METUncertainty(1));
    miset_PUPPI_UnclusEup = met.shiftedPt(pat::MET::UnclusteredEnUp);  //(pat::MET::METUncertainty(10));
    miset_PUPPI_UnclusEdn = met.shiftedPt(pat::MET::UnclusteredEnDown);//(pat::MET::METUncertainty(11));
    
    misphi_PUPPI_JESup = met.shiftedPhi(pat::MET::JetEnUp); //(pat::MET::METUncertainty(2));
    misphi_PUPPI_JESdn = met.shiftedPhi(pat::MET::JetEnDown); //(pat::MET::METUncertainty(3));
    misphi_PUPPI_JERup = met.shiftedPhi(pat::MET::JetResUp); //(pat::MET::METUncertainty(0));
    misphi_PUPPI_JERdn = met.shiftedPhi(pat::MET::JetResDown); //(pat::MET::METUncertainty(1));
    misphi_PUPPI_UnclusEup = met.shiftedPhi(pat::MET::UnclusteredEnUp);  //(pat::MET::METUncertainty(10));
    misphi_PUPPI_UnclusEdn = met.shiftedPhi(pat::MET::UnclusteredEnDown);//(pat::MET::METUncertainty(11));
	
    //See DataFormats/PatCandidates/interface/MET.h for the names of uncertainty sources //
    
  }
  
  }
  
  //MET Filter
  edm::Handle<edm::TriggerResults> METFilterResults;
  e.getByToken(tok_METfilters_, METFilterResults);
  
  if(METFilterResults.isValid()){
    const edm::TriggerNames & metfilterName = e.triggerNames(*METFilterResults);
    //const std::vector<std::string>& names = metfilterName.triggerNames();
        
    //Flag_goodVertices
    unsigned int goodVerticesIndex_ = metfilterName.triggerIndex("Flag_goodVertices");
    Flag_goodVertices_ = METFilterResults.product()->accept(goodVerticesIndex_);
    Flag_globalSuperTightHalo2016Filter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_globalSuperTightHalo2016Filter"));
    Flag_EcalDeadCellTriggerPrimitiveFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter"));
    Flag_BadPFMuonFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_BadPFMuonFilter"));
    Flag_BadPFMuonDzFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_BadPFMuonDzFilter"));
    Flag_hfNoisyHitsFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_hfNoisyHitsFilter"));
    Flag_eeBadScFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_eeBadScFilter"));
  //Flag_ecalBadCalibFilter_ = METFilterResults.product()->accept(metfilterName.triggerIndex("Flag_ecalBadCalibFilter")); // Not applicable for 2022 and 2023

  // End of MET filters //
  }//(METFilterResults.isValid())


}
