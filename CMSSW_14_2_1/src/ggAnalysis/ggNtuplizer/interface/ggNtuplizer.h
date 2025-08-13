#ifndef ggNtuplizer_h
#define ggNtuplizer_h

#include "TTree.h"
#include "TH2D.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include <vector>

//Scale and smearing for Photon                                                                                                                        
#include "EGMCorrectionManager.h"
#include <random>
#include "EGMIDSFManager.h"

using namespace std;
using namespace edm;
using namespace reco;

void setbit(UShort_t& x, UShort_t bit);

class ggNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
 public:
  virtual void beginJob() override;
  virtual void endJob() override;
  
  explicit ggNtuplizer(const edm::ParameterSet&);
  ~ggNtuplizer(); 
 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
 
  void initTriggerFilters(const edm::Event&);
  ULong64_t matchSingleElectronTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchDoubleElectronTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchSinglePhotonTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchDoublePhotonTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchTriplePhotonTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchMuonTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchJetTriggerFilters(double pt, double eta, double phi);
  ULong64_t matchL1TriggerFilters(double pt, double eta, double phi);
  Double_t deltaPhi(Double_t phi1, Double_t phi2);
  Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2);
  Double_t getMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands, const reco::Candidate* ptcl,  
			    double r_iso_min, double r_iso_max, double kt_scale, bool charged_only);

  void branchesGlobalEvent(TTree*);
  void branchesGenInfo    (TTree*, edm::Service<TFileService>&);
  void branchesGenPart    (TTree*);
  void branchesMET        (TTree*);
  void branchesPhotons    (TTree*);
  void branchesElectrons  (TTree*);
  void branchesMuons      (TTree*);
  void fillGlobalEvent(const edm::Event&, const edm::EventSetup&);
  void fillGenInfo    (const edm::Event&);
  void fillGenPart    (const edm::Event&);
  void fillMET        (const edm::Event&, const edm::EventSetup&);
  void fillPhotons    (const edm::Event&, const edm::EventSetup&);
  void fillElectrons  (const edm::Event&, const edm::EventSetup&, math::XYZPoint&);
  void fillMuons      (const edm::Event&, math::XYZPoint&, const reco::Vertex);
  
  void branchesAK4PUPPIJets(TTree* tree);
  void fillAK4PUPPIJets(const edm::Event& e, const edm::EventSetup& es);
 
  bool addFilterInfoMINIAOD_;  
  bool doGenParticles_;
  bool runOnParticleGun_;
  bool runOnSherpa_;
  bool runL1ECALPrefire_;
  bool dumpPFPhotons_;
  bool dumpPDFSystWeight_;
  int  year_;

  double min_pt_AK4jet;
  double max_eta;

  const EcalClusterLazyTools::ESGetTokens ecalClusterToolsESGetTokens_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTop;
  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  edm::EDGetTokenT<double>                         rhoLabel_;
  edm::EDGetTokenT<edm::TriggerResults>            patTrgResultsLabel_;
  edm::EDGetTokenT<GenEventInfoProduct>            generatorLabel_;
  edm::EDGetTokenT<LHEEventProduct>                lheEventLabel_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<vector<reco::GenParticle> >     genParticlesCollection_;
  edm::EDGetTokenT<edm::View<pat::MET> >           pfMETlabel_;
  edm::EDGetTokenT<edm::View<pat::MET> >           puppiMETlabel_;
  edm::EDGetTokenT<edm::View<pat::Electron> >      electronCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> >        photonCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> >          muonCollection_;
  edm::EDGetTokenT<vector<pat::Tau> >              tauCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           esReducedRecHitCollection_; 
  edm::EDGetTokenT<reco::PhotonCollection>         recophotonCollection_;
  edm::EDGetTokenT<edm::View<reco::GsfTrack> >     gsfTracks_;
  edm::EDGetTokenT<reco::PFCandidateCollection>    pfAllParticles_;
  edm::EDGetTokenT<vector<pat::PackedCandidate> >  pckPFCdsLabel_;
  edm::EDGetTokenT<edm::View<reco::Candidate> >    recoCdsLabel_;
  
  edm::EDGetTokenT<reco::JetTagCollection>         boostedDoubleSVLabel_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pckPFCandidateCollection_;

  edm::EDGetTokenT<edm::View<pat::Jet> >           tok_pfjetAK4s_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  
  // for MET filters
  edm::EDGetTokenT<edm::TriggerResults> tok_METfilters_;
  edm::EDGetTokenT<pat::METCollection>tok_mets_, tok_mets_PUPPI_;
  //check
  edm::EDGetToken gsfEle_;

  //Add by me
  edm::EDGetTokenT<double> tok_Rho_;
  edm::EDGetTokenT<reco::GenJetCollection>tok_genjetAK4s_;
  edm::EDGetTokenT<edm::View<pat::Muon>> tok_muons_;

  string year;
  bool isData;
  bool isMC;
  bool isRun3;
  bool isUltraLegacy;
  bool store_electrons, store_muons, store_photons, store_ak4jets, store_CHS_met, store_PUPPI_met;

  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  std::unique_ptr<EGMCorrectionManager> egmCorrectionManager_;
  // ID Scale Factor Manager
  std::unique_ptr<EGMIDSFManager> egmIDSFManager_;
  
  TTree   *tree_;
  TH1F    *hEvents_;
  TH1F    *hPU_;
  TH1F    *hPUTrue_;
  TH1F    *hGenWeight_;
  TH1F    *hSumGenWeight_;
  
  std::mt19937 randomGenerator_;
  std::normal_distribution<double> normalDistribution_;
  int dataYear_;
  std::string dataPeriod_;
  bool applyEGMCorrections_;
  bool isData_;


  //JEC and JER correction file
  std::string mjecL1FastFileAK4, mjecL2RelativeFileAK4, mjecL3AbsoluteFileAK4, mjecL2L3ResidualFileAK4;
  std::string mJECUncFileAK4;
  std::string mPtResoFileAK4, mPtSFFileAK4;
  JetCorrectorParameters *L1FastAK4, *L2RelativeAK4, *L3AbsoluteAK4, *L2L3ResidualAK4;
  vector<JetCorrectorParameters> vecL1FastAK4, vecL2RelativeAK4, vecL3AbsoluteAK4, vecL2L3ResidualAK4;
  FactorizedJetCorrector *jecL1FastAK4, *jecL2RelativeAK4, *jecL3AbsoluteAK4, *jecL2L3ResidualAK4;
  std::vector<JetCorrectionUncertainty*> vsrc;

  //Jet veto flag
  std::string mJetVetoMap;
  TFile* file_jetvetomap;
  TH2D* h_jetvetomap;
  TH2D* h_jetvetomap_eep;

};
#endif
