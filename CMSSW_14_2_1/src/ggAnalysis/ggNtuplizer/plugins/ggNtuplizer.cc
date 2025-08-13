#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" 
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

using namespace std;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) :  
  //hltPrescaleProvider_(ps, consumesCollector(), *this)
  ecalClusterToolsESGetTokens_{consumesCollector()} ,
  caloTop(esConsumes()),
  randomGenerator_(12345),normalDistribution_(0.0, 1.0)
{
  std::cout << "DEBUG: ggNtuplizer constructor called!" << std::endl;

  development_               = ps.getParameter<bool>("development");
  addFilterInfoMINIAOD_      = ps.getParameter<bool>("addFilterInfoMINIAOD");
  doGenParticles_            = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_          = ps.getParameter<bool>("runOnParticleGun");
  runOnSherpa_               = ps.getParameter<bool>("runOnSherpa");
  runL1ECALPrefire_          = ps.getParameter<bool>("runL1ECALPrefire");
  dumpPFPhotons_             = ps.getParameter<bool>("dumpPFPhotons");
  dumpPDFSystWeight_         = ps.getParameter<bool>("dumpPDFSystWeight");
  year_                      = ps.getParameter<int>("year");

  std::cout << "DEBUG: Basic parameters OK" << std::endl;
  
  vtxLabel_                  = consumes<reco::VertexCollection>        (ps.getParameter<InputTag>("VtxLabel"));
  rhoLabel_                  = consumes<double>                        (ps.getParameter<InputTag>("rhoLabel"));

  generatorLabel_            = consumes<GenEventInfoProduct>           (ps.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_             = consumes<LHEEventProduct>               (ps.getParameter<InputTag>("LHEEventLabel"));
  puCollection_              = consumes<vector<PileupSummaryInfo> >    (ps.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (ps.getParameter<InputTag>("genParticleSrc"));
  electronCollection_        = consumes<View<pat::Electron> >          (ps.getParameter<InputTag>("electronSrc"));
  gsfTracks_                 = consumes<View<reco::GsfTrack>>          (ps.getParameter<InputTag>("gsfTrackSrc"));

  photonCollection_          = consumes<View<pat::Photon> >            (ps.getParameter<InputTag>("photonSrc"));
  muonCollection_            = consumes<View<pat::Muon> >              (ps.getParameter<InputTag>("Muons"));
  ebReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("reducedEcalRecHitsEB"));
  eeReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("reducedEcalRecHitsEE"));
  esReducedRecHitCollection_ = consumes<EcalRecHitCollection>          (ps.getParameter<InputTag>("reducedEcalRecHitsES")); 
 
  pckPFCandidateCollection_  = consumes<pat::PackedCandidateCollection>(ps.getParameter<InputTag>("packedPFCands"));
  pckPFCdsLabel_             = consumes<vector<pat::PackedCandidate>>  (ps.getParameter<InputTag>("packedPFCands"));
  recoCdsLabel_              = consumes<View<reco::Candidate>>         (ps.getParameter<InputTag>("packedPFCands"));

  std::cout << "DEBUG: Input tags OK" << std::endl;

  
  // Added by me
  tok_pfjetAK4s_ = consumes<edm::View<pat::Jet>>(ps.getParameter<edm::InputTag>("PFJetsAK4"));
  tok_genjetAK4s_= consumes<reco::GenJetCollection>( ps.getParameter<edm::InputTag>("GENJetAK4"));
  tok_Rho_ = consumes<double>(ps.getParameter<InputTag>("PFRho"));
  min_pt_AK4jet = ps.getUntrackedParameter<double>("minJetPt", 25.);
  max_eta = ps.getUntrackedParameter<double>("maxEta",3.);
  year		=  ps.getUntrackedParameter<string>("YEAR","2018");
  isData    = ps.getUntrackedParameter<bool>("Data",false);
  isMC      = ps.getUntrackedParameter<bool>("MonteCarlo", false);
  isRun3 	= ps.getUntrackedParameter<bool>("isRun3", false);
  isUltraLegacy = ps.getUntrackedParameter<bool>("UltraLegacy", false);
  tok_muons_ = consumes<edm::View<pat::Muon>> ( ps.getParameter<edm::InputTag>("Muons"));
  mJetVetoMap = ps.getParameter<std::string>("JetVetoMap");
  
  store_electron_scalnsmear = ps.getUntrackedParameter<bool>("store_electron_scalnsmear", false);
  store_electrons = ps.getUntrackedParameter<bool>("store_electrons", false);
  store_muons = ps.getUntrackedParameter<bool>("store_muons", false);
  store_photons = ps.getUntrackedParameter<bool>("store_photons", false);
  store_ak4jets = ps.getUntrackedParameter<bool>("store_ak4jets", false);
  store_CHS_met   = ps.getUntrackedParameter<bool>("store_CHS_met", false);
  store_PUPPI_met = ps.getUntrackedParameter<bool>("store_PUPPI_met", false);
  store_electron_idSF      = ps.getUntrackedParameter<bool>("store_electron_idSF", false);
     
  //Scale and Smearing for electron and photon
  dataYear_ = ps.getParameter<int>("dataYear");
  dataPeriod_ = ps.getParameter<std::string>("dataPeriod");
  applyEGMCorrections_ = ps.getParameter<bool>("applyEGMCorrections");

  store_CHS_met   = ps.getUntrackedParameter<bool>("store_CHS_met", false);
  store_PUPPI_met = ps.getUntrackedParameter<bool>("store_PUPPI_met", false);
  tok_mets_= consumes<pat::METCollection> ( ps.getParameter<edm::InputTag>("PFMet"));
  tok_mets_PUPPI_ = consumes<pat::METCollection> ( ps.getParameter<edm::InputTag>("PuppiMet"));

  tok_METfilters_ = consumes<edm::TriggerResults> ( ps.getParameter<edm::InputTag>("MET_Filters"));
  
  isData_ = ps.getParameter<bool>("isData");
  triggerBits_ = consumes<edm::TriggerResults> ( ps.getParameter<edm::InputTag>("bits"));
  triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("TriggerObjects"));
  triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(ps.getParameter<edm::InputTag>("prescales"));

  std::cout << "DEBUG: Generator input tags OK" << std::endl;
  
  //JEC and JER Files

  mjecL1FastFileAK4         = ps.getParameter<std::string>("jecL1FastFileAK4");
  mjecL2RelativeFileAK4     = ps.getParameter<std::string>("jecL2RelativeFileAK4");
  mjecL3AbsoluteFileAK4     = ps.getParameter<std::string>("jecL3AbsoluteFileAK4");
  mjecL2L3ResidualFileAK4   = ps.getParameter<std::string>("jecL2L3ResidualFileAK4");
  mJECUncFileAK4           = ps.getParameter<std::string>("JECUncFileAK4");

  mPtResoFileAK4           = ps.getParameter<std::string>("PtResoFileAK4");
  mPtSFFileAK4             = ps.getParameter<std::string>("PtSFFileAK4");

  std::cout << "JEC and JER Files"<< std::endl;
  
  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  std::cout << "tree_"<< std::endl;
  hEvents_ = fs->make<TH1F>("hEvents","total processed and skimmed events",   2,  0,   2);
  std::cout << "hEvents_"<< std::endl;

  branchesGlobalEvent(tree_);
  std::cout << "branchesGlobalEvent"<< std::endl;

  
  if (doGenParticles_) {
    branchesGenInfo(tree_, fs);
    std::cout << "branchesGenInfo"<< std::endl;
    branchesGenPart(tree_);
    } 
  
  std::cout << "Gen particle "<< std::endl;
 
  branchesPhotons(tree_);
  branchesElectrons(tree_);
  branchesAK4PUPPIJets(tree_);
  branchesMET(tree_);
  branchesMuons(tree_);
  
  std::cout << "leaving ggNtuplizer constructor"<<std::endl;
  
}

ggNtuplizer::~ggNtuplizer() {

}

//Initializing JER and JEC files

void ggNtuplizer::beginJob() {
  std::cout << "DEBUG: ggNtuplizer beginJob() called!" << std::endl;
  
  L1FastAK4       = new JetCorrectorParameters(mjecL1FastFileAK4.c_str());
  L2RelativeAK4   = new JetCorrectorParameters(mjecL2RelativeFileAK4.c_str());
  L3AbsoluteAK4   = new JetCorrectorParameters(mjecL3AbsoluteFileAK4.c_str());
  L2L3ResidualAK4 = new JetCorrectorParameters(mjecL2L3ResidualFileAK4.c_str());

  vecL1FastAK4.push_back(*L1FastAK4);
  vecL2RelativeAK4.push_back(*L2RelativeAK4);
  vecL3AbsoluteAK4.push_back(*L3AbsoluteAK4);
  vecL2L3ResidualAK4.push_back(*L2L3ResidualAK4);

  jecL1FastAK4       = new FactorizedJetCorrector(vecL1FastAK4);
  jecL2RelativeAK4   = new FactorizedJetCorrector(vecL2RelativeAK4);
  jecL3AbsoluteAK4   = new FactorizedJetCorrector(vecL3AbsoluteAK4);
  jecL2L3ResidualAK4 = new FactorizedJetCorrector(vecL2L3ResidualAK4);

  //scale and smearing for electron and photon          
   if (applyEGMCorrections_) {
     std::cout << "DEBUG: Attempting to initialize EGMCorrectionManager with:" << std::endl;
        std::cout << "  dataYear_: " << dataYear_ << std::endl;
        std::cout << "  dataPeriod_: " << dataPeriod_ << std::endl;
    
	egmCorrectionManager_ = std::make_unique<EGMCorrectionManager>(dataYear_, dataPeriod_);
	std::cout<<"egmCorrectionManager_"<<std::endl;
   }

	egmIDSFManager_ = std::make_unique<EGMIDSFManager>(dataYear_, dataPeriod_);

	std::cout << "DEBUG: ggNtuplizer beginJob() completed successfully!" << std::endl;
   
}

void ggNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {
  std::cout << "=== ANALYZE FUNCTION CALLED FOR EVENT "<< " ===" << std::endl;
  std::cout << "DEBUG: ggNtuplizer analyze() called for event " << e.id() << std::endl;

  hEvents_->Fill(0.5);
  
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  edm::Handle<edm::View<pat::Muon>> muons;
  e.getByToken(tok_muons_, muons);

  reco::Vertex vtx;

  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    bool isFake = (v->chi2() == 0 && v->ndof() == 0);
    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v;
      break;
    }
  }

  //  initTriggerFilters(e);
  fillGlobalEvent(e, es);

  if (!e.isRealData()) {
    fillGenInfo(e);
    if (doGenParticles_)
      fillGenPart(e);
  }

  fillPhotons(e, es);
  fillElectrons(e, es, pv);
  fillAK4PUPPIJets(e, es);
  fillMET(e, es);
  fillMuons(e, pv, vtx);
  hEvents_->Fill(1.5);
  
  tree_->Fill();
  std::cout << "DEBUG: Tree filled successfully" << std::endl;
 
}

void ggNtuplizer::endJob() {

  std::cout<<"EndJob ggNtuplizer"<<std::endl;
    // Clean up JEC corrector objects
    if (jecL1FastAK4) delete jecL1FastAK4;
    if (jecL2RelativeAK4) delete jecL2RelativeAK4;
    if (jecL3AbsoluteAK4) delete jecL3AbsoluteAK4;
    if (jecL2L3ResidualAK4) delete jecL2L3ResidualAK4;
    
    // Clean up JEC parameter objects
    if (L1FastAK4) delete L1FastAK4;
    if (L2RelativeAK4) delete L2RelativeAK4;
    if (L3AbsoluteAK4) delete L3AbsoluteAK4;
    if (L2L3ResidualAK4) delete L2L3ResidualAK4;
    
}

DEFINE_FWK_MODULE(ggNtuplizer);
