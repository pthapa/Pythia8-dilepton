
// -*- C++ -*-
//
// Package:    GenStudy/Dimuon
// Class:      Dimuon
// 
/**\class Dimuon Dimuon.cc GenStudy/Dimuon/plugins/Dimuon.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Shawn Gregory Zaleski
//         Created:  Tue, 30 Jun 2015 14:01:36 GMT
//
//
    
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include <FWCore/ServiceRegistry/interface/Service.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <TTree.h>
#include <TVector2.h>
#include <TH1F.h>
#include <TH2F.h>
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include <vector>
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Additional includes
//#include "Pythia8Plugins/HepMC2.h"
//#include "HepMC/IO_GenEvent.h"
//#include <Pythia8Plugins/HepMC2.h> 
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <Pythia8/Pythia.h>
#include <cmath>
#include <complex>

// class declaration
//

class Dimuon : public edm::EDAnalyzer {
public:
  explicit Dimuon(const edm::ParameterSet&);
  ~Dimuon();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

protected:
  
 std::auto_ptr<Pythia8::Pythia> fMasterGen;
  //  Pythia8::Pythia *fMasterGen;
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const reco::Candidate* getDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getLastDaughter(const reco::Candidate* part,int pid);
  const reco::Candidate* getBoson( const reco::GenParticleCollection& genParts);
  const reco::Candidate* getMother(const reco::Candidate* part, int pid);
  //const reco::Candidate* getDYBoson(const reco::Candidate* part int pid)
  bool isBoson(int pid);
  bool isMuon(int pid);
  bool checkBosonStatus(const reco::GenParticleCollection& genParts);


  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;



  struct P4Struct {
    float energy,et,eta,phi,pt,mass,theta;
    void fill(const math::XYZTLorentzVector& p4){
      if(p4.Pt()!=0 && p4.Et()!=0){
	energy = p4.E();
	et = p4.Et();
	eta = p4.Eta();
	phi = p4.Phi();
	pt =p4.Pt();
	mass = p4.mag();
	theta = p4.Theta();
      }else clear();
    }
    void clear(){energy=-999;et=-999;eta=-999;phi=-0;pt=-999;mass=-999;theta=-999;}
    static std::string contents(){return "energy/F:et:eta:phi:pt:mass:theta";}
  }; 

  TTree* tree_;



  TH1F * h_Zmass, *h_Zpt,*h_Zeta,*h_Zphi,*h_Zcharge;
  TH1F *h_muMinusmass,*h_muMinuspt,*h_muMinuseta,*h_muMinusphi,*h_muMinuscharge;
  TH1F *h_muPlusmass,*h_muPluspt,*h_muPluseta,*h_muPlusphi,*h_muPluscharge;
  TH1F *h_dphi,*h_dtheta, *h_dr, *h_thetaMuMinus,*h_thetaMuPlus;
  TH1F *h_massInvar, *h_dimuonPt, *h_dimuonEta, *h_dimuonPhi;
  TH1F *h_cosTheta, *h_tanPhi, *h_csTheta, *h_csPhi;
  TH1F *h_cosThetaPlusInvariantMass, *h_cosThetaMinusInvariantMass;

  TH2F *h2_pt1_vs_pt2,*h2_eta1_vs_eta2,*h2_phi1_vs_phi2;

  TH1F *h_SHat, *h_T_Hat, *h_U_Hat, *h_THat2,*h_UHat2, *h_fractionLR,*h_fractionRL, *h_sumFrac;
  TH1F *h_alphaRunQED;
    
  // Delete it after checking it

TH1F *h_DQuark, *h_UQuark, *h_SQuark, *h_CQuark, *h_BQuark, *h_TQuark;

  P4Struct bosonP4_; // as a sanity check we have the right event...
  P4Struct muMinusP4_;
  P4Struct muPlusP4_;
  int muMinusPID_;
  int muPlusPID_;
  int bosonId_;
  double crossSec, cosTheta, tanPhi, csTheta, csPhi;
  double mCosThetaPlus, mCosThetaMinus;

  // mandelston variables
  double SHat,T_Hat,U_Hat,THat2,UHat2,fractionLR,fractionRL,sumFrac;
  double alphaRunQED;

  //Delete it after checking it 
  int DQuark,UQuark,SQuark,CQuark,BQuark,TQuark;

  int debug_;
  edm::InputTag genPartsTag_;
  int decayParticlePID_;
  edm::InputTag genInfoProduct_;
  edm::EDGetTokenT<GenRunInfoProduct> genInfoProductToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genPartsToken_;


  // Additional declaration for InputTag and Token to access alphaEM

  edm::InputTag genEventInfoProduct_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProductToken_;



  // ----------member data ---------------------------
  
};


void Dimuon::beginJob()
{

  //Acesses pythia parameters in CMSSW environment 
 fMasterGen.reset( new Pythia8::Pythia );
 fMasterGen->init();
 

  edm::Service<TFileService> fs;
  h_Zmass = fs->make<TH1F>("Zmass" , "m", 1000, 0., 600);
  h_Zpt  = fs->make<TH1F>( "Zpt"  , "p_{t}", 500,  0., 2500. );
  h_Zeta = fs->make<TH1F>( "Zeta" , "#eta" , 100, -10., 10.    );
  h_Zphi = fs->make<TH1F>( "Zphi" , "#phi" , 100,  -3.20, 3.20   );
  h_Zcharge = fs->make<TH1F>( "Zcharge" , "Q" ,3,  -1.5, 1.5    );
  h_muMinusmass = fs->make<TH1F>("muMinusmass" , "m", 1000, 0., 500);
  h_muMinuspt  = fs->make<TH1F>( "muMinuspt"  , "p_{t}", 500,  0., 2500. );
  h_muMinuseta = fs->make<TH1F>( "muMinuseta" , "#eta" , 100, -5., 5.    );
  h_muMinusphi = fs->make<TH1F>( "muMinusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muMinuscharge = fs->make<TH1F>( "muMinuscharge" , "Q" ,3,  -1.5, 1.5    );

  h_muPlusmass = fs->make<TH1F>("muPlusmass" , "m", 1000, 0., 500);
  h_muPluspt  = fs->make<TH1F>( "muPluspt"  , "p_{t}", 500,  0., 2500. );
  h_muPluseta = fs->make<TH1F>( "muPluseta" , "#eta" , 100, -5., 5.    );
  h_muPlusphi = fs->make<TH1F>( "muPlusphi" , "#phi" , 100,  -3.15, 3.15   );
  h_muPluscharge = fs->make<TH1F>( "muPluscharge" , "Q" ,3,  -1.5, 1.5    );

  h_dphi = fs->make<TH1F>("delta phi", "#delta #phi", 100, -3.15, 3.15 );       
  h_dtheta = fs->make<TH1F>("delta theta", "#delta #theta", 100, -3.15, 3.15); 
  h_dr = fs->make<TH1F>("delta r", "#delta r", 100, 0, 10);
  h_thetaMuMinus = fs->make<TH1F>("theta muMinus", "#theta", 100, -3.15, 3.15);      
  h_thetaMuPlus = fs->make<TH1F>("theta muPlus", "#theta", 100, -3.15, 3.15); 
  h_massInvar = fs->make<TH1F>("Invariant mass", "Invariant mass", 350, 0., 3500.);
  h_dimuonPt = fs->make<TH1F>("Dimuon Pt", "Dimuon Pt", 500, 0, 2500);
  h_dimuonEta = fs->make<TH1F>("Dimuon eta", "Dimuon #eta", 100, -5, 5);
  h_dimuonPhi = fs->make<TH1F>("Dimuon Phi", "Dimuon #phi", 100, -3.15, 3.15);

  h_cosTheta = fs->make<TH1F>("cosTheta", "cos #theta", 100, -1.01, 1.01);
  h_tanPhi = fs->make<TH1F>("tanPhi", "tan #phi", 100, -1000.0, 1000.0);
  h_csTheta = fs->make<TH1F>("csTheta", "#theta_{CS}", 100, -3.15, 3.15);
  h_csPhi = fs->make<TH1F>("csPhi", "#phi_{CS}", 100, -3.15, 3.15);
  h_cosThetaMinusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaMinus", "InvariantMass_cosThetaMinus", 350, 0., 3500.);
  h_cosThetaPlusInvariantMass = fs->make<TH1F>("InvariantMass_cosThetaPlus", "InvariantMass_cosThetaPlus", 350, 0., 3500.);

  h2_pt1_vs_pt2   = fs->make<TH2F>( "pt1_vs_pt2"   , "p_{t,1} vs. p_{t,2}"   , 500,  0., 2500., 500,  0., 2500.);
  h2_eta1_vs_eta2 = fs->make<TH2F>( "eta1_vs_eta2" , "#eta_{1} vs. #eta_{2}" , 100, -5., 5.   , 100, -5., 5.   );
  h2_phi1_vs_phi2 = fs->make<TH2F>( "phi1_vs_phi2" , "#phi_{1} vs. #phi_{2}" , 100,  -3.15, 3.15  , 100,  -3.15, 3.15  );

  //Book histogram for Mandelston variables
  h_SHat = fs->make<TH1F>("SHat", "SHat", 100, 0.1, 10000000.01);
  h_THat2 = fs->make<TH1F>("THat2", "THat2", 100, 0.1, 10000000000000.01);
  h_UHat2 = fs->make<TH1F>("UHat2", "UHat2", 100, 0.1, 10000000000000.01);
  h_fractionLR= fs->make<TH1F>("fractionLR", "fractionLR", 100, 0.1, 1.01);
  h_fractionRL= fs->make<TH1F>("fractionRL", "fractionRL", 100, 0.1, 1.01);
  h_sumFrac= fs->make<TH1F>("sumFrac", " sumFrac", 100, 0.1, 1.01);

  h_T_Hat = fs->make<TH1F>("T_Hat", "T_Hat", 100, -100000000.1, -1000.01);
  h_U_Hat = fs->make<TH1F>("U_Hat", "U_Hat", 100, -100000000.1, -1000.01);
  h_alphaRunQED = fs->make<TH1F>("alphaRunQED", "alphaRunQED", 100, 0.0001, 0.00001);
 

  //Delete it after checking it

  h_DQuark= fs->make<TH1F>("DQuark", "DQuark", 10, 0.1, 10.01);
  h_UQuark= fs->make<TH1F>("UQuark", "UQuark", 10, 0.1, 10.01);
  h_SQuark= fs->make<TH1F>("SQuark", "SQuark", 10, 0.1, 10.01);
  h_CQuark= fs->make<TH1F>("CQuark", "CQuark", 10, 0.1, 10.01);
  h_BQuark= fs->make<TH1F>("BQuark", "BQuark", 10, 0.1, 10.01);
  h_TQuark= fs->make<TH1F>("TQuark", "TQuark", 10, 0.1, 10.01);
 


  tree_= fs->make<TTree>("pdfTree","PDF Tree");
  // tree_->Branch("evtId",&evtId_,EventId::contents().c_str());
  tree_->Branch("bosonP4",&bosonP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1P4",&muMinusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay2P4",&muPlusP4_,P4Struct::contents().c_str());
  tree_->Branch("decay1PID",&muMinusPID_,"decay1PID/I");
  tree_->Branch("decay2PID",&muPlusPID_,"decay2PID/I");
  tree_->Branch("bosonPID",&bosonId_,"bosonPID/I");
  tree_->Branch("crossSec", &crossSec, "crossSec/D");
  tree_->Branch("cosTheta", &cosTheta, "cosTheta/D");
  tree_->Branch("tanPhi", &tanPhi, "tanPhi/D");
  tree_->Branch("csTheta", &csTheta, "csTheta/D");
  tree_->Branch("csPhi", &csPhi, "csPhi/D");
  tree_->Branch("mCosThetaPlus", &mCosThetaPlus, "mCosThetaPlus/D");
  tree_->Branch("mCosThetaMinus", &mCosThetaMinus, "mCosThetaMinus/D");
  // tree_->Branch("pdfInfo",&pdfInfo_,PDFInfo::contents().c_str());

  // Tree for Mandelston variabales

  tree_->Branch("SHat", &SHat, "SHat/D");
  tree_->Branch("THat2", &THat2, "THat2/D");
  tree_->Branch("UHat2", &UHat2, "UHat2/D");
  tree_->Branch("fractionLR", &fractionLR, "fractionLR/D");
  tree_->Branch("fractionRL", &fractionRL, "fractionRL/D");

  tree_->Branch("sumFrac", &sumFrac, "sumFrac/D");

  tree_->Branch("T_Hat", &T_Hat, "T_Hat/D");
  tree_->Branch("U_Hat", &U_Hat, "U_Hat/D");
  tree_->Branch("alphaRunQED", &alphaRunQED, "alphaRunQED/D");
  

  // Delete it after checking
  tree_->Branch("DQuark", &DQuark, "DQuark/I");
  tree_->Branch("UQuark", &UQuark, "UQuark/I");
  tree_->Branch("SQuark", &SQuark, "SQuark/I");
  tree_->Branch("CQuark", &CQuark, "CQuark/I");
  tree_->Branch("BQuark", &BQuark, "BQuark/I");
  tree_->Branch("TQuark", &TQuark, "TQuark/I");
 

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Dimuon::Dimuon(const edm::ParameterSet& iConfig)

{


  debug_=iConfig.getParameter<int>("debug");
  genPartsTag_=iConfig.getParameter<edm::InputTag>("genPartsTag");
  decayParticlePID_ = iConfig.getParameter<int>("decayParticlePID");
  genInfoProduct_ = iConfig.getParameter<edm::InputTag>("genInfoProduct");
  
  //now do what ever initialization is needed

  genInfoProductToken_ = consumes<GenRunInfoProduct,edm::InRun>(genInfoProduct_);
  genPartsToken_ = consumes<reco::GenParticleCollection>(genPartsTag_);

  // Additional initialization 

  genEventInfoProduct_ = iConfig.getParameter<edm::InputTag>("genEventInfoProduct");
  genEventInfoProductToken_ = consumes<GenEventInfoProduct,edm::InEvent>(genEventInfoProduct_);

}


Dimuon::~Dimuon()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Dimuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genPartsHandle;
  iEvent.getByToken(genPartsToken_,genPartsHandle);
  const reco::GenParticleCollection& genParts = *genPartsHandle;

  // Object declaration using additional infromation
  edm::Handle<GenEventInfoProduct> genEventInfoProductHandle;
  iEvent.getByToken(genEventInfoProductToken_, genEventInfoProductHandle);
  const GenEventInfoProduct& genEventInfoProduct = *genEventInfoProductHandle;


  const reco::Candidate* boson;
  const reco::Candidate* daughter1;
  const reco::Candidate* daughter2;
  const reco::Candidate* mother1;
  //const reco::Candidate* mother2; 
  const reco::Candidate* muMinus;
  const reco::Candidate* muPlus;
  math::XYZTLorentzVectorD dimuon;
  double dimuonPx, dimuonPy, dimuonPz, dimuonPt, pseudorapidity, Phi, mu1Energy, mu2Energy, dimuonQ;
  double thetaCos, thetaCS, phiTan, phiCS;
  double muPlusKPlus, muPlusKMinus, muMinusKPlus, muMinusKMinus, invariantK;

  
  std::cout << "######### IMPLEMENTATION OF STEVE'S CODE:BEFORE FOR LOOP  ################" << std::endl;
  
  
  // Initialize lambdaSquare and etaLR using defalult value 
  //double qCLambda2 = pow(fMasterGen->settings.parm("ContactInteractions:Lambda"),2);
  double qCLambda=1000.0;
  double qCLambda2=qCLambda*qCLambda;
  std::cout << "############### Test_lambdaSquare: defalut value for lambda is 1000 GeV #############################"<< std::endl;
  std::cout <<"Test_lambdaSquare: " <<qCLambda2 << std::endl;
  //int qCetaLR = fMasterGen->settings.mode("ContactInteractions:etaLR"); 
  int qCetaLR =-1;
  std::cout << "############### Test_etaLR: defalut value for etaLR is 0.00  #############################"<< std::endl;
  std::cout <<"Test_etaLR: " << qCetaLR << std::endl;

  // Running coupling constat-> changes event by event and accessible in miniAOD
  double genEventalphaEM =genEventInfoProduct.alphaQED();
  double alpEM=genEventalphaEM;

  std::cout << " ########################################" <<std::endl;

  std::cout << " genEventalphaEM is  alpEM :  " << alpEM <<std::endl;


  std::cout << " ########################################" <<std::endl;


  // May be this one is right
  //  double alpElectromagnetic= fMasterGen->settings.parm("StandardModel:alphaEM0");  
  //std::cout << "alpElectromagnetic: " << alpElectromagnetic <<std::endl; 

  /*
  int qCetaRL = fMasterGen->settings.mode("ContactInteractions:etaRL");
  double alpElectromagnetic= fMasterGen->settings.parm("StandardModel:alphaEM0");
  std::cout << "alpElectromagnetic: " << alpElectromagnetic <<std::endl;
  double qCetaLL =fMasterGen->settings.mode("ContactInteractions:etaLL");
  double qCetaRR= fMasterGen->settings.mode("ContactInteractions:etaRR");;
  //double alpEM_Running = fMasterGen->couplingsPtr->alphaEM(infoPtr.Q2Ren());

  //std::cout <<"alpEM running"<<  alpEM_Running  << std::endl;  
  std::cout << "qCetaLR: " << qCetaLR << std::endl;
  std::cout <<"qCetaRL: " << qCetaRL << std::endl; 
  std::cout <<"qCLambda2: " << qCLambda2 << std::endl;
  std::cout <<"qCetaLL: " << qCetaLL << std::endl;
  std::cout <<"qCetaRR: " << qCetaRR << std::endl;
  */
  
  
  const reco::Candidate* quark;
  const reco::Candidate* antiquark;
  const reco::Candidate* leptonPlus;
  const reco::Candidate* leptonMinus;

  //Define a complex number
  std::complex<double> I(0.0, 1.0);
  // Complex amplitudes.                                                                                                                                                       
  std::complex<double> meLL(0., 0.);
  std::complex<double> meRR(0., 0.);
  std::complex<double> meLR(0., 0.);
  std::complex<double> meRL(0., 0.);
  std::complex<double> meLR_SM(0., 0.);
  std::complex<double> meRL_SM(0., 0.);

  std::cout << "#### IMPLEMENTATION OF STEVE'S CODE::STOP HERE AND GO AFTER THE FOR LOOP ##########" <<std::endl;
  std::cout << "#### ALL THE NECESSARY VARIABLES ARE DEFINED AND IMPLEMENTED AFTER FOR LOOP ##########" <<std::endl;


  

      for(auto &part : genParts){
	if((part.pdgId() == 1 || part.pdgId() == 2 || part.pdgId() == 3 || part.pdgId() == 4 || part.pdgId() == 5 || part.pdgId() == 6) && 
	   (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){ 

	  

	  	  
	  /*std::cout << "##################### INSIDE FOR LOOP  ###################" << std::endl;
	  std::cout << "I guess Z masss width " << qCGZ << std::endl;

	  std::cout << tmPgLl << " " << tmPgRl << std::endl; 
	  
	  int quarkId=part.pdgId();

	  std::cout << "##################### Quark pdgId  ###################" <<  quarkId << std::endl;
	  double tmmPgvl = 0.25 * fMasterGen->couplingsPtr->vf(quarkId);

	  std::cout <<" tmmPgvl" << " " <<  tmmPgvl << std::endl;  
	  */


	  //Modifiacation to grab dileptons form quark
	  quark=&part;
          if(quark->daughter(0)->pdgId() >0 && quark->daughter(1)->pdgId() <0){


            leptonPlus=quark->daughter(1);
	    leptonMinus=quark->daughter(0);

          }
          else {
            leptonPlus=quark->daughter(0);
            leptonMinus=quark->daughter(1);

          }
	 
	  

	  int leptonPlusId=leptonPlus->pdgId();
	  std::cout << "leptonPlusId" << leptonPlusId << std::endl;
	  double tmmPgal = 0.25 * fMasterGen->couplingsPtr->af(leptonPlusId);
	  std::cout << "Double check for tmPgal" << tmmPgal << std::endl;

	  int leptonMinusId=leptonMinus->pdgId();
	  std::cout << "leptonMinusId" << leptonMinusId << std::endl;
	  double tmmmPgal = 0.25 * fMasterGen->couplingsPtr->af(leptonMinusId);
	  std::cout << "Double check for tmPgal" << tmmmPgal << std::endl;
	  

	  
	  std::cout << "Address of leptoPlus: " << &leptonPlus <<std::endl;
	  std::cout << "Address of leptonMinus:" << &leptonMinus <<std::endl;


	  std::cout << "quark Information " << std::endl;

	  std::cout <<"quark Energy: " << quark->energy()<< std::endl;
	  std::cout <<"quark px: " <<quark->px() <<std::endl;
	  std::cout <<"quark py: " <<quark->py() << std::endl;
	  std::cout <<"quark pz: " << quark->pz() <<std::endl;
	  std::cout <<"quark pt: " << quark->pt() <<std::endl;
	  std::cout <<"quark eta: " << quark->eta() <<std::endl;
	  std::cout <<"quark phi: " << quark->phi() <<std::endl;




	  std::cout << "Lepton Information " << std::endl;


	  std::cout <<"leptonPlus Energy: " << leptonPlus->energy()<< std::endl;

	  std::cout <<"leptonMinus Energy: " << leptonMinus->energy()<< std::endl;

	  if(debug_ > 0){ std::cout << "\nFound the quark! " << "\nQuark is: " << part.pdgId() << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
	    //	    mother1 = getMother(part.mother(0), 2212);
	    //std::cout << "\nQuark mother is:" << mother1->pdgId() << std::endl;
	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
	  }
	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
	  boson = nullptr;
	  if(!daughter1 || !daughter2){
	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
	  }
	}


	// Check antiquark
	else if((part.pdgId() == -1 || part.pdgId() == -2 || part.pdgId() == -3 || part.pdgId() == -4 || part.pdgId() == -5 || part.pdgId() == -6) &&
	   (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){

	  antiquark=&part;

	  std::cout <<"antiquark Energy: " << antiquark->energy()<< std::endl;  
	  
	  std::cout << "antiquark Information " << std::endl;

	  std::cout <<"antiquark Energy: " << antiquark->energy()<< std::endl;
	  std::cout <<"antiquark px: " <<antiquark->px() <<std::endl;
	  std::cout <<"antiquark py: " <<antiquark->py() << std::endl;
	  std::cout <<"antiquark pz: " << antiquark->pz() <<std::endl;
	  std::cout <<"antiquark pt: " << antiquark->pt() <<std::endl;
	  std::cout <<"antiquark eta: " << antiquark->eta() <<std::endl;
	  std::cout <<"antiquark phi: " << antiquark->phi() <<std::endl;

	  /*
	  std::cout << "Lepton Information " << std::endl;
	  
	  std::cout <<"leptonPlus Energy: " << leptonPlus->energy()<< std::endl;
	  std::cout <<"leptonPlus px: " <<leptonPlus->px() <<std::endl;
	  std::cout <<"leptonPlus py: " <<leptonPlus->py() << std::endl;
	  std::cout <<"leptonPlus pz: " << leptonPlus->pz() <<std::endl;
	  std::cout <<"leptonPlus pt: " << leptonPlus->pt() <<std::endl;
	  std::cout <<"leptonPlus eta: " << leptonPlus->eta() <<std::endl;
	  std::cout <<"leptonPlus phi: " << leptonPlus->phi() <<std::endl;


	  std::cout <<"leptonMinus Energy: " << leptonMinus->energy()<< std::endl;
	  std::cout <<"leptonMinus px: " <<leptonMinus->px() <<std::endl;
	  std::cout <<"leptonMinus py: " <<leptonMinus->py() << std::endl;
	  std::cout <<"leptonMinus pz: " << leptonMinus->pz() <<std::endl;
	  std::cout <<"leptonMinus pt: " << leptonMinus->pt() <<std::endl;
	  std::cout <<"leptonMinus eta: " << leptonMinus->eta() <<std::endl;

	  std::cout <<"leptonMinus phi: " << leptonMinus->phi() <<std::endl;


	  */


        }
	





	// Check DY process
    
	else if(part.pdgId() == 23 && (abs(part.daughter(0)->pdgId()) == 11 || abs(part.daughter(0)->pdgId()) == 13)){
	  if(debug_ > 0){std::cout << "\nFound the Z boson! " << "\tStatus is: " << part.status() << "\tNumber of daughters are: " <<
	    part.numberOfDaughters() << "\tFirst daughter is:"  << part.daughter(0)->pdgId() << "\tSecond daughter is: " << part.daughter(1)->pdgId() << std::endl;
	  //      if(part.status() < -20 && part.status() > -30){ std::cout << "\nFound the Z boson!";
	  std::cout << "\nkinematic properties of the particles are: " << std::endl;
	  std::cout << "\npT1: " << part.daughter(0)->pt() << "\tpT2: " << part.daughter(1)->pt() << std::endl;
	  std::cout << "\neta1: " << part.daughter(0)->eta() << "\teta2: " << part.daughter(1)->eta() << std::endl;
	  std::cout << "\nphi1: " << part.daughter(0)->phi() << "\tphi2: " << part.daughter(1)->phi() << std::endl;
	  }
	  daughter1 = getLastDaughter(part.daughter(0), part.daughter(0)->pdgId());
	  daughter2 = getLastDaughter(part.daughter(1), part.daughter(1)->pdgId());
	  std::cout << "\nDaughter particle is: " << daughter1->pdgId() << "tStatus is: " << daughter1->status()
		    << "\tDaughter2 is: " << daughter2->pdgId() << "\tStatus is: " << daughter2->status() << std::endl;
	  mother1 = &part;
	  boson = mother1;
	  if(!boson || !daughter1 || !daughter2){
	    std::cout<<"boson::0x"<<std::hex<<boson<<std::dec<<std::endl;
	    std::cout<<"daughter1::0x"<<std::hex<<daughter1<<std::dec<<std::endl;
	    std::cout<<"daughter2::0x"<<std::hex<<daughter2<<std::dec<<std::endl;
	  }

	}



      }


      std::cout << "##### FIGURE OUT HOW TO DEFINE A FOUR VECTOR IN FUTURE ######## " << std::endl;


      //Mandelston varibles goes here

     
      std::vector<double> tempShat; 
      tempShat.push_back(quark->px()+antiquark->px()); 
      tempShat.push_back(quark->py()+antiquark->py()); 
      tempShat.push_back(quark->pz()+antiquark->pz());
      tempShat.push_back(quark->energy()+antiquark->energy());

      double sH= tempShat[3]*tempShat[3] -tempShat[0]*tempShat[0]-tempShat[1]*tempShat[1]-tempShat[2]*tempShat[2];


      // test using final state
      std::vector<double> test_lepton_Shat;
      test_lepton_Shat.push_back(leptonPlus->px()+leptonMinus->px());
      test_lepton_Shat.push_back(leptonPlus->py()+leptonMinus->py());
      test_lepton_Shat.push_back(leptonPlus->pz()+leptonMinus->pz());
      test_lepton_Shat.push_back(leptonPlus->energy()+leptonMinus->energy());

      double test_finalState_sHat=test_lepton_Shat[3]*test_lepton_Shat[3] -test_lepton_Shat[0]*test_lepton_Shat[0]-test_lepton_Shat[1]*test_lepton_Shat[1]-test_lepton_Shat[2]*tempShat[2];

      std::vector<double> tempThat; 
      tempThat.push_back(quark->px()-leptonPlus->px());
      tempThat.push_back(quark->py()-leptonPlus->py()); 
      tempThat.push_back(quark->pz()-leptonPlus->pz());   
      tempThat.push_back(quark->energy()-leptonPlus->energy()); 
      double tH= (tempThat[3]*tempThat[3] -tempThat[0]*tempThat[0]-tempThat[1]*tempThat[1]-tempThat[2]*tempThat[2]);
      double tH2=tH*tH;
      //double tH2=(tempThat[3]*tempThat[3] -tempThat[0]*tempThat[0]-tempThat[1]*tempThat[1]-tempThat[2]*tempThat[2])*(tempThat[3]*tempThat[3] -tempThat[0]*tempThat[0]-tempThat[1]*tempThat[1]-tempThat[2]*tempThat[2]);


      std::vector<double> tempUhat;
      tempUhat.push_back(quark->px()-leptonMinus->px()); 
      tempUhat.push_back(quark->py()-leptonMinus->py());
      tempUhat.push_back(quark->pz()-leptonMinus->pz());
      tempUhat.push_back(quark->energy()-leptonMinus->energy());
     double uH= (tempUhat[3]*tempUhat[3] -tempUhat[0]*tempUhat[0]-tempUhat[1]*tempUhat[1]-tempUhat[2]*tempUhat[2]); 
     double uH2=uH*uH;
     // double uH2= (tempUhat[3]*tempUhat[3] -tempUhat[0]*tempUhat[0]-tempUhat[1]*tempUhat[1]-tempUhat[2]*tempUhat[2])*(tempUhat[3]*tempUhat[3] -tempUhat[0]*tempUhat[0]-tempUhat[1]*tempUhat[1]-tempUhat[2]*tempUhat[2]); 

      std::cout << "########################   MANDELSTON VARIABLES  ##########################################################################" << std::endl;
      std::cout << "sHat I guess:" << sH << std::endl;
      std::cout << "tHat square I guess:" << tH2 << std::endl;
      std::cout << "uHat square I guess:" << uH2 << std::endl;

      std::cout << "tHat I guess:" << tH << std::endl;
      std::cout << "uHat  I guess:" <<uH << std::endl;
      std::cout << " test_finalState_sHat I guess:" << test_finalState_sHat << std::endl;


 std::cout << "############### IMPLEMENT STEVE'S CODE HERE ################" << std::endl;

 
  int quarkId=quark->pdgId();  
 int idAbs=quarkId;
 // int idAbs=5;
 int leptonMinusId=leptonMinus->pdgId(); 
 int idNew=leptonMinusId;


 std::cout << "#########################################################" << std::endl;
 std::cout << "############### TEST AFTER A FOR LOOP:  ################" << std::endl;
 std::cout << "#########################################################" << std::endl;
 std::cout << "############### Test_etaLR: defalut value for etaLR is 0.00  #############################"<< std::endl;
 std::cout <<"Test_etaLR after for loop : " << qCetaLR << std::endl;
 std::cout << "############################## Test_alpEM #############################################"<< std::endl;
 std::cout << "Test_alpEM  after for loop: " << alpEM << std::endl;
 std::cout << "############### Test_lambdaSquare after for loop : #############################"<< std::endl;
 std::cout <<"Test_lambdaSquare: " <<qCLambda2 << std::endl;
 
 


 
 std::cout << "############### TESST HEPMC TO GET FINE STRUCTURE CONSTANT  ################" << std::endl;
// make a pythia event
 //fMasterGen->next();
 // HepMC::Pythia8ToHepMC ToHepMC;
 //HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
 // ToHepMC.fill_next_event(*fMasterGen, hepmcevt );
 //  double alpEM_HEP = hepmcevt->alphaQED();
 // std::cout << "###############  alpEM_HEP ################"<< alpEM_HEP << std::endl;


  
 std::cout << "## ARGUAMENTS (idAbs AND idNew) DO NOT TAKE NEGATIVE VALUE:LETS SEE WHAT HAPPENS !!! ### " << std::endl;
 std::cout << "#fMasterGen->couplingsPtr->vf(-2)&fMasterGen->couplingsPtr->vf(-13)#" <<fMasterGen->couplingsPtr->vf(-2) << " " << fMasterGen->couplingsPtr->vf(-13)<< std::endl;
 std::cout << "#fMasterGen->couplingsPtr->vf(-3)&fMasterGen->couplingsPtr->vf(-11)#" <<fMasterGen->couplingsPtr->vf(-3) << " " << fMasterGen->couplingsPtr->vf(-11)<< std::endl;
 std::cout << "#fMasterGen->couplingsPtr->af(-2)&fMasterGen->couplingsPtr->af(-13)#" <<fMasterGen->couplingsPtr->af(-2) << " " << fMasterGen->couplingsPtr->af(-13) << std::endl;
 std::cout << "#fMasterGen->couplingsPtr->af(-1)&fMasterGen->couplingsPtr->af(-11)#" <<fMasterGen->couplingsPtr->af(-1) << " " << fMasterGen->couplingsPtr->af(-11) << std::endl;
 std::cout << "## fMasterGen->couplingsPtr->vf(13) ### " << fMasterGen->couplingsPtr->vf(13)<< std::endl;
 std::cout << "## fMasterGen->couplingsPtr->af(13) ### " <<fMasterGen->couplingsPtr->af(13) << std::endl;
 std::cout << "## fMasterGen->couplingsPtr->vf(-2) ### " << fMasterGen->couplingsPtr->vf(-2)<< std::endl;
 std::cout << "## fMasterGen->couplingsPtr->af(-2) ### " <<fMasterGen->couplingsPtr->af(-2) << std::endl;
 std::cout << "## fMasterGen->couplingsPtr->vf(2) ### " << fMasterGen->couplingsPtr->vf(2)<< std::endl;
 std::cout << "## fMasterGen->couplingsPtr->af(2) ### " <<fMasterGen->couplingsPtr->af(2) << std::endl;
 //Process name
 // Incoming quarks
 double tmPgvf = 0.25 * fMasterGen->couplingsPtr->vf(idAbs);
 double tmPgaf = 0.25 * fMasterGen->couplingsPtr->af(idAbs);
 //Outgoing fermions
 double tmPgvl = 0.25 * fMasterGen->couplingsPtr->vf(idNew);
 double tmPgal = 0.25 * fMasterGen->couplingsPtr->af(idNew);

 
 double tmPgLf = tmPgvf + tmPgaf;
 double tmPgLl = tmPgvl + tmPgal;
 double tmPgRf = tmPgvf - tmPgaf;
 double tmPgRl = tmPgvl - tmPgal;



 std::cout <<" ############## process name :" << " tmPgLl: " << tmPgLl <<"tmPgRl: " << tmPgRl <<std::endl;
 std::cout << tmPgLl << " " << tmPgRl << std::endl;


 

 // Kinematics
 //double qCmNew  = fMasterGen->particleData.m0(idNew);
 //double qCmNew2 = qCmNew * qCmNew;
 double qCmZ    = fMasterGen->particleData.m0(23);
 double qCmZ2   = qCmZ * qCmZ;
 double qCGZ    = fMasterGen->particleData.mWidth(23);
 double qCGZ2   = qCGZ * qCGZ;


 // Necessary variables to ampitude
 // First term
 // double alpEM =1./137;
 double tmPe2QfQl = 4. * M_PI * alpEM * fMasterGen->couplingsPtr->ef(idAbs) * fMasterGen->couplingsPtr->ef(idNew);
 double qCPropGm   = 1./sH;

 // double alpha= fMasterGen->info.alphaEM();
 std::cout <<"DIRECTLY ACCESSIBLE VARIABLES "<< std::endl;
 std::cout <<" M_P ouputs Value of pi : " <<  M_PI <<std::endl;
 std::cout <<"fMasterGen->couplingsPtr->ef(idAbs) outputs Charge of incoming quark: " << fMasterGen->couplingsPtr->ef(idAbs) <<std::endl;
 std::cout <<" fMasterGen->couplingsPtr->ef(idNew) outputs Charge of outgoing leopton: " << fMasterGen->couplingsPtr->ef(idNew) <<std::endl;
 std::cout <<"alpEM CAN NOT BE ACCESSED DIRECTLY "<<std::endl;
 // std::cout <<" alphaQED() pythia->info.alphaEM();. "<< alpha <<std::endl;
 // std::cout <<"AlphaQED "<< AlphaQED <<std::endl;

 //Second term.Model depended variables are defined using incoming quark and outgoing fermion information
 double tmPe2s2c2 = 4. * M_PI * alpEM;
 double denomPropZ = pow((sH - qCmZ2), 2) + qCmZ2 * qCGZ2;
 double qCrePropZ  = (sH - qCmZ2) / denomPropZ;
 double qCimPropZ  = -qCmZ * qCGZ / denomPropZ;

 //Third term:4. * M_PI * qCetaLR / qCLambda2;

 

 // Amplitudes, M = gamma + Z + CI.                    
 meLL = tmPe2QfQl * qCPropGm
  + tmPe2s2c2 * tmPgLf * tmPgLl * (qCrePropZ + I * qCimPropZ);
 meRR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgRl * (qCrePropZ + I * qCimPropZ);
 meLR = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ)
        + 4. * M_PI * qCetaLR / qCLambda2;
 meRL = tmPe2QfQl * qCPropGm
        + tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ)
        + 4. * M_PI * qCetaLR / qCLambda2;

 // According to Steve's idea, add SM matrix elements for sigmaNew.
 // Define standard model matrix elements of RL and LR model
 
meLR_SM = tmPe2QfQl * qCPropGm
  + tmPe2s2c2 * tmPgLf * tmPgRl * (qCrePropZ + I * qCimPropZ);
        
 meRL_SM = tmPe2QfQl * qCPropGm
  + tmPe2s2c2 * tmPgRf * tmPgLl * (qCrePropZ + I * qCimPropZ);
  
 // Calculate weighting facror
 double sigma0 = 1.0;  
 double sigmaOld = sigma0 * uH2 * std::real(meLL*std::conj(meLL));
 sigmaOld += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
 sigmaOld += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
 sigmaOld += sigma0 * tH2 * std::real(meRL*std::conj(meRL));

 double sigmaNewLR = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
 sigmaNewLR += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
 sigmaNewLR += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
 // sigma += sigma0 * tH2 * std::real(meRL*std::conj(meRL));
 sigmaNewLR += sigma0 * tH2 * std::real(meRL_SM *std::conj(meRL_SM));
 double fracLR = sigmaNewLR / sigmaOld;

 double sigmaNewRL = sigma0 * uH2 *std:: real(meLL*std::conj(meLL));
 sigmaNewRL += sigma0 * uH2 * std::real(meRR*std::conj(meRR));
 //sigmaNew += sigma0 * tH2 * std::real(meLR*std::conj(meLR));
 sigmaNewRL += sigma0 * tH2 * std::real(meRL*std::conj(meRL));  
 sigmaNewRL += sigma0 * tH2 * std::real(meLR_SM*std::conj(meLR_SM));
 double fracRL = sigmaNewRL / sigmaOld;
 // double fracSum += frac;
 //double weight *= frac;

 double sumRL_Plus_LR =fracLR+fracRL;
 std::cout <<"############ FRACTION LR and RL #############" << std::endl;
 std::cout << "fracLR:  " << fracLR << "  " << "fracR:  " <<fracRL << std::endl;
 std::cout << "sumRL_Plus_LR:  " << fracLR+fracRL << std::endl;
 std::cout <<"############ WEIGHT  #############" << std::endl;
 // std::cout << weight << std::endl;



  



  
    if(debug_ > 2){
      std::cout << "Eta of daughter1 is: " << daughter1->eta() << "\n";
      std::cout << "Eta of daughter2 is: " << daughter2->eta() << "\n";
    }

      if(boson){
	bosonId_=boson->pdgId();
	bosonP4_.fill(boson->p4());
	
	h_Zmass->Fill(boson->mass(), fracLR);
	h_Zpt->Fill(boson->pt(), fracLR);
	h_Zeta->Fill(boson->eta(),fracLR);
	h_Zphi->Fill(boson->phi(), fracLR);
	h_Zcharge->Fill(boson->charge(), fracLR);
      }

    if(daughter1->charge() > 0 && daughter2->charge() < 0){
      muMinus = daughter2;
      muPlus = daughter1;
    }
    else if(daughter1->charge() < 0 && daughter2->charge() > 0){
      muMinus = daughter1;
      muPlus = daughter2;
    }

    else return;

    if(debug_ > 0){  
      std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
      std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
	muPlus->phi() << "\tq = " << muPlus->charge();
    }
    
  muMinusP4_.fill(muMinus->p4());
    muMinusPID_=muMinus->pdgId();
    if(debug_ > 0){  
      std::cout<< "\n\nDaughter1: pId = " << muMinus->pdgId() << "\tpT = " << muMinus->pt() << "\teta = " 
	       << muMinus->eta() << "\tphi = " << muMinus->phi() << "\tq = " << muMinus->charge();
    std::cout<< "\nDaughter2: pId = " << muPlus->pdgId() << "\tpT = " << muPlus->pt() << "\teta = " << muPlus->eta() << "\tphi = " << 
      muPlus->phi() << "\tq = " << muPlus->charge();
    }


    h_muMinusmass->Fill(muMinus->mass(),fracLR);
    h_muMinuspt->Fill(muMinus->pt(), fracLR);
    h_muMinuseta->Fill(muMinus->eta(), fracLR);
    h_muMinusphi->Fill(muMinus->phi(), fracLR);
    h_muMinuscharge->Fill(muMinus->charge(), fracLR);
    h_thetaMuMinus->Fill(muMinus->theta(), fracLR);  
  
    muPlusP4_.fill(muPlus->p4());
    muPlusPID_=muPlus->pdgId();

    h_muPlusmass->Fill(muPlus->mass(), fracLR);
    h_muPluspt->Fill(muPlus->pt(), fracLR);
    h_muPluseta->Fill(muPlus->eta(), fracLR);
    h_muPlusphi->Fill(muPlus->phi(), fracLR);
    h_muPluscharge->Fill(muPlus->charge(), fracLR);
    h_thetaMuPlus->Fill(muPlus->theta(), fracLR);    

    muPlusKPlus = (1/sqrt(2))*(muPlus->energy() + muPlus->pz());
    muPlusKMinus = (1/sqrt(2))*(muPlus->energy() - muPlus->pz());
    muMinusKPlus = (1/sqrt(2))*(muMinus->energy() + muMinus->pz());
    muMinusKMinus = (1/sqrt(2))*(muMinus->energy() - muMinus->pz());

    invariantK = (muPlusKPlus*muMinusKMinus - muMinusKPlus*muPlusKMinus);
    std::cout << "\n\nInvariantK is: " << invariantK << std::endl;

    dimuon = muMinus->p4() + muPlus->p4();

    dimuonPt =dimuon.pt();
    dimuonPz = dimuon.pz();
    pseudorapidity = asinh(dimuonPz/dimuonPt);
    dimuonPx = dimuon.px();
    dimuonPy = dimuon.py();
    Phi = acos(dimuonPx/dimuonPt);
    dimuonQ = sqrt(pow(dimuon.energy(),2) - pow(dimuon.pt(),2) - pow(dimuon.pz(),2));
    std::cout << "\n\nDimuon Energy is: " << dimuon.energy() << std::endl << std::endl;
    
    double denominatorTheta, denominatorPhi1, denominatorPhi2, numeratorPhi1, numeratorPhi2;
    double denominatorPhi, numeratorPhi;
    double deltaX, deltaY;
    double invariantMass;

    denominatorTheta = dimuonQ*sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
    thetaCos = (dimuon.pz()/fabs(dimuon.pz()))*(2/denominatorTheta)*invariantK;
    thetaCS = acos(thetaCos);

    denominatorPhi1 = dimuonQ*dimuon.pt();
    numeratorPhi1 = sqrt(pow(dimuonQ, 2) + pow(dimuon.pt(), 2));
    deltaX = muPlus->px() - muMinus->px();
    deltaY = muPlus->py() - muMinus->py();
    denominatorPhi2 = ((deltaX*dimuon.px()) + (deltaY*dimuon.py()));
    numeratorPhi2 = ((deltaX*dimuon.py()) + (deltaY*dimuon.px()));
    numeratorPhi = numeratorPhi1*numeratorPhi2;
    denominatorPhi = denominatorPhi1 * denominatorPhi2;

    phiTan = numeratorPhi/denominatorPhi;
    phiCS = atan(phiTan);


    mu1Energy = muPlus->energy();
    mu2Energy = muMinus->energy();
    std::cout << "\n\nmuon Energies are: " << mu1Energy << "__" << mu2Energy << std::endl << std::endl;
    std::cout << "\ndimuon px_py_pz are: "<< dimuonPx << "_" << dimuonPy << "_" << dimuonPz << std::endl;
 

    //Exaple TTree weighting
    cosTheta=thetaCos*fracLR;

    cosTheta=thetaCos;
    tanPhi=phiTan;
    csTheta=thetaCS;
    csPhi=phiCS;


    h_cosTheta->Fill(thetaCos);
    h_csTheta->Fill(thetaCS);
    h_tanPhi->Fill(phiTan);
    h_csPhi->Fill(phiCS);

    

    std::cout << "\n\n\ncos(Theta_CS) = " << thetaCos << "\tThetaCS = " << thetaCS << std::endl;
    std::cout << "\n\n\nTan(phi_CS) = " << phiTan << "\tPhiCS = " << phiCS << std::endl;

    invariantMass = sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi()))));



    

    //Quark selection: we do not need it, just sanity check
    // Delete it after using it, it is just a garbage

    int dquark,uquark,squark,cquark,bquark,tquark;

    if(idAbs==1){
      dquark=idAbs;
    }
    else if(idAbs==2){
      uquark=idAbs;
    }
    else if(idAbs==3){
      squark=idAbs;
    }
    else if(idAbs==4){
      cquark=idAbs;
      
    }
    else if(idAbs==5){
      bquark=idAbs;
    } 
    else {
      
      tquark=idAbs;
    }


    // Delete it after checking

    DQuark=dquark;
    UQuark=uquark;
    SQuark=squark;
    CQuark=cquark;
    BQuark=bquark;
    TQuark=tquark;

    h_DQuark->Fill(dquark);
    h_UQuark->Fill(uquark);
    h_SQuark->Fill(squark);
    h_CQuark->Fill(cquark);
    h_BQuark->Fill(bquark);
    h_TQuark->Fill(tquark);

    

    // Fill the hostograms for Mandelston variables

    SHat=sH;
    THat2=tH2;
    UHat2=uH2;
     fractionLR=fracLR;
    fractionRL=fracRL;

    sumFrac=sumRL_Plus_LR;

    T_Hat=tH;
    U_Hat=uH;

    alphaRunQED=alpEM;

    h_SHat->Fill(sH);
    h_THat2->Fill(tH2);
    h_UHat2->Fill(uH2);
    h_fractionLR->Fill(fracLR);
    h_fractionRL->Fill(fracRL);
    h_sumFrac->Fill(sumFrac);

    h_T_Hat->Fill(tH);
    h_U_Hat->Fill(uH);
    h_alphaRunQED->Fill(alpEM); 

    if(thetaCos < 0.0){
      h_cosThetaMinusInvariantMass->Fill(invariantMass);
      mCosThetaMinus = invariantMass;
    }
    else{
      h_cosThetaPlusInvariantMass->Fill(invariantMass);
      mCosThetaPlus = invariantMass;
    }


    h_dphi->Fill(TVector2::Phi_mpi_pi(muMinus->phi()- muPlus->phi()));
    h_dtheta->Fill(TVector2::Phi_mpi_pi(muMinus->theta()- muPlus->theta()));
    h_dr->Fill(reco::deltaR(muMinus->p4(),muPlus->p4()));
    h_massInvar->Fill(sqrt(2 * daughter1->pt() * daughter2->pt() *( cosh(daughter1->eta() - daughter2->eta()) - cos(TVector2::Phi_mpi_pi(daughter1->phi() - daughter2->phi())))));
    h_dimuonPt->Fill(dimuonPt);
    h_dimuonEta->Fill(pseudorapidity);
    h_dimuonPhi->Fill(Phi);
    h2_phi1_vs_phi2->Fill(muMinus->phi(),muPlus->phi());  
    h2_eta1_vs_eta2->Fill(muMinus->eta(),muPlus->eta());
    h2_pt1_vs_pt2->Fill(muMinus->pt(),muPlus->pt());


    //  }

  std::cout << "\n\n===========================================================================================================" << std::endl;
  tree_->Fill();  
}

bool Dimuon::isBoson(int pid)
{
  if(pid==23 || abs(pid)==22 || pid==32){
    if(debug_ > 0) std::cout << "\n\nFound Boson\n";
    return true;
  }
  else return false;
}

bool Dimuon::isMuon(int pid){
  if(abs(pid)==11 || abs(pid) ==13){
    if(debug_ > 0) std::cout << "\n\nFound A Muon!\n";
    return true;
  }
  else return false;
}

bool Dimuon::checkBosonStatus( const reco::GenParticleCollection& genParts){
  const reco::Candidate* boson = getBoson(genParts);
  if(boson == nullptr){
    if(debug_ > 0) std::cout << "\nBoson is: "  << boson;
    return false;
  }

  else if( boson->status() != 22){
    if(debug_ > 0)  std::cout <<"\nBoson Status is: "<< boson->status();
    return false;
  }
 
    return true;
 }

const reco::Candidate* Dimuon::getBoson( const reco::GenParticleCollection& genParts)
{
  for(auto &part : genParts){
    if(isBoson(part.pdgId())){
      if(debug_ > 1){
      std::cout << "\npId is: " << part.pdgId();
      std::cout << "\nStatus is: " << part.status();
      }
      return getLastDaughter(&part,part.pdgId());
    }
  }
  return nullptr;
}


const reco::Candidate* Dimuon::getLastDaughter(const reco::Candidate* part,int pid)
{
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return getLastDaughter(part->daughter(partNr),pid);
  }
  return part;
}
       
const reco::Candidate* Dimuon::getDaughter(const reco::Candidate* part,int pid)
{  
  for(size_t partNr =0; part && partNr<part->numberOfDaughters();partNr++){
    if(part->daughter(partNr)->pdgId()==pid) return part->daughter(partNr);
  }
  return nullptr;
}

 const reco::Candidate* Dimuon::getMother(const reco::Candidate* part, int pid)
{
  for(size_t partNr = 0; part && partNr < part->numberOfMothers(); partNr++){
    if(part->mother(partNr)->pdgId() == pid) return getMother(part->mother(partNr),pid);
  
    else if(abs(part->mother(partNr)->pdgId()) == 1 || abs(part->mother(partNr)->pdgId()) == 2 ||
	    abs(part->mother(partNr)->pdgId()) == 3 || abs(part->mother(partNr)->pdgId()) == 4 ||
		 abs(part->mother(partNr)->pdgId()) == 5 || abs(part->mother(partNr)->pdgId()) == 6 ||
	    abs(part->mother(partNr)->pdgId()) == 7 || abs(part->mother(partNr)->pdgId()) == 8 ||
	    abs(part->mother(partNr)->pdgId()) == 23 || abs(part->mother(partNr)->pdgId()) == 32  || 
	    abs(part->mother(partNr)->pdgId()) == 22) return part->mother(partNr);
  }  
  return nullptr;
  
}


// ------------ method called once each job just before starting event loop  ------------

// ------------ method called once each job just after ending the event loop  ------------
void 
Dimuon::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
  
/*void 
Dimuon::beginRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
{
   edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByToken(genInfoProductToken_, genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();
  //  tree_->Fill();
  std::cout<< "Cross Section is: "  << crossSec << std::endl;  
 

}
*/

// ------------ method called when ending the processing of a run  ------------

void 
  Dimuon::endRun(edm::Run const& iRun, edm::EventSetup const& iEventSetup)
  {
   edm::Handle< GenRunInfoProduct > genInfoProduct;
  iRun.getByToken(genInfoProductToken_, genInfoProduct );
  crossSec = genInfoProduct->internalXSec().value();
  std::cout<< "Cross Section is: "  << crossSec << std::endl;  
 
  }

  
// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  Dimuon::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  Dimuon::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Dimuon::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Dimuon);
