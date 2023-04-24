//framework headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

//analysis headers
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"

//ROOT headers
#include <TTree.h>
 
//STL headers 
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
using std::vector;

typedef math::PtEtaPhiELorentzVector LorentzVector;
typedef math::XYZTLorentzVector LorentzVector2;

class DijetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit DijetAnalyzer(const edm::ParameterSet&);
		~DijetAnalyzer() {}
	
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
		struct GenNtuple {
			LorentzVector WideJetNoNu1;
			LorentzVector WideJetNoNu2;
			double MJJNoNu;
			LorentzVector WideJet1;
			LorentzVector WideJet2;
			double MJJ;
		};
	
	private:
		void beginJob() override;
		void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
		void analyze(const edm::Event&, const edm::EventSetup&) override;
		void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
		void endJob() override {}

		//implementation
		bool algorithm(const vector<reco::GenJet>& jets, LorentzVector& wideJet1, LorentzVector& wideJet2, double& mJJ);

		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		TTree* tree;
		//for tree branches
		GenNtuple entry;
		//cuts
		double pt;
		double eta;
		unsigned njet;
		double radius;
		double maxdeta;

		//tokens
		edm::EDGetTokenT<vector<reco::GenJet>> tok_jet, tok_jet_nonu;
};

//
// constructors and destructor
//
DijetAnalyzer::DijetAnalyzer(const edm::ParameterSet& iConfig) :
	tree(nullptr),
	pt(iConfig.getParameter<double>("pt")),
	eta(iConfig.getParameter<double>("eta")),
	njet(iConfig.getParameter<unsigned>("njet")),
	radius(iConfig.getParameter<double>("radius")),
	maxdeta(iConfig.getParameter<double>("maxdeta")),
	tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
	tok_jet_nonu(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetNoNuTag")))
{
	usesResource("TFileService");
}

//
// member functions
//

void DijetAnalyzer::beginJob()
{
	tree = fs->make<TTree>("tree","tree");
	
	tree->Branch("WideJet1", "WideJet1", &entry.WideJet1, 32000, 99);
	tree->Branch("WideJet2", "WideJet2", &entry.WideJet2, 32000, 99);
	tree->Branch("MJJ", &entry.MJJ, "MJJ/D");
	tree->Branch("WideJetNoNu1", "WideJetNoNu1", &entry.WideJetNoNu1, 32000, 99);
	tree->Branch("WideJetNoNu2", "WideJetNoNu2", &entry.WideJetNoNu2, 32000, 99);
	tree->Branch("MJJNoNu", &entry.MJJNoNu, "MJJ/D");
}

bool DijetAnalyzer::algorithm(const vector<reco::GenJet>& jets, LorentzVector& wideJet1, LorentzVector& wideJet2, double& mJJ)
{
	vector<LorentzVector> seedJets;

	if(jets.size()<njet) return false;
	for(unsigned j = 0; j < njet; ++j){
		const auto& jet = jets[j];
		if(jet.pt() <= pt) return false;
		if(abs(jet.eta()) >= eta) return false;
		LorentzVector vjet;
		vjet = jet.p4();
		seedJets.push_back(vjet);
	}

	vector<LorentzVector> wideJets = seedJets;
	for(unsigned j = njet; j < jets.size(); ++j){
		const auto& jet = jets[j];
		if(jet.pt() <= pt) continue;
		if(abs(jet.eta()) >= eta) continue;
		double minDR = 1e10;
		int minDRindex = -1;
		for(unsigned s = 0; s < seedJets.size(); ++s){
			double deltaR = reco::deltaR(jet.p4(),seedJets[s]);
			if(deltaR < minDR){
				minDR = deltaR;
				minDRindex = s;
			}
		}
		if(minDR < radius) wideJets[minDRindex] += jet.p4();
	}

	double deta = abs(wideJets.at(0).eta()-wideJets.at(1).eta());
	if(deta>=maxdeta) return false;

	wideJet1 = wideJets.at(0);
	wideJet2 = wideJets.at(1);
	mJJ = (wideJet1+wideJet2).mass();
	return true;
}

// ------------ method called on each new Event  ------------
void DijetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	entry = GenNtuple();
	entry.MJJ = 0;
	entry.MJJNoNu = 0;

	edm::Handle<vector<reco::GenJet>> h_jet;
	iEvent.getByToken(tok_jet,h_jet);
	bool pass = algorithm(*(h_jet.product()), entry.WideJet1, entry.WideJet2, entry.MJJ);

	edm::Handle<vector<reco::GenJet>> h_jet_nonu;
	iEvent.getByToken(tok_jet_nonu,h_jet_nonu);
	bool pass_nonu = algorithm(*(h_jet_nonu.product()), entry.WideJetNoNu1, entry.WideJetNoNu2, entry.MJJNoNu);

	if(pass or pass_nonu){
		//fill tree
		tree->Fill();
	}
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DijetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;

	desc.add<double>("pt",30.);
	desc.add<double>("eta",2.5);
	desc.add<unsigned>("njet",2);
	desc.add<double>("radius",1.1);
	desc.add<double>("maxdeta",1.1);
	desc.add<edm::InputTag>("JetTag",edm::InputTag("ak4GenJets"));
	desc.add<edm::InputTag>("JetNoNuTag",edm::InputTag("ak4GenJetsNoNu"));

	descriptions.add("DijetAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DijetAnalyzer);
