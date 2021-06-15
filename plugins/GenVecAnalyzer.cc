//framework headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

//analysis headers
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
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
#include <set>
using std::vector;

typedef math::PtEtaPhiELorentzVector LorentzVector;

double TransverseMass(double px1, double py1, double m1, double px2, double py2, double m2){
	double E1 = sqrt(pow(px1,2)+pow(py1,2)+pow(m1,2));
	double E2 = sqrt(pow(px2,2)+pow(py2,2)+pow(m2,2));
	double MTsq = pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2);
	return sqrt(std::max(MTsq,0.0));
}

class GenVecAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
	public:
		explicit GenVecAnalyzer(const edm::ParameterSet&);
		~GenVecAnalyzer() {}
	
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
		
		struct GenNtuple {
			LorentzVector Zprime;
			LorentzVector DarkQuark1;
			LorentzVector DarkQuark2;
			LorentzVector Jet1;
			LorentzVector Jet2;
			LorentzVector Invis1;
			LorentzVector Invis2;
			LorentzVector Met;
			double Mdq;
			double Mmc;
			double Mjj;
			double MT;
		};
	
	private:
		void beginJob() override;
		void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
		void analyze(const edm::Event&, const edm::EventSetup&) override;
		void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
		void endJob() override {}
		
		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		TTree* tree;
		//for tree branches
		GenNtuple entry;

		//tokens
		edm::EDGetTokenT<vector<reco::GenMET>> tok_met;
		edm::EDGetTokenT<vector<reco::GenJet>> tok_jet;
		edm::EDGetTokenT<vector<reco::GenParticle>> tok_part;
};

//
// constructors and destructor
//
GenVecAnalyzer::GenVecAnalyzer(const edm::ParameterSet& iConfig) :
	tree(nullptr),
	tok_met(consumes<vector<reco::GenMET>>(iConfig.getParameter<edm::InputTag>("METTag"))),
	tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
	tok_part(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("PartTag")))
{
	usesResource("TFileService");
}

//
// member functions
//

void GenVecAnalyzer::beginJob()
{
	tree = fs->make<TTree>("tree","tree");
	
	tree->Branch("Zprime", "Zprime", &entry.Zprime, 32000, 99);
	tree->Branch("DarkQuark1", "DarkQuark1", &entry.DarkQuark1, 32000, 99);
	tree->Branch("DarkQuark2", "DarkQuark2", &entry.DarkQuark2, 32000, 99);
	tree->Branch("Jet1", "Jet1", &entry.Jet1, 32000, 99);
	tree->Branch("Jet2", "Jet2", &entry.Jet2, 32000, 99);
	tree->Branch("Invis1", "Invis1", &entry.Invis1, 32000, 99);
	tree->Branch("Invis2", "Invis2", &entry.Invis2, 32000, 99);
	tree->Branch("Met", "Met", &entry.Met, 32000, 99);
	tree->Branch("Mdq", &entry.Mdq, "Mdq/D");
	tree->Branch("Mmc", &entry.Mmc, "Mmc/D");
	tree->Branch("Mjj", &entry.Mjj, "Mjj/D");
	tree->Branch("MT", &entry.MT, "MT/D");
}

// ------------ method called on each new Event  ------------
void GenVecAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	entry = GenNtuple();
	
	edm::Handle<vector<reco::GenMET>> h_met;
	iEvent.getByToken(tok_met,h_met);

	edm::Handle<vector<reco::GenJet>> h_jet;
	iEvent.getByToken(tok_jet,h_jet);

	edm::Handle<vector<reco::GenParticle>> h_part;
	iEvent.getByToken(tok_part,h_part);

	int jet_counter = 0;	
	for(const auto& i_jet : *(h_jet.product())){
		if(jet_counter==0) entry.Jet1 = i_jet.p4();
		else if(jet_counter==1) entry.Jet2 = i_jet.p4();
		else break;
		++jet_counter;
	}

	const double jet_radius = 0.8;
	const std::set<int> stable_ids{51,52,53};
	for(const auto& i_part : *(h_part.product())){
		if(i_part.numberOfDaughters()>0 and stable_ids.find(std::abs(i_part.daughter(0)->pdgId()))!=stable_ids.end()){
			double dr1 = 1e10, dr2 = 1e10;
			if(jet_counter>0) dr1 = reco::deltaR(entry.Jet1,i_part.p4());
			if(jet_counter>1) dr2 = reco::deltaR(entry.Jet2,i_part.p4());
			if(dr1 < dr2 and dr1 < jet_radius) entry.Invis1 += i_part.p4();
			else if(dr2 < dr1 and dr2 < jet_radius) entry.Invis2 += i_part.p4();
		}
	}

	const auto& i_met = h_met->front();
	entry.Met = i_met.p4();

	//get truth: last zprime, first dark quarks
	const int zprime_id = 4900023;
	for (const auto& part_i : *(h_part.product())) {
		if(abs(part_i.pdgId())==zprime_id and part_i.isLastCopy()){
			entry.Zprime = part_i.p4();
			if(part_i.numberOfDaughters()>0) entry.DarkQuark1 = part_i.daughter(0)->p4();
			if(part_i.numberOfDaughters()>1) entry.DarkQuark2 = part_i.daughter(1)->p4();
			break;
		}
	}

	//compute mass variables
	auto vjj = entry.Jet1 + entry.Jet2;
	entry.Mdq = (entry.DarkQuark1 + entry.DarkQuark2).mass();
	entry.Mmc = (vjj + entry.Invis1 + entry.Invis2).mass();
	entry.Mjj = vjj.mass();
	entry.MT = TransverseMass(vjj.px(),vjj.py(),vjj.mass(),entry.Met.px(),entry.Met.py(),0.0);
	
	//fill tree
	tree->Fill();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenVecAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.add<edm::InputTag>("METTag",edm::InputTag("genMetTrue"));
	desc.add<edm::InputTag>("JetTag",edm::InputTag("ak8GenJetsNoNu"));
	desc.add<edm::InputTag>("PartTag",edm::InputTag("genParticles"));
	
	descriptions.add("GenVecAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenVecAnalyzer);
