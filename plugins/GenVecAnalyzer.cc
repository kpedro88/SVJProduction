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
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

//ROOT headers
#include <TTree.h>
 
//STL headers 
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
using std::vector;

//user headers
#include "SVJ/Production/interface/lester_mt2_bisect.h"

typedef math::PtEtaPhiELorentzVector LorentzVector;
typedef math::XYZTLorentzVector LorentzVector2;

namespace parse {
	//generalization for processing a line
	inline void process(const std::string& line, char delim, std::vector<std::string>& fields){
		std::stringstream ss(line);
		std::string field;
		while(getline(ss,field,delim)){
			fields.push_back(field);
		}
	}
}

double TransverseMass(double px1, double py1, double m1, double px2, double py2, double m2){
	double E1 = sqrt(pow(px1,2)+pow(py1,2)+pow(m1,2));
	double E2 = sqrt(pow(px2,2)+pow(py2,2)+pow(m2,2));
	double MTsq = pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2);
	return sqrt(std::max(MTsq,0.0));
}

class GenVecAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchLuminosityBlocks> {
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
			int Ninv;
			int Ninv1;
			int Ninv2;
			int Nvis;
			int Nvis1;
			int Nvis2;
			int Nconst1;
			int Nconst2;
			double Mdq;
			double Mmc;
			double Mjj;
			double MT;
			double MT2;
			double MAOS;
			LorentzVector Met1;
			LorentzVector Met2;
			double mZprime;
			double mDark;
			double rinv;
			double alpha;
		};
	
	private:
		void beginJob() override;
		void doBeginRun_(const edm::Run&, const edm::EventSetup&) override {}
		void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		void analyze(const edm::Event&, const edm::EventSetup&) override;
		void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override {}
		void doEndRun_(const edm::Run&, const edm::EventSetup&) override {}
		void endJob() override {}

		//helper
		void getSVJComment(const GenLumiInfoHeader& gen);
		
		// ----------member data ---------------------------
		edm::Service<TFileService> fs;
		TTree* tree;
		//for tree branches
		GenNtuple entry;

		//tokens
		edm::EDGetTokenT<vector<reco::GenMET>> tok_met;
		edm::EDGetTokenT<vector<reco::GenJet>> tok_jet;
		edm::EDGetTokenT<vector<reco::GenParticle>> tok_part;
		edm::EDGetTokenT<GenLumiInfoHeader> tok_scan;
		std::vector<double> signalParameters_;
};

//
// constructors and destructor
//
GenVecAnalyzer::GenVecAnalyzer(const edm::ParameterSet& iConfig) :
	tree(nullptr),
	tok_met(consumes<vector<reco::GenMET>>(iConfig.getParameter<edm::InputTag>("METTag"))),
	tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
	tok_part(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("PartTag"))),
	tok_scan(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator")))
{
	usesResource("TFileService");
}

//
// member functions
//

void GenVecAnalyzer::beginJob()
{
	asymm_mt2_lester_bisect::disableCopyrightMessage();

	tree = fs->make<TTree>("tree","tree");
	
	tree->Branch("Zprime", "Zprime", &entry.Zprime, 32000, 99);
	tree->Branch("DarkQuark1", "DarkQuark1", &entry.DarkQuark1, 32000, 99);
	tree->Branch("DarkQuark2", "DarkQuark2", &entry.DarkQuark2, 32000, 99);
	tree->Branch("Jet1", "Jet1", &entry.Jet1, 32000, 99);
	tree->Branch("Jet2", "Jet2", &entry.Jet2, 32000, 99);
	tree->Branch("Invis1", "Invis1", &entry.Invis1, 32000, 99);
	tree->Branch("Invis2", "Invis2", &entry.Invis2, 32000, 99);
	tree->Branch("Met", "Met", &entry.Met, 32000, 99);
	tree->Branch("Ninv", &entry.Ninv, "Ninv/I");
	tree->Branch("Ninv1", &entry.Ninv1, "Ninv1/I");
	tree->Branch("Ninv2", &entry.Ninv2, "Ninv2/I");
	tree->Branch("Nvis", &entry.Nvis, "Nvis/I");
	tree->Branch("Nvis1", &entry.Nvis1, "Nvis1/I");
	tree->Branch("Nvis2", &entry.Nvis2, "Nvis2/I");
	tree->Branch("Nconst1", &entry.Nconst1, "Nconst1/I");
	tree->Branch("Nconst2", &entry.Nconst2, "Nconst2/I");
	tree->Branch("Mdq", &entry.Mdq, "Mdq/D");
	tree->Branch("Mmc", &entry.Mmc, "Mmc/D");
	tree->Branch("Mjj", &entry.Mjj, "Mjj/D");
	tree->Branch("MT", &entry.MT, "MT/D");
	tree->Branch("MT2", &entry.MT2, "MT2/D");
	tree->Branch("MAOS", &entry.MAOS, "MAOS/D");
	tree->Branch("Met1", "Met1", &entry.Met1, 32000, 99);
	tree->Branch("Met2", "Met2", &entry.Met2, 32000, 99);
	tree->Branch("mZprime", &entry.mZprime, "mZprime/D");
	tree->Branch("mDark", &entry.mDark, "mDark/D");
	tree->Branch("rinv", &entry.rinv, "rinv/D");
	tree->Branch("alpha", &entry.alpha, "alpha/D");
}

//parse GenLumiInfo for SVJ
//from https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/Utils/src/SignalScanProducer.cc
void GenVecAnalyzer::getSVJComment(const GenLumiInfoHeader& gen){
	signalParameters_.clear();

	const std::map<std::string,double> alpha_vals{
		{"peak",-2.},
		{"high",-1.},
		{"low",-3.},
	};
	const std::string& model = gen.configDescription();
	if(model.empty()) return;
	std::vector<std::string> fields;
	parse::process(model,'_',fields);

	//format: SVJ_s-channel_mZprime-X_mDark-Y_rinv-Z_alpha-W
	for(const auto& f : fields){
		std::vector<std::string> subfields;
		parse::process(f,'-',subfields);
		if(subfields.size()!=2 or subfields[0]=="s") continue;
		double val = 0.;
		if(subfields[0]=="alpha") val = alpha_vals.at(subfields[1]);
		else {
			std::stringstream sval(subfields[1]);
			sval >> val;
		}
		signalParameters_.push_back(val);
	}
}

void GenVecAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
	edm::Handle<GenLumiInfoHeader> gen_header;
	iLumi.getByToken(tok_scan, gen_header);
	getSVJComment(*gen_header);
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
	entry.Nconst1 = 0;
	entry.Nconst2 = 0;
	for(const auto& i_jet : *(h_jet.product())){
		if(jet_counter==0){
			entry.Jet1 = i_jet.p4();
			entry.Nconst1 = i_jet.numberOfDaughters();
		}
		else if(jet_counter==1){
			entry.Jet2 = i_jet.p4();
			entry.Nconst2 = i_jet.numberOfDaughters();
		}
		else break;
		++jet_counter;
	}

	const double jet_radius = 0.8;
	const std::set<int> darkhadron_ids{4900111,4900113,4900211,4900213};
	const std::set<int> stable_ids{51,52,53};
	entry.Ninv = 0;
	entry.Ninv1 = 0;
	entry.Ninv2 = 0;
	entry.Nvis = 0;
	entry.Nvis1 = 0;
	entry.Nvis2 = 0;
	for(const auto& i_part : *(h_part.product())){
		//todo: for unstable dark hadrons, can actually make list of all that produced genjet constituents
		if(darkhadron_ids.find(std::abs(i_part.pdgId()))!=darkhadron_ids.end() and i_part.numberOfDaughters()>0){
			bool invis = stable_ids.find(std::abs(i_part.daughter(0)->pdgId()))!=stable_ids.end();
			if(invis) ++entry.Ninv;
			else ++entry.Nvis;
			double dr1 = 1e10, dr2 = 1e10;
			if(jet_counter>0) dr1 = reco::deltaR(entry.Jet1,i_part.p4());
			if(jet_counter>1) dr2 = reco::deltaR(entry.Jet2,i_part.p4());
			if(dr1 < dr2 and dr1 < jet_radius){
				if(invis){
					entry.Invis1 += i_part.p4();
					++entry.Ninv1;
				}
				else ++entry.Nvis1;
			}
			else if(dr2 < dr1 and dr2 < jet_radius){
				if(invis){
					entry.Invis2 += i_part.p4();
					++entry.Ninv2;
				}
				else ++entry.Nvis2;
			}
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

	//compute MT2 and MAOS
	double METx = i_met.px();
	double METy = i_met.py();
	entry.MT2 = asymm_mt2_lester_bisect::get_mT2(
		entry.Jet1.mass(), entry.Jet1.px(), entry.Jet1.py(),
		entry.Jet2.mass(), entry.Jet2.px(), entry.Jet2.py(),
		METx, METy, 0.0, 0.0, 0
	);
	//get invisible systems from MT2
	auto MET1 = asymm_mt2_lester_bisect::ben_findsols(entry.MT2,
		entry.Jet1.px(), entry.Jet1.py(), entry.Jet1.mass(), 0.0,
		entry.Jet2.px(), entry.Jet2.py(),
		METx, METy, entry.Jet2.mass(), 0.0
	);
	double MET1x = MET1.first;
	double MET1y = MET1.second;
	double MET1t = std::sqrt(std::pow(MET1x,2)+std::pow(MET1y,2));
	double MET2x = METx - MET1x;
	double MET2y = METy - MET1y;
	double MET2t = std::sqrt(std::pow(MET2x,2)+std::pow(MET2y,2));

	//use MAOS scheme 2 ("modified") to estimate longitudinal momenta of invisible systems
	double MET1z = MET1t*entry.Jet1.Pz()/entry.Jet1.Pt();
	double MET2z = MET2t*entry.Jet2.Pz()/entry.Jet2.Pt();

	//MAOS scheme 2 assumes met is massless: E = p
	LorentzVector2 tmp1(MET1x,MET1y,MET1z,std::sqrt(MET1t*MET1t+MET1z*MET1z));
	entry.Met1 = LorentzVector(tmp1.pt(),tmp1.eta(),tmp1.phi(),tmp1.energy());
	LorentzVector2 tmp2(MET2x,MET2y,MET2z,std::sqrt(MET2t*MET2t+MET2z*MET2z));
	entry.Met2 = LorentzVector(tmp2.pt(),tmp2.eta(),tmp2.phi(),tmp2.energy());

	//construct invariant mass of parent
	entry.MAOS = (entry.Jet1+entry.Jet2+entry.Met1+entry.Met2).mass();

	//check signal scan info
	if(!signalParameters_.empty()){
		entry.mZprime = signalParameters_[0];
		entry.mDark = signalParameters_[1];
		entry.rinv = signalParameters_[2];
		entry.alpha = signalParameters_[3];
	}
	else {
		entry.mZprime = 0.;
		entry.mDark = 0.;
		entry.rinv = 0.;
		entry.alpha = 0.;
	}
	
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
