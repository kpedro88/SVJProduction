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
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//ROOT headers
#include <TTree.h>
 
//STL headers 
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include <map>
#include <string>
#include <sstream>
#include <iomanip>
using std::vector;

//user headers
#include "SVJ/Production/interface/lester_mt2_bisect.h"

typedef math::PtEtaPhiELorentzVector LorentzVector;
typedef math::XYZTLorentzVector LorentzVector2;

typedef std::unordered_set<unsigned> PidSet;
typedef const reco::Candidate* CandPtr;
typedef std::unordered_set<CandPtr> CandSet;

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
			LorentzVector Mediator1;
			LorentzVector Mediator2;
			LorentzVector SMQuark1;
			LorentzVector SMQuark2;
			LorentzVector Jet1;
			LorentzVector Jet2;
			LorentzVector Jet3;
			LorentzVector Jet4;
			LorentzVector Invis;
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
			double MTx;
			double MT2;
			double MAOS;
			double Meff;
			double Girth1;
			double Girth2;
			double Girth3;
			LorentzVector Met1;
			LorentzVector Met2;
			double PairMT2;
			vector<int> PairMT2ID;
			double PairMT2Reco;
			vector<int> PairMT2RecoID;
			int nMediator;
			double mMediator;
			double mDark;
			double rinv;
			double alpha;
			double yukawa;
			double Weight;
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
		void getSVJComment(const std::string& model);

		//helpers for t-channel ID/categorization
		//extracted from: https://github.com/TreeMaker/TreeMaker/blob/Run2_UL/Utils/src/HiddenSectorProducer.cc
		bool isParticle(const PidSet& darkList, const reco::GenParticle& part) const;
		bool isParticle(const PidSet& darkList, CandPtr part) const;
		bool isParticle(const PidSet& darkList, int pid) const;
		void firstDark(CandPtr part, CandSet& firstMd, CandPtr& firstQdM1, CandPtr& firstQdM2, CandPtr& firstQsM1, CandPtr& firstQsM2, bool& secondDM, bool& secondSM);
		void medDecay(CandPtr part, CandPtr& firstQdM1, CandPtr& firstQdM2, CandPtr& firstQsM1, CandPtr& firstQsM2, bool& secondDM, bool& secondSM);
		std::vector<int> matchJetPart(const vector<reco::GenJet>& jets, const vector<CandPtr>& parts) const;
		double girth(const reco::GenJet& jet) const;

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
		edm::EDGetTokenT<GenEventInfoProduct> tok_weight;
		std::map<std::string,double> signalParameters_;
		PidSet DarkSMediatorIDs_, DarkTMediatorIDs_, DarkQuarkIDs_, DarkHadronIDs_, DarkGluonIDs_, DarkStableIDs_, DarkFirstIDs_, SMQuarkIDs_;
		double coneSize_{0.8};
		bool signal;
};

//
// constructors and destructor
//
GenVecAnalyzer::GenVecAnalyzer(const edm::ParameterSet& iConfig) :
	tree(nullptr),
	tok_met(consumes<vector<reco::GenMET>>(iConfig.getParameter<edm::InputTag>("METTag"))),
	tok_jet(consumes<vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("JetTag"))),
	tok_part(consumes<vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("PartTag"))),
	tok_scan(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),
	tok_weight(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
	DarkSMediatorIDs_{4900023},
	DarkTMediatorIDs_{4900001,4900002,4900003,4900004,4900005,4900006},
	DarkQuarkIDs_{4900101,4900102},
	DarkHadronIDs_{4900111,4900113,4900211,4900213},
	DarkGluonIDs_{4900021},
	DarkStableIDs_{51,52,53},
	SMQuarkIDs_{1,2,3,4,5,6},
    signal(iConfig.getParameter<bool>("signal"))
{
	usesResource("TFileService");
	const auto& model = iConfig.getParameter<std::string>("model");
	if(!model.empty()) getSVJComment(model);

	//combined list of possible first particles
	DarkFirstIDs_.insert(DarkTMediatorIDs_.begin(),DarkTMediatorIDs_.end());
	DarkFirstIDs_.insert(DarkQuarkIDs_.begin(),DarkQuarkIDs_.end());
	DarkFirstIDs_.insert(DarkGluonIDs_.begin(),DarkGluonIDs_.end());
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
	tree->Branch("Mediator1", "Mediator1", &entry.Mediator1, 32000, 99);
	tree->Branch("Mediator2", "Mediator2", &entry.Mediator2, 32000, 99);
	tree->Branch("SMQuark1", "SMQuark1", &entry.SMQuark1, 32000, 99);
	tree->Branch("SMQuark2", "SMQuark2", &entry.SMQuark2, 32000, 99);
	tree->Branch("Jet1", "Jet1", &entry.Jet1, 32000, 99);
	tree->Branch("Jet2", "Jet2", &entry.Jet2, 32000, 99);
	tree->Branch("Jet3", "Jet3", &entry.Jet3, 32000, 99);
	tree->Branch("Jet4", "Jet4", &entry.Jet4, 32000, 99);
	tree->Branch("Invis", "Invis", &entry.Invis, 32000, 99);
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
	tree->Branch("MTx", &entry.MTx, "MTx/D");
	tree->Branch("MT2", &entry.MT2, "MT2/D");
	tree->Branch("MAOS", &entry.MAOS, "MAOS/D");
	tree->Branch("Meff", &entry.Meff, "Meff/D");
	tree->Branch("Girth1", &entry.Girth1, "Girth1/D");
	tree->Branch("Girth2", &entry.Girth2, "Girth2/D");
	tree->Branch("Girth3", &entry.Girth3, "Girth3/D");
	tree->Branch("Met1", "Met1", &entry.Met1, 32000, 99);
	tree->Branch("Met2", "Met2", &entry.Met2, 32000, 99);
	tree->Branch("PairMT2", &entry.PairMT2, "PairMT2/D");
	tree->Branch("PairMT2ID", "PairMT2ID", &entry.PairMT2ID, 32000, 99);
	tree->Branch("PairMT2Reco", &entry.PairMT2Reco, "PairMT2Reco/D");
	tree->Branch("PairMT2RecoID", "PairMT2RecoID", &entry.PairMT2RecoID, 32000, 99);
	tree->Branch("nMediator", &entry.nMediator, "nMediator/I");
	tree->Branch("mMediator", &entry.mMediator, "mMediator/D");
	tree->Branch("mDark", &entry.mDark, "mDark/D");
	tree->Branch("rinv", &entry.rinv, "rinv/D");
	tree->Branch("alpha", &entry.alpha, "alpha/D");
	tree->Branch("yukawa", &entry.yukawa, "yukawa/D");
	tree->Branch("Weight", &entry.Weight, "Weight/D");
}

//parse GenLumiInfo for SVJ
//from https://github.com/TreeMaker/TreeMaker/blob/Run2_2017/Utils/src/SignalScanProducer.cc
void GenVecAnalyzer::getSVJComment(const GenLumiInfoHeader& gen){
	const std::string& model = gen.configDescription();
	if(model.empty()) return;
	getSVJComment(model);
}

void GenVecAnalyzer::getSVJComment(const std::string& model){
	signalParameters_.clear();

	const std::map<std::string,double> alpha_vals{
		{"peak",-2.},
		{"high",-1.},
		{"low",-3.},
	};
	std::vector<std::string> fields;
	parse::process(model,'_',fields);

	//format 1: SVJ_s-channel_mZprime-X_mDark-Y_rinv-Z_alpha-W
	//format 2: SVJ_t-channel_[nMed-N]_mMed-X_mDark-Y_rinv-Z_alpha-W_yukawa-K
	for(const auto& f : fields){
		std::vector<std::string> subfields;
		parse::process(f,'-',subfields);
		if(subfields.size()!=2 or subfields[0].size()==1) continue;
		double val = 0.;
		if(subfields[0]=="alpha") val = alpha_vals.at(subfields[1]);
		else {
			std::stringstream sval(subfields[1]);
			sval >> val;
		}
		signalParameters_[subfields[0]] = val;
	}
}

void GenVecAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
	if(!signal) return;
	edm::Handle<GenLumiInfoHeader> gen_header;
	iLumi.getByToken(tok_scan, gen_header);
	getSVJComment(*gen_header);
}

bool GenVecAnalyzer::isParticle(const PidSet& darkList, const reco::GenParticle& part) const {
  return isParticle(darkList, part.pdgId());
}

bool GenVecAnalyzer::isParticle(const PidSet& darkList, CandPtr part) const {
  return isParticle(darkList, part->pdgId());
}

bool GenVecAnalyzer::isParticle(const PidSet& darkList, int pid) const {
  return darkList.find(std::abs(pid)) != darkList.end();
}

// this function intends to collect immediate non-mediator daughters of the mediators. These mediator daughters can then be used to reconstruct the mass of the mediator.
void GenVecAnalyzer::medDecay(CandPtr part, CandPtr& firstQdM1, CandPtr& firstQdM2, CandPtr& firstQsM1, CandPtr& firstQsM2, bool& secondDM, bool& secondSM) {
	for(unsigned i = 0; i < part->numberOfDaughters(); i++){
		CandPtr dau = part->daughter(i);
		// if the first dark mediator's daughter is still a dark mediator, then check the daughters of this daughter dark mediator until we get daughters that are not dark mediator
		if(isParticle(DarkTMediatorIDs_,dau)) medDecay(dau,firstQdM1,firstQdM2,firstQsM1,firstQsM2,secondDM,secondSM);
		else{
			// a mediator decays into a dark and an SM quark. Here we are collecting the dark quarks from the mediators while labeling them firstQdM1 and firstQdM2 depending on which mediator the quarks came from.
			if(isParticle(DarkQuarkIDs_,dau)){
				// this condition makes sure that the firstQdM1 (firstQdM2) and firstQsM1 (firstQsM2) came from the same mediator.
				// The labels 1 and 2 have no significance other than making sure that we get the correct pairings of dark and SM quarks from the mediators.
				if(secondDM == false){
					entry.Mediator1 = part->p4();
					entry.DarkQuark1 = dau->p4();
					firstQdM1 = dau;
					secondDM = true;
				}
				else {
					entry.Mediator2 = part->p4();
					entry.DarkQuark2 = dau->p4();
					firstQdM2 = dau;
				}
			}
			// Here we are collecting the SM quarks from the mediators which assigning them to firstQsM1 and firstQsM2 depending on which mediator the quarks came from.
			else if(isParticle(SMQuarkIDs_,dau)){
				if(secondSM == false){
					entry.SMQuark1 = dau->p4();
					firstQsM1 = dau;
					secondSM = true;
				}
				else {
					entry.SMQuark2 = dau->p4();
					firstQsM2 = dau;
				}
			}
		}
	}
}

void GenVecAnalyzer::firstDark(CandPtr part, CandSet& firstMd, CandPtr& firstQdM1, CandPtr& firstQdM2, CandPtr& firstQsM1, CandPtr& firstQsM2, bool& secondDM, bool& secondSM) {
	if(isParticle(DarkFirstIDs_,part)){
		CandPtr parent = part->mother(0);
		if(isParticle(DarkFirstIDs_,parent)) firstDark(parent, firstMd, firstQdM1, firstQdM2, firstQsM1, firstQsM2, secondDM, secondSM);
		else{
			// SM parent of first dark particles
			// looping through the daughters of this SM parent
			for(unsigned i = 0; i < part->numberOfDaughters(); i++){
				CandPtr dau = part->daughter(i);
				if (isParticle(DarkTMediatorIDs_,dau)){
					if(firstMd.find(dau)==firstMd.end()){
						firstMd.insert(dau);
						// once a mediator daughter is found, we look for the descendants of the mediator
						medDecay(dau,firstQdM1,firstQdM2,firstQsM1,firstQsM2,secondDM,secondSM);
					}
				}
			}
		}
	}
}

// match particles from mediator to jets
//based on official JetMET matching procedure: equiv to set<DR,jet_index,gen_index> sorted by DR
//from https://github.com/cms-jet/JetMETAnalysis/blob/master/JetUtilities/plugins/MatchRecToGen.cc
std::vector<int> GenVecAnalyzer::matchJetPart(const vector<reco::GenJet>& jets, const vector<CandPtr>& parts) const {
  std::map<double,std::pair<unsigned,unsigned>> matchMap;
  for(unsigned j = 0; j < jets.size(); ++j){
    for(unsigned p = 0; p < parts.size(); ++p){
      matchMap.emplace(std::piecewise_construct, std::forward_as_tuple(reco::deltaR(jets[j],*parts[p])), std::forward_as_tuple(j,p));
    }
  }

  std::vector<int> jetIndex(parts.size(),-1);
  std::unordered_set<unsigned> j_used, p_used;
  for(const auto& matchItem: matchMap){
    unsigned j = matchItem.second.first;
    unsigned p = matchItem.second.second;
    if(j_used.find(j)==j_used.end() and p_used.find(p)==p_used.end()){
      jetIndex[p] = j;
      j_used.insert(j);
      p_used.insert(p);
    }
  }

  return jetIndex;
}

double GenVecAnalyzer::girth(const reco::GenJet& jet) const {
	double girth_val = 0;
	for(unsigned k = 0; k < jet.numberOfDaughters(); ++k){
		const reco::Candidate* part = jet.daughter(k);
		float dR = reco::deltaR(jet.p4(),part->p4());
		float pT = part->pt();
		girth_val += pT*dR;
	}
	girth_val /= jet.pt();
	return girth_val;
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
	if(signal) iEvent.getByToken(tok_part,h_part);

	int jet_counter = 0;
	entry.Nconst1 = 0;
	entry.Nconst2 = 0;
	entry.Girth1 = 0;
	entry.Girth2 = 0;
	entry.Girth3 = 0;
	for(const auto& i_jet : *(h_jet.product())){
		if(jet_counter==0){
			entry.Jet1 = i_jet.p4();
			entry.Nconst1 = i_jet.numberOfDaughters();
			entry.Girth1 = girth(i_jet);
		}
		else if(jet_counter==1){
			entry.Jet2 = i_jet.p4();
			entry.Nconst2 = i_jet.numberOfDaughters();
			entry.Girth2 = girth(i_jet);
		}
		else if(jet_counter==2){
			entry.Jet3 = i_jet.p4();
			entry.Girth3 = girth(i_jet);
		}
		else if(jet_counter==3){
			entry.Jet4 = i_jet.p4();
		}
		else break;
		++jet_counter;
	}

	const auto& i_met = h_met->front();
	entry.Met = i_met.p4();
	double METx = i_met.px();
	double METy = i_met.py();

	const double jet_radius = 0.8;

	entry.Ninv = 0;
	entry.Ninv1 = 0;
	entry.Ninv2 = 0;
	entry.Nvis = 0;
	entry.Nvis1 = 0;
	entry.Nvis2 = 0;

	if(signal){
		//t-channel jet categorization
		//naming scheme: dark particles Pd, SM particles Ps; P = D (generic daughters), Q (quarks), G (gluons), M (mediators)
		CandSet firstMd;
		CandPtr firstQdM1, firstQdM2, firstQsM1, firstQsM2;
		bool secondDM = false, secondSM = false;
		std::vector<const reco::GenJet*> dQM1Js,dQM2Js,SMM1Js,SMM2Js;

		for(const auto& i_part : *(h_part.product())){
			//todo: for unstable dark hadrons, can actually make list of all that produced genjet constituents
			if(isParticle(DarkHadronIDs_,i_part) and i_part.numberOfDaughters()>0){
				bool invis = isParticle(DarkStableIDs_,i_part.daughter(0));
				if(invis) {
					++entry.Ninv;
					entry.Invis += i_part.p4();
				}
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

			//get truth: last zprime, first dark quarks
			if(isParticle(DarkSMediatorIDs_,i_part) and i_part.isLastCopy()){
				entry.Zprime = i_part.p4();
				if(i_part.numberOfDaughters()>0) entry.DarkQuark1 = i_part.daughter(0)->p4();
				if(i_part.numberOfDaughters()>1) entry.DarkQuark2 = i_part.daughter(1)->p4();
			}

			//t-channel loop over gen particles
			firstDark(&i_part, firstMd, firstQdM1, firstQdM2, firstQsM1, firstQsM2, secondDM, secondSM);
		}

		//t-channel loop over gen jets
		vector<CandPtr> firsts{firstQdM1,firstQdM2,firstQsM1,firstQsM2};
		vector<int> jetIndex(firsts.size(),-1);
		if(firstMd.size()==2){
			jetIndex = matchJetPart(*(h_jet.product()),firsts);
			entry.PairMT2ID = vector<int>(h_jet->size(),0);
			for(unsigned p = 0; p < jetIndex.size(); ++p){
				if(jetIndex[p]>-1) entry.PairMT2ID[jetIndex[p]] = p+1;
			}
		}
		else {
			entry.PairMT2ID = std::vector<int>(h_jet->size(),0);
		}

		//compute MT2 for t-channel (paired jets)
		bool geqOneJet1 = jetIndex[0] >= 0 or jetIndex[2] >= 0;
		bool geqOneJet2 = jetIndex[1] >= 0 or jetIndex[3] >= 0;
		if(geqOneJet1 and geqOneJet2){
			// using the highest pT jet if more than one jet contains the same particle
			LorentzVector FJet0;
			if(jetIndex[0] >= 0) FJet0 += h_jet->at(jetIndex[0]).p4();
			if(jetIndex[2] >= 0) FJet0 += h_jet->at(jetIndex[2]).p4();
			LorentzVector FJet1;
			if(jetIndex[1] >= 0) FJet1 += h_jet->at(jetIndex[1]).p4();
			if(jetIndex[3] >= 0) FJet1 += h_jet->at(jetIndex[3]).p4();
			entry.PairMT2 = asymm_mt2_lester_bisect::get_mT2(
				FJet0.M(), FJet0.Px(), FJet0.Py(),
				FJet1.M(), FJet1.Px(), FJet1.Py(),
				METx, METy, 0.0, 0.0, 0
			);
		}
		else
		{
			entry.PairMT2 = (!geqOneJet1 and !geqOneJet2) ? 0 : (!geqOneJet1 ? -1 : -2);
		}
	}

	//compute mass variables
	auto vjj = entry.Jet1 + entry.Jet2;
	entry.Mdq = (entry.DarkQuark1 + entry.DarkQuark2).mass();
	entry.Mmc = (vjj + entry.Invis1 + entry.Invis2).mass();
	entry.Mjj = vjj.mass();
	entry.MT = TransverseMass(vjj.px(),vjj.py(),vjj.mass(),entry.Met.px(),entry.Met.py(),0.0);
	//for boosted case
	entry.MTx = TransverseMass(entry.Jet2.px(),entry.Jet2.py(),entry.Jet2.mass(),entry.Met.px(),entry.Met.py(),0.0);
	entry.Meff = entry.Jet2.pt()+entry.Met.pt();

	//compute MT2 and MAOS
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

	//"reco-level" MT2: find pairs of jets with most similar invariant mass values
	//based on: https://github.com/cms-svj/t-channel_Analysis/blob/2592bd5a5b313747f51d59bd1dc10842b28d6f64/utils/utility.py#L99
	if(jet_counter>=4){
		vector<vector<int>> combos{{0,1,2,3},{0,2,1,3},{0,3,1,2}};
		vector<LorentzVector*> jets{&entry.Jet1,&entry.Jet2,&entry.Jet3,&entry.Jet4};
		int min_index = -1;
		double min_val = 1e10;
		for(unsigned c = 0; c < combos.size(); ++c){
			const auto& combo = combos[c];
			vector<LorentzVector*> _jets{jets[combo[0]],jets[combo[1]],jets[combo[2]],jets[combo[3]]};
			double diff = std::abs((*_jets[0]+*_jets[1]).mass() - (*_jets[2]+*_jets[3]).mass());
			if(diff<min_val){
				min_index = c;
				min_val = diff;
			}
		}
		const auto& combo = combos[min_index];
		auto FJet0 = *jets[combo[0]] + *jets[combo[1]];
		auto FJet1 = *jets[combo[2]] + *jets[combo[3]];
		entry.PairMT2Reco = asymm_mt2_lester_bisect::get_mT2(
			FJet0.M(), FJet0.Px(), FJet0.Py(),
			FJet1.M(), FJet1.Px(), FJet1.Py(),
			METx, METy, 0.0, 0.0, 0
		);
		entry.PairMT2RecoID = combo;
	}
	else {
		entry.PairMT2RecoID = std::vector<int>(jet_counter,0);
		entry.PairMT2Reco = 0;
	}

	//check signal scan info
	if(!signalParameters_.empty()){
		entry.nMediator = signalParameters_["nMed"];
		entry.mMediator = signalParameters_["mMed"];
		entry.mDark = signalParameters_["mDark"];
		entry.rinv = signalParameters_["rinv"];
		entry.alpha = signalParameters_["alpha"];
		entry.yukawa = signalParameters_["yukawa"];
	}
	else {
		entry.nMediator = 0;
		entry.mMediator = 0.;
		entry.mDark = 0.;
		entry.rinv = 0.;
		entry.alpha = 0.;
		entry.yukawa = 0.;
	}

	entry.Weight = 1;
	if(!signal){
		edm::Handle<GenEventInfoProduct> h_weight;
		iEvent.getByToken(tok_weight,h_weight);
		entry.Weight = h_weight->weight();
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
	desc.add<std::string>("model","");
	desc.add<bool>("signal",true);
	
	descriptions.add("GenVecAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenVecAnalyzer);
