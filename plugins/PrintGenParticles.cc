//framework headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Utilities/interface/Exception.h"

//analysis headers
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//STL headers 
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace {
	template <class T, class O>
	void printvec(const std::vector<T>& vec, O& out, const std::string& delim){
		if(!vec.empty()){
			//avoid trailing delim
			std::copy(vec.begin(),vec.end()-1,std::ostream_iterator<T>(out,delim.c_str()));
			//last element
			out << vec.back();
		}
	}
}

class PrintGenParticles : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
	explicit PrintGenParticles(const edm::ParameterSet&);
	~PrintGenParticles() {}

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	void analyze(const edm::Event&, const edm::EventSetup&) override;

	//members
	std::string fileName;
	edm::EDGetTokenT<reco::CandidatePtrVector> tok_part;
};

PrintGenParticles::PrintGenParticles(const edm::ParameterSet& iConfig) :
	fileName(iConfig.getParameter<std::string>("fileName")),
	tok_part(consumes<reco::CandidatePtrVector>(iConfig.getParameter<edm::InputTag>("PartTag")))
{
	//remove unnecessary prefix
	std::string pref("file:");
	if(fileName.rfind(pref,0)==0){
		fileName = fileName.substr(pref.size());
	}
	std::string suff(".root");
	auto pos = fileName.find(suff);
	if(pos!=std::string::npos){
		fileName.erase(pos, suff.size());
	}
}

void PrintGenParticles::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
	edm::Handle<reco::CandidatePtrVector> h_part;
	iEvent.getByToken(tok_part,h_part);

	std::stringstream fname;
	fname << fileName << "_" << iEvent.id().event() << ".dat";
	std::ofstream file(fname.str());
	if(!file.is_open()){
		throw cms::Exception("Could not open file: "+fname.str());
	}
	for(const auto& i_part : *(h_part.product())){
		printvec<double>({i_part->energy(),i_part->px(),i_part->py(),i_part->pz()},file,"\t");
		file << std::endl;
	}
	file.close();
}

void PrintGenParticles::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.add<edm::InputTag>("PartTag",edm::InputTag("genParticlesForJetsNoNu"));
	desc.add<std::string>("fileName","");

	descriptions.add("PrintGenParticles",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrintGenParticles);
