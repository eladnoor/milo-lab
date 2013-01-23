#include <sstream>
#include "holders.h"
#include "dt_smarts.h"


/* Smiles Holder */

void SmilesHolder::addSmile (std::string smile) {
	_smilesList.insert(smile);
}

void SmilesHolder::addSmiles (SmilesHolder smiles) {
	auto smilesList = smiles.getList();
	for (auto smile = smilesList.begin(); smile != smilesList.end(); ++smile) {
		this->addSmile(*smile);
	}
}

bool SmilesHolder::exists (std::string smile) {
	return _smilesList.count(smile) > 0;
}

const std::set<std::string> SmilesHolder::getList () {
	return _smilesList;
}
std::string SmilesHolder::getListAsString (char delimiter) {
	std::stringstream s;
	auto list = getList();
	auto iterator = list.begin();
	if (iterator != list.end()) {
		s << *iterator;
		for (++iterator; iterator != list.end(); ++iterator) {
			s << delimiter << *iterator;
		}
	}
	return s.str();
}

SmilesHolder::~SmilesHolder () {}



std::string getReactionProduct (std::string reaction) {
	return reaction.substr(reaction.find_last_of('>')+1);
}


/* Reaction Holder */

void ReactionHolder::addReaction (const dt_String smirk, dt_Integer length) {
	dt_Integer slen;
	boost::shared_ptr<DaylightHandle> ptr (new DaylightHandle(dt_smirkin(slen,dt_smarts_opt(&slen, length, smirk, false))));
	_reactionsList.insert(_reactionsList.end(), ptr);
}


const std::vector<SmilesHolder> ReactionHolder::runReactionsOverReactant (DaylightHandle reactant) {
	std::vector<SmilesHolder> productVector;
	dt_Integer slen;
	dt_Handle productList, rxn, atom, atoms;
	dt_String out;
	for (auto reaction = _reactionsList.begin(); reaction != _reactionsList.end(); ++reaction) {
		SmilesHolder smi;
		productList = dt_transform(reaction->get()->get(), reactant.get(),  DX_FORWARD, 0);

		while (NULL_OB != (rxn = dt_next(productList))) {
			atoms = dt_stream(rxn, TYP_ATOM);
			while (NULL_OB != (atom = dt_next(atoms))) {
				dt_setmap(atom, NULL_OB);
			}
			out = dt_cansmiles(&slen, rxn, true);
			_collector.addSmile(out);
			smi.addSmile(getReactionProduct(out));
		}
	productVector.insert(productVector.end(), smi);
	}
	dt_dealloc(productList);
	dt_dealloc(rxn);
	return productVector;
}

const std::vector<SmilesHolder> ReactionHolder::runReactionsOverReactant (const dt_String reactant, dt_Integer length) {
	DaylightHandle reactantHandle = DaylightHandle(dt_smilin(length, reactant));
	return runReactionsOverReactant(reactantHandle);
}

std::vector<boost::shared_ptr<DaylightHandle>> ReactionHolder::getList () {
	return _reactionsList;
}

SmilesHolder ReactionHolder::collectAllReactions () {
	return _collector;
}


ReactionHolder::~ReactionHolder () {}
