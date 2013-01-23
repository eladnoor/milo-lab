#pragma once

#include <vector>
#include <set>
#include <string>
#include <memory>
#include <boost/shared_ptr.hpp>
#include "dt_smiles.h"

/**
 * This class is just a wrapper for the dt_Handle type for automatic deallocation.
 */
class DaylightHandle {
private:
	dt_Handle _handle;
public:
	/**
	 * Initiation using dt_handle.
	 */
	DaylightHandle (dt_Handle handle) { _handle = handle; }
	/**
	 * Returns an dt_handle.
	 */
	dt_Handle get () { return _handle; }
	virtual ~DaylightHandle () { dt_dealloc(_handle); }
};




/**
 * This class holds list of smiles. The main importance of this class is that it does
 * not hold the same smiles twice.
 * For example, CCO and OCC will appear only once in the holder.
 */
class SmilesHolder {
private:
	std::set<std::string> _smilesList;
public:
	/**
	 * Adds a single smile.
	 */
	void addSmile (std::string smile);
	/**
	 * Adds a list of smiles, as SmilesHolder list
	 */
	void addSmiles (SmilesHolder smiles);
	/**
	 * Check whether a given smile exists.
	 */
	bool exists (std::string smile);
	/**
	 * Returns the list of smiles as a vector.
	 */
	const std::set<std::string> getList ();
	/**
	 * Returns the list of smiles as a string.
	 * Optional: add a separator for the list.
	 */
	std::string getListAsString (char delimiter = '.');
	virtual ~SmilesHolder ();
};




/**
 * This class holds a list of reaction and can perform several tasks above them,
 * like execute them.
 */
class ReactionHolder {
private:
	std::vector<boost::shared_ptr<DaylightHandle>> _reactionsList;
	SmilesHolder _collector;
public:
	/**
	 * Adds a new reaction to the list of reactions, where reaction is an smirk string.
	 */
	void addReaction (const dt_String smirk, dt_Integer length);
	/**
	 * For each reaction in the holder, returns a sequence holder of all the products
	 * of the reaction (the reactant given as handle).
	 */
	const std::vector<SmilesHolder> runReactionsOverReactant (DaylightHandle reactant);
	/**
	 * For each reaction in the holder, returns a sequence holder of all the products
	 * of the reaction (the reactant given as smiles).
	 */
	const std::vector<SmilesHolder> runReactionsOverReactant (const dt_String reactant, dt_Integer length);
	/**
	 * Returns a list of all the reactions.
	 */
	std::vector<boost::shared_ptr<DaylightHandle>> getList ();
	/**
	 * collectAllReactions returns a SmilesHolder object which contains all the reactions
	 * that already ran in using this object.
	 */
	SmilesHolder collectAllReactions ();
	virtual ~ReactionHolder ();
};
