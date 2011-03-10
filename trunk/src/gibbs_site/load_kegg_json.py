#!/usr/bin/env python

import json
import logging

from util import django_utils
from util import kegg

django_utils.SetupDjango()

from gibbs import models


def GetOrCreateNames(names_list):
    """Find all the names in the database.
    
    Create them if they are not present.
    """
    return [models.CommonName.GetOrCreate(n)
            for n in names_list] 
    
    
def GetReactions(rids_list):
    """Find all the given reactions in the database.
    
    Skip those that are not present.
    """
    rxns = []
    for rid in rids_list:
        try:
            rxns.append(models.StoredReaction.objects.get(kegg_id=rid))
        except Exception:
            logging.warning('Failed to retrieve reaction %s', rid)
            continue
    return rxns
        
         
def LoadKeggCompounds(kegg_json_filename='data/kegg_compounds.json'):
    parsed_json = json.load(open(kegg_json_filename))
    
    for cd in parsed_json:
        try:
            cid = cd['CID']
            formula = cd['formula']
            mass = float(cd['mass'])
            inchi = cd['InChI']
            num_electrons = cd.get('num_electrons')
            c = models.Compound(kegg_id=cid,
                                formula=formula,
                                inchi=inchi,
                                mass=mass)
            
            if num_electrons != None:
                c.num_electrons = int(num_electrons)
                
            # Need to save before setting up many-to-many relationships.
            c.save()
                        
            names = GetOrCreateNames(cd['names'])            
            for n in names:
                c.common_names.add(n)
            c.save()        
        except Exception, e:
            logging.warning('Error parsing cid %s', cid)
            logging.error(e)
            continue
        

def LoadKeggReactions(reactions_json_filename='data/kegg_reactions.json'):
    parsed_json = json.load(open(reactions_json_filename))

    for rd in parsed_json:
        try:
            rid = rd['RID']
            rxn = rd['reaction']
            reactants = []
            products = []
            for coeff, cid in rxn:
                reactant = models.Reactant.GetOrCreate(cid, abs(coeff))
                if coeff < 0:
                    reactants.append(reactant)
                else:
                    products.append(reactant)
                
            # Need to save once.
            rxn = models.StoredReaction(kegg_id=rid)
            rxn.save()
            
            for reactant in reactants:
                rxn.reactants.add(reactant)
            for product in products:
                rxn.products.add(product)
            rxn.save()
        except Exception, e:
            logging.warning('Error parsing rid %s', rid)
            logging.error(e)
            continue


def LoadKeggEnzymes(enzymes_json_filename='data/kegg_enzymes.json'):
    parsed_json = json.load(open(enzymes_json_filename))

    for ed in parsed_json:
        try:
            ec = ed['EC']
            names = GetOrCreateNames(ed['names'])
            reactions = GetReactions(ed['reaction_ids'])
            
            enz = models.Enzyme(ec=ec)
            enz.save()
            
            for name in names:
                enz.common_names.add(name)
            for rxn in reactions:
                enz.reactions.add(rxn)
            enz.save()
        except Exception, e:
            logging.warning('Error parsing ec %s', ec)
            logging.error(e)
            continue


def LoadAllKeggData():
    LoadKeggCompounds()
    LoadKeggReactions()
    LoadKeggEnzymes()    


if __name__ == '__main__':
    LoadAllKeggData()
