#!/usr/bin/python

import pybel

def FindSmarts(mol, smarts_str):
    """
    Corrects the pyBel version of Smarts.findall() which returns results as tuples,
    with 1-based indices even though Molecule.atoms is 0-based.

    Args:
        mol: the molecule to search in.
        smarts_str: the SMARTS query to search for.
    
    Returns:
        The re-mapped list of SMARTS matches.
    """
    shift_left = lambda m: [(n - 1) for n in m] 
    return map(shift_left, pybel.Smarts(smarts_str).findall(mol))