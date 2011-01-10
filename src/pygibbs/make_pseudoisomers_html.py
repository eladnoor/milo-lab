#!/usr/bin/python

import logging

from pygibbs import group_decomposition
from pygibbs import pseudoisomers_data
from pygibbs import templates


def main():
    pdata = pseudoisomers_data.PseudoisomersData.FromFile(
        '../data/thermodynamics/dG0.csv')
    decomposer = group_decomposition.GroupDecomposer.FromGroupsFile(
        '../data/thermodynamics/groups_species.csv')
    
    train_pseudoisomers = []
    test_pseudoisomers = []
    skip_pseudoisomers = []
    for pisomer in pdata:
        mol = pisomer.Mol()
        decomposition = None
        if mol:
            decomposition = decomposer.Decompose(pisomer.Mol())
            
            if decomposition.unassigned_nodes:
                logging.warning('%s didn\'t decompose', pisomer)
            
        pisomer_dict = {'data': pisomer,
                        'decomposition': decomposition}
        if pisomer.Train():
            train_pseudoisomers.append(pisomer_dict)
        elif pisomer.Test():
            test_pseudoisomers.append(pisomer_dict)
        elif pisomer.Skip():
            skip_pseudoisomers.append(pisomer_dict)
        else:
            logging.warning('Found a pseudoisomer that is not marked as'
                            ' test, train, or skip.')

    template_data = {'pseudoisomer_collections': 
                     [{'name': 'Train',
                       'pseudoisomers': train_pseudoisomers},
                      {'name': 'Test',
                       'pseudoisomers': test_pseudoisomers},
                      {'name': 'Skip',
                       'pseudoisomers': skip_pseudoisomers}
                     ]}
    templates.render_to_file('pseudoisomers_ground_truth.html',
                             template_data,
                             '../res/pseudoisomers_ground_truth.html')
    
if __name__ == '__main__':
    main()
    logging.info('Done')