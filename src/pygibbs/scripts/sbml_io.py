import sys
import libsbml
import semanticSBML.annotate
import csv
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from optparse import OptionParser
import logging
from pygibbs.kegg import Kegg

def MakeOpts():
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--sbml_model_filename", 
                          dest="sbml_model_filename",
                          default=None,
                          help="The SBML model filename")
    opt_parser.add_option("-o", "--csv_output_filename", 
                          dest="csv_output_filename",
                          default=None,
                          help="The output filename (CSV format)")
    opt_parser.add_option("-k", "--kegg_database_location", 
                          dest="kegg_db_filename",
                          default="../data/public_data.sqlite",
                          help="The KEGG database location")
    opt_parser.add_option("-c", "--compounds_out_filename",
                          dest="compounds_out_filename",
                          default="../res/kegg_compounds.json",
                          help="Compounds output filename.")
    opt_parser.add_option("-t", "--thermodynamics_filename",
                          dest="thermo_filename",
                          default='../data/thermodynamics/dG0.csv',
                          help="Thermodynamics filename")
    opt_parser.add_option("-d", "--database_location", 
                          dest="db_filename",
                          default="../res/gibbs.sqlite",
                          help="The Thermodynamic database location")
    opt_parser.add_option("-g", "--gc_table_name",
                          dest="gc_table_name",
                          default='gc_pseudoisomers',
                          help="Group Contribution Table Name")    
    return opt_parser

options, _ = MakeOpts().parse_args(sys.argv)
if options.sbml_model_filename == None:
    raise ValueError("Must provide a SBML model")

print 'SBML model filename:', options.sbml_model_filename
print 'CSV output filename:', options.csv_output_filename
print 'KEGG Database filename:', options.kegg_db_filename
print 'Observed Thermodynamics filename:', options.thermo_filename
print 'Thermodynamic Database filename:', options.db_filename
print 'Group Contribution Table Name:', options.gc_table_name

db = SqliteDatabase(options.db_filename)
observed_thermo = PsuedoisomerTableThermodynamics.FromCsvFile(
    options.thermo_filename)
if not db.DoesTableExist(options.gc_table_name):
    raise ValueError('The table %s does not exist in the database. '
                     'Please run the groups.py script and try again.'
                     % options.gc_table_name)
thermo = PsuedoisomerTableThermodynamics.FromDatabase(
    db, options.gc_table_name)
thermo.override_data(observed_thermo)
kegg = Kegg.getInstance()

document = libsbml.readSBML(options.sbml_model_filename)
if document.getNumErrors():
    raise Exception('cannot read SBML model from file %s due to error: %s' % 
                    (options.sbml_model_filename, document.getError(0).getMessage()))
model = document.getModel()
logging.info('Done parsing the model: ' + model.getName())
model_annotation = semanticSBML.annotate.ModelElementsAnnotations(model, suppress_errors=True)

species_ids_to_names = dict([(s.getId(), s.getName()) for s in model.getListOfSpecies()])
reaction_ids_to_names = dict([(r.getId(), r.getName()) for r in model.getListOfReactions()])

rowdicts = []
# go through all the elements
for element_annotation in model_annotation.getElementAnnotations():
    # and through all their annotations
    for annotation in element_annotation.getAnnotations():
        if annotation.db == 'KEGG Compound':
            rowdict = {}
            cid = int(annotation.id[1:])
            rowdict['Quantity Type'] = 'standard chemical potential'
            rowdict['Miriam ID'] = 'urn:miriam:kegg.compound:' + annotation.id
            rowdict['SBML Element ID'] = element_annotation.id
            rowdict['SBML Element Name'] = species_ids_to_names.get(element_annotation.id, '')
            rowdict['KEGG name'] = kegg.cid2name(cid)
            #rowdict['kegg id'] = annotation.id
            rowdict['Value'] = '%.1f' % thermo.cid2dG0_tag(cid)
            rowdict['Unit'] = 'kJ/mol'
            rowdicts.append(rowdict)
        if annotation.db == 'KEGG Reaction':
            rowdict = {}
            rid = int(annotation.id[1:])
            reaction = kegg.rid2reaction(rid)
            rowdict['Quantity Type'] = 'standard Gibbs energy of reaction'
            rowdict['Miriam ID'] = 'urn:miriam:kegg.reaction:' + annotation.id
            rowdict['SBML Element ID'] = element_annotation.id
            rowdict['SBML Element Name'] = reaction_ids_to_names.get(element_annotation.id, '')
            rowdict['KEGG name'] = reaction.name
            #rowdict['kegg id'] = annotation.id
            rowdict['Value'] = '%.1f' % reaction.PredictReactionEnergy(thermo)
            rowdict['Unit'] = 'kJ/mol'
            rowdicts.append(rowdict)
            
fieldnames=['Quantity Type', 'SBML Element Name', 'SBML Element ID', 
            'Miriam ID', 'KEGG name', 'Value', 'Unit']            
csv_writer = csv.DictWriter(open(options.csv_output_filename, 'w'), fieldnames)
csv_writer.writer.writerow(fieldnames)
csv_writer.writerows(rowdicts)