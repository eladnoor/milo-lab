import sys
import libsbml
import semanticSBML.annotate
import csv
from toolbox.database import SqliteDatabase
from pygibbs.thermodynamics import PsuedoisomerTableThermodynamics
from optparse import OptionParser

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

document = libsbml.readSBML(options.sbml_model_filename)
if document.getNumErrors():
    raise Exception('cannot read SBML model from file %s due to error: %s' % 
                    (options.sbml_model_filename, document.getError(0).getMessage()))

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

model_annotation = semanticSBML.annotate.ModelElementsAnnotations(document.getModel(), suppress_errors=True)

rowdicts = []
# go through all the elements
for element_annotation in model_annotation.getElementAnnotations():
    # and through all their annotations
    for annotation in element_annotation.getAnnotations():
        if annotation.db == 'KEGG Compound':
            rowdict = {}
            cid = int(annotation.id[1:])
            rowdict['cid'] = annotation.id
            rowdict['element_id'] = element_annotation.id
            rowdict['dG0\'_f'] = thermo.cid2dG0_tag(cid)
            rowdicts.append(rowdict)
            
csv_writer = csv.DictWriter(open(options.csv_output_filename, 'w'), 
    fieldnames=['element_id', 'cid', 'dG0\'_f'])
csv_writer.writer.writerow(['element_id', 'cid', 'dG0\'_f'])
csv_writer.writerows(rowdicts)