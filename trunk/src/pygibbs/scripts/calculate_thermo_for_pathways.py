import sys
import csv
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg,\
    default_T
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-c", "--input_filename",
                          dest="input_filename",
                          default="../data/thermodynamics/pathways.txt",
                          help="The file to read for pathways to analyze.")
    opt_parser.add_option("-o", "--csv_output_filename",
                          dest="csv_output_filename",
                          default=None,
                          help="output CSV file with dGs and groups")
    opt_parser.add_option("-p", "--ph", action="store", type="float",
                          dest="pH", default=default_pH,
                          help="The width of the range of pH values to plot")
    opt_parser.add_option("-i", "--ionic_strength", action="store", type="float",
                          dest="I", default=default_I,
                          help="The Ionic strength of the solution (in M)")
    opt_parser.add_option("-m", "--pmg", action="store", type="float",
                          dest="pMg", default=default_pMg,
                          help="The pMg of the solution")
    opt_parser.add_option("-t", "--temperature", action="store", type="float",
                          dest="T", default=default_T,
                          help="The T of the solution (in K)")
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="merged",
                          help="The thermodynamic data to use")
    return opt_parser

def CalculateThermo():
    estimators = LoadAllEstimators()
    parser = MakeOpts(estimators)
    options, _ = parser.parse_args(sys.argv)
    if options.input_filename is None:
        sys.stderr.write(parser.get_usage())
        sys.exit(-1)
    
    estimator = estimators[options.thermodynamics_source]
    
    pH, I, pMg, T = options.pH, options.I, options.pMg, options.T
    kegg = Kegg.getInstance()

    if options.csv_output_filename is not None:
        out_fp = open(options.csv_output_filename, 'w')
        print "writing results to %s ... " % options.csv_output_filename
    else:
        out_fp = sys.stdout

    entry2fields_map = ParsedKeggFile.FromKeggFile(options.input_filename)
    all_reactions = []
    for key in sorted(entry2fields_map.keys()):
        field_map = entry2fields_map[key]
        p_data = PathwayData.FromFieldMap(field_map)
        if p_data.skip:
            continue

        cid_mapping = p_data.cid_mapping
        field_map = p_data.field_map
        _, _, _, reactions = kegg.parse_explicit_module_to_reactions(field_map, cid_mapping)
        all_reactions += reactions
    S, cids = kegg.reaction_list_to_S(all_reactions)
    dG0_r = estimator.GetTransfromedReactionEnergies(S, cids)
    
    csv_writer = csv.writer(out_fp)
    csv_writer.writerow(['reaction', 'dG0\'', 'pH', 'I', 'pMg', 'T'])

    for r, reaction in enumerate(all_reactions):
        csv_writer.writerow([reaction.FullReactionString(), dG0_r[r, 0], pH, I, pMg, T])

if __name__ == '__main__':
    CalculateThermo()