import sys
import csv
from optparse import OptionParser
from pygibbs.kegg import Kegg
from pygibbs.thermodynamic_constants import default_I, default_pH, default_pMg,\
    default_T
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_errors import MissingCompoundFormationEnergy

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-c", "--csv_input_filename",
                          dest="csv_input_filename",
                          default=None,
                          help="input CSV file with InChIs")
    opt_parser.add_option("-o", "--csv_output_filename",
                          dest="csv_output_filename",
                          default=None,
                          help="output CSV file with dGs and groups")
    opt_parser.add_option("-b", "--biochemical",
                          dest="biochemical",
                          default=True,
                          action="store_true",
                          help="calculate the biochemical (transformed) energy at standard conditions")
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
                          default="PGC",
                          help="The thermodynamic data to use")
    return opt_parser

def CalculateThermo():
    estimators = LoadAllEstimators()
    parser = MakeOpts(estimators)
    options, _ = parser.parse_args(sys.argv)
    if options.csv_input_filename is None:
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
    
    csv_reader = csv.reader(open(options.csv_input_filename, 'r'))
    headers = csv_reader.next()
    cid_index = headers.index('kegg id')
    
    csv_writer = csv.writer(out_fp)
    if options.biochemical:
        csv_writer.writerow(headers + ['dG0\'', 'pH', 'I', 'pMg', 'T'])
    else:
        csv_writer.writerow(headers + ['dG0', 'nH', 'z', 'nMg'])

    for row in csv_reader:
        try:
            cid = int(row[cid_index][1:])
        except ValueError:
            csv_writer.writerow(row + ['NaN', 'cannot parse KEGG ID'])
            continue
        try:
            if options.biochemical:
                dG0_prime = estimator.cid2dG0_tag(cid, pH=pH, I=I, pMg=pMg, T=T)
                csv_writer.writerow(row + [dG0_prime, pH, I, pMg, T])
            else:
                for nH, z, nMg, dG0 in estimator.cid2PseudoisomerMap(cid).ToMatrix():
                    csv_writer.writerow(row + [dG0, nH, z, nMg])
        except MissingCompoundFormationEnergy as e:
            csv_writer.writerow(row + ['NaN', str(e)])

if __name__ == '__main__':
    CalculateThermo()