import sys
from optparse import OptionParser
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_constants import default_pH, default_I, default_T,\
    default_pMg

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    opt_parser = OptionParser()
    opt_parser.add_option("-s", "--thermodynamics_source",
                          dest="thermodynamics_source",
                          type="choice",
                          choices=estimators.keys(),
                          default="UGC",
                          help="The thermodynamic data to use")
    opt_parser.add_option("-c", "--csv_out_filename",
                          dest="csv_out_filename",
                          default="../res/thermo.csv",
                          help="CSV output filename")
    opt_parser.add_option("-b", "--biochemical",
                          dest="biochemical",
                          default=False,
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
    return opt_parser

def ExportThermo():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    thermo = estimators[options.thermodynamics_source]

    print 'Thermodynamic source:', thermo.name
    print 'CSV output filename:', options.csv_out_filename

    if not options.biochemical:
        thermo.write_chemical_data_to_csv(options.csv_out_filename)
    else:
        thermo.SetConditions(pH=options.pH, I=options.I,
                             pMg=options.pMg, T=options.T)
        thermo.write_biochemical_data_to_csv(options.csv_out_filename)
    
if __name__ == "__main__":
    ExportThermo()
