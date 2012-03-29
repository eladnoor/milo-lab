import sys
from optparse import OptionParser
from pygibbs.nist_verify import LoadAllEstimators

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
    return opt_parser

def ExportThermo():
    estimators = LoadAllEstimators()
    options, _ = MakeOpts(estimators).parse_args(sys.argv)
    thermo = estimators[options.thermodynamics_source]

    print 'Thermodynamic source:', thermo.name
    print 'CSV output filename:', options.csv_out_filename

    thermo.write_data_to_csv(options.csv_out_filename)
    
if __name__ == "__main__":
    ExportThermo()
