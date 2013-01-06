import sys
from argparse import ArgumentParser
from pygibbs.nist_verify import LoadAllEstimators
from pygibbs.thermodynamic_constants import default_pH, default_I, default_T,\
    default_pMg

def MakeOpts(estimators):
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser()
    parser.add_argument("-s", "--thermodynamics_source",
                        choices=estimators.keys(),
                        default="UGC",
                        help="The thermodynamic data to use")
    parser.add_argument("-c", "--csv_out_filename",
                        default="../res/thermo.csv",
                        help="CSV output filename")
    parser.add_argument("-o", "--output_type",
                        choices=['reaction', 'formation', 'chemical'],
                        default="reaction",
                        help="what kind of Gibbs energies to print")
    parser.add_argument("-p", "--ph", type=float,
                        default=default_pH,
                        help="The width of the range of pH values to plot")
    parser.add_argument("-i", "--ionic_strength", type=float,
                        default=default_I,
                        help="The Ionic strength of the solution (in M)")
    parser.add_argument("-m", "--pmg", type=float,
                        default=default_pMg,
                        help="The pMg of the solution")
    parser.add_argument("-t", "--temperature", type=float,
                        default=default_T,
                        help="The T of the solution (in K)")
    return parser

def ExportThermo():
    estimators = LoadAllEstimators()
    options = MakeOpts(estimators).parse_args()
    thermo = estimators[options.thermodynamics_source]

    print 'Thermodynamic source:', thermo.name
    print 'CSV output filename:', options.csv_out_filename

    thermo.SetConditions(pH=options.ph, I=options.ionic_strength,
                         pMg=options.pmg, T=options.temperature)
    if options.output_type == 'reaction':
        thermo.WriteBiochemicalReactionEnergiesToCsv(options.csv_out_filename)
    elif options.output_type == 'formation':
        thermo.WriteBiochemicalFormationEnergiesToCsv(options.csv_out_filename)
    else:
        thermo.WriteChemicalFormationEnergiesToCsv(options.csv_out_filename)
    
if __name__ == "__main__":
    ExportThermo()
