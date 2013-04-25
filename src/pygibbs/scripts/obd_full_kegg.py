import numpy as np
import logging
from pygibbs.kegg import Kegg
from toolbox.html_writer import HtmlWriter
from pygibbs.obd_dual import KeggPathway
from pygibbs.thermodynamic_estimators import LoadAllEstimators
from argparse import ArgumentParser
from pygibbs.kegg_errors import KeggReactionNotBalancedException
from pygibbs.kegg_parser import ParsedKeggFile
from pygibbs.pathway import PathwayData

def AnalyzeKeggModule(thermo, mid, bounds, html_writer):
    d = {}
    d['OBD [kJ/mol]'] = "N/A"
    
    kegg = Kegg.getInstance()
    S, rids, fluxes, cids = kegg.get_module(mid)
    
    thermo.bounds = bounds.GetOldStyleBounds(cids)
    d['Length'] = len(rids)

    # the S matrix already has the coefficients in the correct direction
    fluxes = [abs(f) for f in fluxes]
    
    for rid in rids:
        r = kegg.rid2reaction(rid)
        try:
            r.Balance(balance_water=True, exception_if_unknown=True)
        except KeggReactionNotBalancedException:
            logging.warning('R%05d is not a balanced reaction, skipping module'
                            % rid)
            return d
        
    
    dG0_r_prime = thermo.GetTransfromedReactionEnergies(S, cids)
    if np.any(np.isnan(dG0_r_prime)):
        logging.warning("Cannot analyze module M%05d because some of the "
                        "Gibbs energies cannot be calculated." % mid)
        return d

    keggpath = KeggPathway(S, rids, fluxes, cids, reaction_energies=dG0_r_prime,
                           cid2bounds=thermo.bounds, c_range=thermo.c_range)
    obd, params = keggpath.FindOBD()
    keggpath.WriteResultsToHtmlTables(html_writer, params['concentrations'],
            params['reaction prices'], params['compound prices'])
    d['OBD [kJ/mol]'] = obd
    return d

def MakeArgParser(estimators):
    """Returns an OptionParser object with all the default options."""
    parser = ArgumentParser(description='A script for running OBD analysis for a gradient of concentrations')
    parser.add_argument('-c', '--config_fname', action='store',
                        required=False, default='../data/obd/full_kegg_config.txt',
                        help='the configuration file for the OBD analysis')
    parser.add_argument('-o', '--output_prefix', action='store',
                        required=False, default='../res/obd_full_kegg',
                        help='the prefix for the output files (*.html and *.csv)')
    parser.add_argument('-s', '--thermodynamics_source', action='store',
                        required=False, default="UGC",
                        choices=estimators.keys(),
                        help="The thermodynamic data to use")
    return parser
    
def main():
    estimators = LoadAllEstimators()
    parser = MakeArgParser(estimators)
    args = parser.parse_args()

    thermo = estimators[args.thermodynamics_source]

    kegg_file = ParsedKeggFile.FromKeggFile(args.config_fname)
    entries = kegg_file.entries()
    if len(entries) == 0:
        raise ValueError('No entries in configuration file')
    entry = 'CONFIGURATION'
    if entry not in entries:
        logging.warning('Configuration file does not contain the entry "CONFIGURATION". '
                        'Using the first entry by default: %s' % entries[0])
        entry = entries[0]
    p_data = PathwayData.FromFieldMap(kegg_file[entry])
    thermo.SetConditions(pH=p_data.pH, I=p_data.I, T=p_data.T, pMg=p_data.pMg)
    thermo.c_range = p_data.c_range
    bounds = p_data.GetBounds()
    
    html_writer = HtmlWriter(args.output_prefix + ".html")

    rowdicts = []
    headers = ['Module', 'Name', 'OBD [kJ/mol]', 'Length']
    kegg = Kegg.getInstance()
    for mid in kegg.get_all_mids():
        html_writer.write('<h2 id=M%05d>M%05d: %s</h2>' %
                          (mid, mid, kegg.get_module_name(mid)))
        try:
            d = AnalyzeKeggModule(thermo, mid, bounds, html_writer)
        except KeyError:
            continue
        d['Module'] = '<a href="#M%05d">M%05d</a>' % (mid, mid)
        d['Name'] = kegg.get_module_name(mid)
        rowdicts.append(d)
    
    rowdicts.sort(key=lambda x:x['OBD [kJ/mol]'])
    html_writer.write_table(rowdicts, headers, decimal=1)
    html_writer.close()
    
if __name__ == "__main__":
    main()