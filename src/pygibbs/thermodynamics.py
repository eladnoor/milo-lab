import csv
import pseudoisomer
import pylab
import logging
import json

from thermodynamic_constants import default_T, default_pH, default_I, default_pMg
from pygibbs import thermodynamic_constants

class MissingCompoundFormationEnergy(Exception):
    def __init__(self, value, cid=0):
        self.value = value
        self.cid = cid
    def __str__(self):
        return repr(self.value)    

class Thermodynamics(object):
    def __init__(self):
        self.pH = default_pH
        self.I = default_I
        self.T = default_T
        self.pMg = default_pMg
        
        self.c_mid = 1e-3
        self.c_range = (1e-6, 1e-2)
        self.bounds = {}
        self.cid2source_string = {}
        self.anchors = set()

    def cid2PseudoisomerMap(self, cid):
        """
            Given a CID, returns the entire set of pseudoisomers (using class
            PseudoisomerMap).
        """
        raise NotImplementedError

    def get_all_cids(self):
        raise NotImplementedError

    def cid2dG0(self, cid, nH, nMg=0):
        """
            Given a CID and protonation state (number of hydrogens),
            returns the chemical dG0_f (i.e. untransformed Gibbs free energy
            of formation)
            
            One can also provide nMg (number of bound Mg2+ ions),
            but the default is 0
        """
        pmap = self.cid2PseudoisomerMap(cid)
        for p_nH, p_z, p_nMg, dG0 in pmap.ToMatrix():
            if nH == p_nH and nMg == p_nMg:
                return dG0
        raise MissingCompoundFormationEnergy("C%05d doesn't have a formation "
            "energy for the species (nH=%d, nMg=%d)" % (cid, nH, nMg))

    def cid2dG0_tag(self, cid, pH=None, pMg=None, I=None, T=None):
        """
            Input:
                A CID of a compound and the aqueous conditions 
                (pH, I, pMg)
            
            Returns:
                The biochemical dG'0_f (i.e. transformed Gibbs free energy
                of formation)
        """
        pH = pH or self.pH
        I = I or self.I
        T = T or self.T
        pMg = pMg or self.pMg
        return self.cid2PseudoisomerMap(cid).Transform(pH, pMg, I, T, most_abundant=False)
    
    def reaction_to_dG0(self, sparse_reaction, pH=None, pMg=None, I=None, T=None):
        """
            Input:
                A reaction in sparse representation and the aqueous conditions 
                (pH, I, pMg)
            
            Returns:
                The biochemical dG'0_r (i.e. transformed changed in Gibbs free 
                energy of reaction)
        """
        return sum([coeff * self.cid2dG0_tag(cid, pH, pMg, I, T) for 
                    (cid, coeff) in sparse_reaction.iteritems()])
    
    def cid_to_bounds(self, cid, use_default=True):
        (curr_c_min, curr_c_max) = self.bounds.get(cid, (None, None))
        if (curr_c_min == None and use_default):
            curr_c_min = self.c_range[0]
        if (curr_c_max == None and use_default):
            curr_c_max = self.c_range[1]
        return (curr_c_min, curr_c_max)

    @staticmethod
    def pmap_to_table(pmap, pH=default_pH, pMg=default_pMg, I=default_I, T=default_T):
        s = ""
        s += "%2s | %2s | %3s | %7s | %7s\n" % ("nH", "z", "nMg", "dG0_f", "dG0'_f")
        s += "-" * 35 + "\n"
        for nH, z, nMg, dG0 in pmap.ToMatrix():
            dG0_tag = thermodynamic_constants.transform(dG0, nH, nMg, z, pH, pMg, I, T)
            s += "%2d | %2d | %3d | %7.1f | %7.1f\n" % (nH, z, nMg, dG0, dG0_tag)
        return s     

    def display_pmap(self, cid):
        for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
            print "C%05d | %2d | %2d | %3d | %6.2f" % (cid, nH, z, nMg, dG0)
    
    def WriteDataToHtml(self, html_writer, kegg):
        dict_list = []
        for cid in self.get_all_cids():
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                dict = {}
                dict['cid'] = 'C%05d' % cid
                dict['name'] = kegg.cid2name(cid)
                dict['nH'] = '%d' % nH
                dict['z'] = '%d' % z
                dict['nMg'] = '%d' % nMg
                dict['dG0_f'] = '%.2f' % dG0
                if cid in self.anchors:
                    dict['anchor'] = 'yes'
                else:
                    dict['anchor'] = 'no'
                dict_list.append(dict)
        
        html_writer.write_table(dict_list, ['cid', 'name', 'nH', 'z', 
                                            'nMg', 'dG0_f', 'anchor'])
    
    def write_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['cid', 'nH', 'z', 'nMg', 'dG0'])
        for cid in self.get_all_cids():
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                writer.writerow([cid, nH, z, nMg, dG0])

    def write_data_to_json(self, json_fname, kegg):
        formations = []
        for cid in self.get_all_cids():
            h = {}
            h['cid'] = cid
            try:
                h['name'] = kegg.cid2name(h['cid'])
            except KeyError:
                h['name'] = None
            try:
                h['inchi'] = kegg.cid2inchi(h['cid'])
            except KeyError:
                h['inchi'] = None
            h['source'] = self.cid2source_string.get(cid, None)
            h['species'] = []
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                h['species'].append({"nH":nH, "z":z, "nMg":nMg, "dG0_f":dG0})
            formations.append(h)

        json_file = open(json_fname, 'w')
        json_file.write(json.dumps(formations, indent=4))
        json_file.close()
        
    def override_data(self, other):
        for cid in other.get_all_cids():
            self.cid2source_string[cid] = other.cid2source_string.get(cid, None)
            self.cid2pmap_dict[cid] = other.cid2pmap_dict[cid]
            
    def write_transformed_data_to_csv(self, csv_fname):
        writer = csv.writer(open(csv_fname, 'w'))
        writer.writerow(['cid', 'pH', 'pMg', 'I', 'T', 'dG0_tag'])
        for cid in self.get_all_cids():
            dG0_tag = self.cid2dG0_tag(cid)
            writer.writerow([cid, self.pH, self.pMg, self.I, self.T, dG0_tag])
            
    def ToDatabase(self, db, table_name):
        db.CreateTable(table_name, "cid INT, nH INT, z INT, nMg INT, "
                       "dG0_f REAL, anchor BOOL")
        for cid in self.get_all_cids():
            for nH, z, nMg, dG0 in self.cid2PseudoisomerMap(cid).ToMatrix():
                db.Insert(table_name, [cid, nH, z, nMg, dG0, cid in self.anchors])
        db.Commit()

    def FromDatabase(self, db, table_name):
        self.cid2pmap_dict = {}
        self.anchors = set()
        for row in db.DictReader(table_name):
            self.cid2pmap_dict.setdefault(row['cid'], pseudoisomer.PseudoisomerMap())
            self.cid2pmap_dict[row['cid']].Add(row['nH'], row['z'], row['nMg'], 
                                               row['dG0_f'])
            if row['anchor']:
                self.anchors.add(row['cid'])
    
    def write_pseudoisomers_to_html(self, html_writer, kegg, cids):
        # calculate the dG0_f of each compound
        dG0_f = pylab.zeros((len(cids), 1))
        html_writer.write('<table border="1">\n')
        html_writer.write('  ' + '<td>%s</td>'*6 % ("KEGG CID", "Compound Name", 
                                                    "dG0_f [kJ/mol]", "nH", "z", 
                                                    "nMg") + '\n')
        for c, cid in enumerate(cids):
            name = kegg.cid2name(cid)
            try:
                for (nH, z, nMg, dG0) in self.cid2PseudoisomerMap(cid).ToMatrix():
                    html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>%.2f</td><td>%d</td><td>%d</td><td>%d</td></tr>\n' % \
                                      (kegg.cid2link(cid), cid, name, dG0, nH, z, nMg))
                dG0_f[c] = self.cid2dG0_tag(cid)
            
            except MissingCompoundFormationEnergy:
                # this is okay, since it means this compound's dG_f will be unbound, but only if it doesn't appear in the total reaction
                dG0_f[c] = pylab.nan
                html_writer.write('<tr><td><a href="%s">C%05d</a></td><td>%s</td><td>N/A</td><td>N/A</td><td>N/A</td></tr>\n' % \
                                  (kegg.cid2link(cid), cid, name))
        html_writer.write('</table>\n')
        return dG0_f
    
class PsuedoisomerTableThermodynamics(Thermodynamics):
    def __init__(self):
        Thermodynamics.__init__(self)
        self.cid2pmap_dict = {}
        
    @staticmethod
    def FromDatabase(db, table_name):
        """
            Imports the pseudoisomer maps from a CSV file, with these headers:
            'cid', 'nH', 'z', 'nMg', 'dG0'
        """
        thermo = PsuedoisomerTableThermodynamics()
        for row in db.DictReader(table_name):
            if row.get('use for', None) == 'skip':
                continue
            if not row['cid']:
                continue
            thermo.cid2pmap_dict.setdefault(row['cid'], pseudoisomer.PseudoisomerMap())
            thermo.cid2pmap_dict[row['cid']].Add(row['nH'], row['z'], 
                                                 row['nMg'], row['dG0'])
            ref = row.get('ref', '')
            if row['cid'] in thermo.cid2source_string and \
                    thermo.cid2source_string[row['cid']] != ref:
                logging.warning('There are conflicting references for C%05d in '
                                'table %s' % (row['cid'], table_name))
            else:
                thermo.cid2source_string[row['cid']] = ref
            
        return thermo
        
    @staticmethod
    def FromCsvFile(filename):
        """
            Imports the pseudoisomer maps from a CSV file, with these headers:
            'cid', 'nH', 'z', 'nMg', 'dG0'
        """
        thermo = PsuedoisomerTableThermodynamics()
        for row in csv.DictReader(open(filename, 'r')):
            if 'use for' in row and row['use for'] == 'skip':
                continue
            
            cid = int(row['cid'])
            if not cid:
                continue
            
            nH = int(row['nH'])
            z = int(row['z'])
            nMg = int(row['nMg'])
            dG0 = float(row['dG0'])
            thermo.cid2pmap_dict.setdefault(cid, pseudoisomer.PseudoisomerMap())
            thermo.cid2pmap_dict[cid].Add(nH, z, nMg, dG0)
            
            ref = row.get('ref', '')
            if cid in thermo.cid2source_string and thermo.cid2source_string[cid] != ref:
                logging.warning('There are conflicting references for C%05d in '
                                '%s' % (cid, filename))
            else:
                thermo.cid2source_string[cid] = ref
            
        return thermo

    def cid2PseudoisomerMap(self, cid):
        if (cid in self.cid2pmap_dict):
            return self.cid2pmap_dict[cid]
        else:
            raise MissingCompoundFormationEnergy("The compound C%05d does not "
                "have a value for its formation energy of any of its "
                "pseudoisomers" % cid, cid)

    def get_all_cids(self):
        return sorted(self.cid2pmap_dict.keys())

    def test(self):
        pMg_vec = pylab.arange(0, 10, 0.1)
        dG_vec = []
        dG_f_mat = []
        sparse = {20:-1, 13:-1, 147:1, 119:1}
        #sparse = {2:-1, 1:-1, 20:1, 13:1}
        #sparse = {2:-1, 1:-1, 8:1, 9:1}
        for pMg in pMg_vec:
            dG = 0
            dG_f_vec = []
            for cid, coeff in sparse.iteritems():
                dG_f = self.cid2PseudoisomerMap(cid).Transform(pH=7.4, 
                    pMg=pMg, I=0.0, T=303.1)
                dG_f_vec.append(dG_f)
                dG += coeff * dG_f
            dG_f_mat.append(dG_f_vec)
            dG_vec.append(dG)
        pylab.plot(pMg_vec, dG_vec)
        #pylab.plot(pMg_vec, dG_f_mat)
        pylab.show()

if __name__ == "__main__":
    T = PsuedoisomerTableThermodynamics.FromCsvFile(
        '../data/thermodynamics/alberty_pseudoisomers.csv')
    T.test()