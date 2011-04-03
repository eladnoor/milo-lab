import gzip, sqlite3, os, sys, csv, re, urllib, zipfile
import libsbml # obtained from http://sbml.org/Software/libSBML (version 4.1.0)
from pygibbs.groups import GroupContribution
from toolbox.sql import write_csv_table
from toolbox.molecule import Molecule

def connect_db(fname='../res/gibbs.sqlite'):
    try:
        os.mkdir('../res')
    except OSError:
        pass
    
    comm = sqlite3.connect(fname)
    return comm

def parse_ChEBI(cursor):
    sys.stderr.write("Parsing the chebiId_inchi file ... ")
    cursor.execute("DROP TABLE IF EXISTS chebiId_inchi")
    cursor.execute("CREATE TABLE chebiId_inchi (chebi_id INT, inchi TEXT)")
    cursor.execute("DROP INDEX IF EXISTS chebiId_inchi_idx")
    cursor.execute("CREATE UNIQUE INDEX chebiId_inchi_idx ON chebiId_inchi (chebi_id);")
    
    fname = '../chebi/chebiId_inchi.tsv.zip'
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chebiId_inchi.tsv.zip'
    if (not os.path.exists(fname)):
        sys.stderr.write("Downloading chebiId_inchi from the ChEBI FTP site ... ")
        urllib.urlretrieve(url, fname)
        
    chebi_file = zipfile.ZipFile(fname, 'r')
    csv_reader = csv.reader(chebi_file.open('chebiId_inchi.tsv'), delimiter='\t')
    csv_reader.next() # skip title row
    for row in csv_reader:
        (chebi_id, inchi) = row
        cursor.execute("INSERT INTO chebiId_inchi VALUES(?,?)", (int(chebi_id), inchi))
    sys.stderr.write("[DONE]\n")
    chebi_file.close()

def parse_MIRIAM_annotation(xml):
    if (xml == None):
        return []
    return re.findall('rdf:resource=\"urn:miriam:([^:^\"]+):([^:^\"]+)\"', xml.toXMLString())

def parse_SBML(cursor):
    sys.stderr.write("Parsing the SBML model file ... ")
    fname = '../data/thermodynamics/yeast_4.02.xml'
    if (not os.path.exists(fname)):
        sys.stderr.write("Please download the yeast SBML from http://www.comp-sys-bio.org/yeastnet/ and place in at: %s\n" % fname)
        sys.exit(-1)
    yeast_SBML = libsbml.readSBML(fname)
    yeast_model = yeast_SBML.getModel()
    
    cursor.execute("DROP TABLE IF EXISTS yeast_reaction")
    cursor.execute("CREATE TABLE yeast_reaction (reaction_id TEXT, name TEXT)")
    
    cursor.execute("DROP INDEX IF EXISTS yeast_reaction_idx")
    cursor.execute("CREATE UNIQUE INDEX yeast_reaction_idx ON yeast_reaction (reaction_id);")

    cursor.execute("DROP TABLE IF EXISTS yeast_reaction2species")
    cursor.execute("CREATE TABLE yeast_reaction2species (reaction_id TEXT, species_id TEXT, stoichiometry REAL)")

    cursor.execute("DROP INDEX IF EXISTS yeast_reaction2species_idx1")
    cursor.execute("CREATE INDEX yeast_reaction2species_idx1 ON yeast_reaction2species (reaction_id);")
    
    cursor.execute("DROP INDEX IF EXISTS yeast_reaction2species_idx2")
    cursor.execute("CREATE INDEX yeast_reaction2species_idx2 ON yeast_reaction2species (species_id);")
    for reaction in yeast_model.getListOfReactions():
        cursor.execute("INSERT INTO yeast_reaction VALUES(?,?)", [reaction.getId(), reaction.getName()])
        for sp_ref in reaction.getListOfReactants():
            cursor.execute("INSERT INTO yeast_reaction2species VALUES(?,?,?)", \
                           [reaction.getId(), sp_ref.getSpecies(), -sp_ref.getStoichiometry()])
        for sp_ref in reaction.getListOfProducts():
            cursor.execute("INSERT INTO yeast_reaction2species VALUES(?,?,?)", \
                           [reaction.getId(), sp_ref.getSpecies(), sp_ref.getStoichiometry()])
    
    cursor.execute("DROP TABLE IF EXISTS yeast_species")
    cursor.execute("CREATE TABLE yeast_species (species_id TEXT, name TEXT, compartment_id TEXT, species_type_id TEXT)")
    cursor.execute("DROP INDEX IF EXISTS yeast_species_idx")
    cursor.execute("CREATE UNIQUE INDEX yeast_species_idx ON yeast_species (species_id);")
    for species in yeast_model.getListOfSpecies():
        cursor.execute("INSERT INTO yeast_species VALUES(?,?,?,?)", \
                       [species.getId(), species.getName(), species.getCompartment(), species.getSpeciesType()])

    cursor.execute("DROP TABLE IF EXISTS yeast_species_type")
    cursor.execute("CREATE TABLE yeast_species_type (species_type_id TEXT, name TEXT, ann_name TEXT, ann_value TEXT)")
    cursor.execute("DROP INDEX IF EXISTS yeast_species_type_idx")
    cursor.execute("CREATE INDEX yeast_species_type_idx ON yeast_species_type (species_type_id);")

    cursor.execute("DROP TABLE IF EXISTS yeast_species_type2chebi")
    cursor.execute("CREATE TABLE yeast_species_type2chebi (species_type_id TEXT, chebi_id INT)")
    cursor.execute("DROP INDEX IF EXISTS yeast_species_type2chebi_idx1")
    cursor.execute("CREATE UNIQUE INDEX yeast_species_type2chebi_idx1 ON yeast_species_type2chebi (species_type_id);")
    cursor.execute("DROP INDEX IF EXISTS yeast_species_type2chebi_idx2")
    cursor.execute("CREATE INDEX yeast_species_type2chebi_idx2 ON yeast_species_type2chebi (chebi_id);")
    for type in yeast_model.getListOfSpeciesTypes():
        for (name, value) in parse_MIRIAM_annotation(type.getAnnotation()):
            cursor.execute("INSERT INTO yeast_species_type VALUES(?,?,?,?)", \
                           [type.getId(), type.getName(), name, value])
            if (name == 'obo.chebi'):
                for chebi_id in re.findall('CHEBI%3A(\d+)', value):
                    cursor.execute("INSERT INTO yeast_species_type2chebi VALUES(?,?)", \
                                   [type.getId(), chebi_id])
    
    cursor.execute("DROP TABLE IF EXISTS yeast_compartment")
    cursor.execute("CREATE TABLE yeast_compartment (compartment_id TEXT, name TEXT)")
    cursor.execute("DROP INDEX IF EXISTS yeast_compartment_idx")
    cursor.execute("CREATE UNIQUE INDEX yeast_compartment_idx ON yeast_compartment (compartment_id);")
    for compartment in yeast_model.getListOfCompartments():
        cursor.execute("INSERT INTO yeast_compartment VALUES(?,?)", \
                       [compartment.getId(), compartment.getName()])
    
            
    sys.stderr.write("[DONE]\n")
    
def create_species2inchi(cursor):
    """
        This table is actually only used for determining which of the Chebi compounds need to be used.
        Since the Chebi database is much larger than the yeast model, it would take a very long time to calculate
        the formation energies for all its compounds, and we really need only the ones that are in the yeast DB.
    """
    
    cursor.execute("DROP TABLE IF EXISTS yeast_species2inchi")
    cursor.execute("CREATE TABLE yeast_species2inchi (species_id TEXT, inchi TEXT)")
    cursor.execute("DROP INDEX IF EXISTS yeast_species2inchi_idx")
    cursor.execute("CREATE INDEX yeast_species2inchi_idx ON yeast_species2inchi (species_id);")
    cursor.execute("INSERT INTO yeast_species2inchi SELECT S.species_id, C.inchi " \
                   "FROM yeast_species S, yeast_species_type2chebi T, chebiId_inchi C " \
                   "WHERE S.species_type_id = T.species_type_id AND T.chebi_id = C.chebi_id")

def add_thermodynamics(cursor):
    from groups import GroupMissingTrainDataError, GroupDecompositionError

    gc = GroupContribution(sqlite_name="gibbs.sqlite", html_name="pathologic")
    gc.init()
    
    cursor.execute("DROP TABLE IF EXISTS yeast_inchi2thermo")
    cursor.execute("CREATE TABLE yeast_inchi2thermo (inchi TEXT, charge INT, nH INT, dG0_f REAL)")
    cursor.execute("DROP INDEX IF EXISTS yeast_inchi2thermo_idx")
    cursor.execute("CREATE INDEX yeast_inchi2thermo_idx ON yeast_inchi2thermo (inchi);")
    
    inchi_list = []
    for row in cursor.execute("SELECT distinct(inchi) " \
                              "FROM yeast_species2inchi WHERE inchi IS NOT NULL"):
        inchi = row[0]
        inchi_list.append(str(inchi))
    
    for inchi in inchi_list:
        try:
            mol = Molecule.FromInChI(str(inchi))
            pmap = gc.Mol2PseudoisomerMap(mol)
            for ((z, nH), dG0) in pmap.iteritems():
                cursor.execute("INSERT INTO yeast_inchi2thermo VALUES(?,?,?,?)", [inchi, z, nH, dG0])
        except (IOError, GroupMissingTrainDataError, GroupDecompositionError):
            sys.stderr.write("Cannot convert the following InChI to a pybel Molecule")
            #cursor.execute("INSERT INTO yeast_inchi2thermo VALUES(?,?,?,?)", [inchi, None, None, None])    

def get_model(cursor, pH=7.0, I=0.2, T=300):
    """
        returns a tuple (species, reactions)
        species is a list of 3-tuples, which contains (species ID, compartment ID, dG0 formation)
        reactions is a list of 2-tuples, containing (reaction ID, reaction name)
    """
    inchi2pmap = {}
    for row in cursor.execute("SELECT inchi, charge, nH, dG0_f FROM yeast_inchi2thermo"):
        (inchi, z, nH, dG0) = row
        inchi2pmap.setdefault(inchi, {})[(z, nH)] = dG0
    
    species = []
    sp2column = {}
    for row in cursor.execute("SELECT S.species_id, S.compartment_id, I.inchi " \
                              "FROM yeast_species S LEFT OUTER JOIN yeast_species2inchi I " \
                              "ON S.species_id = I.species_id " \
                              "ORDER BY S.species_id"):
        (species_id, compartment_id, inchi) = row
        species_id = str(species_id)
        if (inchi != None and inchi in inchi2pmap):
            dG0 = inchi2pmap[inchi].Transform(pH=pH, I=I, T=T)
        else:
            dG0 = None
        species.append((species_id, compartment_id, dG0))
        sp2column[species_id] = len(species) - 1

    reaction2sparse = {}
    for row in cursor.execute("SELECT reaction_id, species_id, stoichiometry FROM yeast_reaction2species"):
        (reaction_id, species_id, coeff) = row
        reaction_id = str(reaction_id)
        species_id = str(species_id)
        
        if (reaction_id not in reaction2sparse):
            reaction2sparse[reaction_id] = {}
        reaction2sparse[reaction_id][species_id] = coeff
    
    reactions = sorted(reaction2sparse.iteritems())

    return (species, reactions)
    
def get_name2rid(cursor):
    name2rid = {}
    for row in cursor.execute("SELECT reaction_id, name FROM yeast_reaction"):
        (rid, name) = row
        name2rid[name] = rid
    return name2rid

def export_tables(cursor):
    if (not os.path.exists('../res')):
        os.mkdir('../res')
    if (not os.path.exists('../res/yeast')):
        os.mkdir('../res/yeast')
        
    write_csv_table(cursor, '../res/yeast/chebiId_inchi.csv', 'chebiId_inchi')

    tables_to_write = ['inchi2thermo', 'reaction', 'reaction2species', 'species', 'species2inchi', \
                       'compartment', 'species_type', 'species_type2chebi', 'reaction2species']
    for tablename in tables_to_write:
        write_csv_table(cursor, '../res/yeast/%s.csv' % tablename, 'yeast_%s' % tablename)
        
if (__name__ == "__main__"):
    comm = connect_db()
    cursor = comm.cursor()
    parse_ChEBI(cursor); comm.commit()
    parse_SBML(cursor); comm.commit()
    create_species2inchi(cursor); comm.commit()
    add_thermodynamics(cursor); comm.commit()
    export_tables(cursor)