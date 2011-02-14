from pygibbs.thermodynamics import CsvFileThermodynamics
from pygibbs.kegg import Kegg
from toolbox.database import SqliteDatabase
from pygibbs.groups import GroupContribution

def main():
    # References:
    # Alberty 2006: Alberty R. A. - Biochemical Thermodynamics: Applications of Mathematica (Methods of Biochemical Analysis), Wiley 2006
    # 
    
    db = SqliteDatabase('../res/gibbs.sqlite')
    kegg = Kegg(db)
    alberty = CsvFileThermodynamics('../data/thermodynamics/dG0.csv')
    
    G = GroupContribution(db=db, kegg=kegg)
    G.init()
    G.override_data(alberty)
    G.write_data_to_json('../res/pseudoisomers.json', kegg)
    
if __name__ == "__main__":
    main()
