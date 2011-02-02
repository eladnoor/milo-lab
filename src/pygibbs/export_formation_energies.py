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
    alberty.write_data_to_json('../res/alberty.json', kegg)
    
    G = GroupContribution(db)
    G.init()
    G.write_data_to_json('../res/group_contribution.json', kegg)
    
if __name__ == "__main__":
    main()
