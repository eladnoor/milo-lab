from pygibbs.thermodynamic_estimators import LoadAllEstimators
from pygibbs.nist import Nist
from toolbox.html_writer import HtmlWriter

def main():
    html_writer = HtmlWriter("../res/formation_resolve.html")
    estimators = LoadAllEstimators()
    for name in ['alberty']:
        thermo = estimators[name]
        nist = Nist()
        nist.verify_formation(html_writer=html_writer, 
                              thermodynamics=thermo,
                              name=name)
    html_writer.close()
    
if __name__ == "__main__":
    main()