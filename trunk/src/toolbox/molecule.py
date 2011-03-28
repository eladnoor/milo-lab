import pybel, indigo, indigo_renderer, uuid, rsvg, gtk
import openbabel
import types

class Molecule(object):

    # for more rendering options visit:
    # http://www.ggasoftware.com/opensource/indigo/api/options#rendering
    _indigo = indigo.Indigo()
    _renderer = indigo_renderer.IndigoRenderer(_indigo)
    _indigo.setOption('render-output-format', 'svg')
    _indigo.setOption('render-margins', 10, 10)
    _indigo.setOption('render-stereo-style', 'none')
    _indigo.setOption('render-implicit-hydrogens-visible', False)
    _indigo.setOption('render-coloring', True)
    _indigo.setOption('render-bond-length', 20.0)
    _indigo.setOption('render-label-mode', 'hetero')
    
    _obConversion = openbabel.OBConversion()
    
    @staticmethod
    def SetBondLength(l):
        Molecule._indigo.setOption('render-bond-length', l)
    
    def __init__(self):
        self.title = None
        self.obmol = openbabel.OBMol()
        self.pybel_mol = None
        self.indigo_mol = None
        self.smiles = ""
        self.inchi = ""
    
    def __str__(self):
        return self.title or self.smiles or self.inchi or ""
        
    def __len__(self):
        return self.GetNumAtoms()
    
    def SetTitle(self, title):
        self.title = title 
    
    @staticmethod
    def FromSmiles(smiles):
        m = Molecule()
        m.smiles = smiles
        Molecule._obConversion.SetInAndOutFormats("smiles", "mol")
        Molecule._obConversion.ReadString(m.obmol, m.smiles)
        m.pybel_mol = pybel.Molecule(m.obmol)
        Molecule._obConversion.SetInAndOutFormats("mol", "inchi")
        m.inchi = Molecule._obConversion.WriteString(m.obmol).strip()
        
        m.indigo_mol = Molecule._indigo.loadMolecule(smiles)
        m.SetTitle(smiles)
        return m
        
    @staticmethod
    def FromInChI(inchi):
        m = Molecule()
        m.inchi = inchi
        Molecule._obConversion.SetInAndOutFormats("inchi", "mol")
        Molecule._obConversion.ReadString(m.obmol, m.inchi)
        m.pybel_mol = pybel.Molecule(m.obmol)
        Molecule._obConversion.SetInAndOutFormats("mol", "smiles")
        m.smiles = Molecule._obConversion.WriteString(m.obmol).strip()
        
        m.indigo_mol = Molecule._indigo.loadMolecule(m.smiles)
        m.SetTitle(inchi)
        return m

    def RemoveHydrogens(self):
        self.pybel_mol.removeh()
    
    def ToOBMol(self):
        return self.obmol
    
    def ToPybelMol(self):
        return self.pybel_mol

    def ToInChI(self):
        return self.inchi
    
    def ToSmiles(self):
        return self.smiles
    
    def GetNumHydrogens(self):
        """
            Returns the number of hydrogen atoms in a compound.
            It is calculated by subtracting the number of heavy atoms (anything bigger than H)
            from the total number of atoms.
        """
        return self.obmol.NumAtoms() - self.obmol.NumHvyAtoms()
    
    def GetTotalCharge(self):
        return self.obmol.GetTotalCharge()
    
    def GetNumElectrons(self):
        """Calculates the number of electrons in a given molecule."""
        atom_bag = {}
        for i in xrange(self.obmol.NumAtoms()):
            atom = self.obmol.GetAtom(i+1)
            atom_bag.setdefault(atom.GetAtomicNum(), 0)
            atom_bag[atom.GetAtomicNum()] += 1
        n_protons = sum([an*cnt for (an, cnt) in atom_bag.iteritems()])
        return n_protons - self.obmol.GetTotalCharge()
    
    def GetNumAtoms(self):
        return self.obmol.NumAtoms()

    def GetAtoms(self):
        return self.pybel_mol.atoms
    
    def FindSmarts(self, smarts):
        """
        Corrects the pyBel version of Smarts.findall() which returns results as tuples,
        with 1-based indices even though Molecule.atoms is 0-based.
    
        Args:
            mol: the molecule to search in.
            smarts_str: the SMARTS query to search for.
        
        Returns:
            The re-mapped list of SMARTS matches.
        """
        if type(smarts) == types.StringType:
            smarts = pybel.Smarts(smarts)
        shift_left = lambda m: [(n - 1) for n in m] 
        return map(shift_left, smarts.findall(self.pybel_mol))

    def ToSVG(self, comment=None):
        if comment:
            Molecule._indigo.setOption('render-comment', comment)
        else:
            Molecule._indigo.setOption('render-comment', '')
        self.indigo_mol.aromatize()
        self.indigo_mol.layout()
        svg_str = Molecule._renderer.renderToBuffer(self.indigo_mol).tostring()
        id = str(uuid.uuid4())
        i = 0
        while True:
            symbol = 'glyph0-%d' % i
            if svg_str.find('id="' + symbol + '"') != -1:
                svg_str = svg_str.replace('id="' + symbol + '"', 
                                          'id="' + id + "_" + symbol + '"')
                svg_str = svg_str.replace('href="#' + symbol + '"', 
                                          'href="#' + id + "_" + symbol + '"')
            else:
                break
            i += 1
        return svg_str
        
    def Draw(self, show_title=False):
        def expose_cairo(win, event, svg):
            cr = win.window.cairo_create()
            svg.render_cairo(cr)
            return True
        
        if show_title:
            svg = rsvg.Handle(data=self.ToSVG(self.title))
        else:
            svg = rsvg.Handle(data=self.ToSVG())
        _x, _y, w, h = svg.get_dimension_data()
        win = gtk.Window()
        win.resize(int(w), int(h))
        win.connect("delete-event", lambda w, e: gtk.main_quit())
        win.connect("expose-event", expose_cairo, svg)
        win.show_all()
        win.connect("destroy", lambda w: gtk.main_quit())
        gtk.main()
        
if __name__ == "__main__":
    Molecule.SetBondLength(50.0)
    #m = Molecule.FromSmiles('CC(=O)O')
    #m.SetTitle('acetate')
    m = Molecule.FromInChI('InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1')
    m.SetTitle('ATP')
    print m.ToSmiles()
    print m.ToInChI()
    obmol = m.ToOBMol()
    print m.GetNumElectrons()
    m.Draw(True)
