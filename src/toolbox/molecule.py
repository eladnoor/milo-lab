import pybel
import indigo, indigo_renderer
import uuid, rsvg, gtk
import openbabel
import types
import re
import glib

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
    _obElements = openbabel.OBElementTable()
    
    @staticmethod
    def GetNumberOfElements():
        return Molecule._obElements.GetNumberOfElements()
    
    @staticmethod
    def GetAllElements():
        return [Molecule._obElements.GetSymbol(i) for i in 
                xrange(Molecule.GetNumberOfElements())]

    @staticmethod
    def GetSymbol(atomic_num):
        return Molecule._obElements.GetSymbol(atomic_num)
            
    @staticmethod
    def GetAtomicNum(elem):
        if type(elem) == types.UnicodeType:
            elem = str(elem)
        return Molecule._obElements.GetAtomicNum(elem)
    
    @staticmethod
    def SetBondLength(l):
        Molecule._indigo.setOption('render-bond-length', l)
    
    @staticmethod
    def VerifySmarts(smarts):
        try:
            pybel.Smarts(smarts)
            return True
        except IOError:
            return False
    
    def __init__(self):
        self.title = None
        self.obmol = openbabel.OBMol()
        self.pybel_mol = None
        self.smiles = None
        self.inchi = None
    
    def __str__(self):
        return self.title or self.smiles or self.inchi or ""
        
    def __len__(self):
        return self.GetNumAtoms()
    
    def Clone(self):
        tmp = Molecule()
        tmp.title = self.title
        tmp.obmol = openbabel.OBMol(self.obmol)
        tmp.pybel_mol = pybel.Molecule(tmp.obmol)
        tmp.smiles = self.smiles
        tmp.inchi = self.inchi
        return tmp
    
    def SetTitle(self, title):
        self.title = title 
    
    @staticmethod
    def FromSmiles(smiles):
        m = Molecule()
        m.smiles = smiles
        Molecule._obConversion.SetInAndOutFormats("smiles", "mol")
        Molecule._obConversion.ReadString(m.obmol, m.smiles)
        m.pybel_mol = pybel.Molecule(m.obmol)
        m.SetTitle(smiles)
        return m
        
    @staticmethod
    def FromInChI(inchi):
        m = Molecule()
        m.inchi = inchi
        Molecule._obConversion.SetInAndOutFormats("inchi", "mol")
        Molecule._obConversion.ReadString(m.obmol, m.inchi)
        m.pybel_mol = pybel.Molecule(m.obmol)
        m.SetTitle(inchi)
        return m
    
    @staticmethod
    def Smiles2InChI(smiles):
        Molecule._obConversion.SetInAndOutFormats("smiles", "inchi")
        obmol = openbabel.OBMol()
        Molecule._obConversion.ReadString(obmol, str(smiles))
        return Molecule._obConversion.WriteString(obmol).strip()

    @staticmethod
    def InChI2Smiles(inchi):
        Molecule._obConversion.SetInAndOutFormats("inchi", "smiles")
        obmol = openbabel.OBMol()
        Molecule._obConversion.ReadString(obmol, inchi)
        return Molecule._obConversion.WriteString(obmol).split()[0]

    def RemoveHydrogens(self):
        self.pybel_mol.removeh()
    
    def RemoveAtoms(self, indices):
        self.obmol.BeginModify()
        for i in sorted(indices, reverse=True):
            self.obmol.DeleteAtom(self.obmol.GetAtom(i+1))
        self.obmol.EndModify()
        self.smiles = None
        self.inchi = None
        
    def SetAtomicNum(self, index, new_atomic_num):
        self.obmol.GetAtom(index+1).SetAtomicNum(new_atomic_num)
        self.smiles = None
        self.inchi = None
        
    def ToOBMol(self):
        return self.obmol
    
    def ToPybelMol(self):
        return self.pybel_mol

    def UpdateInChI(self):
        Molecule._obConversion.SetInAndOutFormats("mol", "inchi")
        self.inchi = Molecule._obConversion.WriteString(self.obmol).strip()

    def ToInChI(self):
        """ 
            Lazy storage of the InChI identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.inchi:
            self.UpdateInChI()
        return self.inchi
    
    def UpdateSmiles(self):
        Molecule._obConversion.SetInAndOutFormats("mol", "smiles")
        self.smiles = Molecule._obConversion.WriteString(self.obmol).strip()
    
    def ToSmiles(self):
        """ 
            Lazy storage of the SMILES identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.smiles:
            self.UpdateSmiles()
        return self.smiles
    
    def GetFormula(self):
        tokens = re.findall('InChI=1S?/([0-9A-Za-z\.]+)/', self.ToInChI())
        if len(tokens) == 1:
            return tokens[0]
        elif len(tokens) > 1:
            raise ValueError('Bad InChI: ' + self.ToInChI())
        else:
            return ''
    
    def GetExactMass(self):
        return self.obmol.GetExactMass()
    
    def GetAtomBagAndCharge(self):
        inchi = self.ToInChI()

        fixed_charge = 0
        for s in re.findall('/q([0-9\+\-]+)', inchi):
            fixed_charge += int(s)

        fixed_protons = 0
        for s in re.findall('/p([0-9\+\-]+)', inchi):
            fixed_protons += int(s)
        
        formula = self.GetFormula()

        atom_bag = {}
        for mol_formula_times in formula.split('.'):
            for times, mol_formula in re.findall('^(\d+)?(\w+)', mol_formula_times):
                if not times:
                    times = 1
                else:
                    times = int(times)
                for atom, count in re.findall("([A-Z][a-z]*)([0-9]*)", mol_formula):
                    if count == '':
                        count = 1
                    else:
                        count = int(count)
                    atom_bag[atom] = atom_bag.get(atom, 0) + count * times
        
        if fixed_protons:
            atom_bag['H'] = atom_bag.get('H', 0) + fixed_protons
            fixed_charge += fixed_protons
        return atom_bag, fixed_charge
        
    def GetHydrogensAndCharge(self):
        atom_bag, charge = self.GetAtomBagAndCharge()
        return atom_bag.get('H', 0), charge
        
    def GetNumElectrons(self):
        """Calculates the number of electrons in a given molecule."""
        atom_bag, fixed_charge = self.GetAtomBagAndCharge()
        n_protons = 0
        for elem, count in atom_bag.iteritems():
            n_protons += count * self._obElements.GetAtomicNum(elem)
        return n_protons - fixed_charge
    
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
        try:
            indigo_mol = Molecule._indigo.loadMolecule(self.ToSmiles())
            indigo_mol.aromatize()
            indigo_mol.layout()
            svg_str = Molecule._renderer.renderToBuffer(indigo_mol).tostring()
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
        except indigo.IndigoException as e:
            return "<b>Indigo error</b>: %s</br>\n" % str(e)
        
    def Draw(self, show_title=False):
        def expose_cairo(win, event, svg):
            cr = win.window.cairo_create()
            svg.render_cairo(cr)
            return True
        
        try:
            if show_title:
                svg = rsvg.Handle(data=self.ToSVG(self.title))
            else:
                svg = rsvg.Handle(data=self.ToSVG())
        except glib.GError: #@UndefinedVariable
            return
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
    #m = Molecule.FromSmiles('CC(=O)[O-]'); m.SetTitle('acetate')
    #m = Molecule.FromSmiles('S[Fe+3]1(S)S[Fe+3](S1)(S)S'); m.SetTitle('oxidized ferredoxin')
    #m = Molecule.FromInChI('InChI=1S/p+1'); m.SetTitle('proton')
    m = Molecule.FromInChI('InChI=1/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1'); m.SetTitle('NAD+')
    #m = Molecule.FromInChI('InChI=1/C5H14NO/c1-6(2,3)4-5-7/h7H,4-5H2,1-3H3/q+1'); m.SetTitle('choline')
    #m = Molecule.FromInChI('InChI=1/CH2O3/c2-1(3)4/h(H2,2,3,4)/p-1'); m.SetTitle('carbonate')
    #m = Molecule.FromInChI('InChI=1/CO2/c2-1-3'); m.SetTitle('CO2')
    #m = Molecule.FromInChI('InChI=1/CO/c1-2'); m.SetTitle('CO')
    #m = Molecule.FromInChI('InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'); m.SetTitle('ATP')
    
    print m.ToSmiles()
    print m.ToInChI()
    obmol = m.ToOBMol()
    print 'atom bag = %s, charge = %d' % m.GetAtomBagAndCharge()
    print 'no. e- =', m.GetNumElectrons()
    print 'nH = %d, charge = %d' % m.GetHydrogensAndCharge()
    print 'no. atoms =', m.GetNumAtoms()
    m.Draw(True)
    
