import unittest
import pybel
import indigo, indigo_renderer
import uuid, rsvg, gtk #@UnresolvedImport
import openbabel
import types
import re
import glib
from pygibbs.thermodynamic_constants import default_T, default_pH
import sys
import csv

class OpenBabelError(Exception):
    pass

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
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("smiles")
        if not obConversion.ReadString(m.obmol, m.smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        try:
            m.UpdateInChI()
            m.UpdatePybelMol()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from SMILES: " + smiles)
        m.SetTitle(smiles)
        return m
        
    @staticmethod
    def FromInChI(inchi):
        m = Molecule()
        m.inchi = inchi
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("inchi")
        obConversion.ReadString(m.obmol, m.inchi)
        try:
            m.UpdateSmiles()
            m.UpdatePybelMol()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from InChI: " + inchi)
        m.SetTitle(inchi)
        return m
    
    @staticmethod
    def FromMol(mol):
        m = Molecule()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("mol")
        obConversion.ReadString(m.obmol, mol)
        try:
            m.UpdateInChI()
            m.UpdateSmiles()
            m.UpdatePybelMol()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from MOL file:\n" + mol)
        m.SetTitle("")
        return m

    @staticmethod
    def FromOBMol(obmol):
        m = Molecule()
        m.obmol = obmol
        try:
            m.UpdateInChI()
            m.UpdateSmiles()
            m.UpdatePybelMol()
        except OpenBabelError:
            raise OpenBabelError("Failed to create Molecule from OBMol")
        m.SetTitle("")
        return m
    
    @staticmethod
    def _FromFormat(s, fmt='inchi'):
        if fmt == 'smiles' or fmt == 'smi':
            return Molecule.FromSmiles(s)
        if fmt == 'inchi':
            return Molecule.FromInChI(s)
        if fmt == 'mol':
            return Molecule.FromMol(s)
        if fmt == 'obmol':
            return Molecule.FromOBMol(s)
    
    @staticmethod
    def _ToFormat(obmol, fmt='inchi'):
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat(fmt)
        res = obConversion.WriteString(obmol)
        if not res:
            raise OpenBabelError("Cannot convert OBMol to %s" % fmt)
        if fmt == 'smiles' or fmt == 'smi':
            res = res.split()
            if res == []:
                raise OpenBabelError("Cannot convert OBMol to %s" % fmt)
            else:
                return res[0]
        elif fmt == 'inchi':
            return res.strip()
        else:
            return res
        
    @staticmethod
    def Smiles2InChI(smiles):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smiles", "inchi")
        obmol = openbabel.OBMol()
        if not obConversion.ReadString(obmol, smiles):
            raise OpenBabelError("Cannot read the SMILES string: " + smiles)
        return obConversion.WriteString(obmol).strip()

    @staticmethod
    def InChI2Smiles(inchi):
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("inchi", "smiles")
        obmol = openbabel.OBMol()
        if not obConversion.ReadString(obmol, inchi):
            raise OpenBabelError("Cannot read the InChI string: " + inchi)
        return obConversion.WriteString(obmol).split()[0]
        
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

    def ToFormat(self, fmt='inchi'):
        return Molecule._ToFormat(self.obmol, fmt=fmt)
    
    def ToMolfile(self):
        return self.ToFormat('mol')

    def UpdateInChI(self):
        self.inchi = Molecule._ToFormat(self.obmol, 'inchi')

    def ToInChI(self):
        """ 
            Lazy storage of the InChI identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.inchi:
            self.UpdateInChI()
        return self.inchi
    
    def UpdateSmiles(self):
        self.smiles = Molecule._ToFormat(self.obmol, 'smiles')
    
    def ToSmiles(self):
        """ 
            Lazy storage of the SMILES identifier (calculate once only when 
            asked for and store for later use).
        """
        if not self.smiles:
            self.UpdateSmiles()
        return self.smiles
    
    def UpdatePybelMol(self):
        self.pybel_mol = pybel.Molecule(self.obmol)
    
    @staticmethod
    def _GetFormulaFromInChI(inchi):
        tokens = re.findall('/f([0-9A-Za-z\.]+/)', inchi)
        if len(tokens) == 0:
            tokens = re.findall('InChI=1S?/([0-9A-Za-z\.]+)', inchi)

        if len(tokens) == 1:
            return tokens[0]
        elif len(tokens) > 1:
            raise ValueError('Bad InChI: ' + inchi)
        else:
            return ''

    @staticmethod
    def _GetAtomBagAndChargeFromInChI(inchi):
        fixed_charge = 0
        for q in re.findall('/q([0-9\+\-\;]+)', inchi):
            for s in q.split(';'): 
                if s:
                    fixed_charge += int(s)

        fixed_protons = 0
        for p in re.findall('/p([0-9\+\-\;]+)', inchi):
            for s in p.split(';'):
                if s:
                    fixed_protons += int(s)
        
        formula = Molecule._GetFormulaFromInChI(inchi)

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

    @staticmethod
    def _GetNumElectronsFromInChI(inchi):
        """Calculates the number of electrons in a given molecule."""
        atom_bag, fixed_charge = Molecule._GetAtomBagAndChargeFromInChI(inchi)
        n_protons = 0
        for elem, count in atom_bag.iteritems():
            n_protons += count * Molecule._obElements.GetAtomicNum(elem)
        return n_protons - fixed_charge

    def GetFormula(self):
        return Molecule._GetFormulaFromInChI(self.ToInChI())
    
    def GetExactMass(self):
        return self.obmol.GetExactMass()
    
    def GetAtomBagAndCharge(self):
        return Molecule.__GetAtomBagAndChargeFromInChI(self.ToInChI())
        
    def GetHydrogensAndCharge(self):
        atom_bag, charge = self.GetAtomBagAndCharge()
        return atom_bag.get('H', 0), charge
        
    def GetNumElectrons(self):
        return Molecule._GetNumElectronsFromInChI(self.ToInChI())
    
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

    def GetAtomCharges(self):
        """
            Returns:
                A list of charges, according to the number of atoms
                in the molecule
        """
        return [atom.formalcharge for atom in self.pybel_mol.atoms]

    @staticmethod
    def _GetDissociationTable(molstring, fmt='inchi', mid_pH=default_pH, 
                              min_pKa=0, max_pKa=14, T=default_T):
        """
            Returns the relative potentials of pseudoisomers,
            relative to the most abundant one at pH 7.
        """
        from pygibbs.dissociation_constants import DissociationTable
        from toolbox import chemaxon

        diss_table = DissociationTable()
        try:
            pKa_table, major_ms = chemaxon.GetDissociationConstants(molstring, 
                                                                    mid_pH=mid_pH)

            mol = Molecule.FromSmiles(major_ms)
            nH, z = mol.GetHydrogensAndCharge()
            diss_table.SetMolString(nH, nMg=0, s=major_ms)
            diss_table.SetCharge(nH, z, nMg=0)
            
            pKa_higher = [x for x in pKa_table if mid_pH < x[0] < max_pKa]
            pKa_lower = [x for x in pKa_table if mid_pH > x[0] > min_pKa]
            for i, (pKa, _, smiles_above) in enumerate(sorted(pKa_higher)):
                diss_table.AddpKa(pKa, nH_below=(nH-i), nH_above=(nH-i-1),
                                  nMg=0, ref='ChemAxon', T=T)
                diss_table.SetMolString((nH-i-1), nMg=0, s=smiles_above)
    
            for i, (pKa, smiles_below, _) in enumerate(sorted(pKa_lower, reverse=True)):
                diss_table.AddpKa(pKa, nH_below=(nH+i+1), nH_above=(nH+i),
                                  nMg=0, ref='ChemAxon', T=T)
                diss_table.SetMolString((nH+i+1), nMg=0, s=smiles_below)
        except chemaxon.ChemAxonError:
            mol = Molecule._FromFormat(molstring, fmt)
            diss_table.SetOnlyPseudoisomerMolecule(mol)
            
        return diss_table

    def GetDissociationTable(self, fmt='inchi', mid_pH=default_pH, 
                           min_pKa=0, max_pKa=14, T=default_T):
        """
            Returns the relative potentials of pseudoisomers,
            relative to the most abundant one at pH 7.
        """
        
        return Molecule._GetDissociationTable(self.ToInChI(), 'inchi',
                                            mid_pH, min_pKa, max_pKa, T)

class MoleculeTest(unittest.TestCase):
    
    def setUp(self):
        # pairs of (inchi, number of electrons)
        if False:
            self.inchis = [
                  {'inchi':"InChI=1S/C34H34N4O4.Fe/c1-7-21-17(3)25-13-26-19(5)23(9-11-33(39)40)31(37-26)16-32-24(10-12-34(41)42)20(6)28(38-32)15-30-22(8-2)18(4)27(36-30)14-29(21)35-25;/h7-8,13-16H,1-2,9-12H2,3-6H3,(H4,35,36,37,38,39,40,41,42);/q;+2/p-2/b25-13-,26-13-,27-14-,28-15-,29-14-,30-15-,31-16-,32-16-;", 'n_e':322},
                  {'inchi':"InChI=1/Mn/q+2", 'n_e':23},
                  {'inchi':"InChI=1S/Cu/q+2", 'n_e':27},
                  {'inchi':"InChI=1S/C62H90N13O14P.C10H12N5O3.Co/c1-29-20-39-40(21-30(29)2)75(28-70-39)57-52(84)53(41(27-76)87-57)89-90(85,86)88-31(3)26-69-49(83)18-19-59(8)37(22-46(66)80)56-62(11)61(10,25-48(68)82)36(14-17-45(65)79)51(74-62)33(5)55-60(9,24-47(67)81)34(12-15-43(63)77)38(71-55)23-42-58(6,7)35(13-16-44(64)78)50(72-42)32(4)54(59)73-56;1-4-6(16)7(17)10(18-4)15-3-14-5-8(11)12-2-13-9(5)15;/h20-21,23,28,31,34-37,41,52-53,56-57,76,84H,12-19,22,24-27H2,1-11H3,(H15,63,64,65,66,67,68,69,71,72,73,74,77,78,79,80,81,82,83,85,86);2-4,6-7,10,16-17H,1H2,(H2,11,12,13);/q;;+2/p-2/t31-,34-,35-,36-,37+,41-,52?,53-,56?,57+,59-,60+,61+,62+;4-,6-,7-,10-;/m11./s1", 'n_e':836},
                  {'inchi':"InChI=1S/Ni/q+2", 'n_e':26},
                  {'inchi':"InChI=1S/C62H90N13O14P.Co/c1-29-20-39-40(21-30(29)2)75(28-70-39)57-52(84)53(41(27-76)87-57)89-90(85,86)88-31(3)26-69-49(83)18-19-59(8)37(22-46(66)80)56-62(11)61(10,25-48(68)82)36(14-17-45(65)79)51(74-62)33(5)55-60(9,24-47(67)81)34(12-15-43(63)77)38(71-55)23-42-58(6,7)35(13-16-44(64)78)50(72-42)32(4)54(59)73-56;/h20-21,23,28,31,34-37,41,52-53,56-57,76,84H,12-19,22,24-27H2,1-11H3,(H15,63,64,65,66,67,68,69,71,72,73,74,77,78,79,80,81,82,83,85,86);/q;+2/p-2/t31-,34-,35-,36-,37+,41-,52?,53-,56?,57+,59-,60+,61+,62+;/m1./s1", 'n_e':705},
                  {'inchi':"InChI=1S/C42H46N4O16.Fe/c1-41(17-39(59)60)23(5-9-35(51)52)29-14-27-21(11-37(55)56)19(3-7-33(47)48)25(43-27)13-26-20(4-8-34(49)50)22(12-38(57)58)28(44-26)15-31-42(2,18-40(61)62)24(6-10-36(53)54)30(46-31)16-32(41)45-29;/h13-16,23-24H,3-12,17-18H2,1-2H3,(H10,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62);/q;+2/p-2/t23-,24-,41+,42+;/m1./s1", 'n_e':478},
                  {'inchi':"InChI=1S/C62H90N13O14P.Co/c1-29-20-39-40(21-30(29)2)75(28-70-39)57-52(84)53(41(27-76)87-57)89-90(85,86)88-31(3)26-69-49(83)18-19-59(8)37(22-46(66)80)56-62(11)61(10,25-48(68)82)36(14-17-45(65)79)51(74-62)33(5)55-60(9,24-47(67)81)34(12-15-43(63)77)38(71-55)23-42-58(6,7)35(13-16-44(64)78)50(72-42)32(4)54(59)73-56;/h20-21,23,28,31,34-37,41,52-53,56-57,76,84H,12-19,22,24-27H2,1-11H3,(H15,63,64,65,66,67,68,69,71,72,73,74,77,78,79,80,81,82,83,85,86);/q;+1/p-1/t31-,34-,35-,36-,37+,41-,52?,53-,56?,57+,59-,60+,61+,62+;/m1./s1", 'n_e':706},
                  {'inchi':"InChI=1S/Cd/q+2", 'n_e':46},
                  {'inchi':"InChI=1S/C35H36N4O5.Mg/c1-8-19-15(3)22-12-24-17(5)21(10-11-28(40)41)32(38-24)30-31(35(43)44-7)34(42)29-18(6)25(39-33(29)30)14-27-20(9-2)16(4)23(37-27)13-26(19)36-22;/h8,12-14,17,21,31H,1,9-11H2,2-7H3,(H3,36,37,38,39,40,41,42);/q;+2/p-2/t17-,21-,31+;/m0./s1", 'n_e':324},
                  {'inchi':"InChI=1S/C27H42N7O20P3S/c1-27(2,22(41)25(42)30-6-5-16(36)29-7-8-58-18(39)9-14(35)3-4-17(37)38)11-51-57(48,49)54-56(46,47)50-10-15-21(53-55(43,44)45)20(40)26(52-15)34-13-33-19-23(28)31-12-32-24(19)34/h12-13,15,20-22,26,40-41H,3-11H2,1-2H3,(H,29,36)(H,30,42)(H,37,38)(H,46,47)(H,48,49)(H2,28,31,32)(H2,43,44,45)/t15-,20-,21-,22+,26-/m1/s1", 'n_e':474},
                  {'inchi':"InChI=1S/C35H34N4O5.Mg/c1-8-19-15(3)22-12-24-17(5)21(10-11-28(40)41)32(38-24)30-31(35(43)44-7)34(42)29-18(6)25(39-33(29)30)14-27-20(9-2)16(4)23(37-27)13-26(19)36-22;/h8,12-14,31H,1,9-11H2,2-7H3,(H3,36,37,38,39,40,41,42);/q;+2/p-2/b22-12-,23-13-,24-12-,25-14-,26-13-,27-14-,32-30-;/t31-;/m1./s1", 'n_e':322},
                  {'inchi':"InChI=1S/C34H34N4O4.Mg/c1-7-21-17(3)25-13-26-19(5)23(9-11-33(39)40)31(37-26)16-32-24(10-12-34(41)42)20(6)28(38-32)15-30-22(8-2)18(4)27(36-30)14-29(21)35-25;/h7-8,13-16H,1-2,9-12H2,3-6H3,(H4,35,36,37,38,39,40,41,42);/q;+2/p-2/b25-13-,26-13-,27-14-,28-15-,29-14-,30-15-,31-16-,32-16-;", 'n_e':308},
                  {'inchi':"InChI=1/C5H10O4/c1-3(7)5(9)4(8)2-6/h2-5,7-9H,1H3/t3-,4+,5-/m1/s1", 'n_e':72},
                  {'inchi':"InChI=1S/C39H60O4/Dummy", "n_e":326},
                  {'inchi':"InChI=1/C10H16N4O4/c11-6(9(15)16)1-2-8-13-4-5(14-8)3-7(12)10(17)18/h4,6-7H,1-3,11-12H2,(H,13,14)(H,15,16)(H,17,18)/t6?,7-/m0/s1/f/h13,15,17H", 'n_e':None}]
        else:
            self.inchis = csv.DictReader(open('../data/examples/inchi.csv'))

    def testNumElectrons(self):    
        for line_num, d in enumerate(self.inchis):
            n_e_inchi = Molecule._GetNumElectronsFromInChI(d['inchi'])
            try:
                if d['n_e']:
                    self.assertEquals(n_e_inchi, int(d['n_e']))
            except AssertionError as e:
                sys.stderr.write("mismatch in line %d: %s\n" % (line_num+2, d['inchi']))
                raise e

def Suite():
    suites = (unittest.makeSuite(MoleculeTest, 'test'),)
    return unittest.TestSuite(suites)

if __name__ == "__main__":
    unittest.main()
    sys.exit(0)
    from toolbox.html_writer import HtmlWriter

    Molecule.SetBondLength(50.0)

    #m = Molecule.FromInChI('InChI=1S/Fe') # Iron
    #m = Molecule.FromSmiles('CC(O)=O'); m.SetTitle('acetate')
    #m = Molecule.FromSmiles('S[Fe+3]1(S)S[Fe+3](S1)(S)S'); m.SetTitle('oxidized ferredoxin')
    #m = Molecule.FromInChI('InChI=1S/p+1'); m.SetTitle('proton')
    #m = Molecule.FromInChI('InChI=1/C21H27N7O14P2/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(32)14(30)11(41-21)6-39-44(36,37)42-43(34,35)38-5-10-13(29)15(31)20(40-10)27-3-1-2-9(4-27)18(23)33/h1-4,7-8,10-11,13-16,20-21,29-32H,5-6H2,(H5-,22,23,24,25,33,34,35,36,37)/p+1/t10-,11-,13-,14-,15-,16-,20-,21-/m1/s1'); m.SetTitle('NAD+')
    #m = Molecule.FromInChI('InChI=1/C5H14NO/c1-6(2,3)4-5-7/h7H,4-5H2,1-3H3/q+1'); m.SetTitle('choline')
    #m = Molecule.FromInChI('InChI=1/CH2O3/c2-1(3)4/h(H2,2,3,4)/p-1'); m.SetTitle('carbonate')
    #m = Molecule.FromInChI('InChI=1/CO2/c2-1-3'); m.SetTitle('CO2')
    #m = Molecule.FromInChI('InChI=1/CO/c1-2'); m.SetTitle('CO')
    #m = Molecule.FromInChI('InChI=1/C10H16N5O13P3/c11-8-5-9(13-2-12-8)15(3-14-5)10-7(17)6(16)4(26-10)1-25-30(21,22)28-31(23,24)27-29(18,19)20/h2-4,6-7,10,16-17H,1H2,(H,21,22)(H,23,24)(H2,11,12,13)(H2,18,19,20)/t4-,6-,7-,10-/m1/s1'); m.SetTitle('ATP')
    #m = Molecule.FromSmiles("P(=O)(O)(O)O")
    
    #print m.ToFormat('mol')
    #print m.ToFormat('mol2')
    #print m.ToFormat('smi')
    #print m.ToFormat('inchi')
    #print m.ToFormat('sdf')

    diss_table = Molecule._GetDissociationTable('C(=O)(O)CN', fmt='smiles',
                 mid_pH=default_pH, min_pKa=0, max_pKa=14, T=default_T)
    print "glycine\n", diss_table
    
    html_writer = HtmlWriter('../res/molecule.html')
    from pygibbs.kegg import Kegg
    kegg = Kegg.getInstance()
    html_writer.write('<h1>pKa estimation using ChemAxon</h1>\n')
    for cid in [41]:
        m = kegg.cid2mol(cid)
        html_writer.write("<h2>C%05d : %s</h2>\n" % (cid, str(m)))
        diss_table = m.GetDissociationTable()
        pmap = diss_table.GetPseudoisomerMap()
        diss_table.WriteToHTML(html_writer)
        pmap.WriteToHTML(html_writer)
        html_writer.write("</p>\n")
        #print m.GetDissociationConstants()
        #print m.GetMacrospecies()

    #obmol = m.ToOBMol()
    #print 'atom bag = %s, charge = %d' % m.GetAtomBagAndCharge()
    #print 'no. e- =', m.GetNumElectrons()
    #print 'nH = %d, charge = %d' % m.GetHydrogensAndCharge()
    #print 'no. atoms =', m.GetNumAtoms()
    #m.Draw(True)
    
