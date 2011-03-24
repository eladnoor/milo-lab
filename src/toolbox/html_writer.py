#!/usr/bin/python

"""\
html_writer.py - Construct HTML pages

"""

import os, math, types, openbabel, oasa
import xml.dom.minidom
from toolbox.util import _mkdir, get_current_svn_revision

class BaseHtmlWriter:
    def __init__(self):
        self.div_counter = 0
        pass

    def relative_to_full_path(self, relpath):
        raise Exception("class not implemented")

    def write(self):
        raise Exception("class not implemented")

    def write_header(self):
        self.write('<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n')
        self.write('<head>\n')
        self.write('<script type="text/javascript" src="expandCollapse.js"></script>\n')
        self.write('</head>\n')
        self.write('<html>\n<body>\n')
        r = get_current_svn_revision()
        if r:
            self.write('<a href=https://code.google.com/p/milo-lab/source/browse/trunk/?r=%d>milo-lab SVN r%d</a></br>' % (r,r))

    def write_js(self, path):
        if (os.path.exists(path + '/expandCollapse.js')):
            return
        file = open(path + '/expandCollapse.js', 'w')
        file.write("""function toggleMe(a){
  var e=document.getElementById(a);
  if(!e)return true;
  if(e.style.display=="none"){
    e.style.display="block"
  } else {
    e.style.display="none"
  }
  return true;
}
""")
        file.close()
        
    def write_ol(self, l):
        self.write("<ol>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % str(mem))
        self.write("</ol>\n")
            
    def write_ul(self, l):
        self.write("<ul>\n")
        for mem in l:
            self.write("  <li>%s</li>\n" % str(mem))
        self.write("</ul>\n")
        
    def write_table(self, dict_list, headers=None, border=1):
        def to_string(x):
            if type(x) == types.StringType:
                return x
            else:
                return str(x)
        
        if not headers:
            headers = set()
            for dict in dict_list:
                for key in dict.keys():
                    headers.add(to_string(key))
            headers = sorted(headers)
        
        self.write('<table border=%d>\n' % border)
        self.write('<tr><td><b>' + '</b></td><td><b>'.join(headers) + '</b></td></tr>\n')
        for dict in dict_list:
            values = [to_string(dict.get(key, "")) for key in headers]
            self.write('<tr><td>' + '</td><td>'.join(values) + '</td></tr>\n')
        self.write("</table>\n")
        
    def insert_toggle(self, div_id=None):
        if not div_id:
            div_id = "DIV%05d" % self.div_counter
            self.div_counter += 1
        elif type(div_id) != types.StringType:
            raise ValueError("HTML div ID must be a string")
        self.write('<input type="button" class="button" onclick="return toggleMe(\'%s\')" value="Show">\n' % div_id)
        return div_id
    
    def start_div(self, div_id):
        self.write('<div id="%s" style="display:none">' % div_id)
        
    def end_div(self):
        self.write('</div>\n')
    
    def embed_img(self, fig_fname, alternative_string=""):
        self.write('<img src="' + fig_fname + '" atl="' + alternative_string + '" />')
    
    def embed_svg(self, fig_fname, width=320, height=240, name=''):
        self.write('<object data="%s" type="image/svg+xml" width="%d" height="%d" name="%s"/></object>' % (fig_fname, width, height, name))
    
    def embed_matplotlib_figure(self, fig, width=320, height=240):
        """
            Converts the figure to an SVG DOM and uses the inline SVG option to 
            add it directly into the HTML (without creating a separate SVG file).
        """
        fig.savefig('.svg', format='svg')
        self.extract_svg_from_file('.svg', width, height)
        os.remove('.svg')

    def embed_dot_inline(self, Gdot, width=320, height=240):
        """
            Converts the DOT graph to an SVG DOM and uses the inline SVG option to 
            add it directly into the HTML (without creating a separate SVG file).
        """
        Gdot.write('.svg', prog='dot', format='svg')
        self.extract_svg_from_file('.svg', width, height)
        os.remove('.svg')

    def embed_dot(self, Gdot, name, width=320, height=240):
        """
            Converts the DOT graph to an SVG DOM and uses the inline SVG option to 
            add it directly into the HTML (without creating a separate SVG file).
        """
        svg_filename = self.relative_to_full_path(name + '.svg')
        Gdot.write(svg_filename, prog='dot', format='svg')
        self.embed_svg(svg_filename, width=width, height=height, name=name)

    @staticmethod
    def pybel_mol_to_oasa_mol(mol):
        etab = openbabel.OBElementTable()
        oasa_mol = oasa.molecule()
        for atom in mol.atoms:
            v = oasa_mol.create_vertex()
            v.symbol = etab.GetSymbol(atom.atomicnum)
            v.charge = atom.formalcharge
            oasa_mol.add_vertex(v)

        for bond in openbabel.OBMolBondIter(mol.OBMol):
            e = oasa_mol.create_edge()
            e.order = bond.GetBO()
            if bond.IsHash():
                e.type = "h"
            elif bond.IsWedge():
                e.type = "w"
            oasa_mol.add_edge(bond.GetBeginAtomIdx() - 1,
                         bond.GetEndAtomIdx() - 1,
                         e)
        # I'm sure there's a more elegant way to do the following, but here goes...
        # let's set the stereochemistry around double bonds
        mol.write("can") # Perceive UP/DOWNness
        for bond in openbabel.OBMolBondIter(mol.OBMol):
            ends = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if bond.GetBO() == 2:
                stereobonds = [[b for b in openbabel.OBAtomBondIter(mol.OBMol.GetAtom(x)) if b.GetIdx() != bond.GetIdx() and (b.IsUp() or b.IsDown())]
                               for x in ends]
                if stereobonds[0] and stereobonds[1]: # Needs to be defined at either end
                    if stereobonds[0][0].IsUp() == stereobonds[1][0].IsUp():
                        # Either both up or both down
                        stereo = oasa.stereochemistry.cis_trans_stereochemistry.SAME_SIDE
                    else:
                        stereo = oasa.stereochemistry.cis_trans_stereochemistry.OPPOSITE_SIDE
                    atomids = [(b[0].GetBeginAtomIdx(), b[0].GetEndAtomIdx()) for b in stereobonds]
                    extremes = []
                    for id, end in zip(ends, atomids):
                        if end[0] == id:
                            extremes.append(end[1])
                        else:
                            extremes.append(end[0])
                    center = oasa_mol.get_edge_between(oasa_mol.atoms[ends[0] - 1], oasa_mol.atoms[ends[1] - 1])
                    st = oasa.stereochemistry.cis_trans_stereochemistry(
                              center = center, value = stereo,
                              references = (oasa_mol.atoms[extremes[0] - 1], oasa_mol.atoms[ends[0] - 1],
                                            oasa_mol.atoms[ends[1] - 1], oasa_mol.atoms[extremes[1] - 1]))
                    oasa_mol.add_stereochemistry(st)
        
        oasa_mol.remove_unimportant_hydrogens()
        oasa.coords_generator.calculate_coords(oasa_mol, bond_length=30)
        return oasa_mol
    
    @staticmethod
    def oasa_mol_to_canvas(oasa_mol):
        maxx = max([(v.x or 0) for v in oasa_mol.vertices])
        minx = min([(v.x or 0) for v in oasa_mol.vertices])
        maxy = max([(v.y or 0) for v in oasa_mol.vertices])
        miny = min([(v.y or 0) for v in oasa_mol.vertices])
        maxcoord = max(maxx - minx, maxy - miny)
        fontsize = 16
        bondwidth = 6
        linewidth = 2
        if maxcoord > 270: # 300  - margin * 2
            for v in oasa_mol.vertices:
                v.x *= 270. / maxcoord
                v.y *= 270. / maxcoord
            fontsize *= math.sqrt(270. / maxcoord)
            bondwidth *= math.sqrt(270. / maxcoord)
            linewidth *= math.sqrt(270. / maxcoord)
        canvas = oasa.cairo_out.cairo_out()
        canvas.show_hydrogens_on_hetero = True
        canvas.font_size = fontsize
        canvas.bond_width = bondwidth
        canvas.line_width = linewidth
        return canvas
            
    def embed_molecule_as_png(self, mol, png_filename, width=50, height=50):
        """
            Create a 2D depiction of the molecule and adds it as a PNG into the HTML
            OASA is used for 2D coordinate generation and depiction.
        """
        full_filename = self.relative_to_full_path(png_filename)
        
        dir_path = os.path.abspath(os.path.dirname(full_filename))
        if not os.path.exists(dir_path):
            _mkdir(dir_path)
        
        oasa_mol = HtmlWriter.pybel_mol_to_oasa_mol(mol)
        canvas = HtmlWriter.oasa_mol_to_canvas(oasa_mol)
        canvas.mol_to_cairo(oasa_mol, full_filename, format='png')
        self.embed_img(png_filename, mol.title)

    def embed_molecule_as_svg(self, mol, width=50, height=50):
        """
            Create a 2D depiction of the molecule and adds it as an SVG into the HTML
            OASA is used for 2D coordinate generation and depiction.
        """
        
        oasa_mol = HtmlWriter.pybel_mol_to_oasa_mol(mol)
        canvas = HtmlWriter.oasa_mol_to_canvas(oasa_mol)
        canvas.mol_to_cairo(oasa_mol, '.svg', format='svg')
        self.extract_svg_from_file('.svg', width, height)
        os.remove('.svg')

    def extract_svg_from_xmldom(self, dom, width=320, height=240):
        svg = dom.getElementsByTagName("svg")[0]
        svg.setAttribute('width', '%dpt' % width)
        svg.setAttribute('height', '%dpt' % height)
        self.write(svg.toprettyxml(indent='  ', newl=''))
                
    def extract_svg_from_file(self, fname, width=320, height=240):
        xmldom = xml.dom.minidom.parse(fname)
        self.extract_svg_from_xmldom(xmldom, width, height)
    
    def branch(self, relative_path, link_text=None):
        """
            Branches the HTML file by creating a new HTML and adding a link to it with the desired text
        """
        if (link_text == None):
            link_text = relative_path
            
        self.write("<a href=\"" + relative_path + ".html\">" + link_text + "</a>")
        return HtmlWriter(os.path.join(self.filepath, relative_path + ".html"))
    
    def close(self):
        self.write("</body>\n</html>\n")

class NullHtmlWriter(BaseHtmlWriter):
    def __init__(self):
        BaseHtmlWriter.__init__(self)
        self.filename = None
    
    def write(self, str):
        pass
    
    def relative_to_full_path(self, relpath):
        pass


class HtmlWriter(BaseHtmlWriter):
    def __init__(self, filename, force_path_creation=True, flush_always=True):
        BaseHtmlWriter.__init__(self)
        self.filename = filename
        self.filepath = os.path.dirname(filename)
        self.flush_always = flush_always
        if (not os.path.exists(self.filepath)):
            if (force_path_creation and not os.path.exists(self.filepath)):
                _mkdir(self.filepath)
            else:
                raise Exception("cannot write to HTML file %s since the directory doesn't exist" % filename)
        
        self.file = open(self.filename, "w")
        self.write_header()
        self.write_js(self.filepath)
    
    def relative_to_full_path(self, relpath):
        return self.filepath + "/" + relpath

    def write(self, str):
        if (self.file == None):
            raise Exception("cannot write to this HTML since it is already closed")
        self.file.write(str)
        if (self.flush_always):
            self.file.flush()
    
    def __del__(self):
        if self.file:
            self.close()
        
    def close(self):
        BaseHtmlWriter.close(self)
        self.file.flush()
        self.file.close()
        self.file = None

def test():
    html_write = HtmlWriter("../res/test.html")
    html_write.write("<h1>hello world</h1>\n")
    import pylab
    fig = pylab.figure()
    pylab.plot([1, 2, 3, 4], [4, 3, 2, 1], 'g--')
    html_write.embed_matplotlib_figure(fig, 1000, 1000)

if __name__ == '__main__':
    test()