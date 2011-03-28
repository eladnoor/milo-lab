# Indigo code acquired from http://www.ggasoftware.com/
from indigo import Indigo
from indigo_renderer import IndigoRenderer
import uuid
import rsvg
import gtk

BORDER_WIDTH = 0

# for more rendering options visit:
# http://www.ggasoftware.com/opensource/indigo/api/options#rendering

_indigo = Indigo()
_renderer = IndigoRenderer(_indigo)
_indigo.setOption('render-output-format', 'svg')
_indigo.setOption('render-margins', 10, 10)
_indigo.setOption('render-stereo-style', 'none')
_indigo.setOption('render-implicit-hydrogens-visible', False)
_indigo.setOption('render-coloring', True)
_indigo.setOption('render-bond-length', 20.0)
_indigo.setOption('render-label-mode', 'hetero')

def smiles2svg(smiles, comment=''):
    _indigo.setOption('render-comment', comment)
    m = _indigo.loadMolecule(smiles)
    m.aromatize()
    m.layout()
    s = _renderer.renderToBuffer(m).tostring()
    id = str(uuid.uuid4())
    i = 0
    while True:
        symbol = 'glyph0-%d' % i
        if s.find('id="' + symbol + '"') != -1:
            s = s.replace('id="' + symbol + '"', 'id="' + id + "_" + symbol + '"')
            s = s.replace('href="#' + symbol + '"', 'href="#' + id + "_" + symbol + '"')
        else:
            break
        i += 1
    return s

def display_svg(svg_string):
    def expose_cairo(win, event, svg):
        cr = win.window.cairo_create()
        svg.render_cairo(cr)
        return True
    
    svg = rsvg.Handle(data=svg_string)
    _x, _y, w, h = svg.get_dimension_data()
    win = gtk.Window()
    win.resize(int(w), int(h))
    win.connect("delete-event", lambda w, e: gtk.main_quit())
    win.connect("expose-event", expose_cairo, svg)
    win.show_all()
    win.connect("destroy", lambda w: gtk.main_quit())
    gtk.main()

if __name__ == "__main__":
    s1 = smiles2svg('C(=O)O')
    display_svg(s1)
