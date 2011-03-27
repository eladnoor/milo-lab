# Indigo code acquired from http://www.ggasoftware.com/
from indigo import Indigo
from indigo_renderer import IndigoRenderer
from toolbox import html_writer
import random

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

def smiles2svg(smiles, comment=''):
    _indigo.setOption('render-comment', comment)
    m = _indigo.loadMolecule(smiles)
    m.aromatize()
    m.layout()
    #m.saveMolfile("/home/eladn/Desktop/mol2.mol")
    s = _renderer.renderToBuffer(m).tostring()
    hash = "%09d" % random.randrange(0, 1e9)
    i = 0
    while True:
        symbol = 'glyph0-%d' % i
        if s.find(symbol) != -1:
            s = s.replace(symbol, hash + "_" + symbol)
        else:
            break
        i += 1
    return s

if __name__ == "__main__":
    if False:
        print smiles2svg('C#N')
    else:
        f = html_writer.HtmlWriter('../res/tmp.html')
        #s1 = smiles2svg('C(C1C(C(C(n2cnc3c(N)ncnc23)O1)O)O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O')
        #s1 = smiles2svg('n2cnc3c(N)ncnc23')
        s1 = smiles2svg('N#N')
        
        f.write('<div><p>1:</br>%s</p></div>\n' % (s1))
        s2 = smiles2svg('O=O')
        f.write('<div><p>2:</br>%s</p></div>\n' % (s2))
        f.close()
