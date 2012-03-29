import indigo, indigo_renderer
        
def draw(smiles, fname):
    _indigo = indigo.Indigo()
    _renderer = indigo_renderer.IndigoRenderer(_indigo)
    _indigo.setOption('render-output-format', 'png')
    _indigo.setOption('render-image-size', 250, 200)
    _indigo.setOption('render-margins', 10, 10)
    _indigo.setOption('render-stereo-style', 'none')
    _indigo.setOption('render-implicit-hydrogens-visible', False)
    _indigo.setOption('render-coloring', True)
    _indigo.setOption('render-bond-length', 50.0)
    _indigo.setOption('render-label-mode', 'hetero')
    
    indigo_mol = _indigo.loadMolecule(smiles)
    indigo_mol.aromatize()
    indigo_mol.layout()
    
    _renderer.renderToFile(indigo_mol, fname)


draw("CC(=O)O", "/home/eladn/Desktop/acetate.png")
draw("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "/home/eladn/Desktop/glucose.png")
draw("O=C(NCCS)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3OP(=O)(O)O", "/home/eladn/Desktop/coa.png")