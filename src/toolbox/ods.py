import ooolib

def ODSDictReader(f, sheet_index=0):
    row_dicts = []
    doc = ooolib.Calc()
    doc.load(f)
    doc.set_sheet_index(sheet_index)
    Nc, Nr = doc.get_sheet_dimensions()
    titles = [doc.get_cell_value(c+1, 1)[1] for c in xrange(Nc)]
    for r in xrange(1, Nr):
        dict = {}
        for c in xrange(Nc):
            v = doc.get_cell_value(c+1, r+1)
            if v:
                dict[titles[c]] = v[1]
            else:
                dict[titles[c]] = None
        row_dicts.append(dict)
    return row_dicts

