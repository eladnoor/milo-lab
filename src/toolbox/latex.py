import types

def table2LaTeX(dict_list, headers=None):
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

    s = ' & '.join(headers) + "\\\\\n"
    s += "\\hline\n"
    for i, d in enumerate(dict_list):
        d['#'] = '%d' % i
        values = [to_string(d.get(key, "")) for key in headers]
        s += ' & '.join(values) + "\\\\\n"
    return s