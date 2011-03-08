from HTMLParser import HTMLParser, HTMLParseError
from urllib2 import urlopen
from copy import deepcopy
from toolbox.database import SqliteDatabase
import logging

class ShowAll(HTMLParser):
    def __init__(self, url):
        HTMLParser.__init__(self)
        self._data = {}
        self._current_id = None
        self._current_textfield = ""
        req = urlopen(url)
        self.feed(req.read())
        
    def handle_starttag(self, tag, attrs):
        if (tag == "a"):
            for (attr, value) in attrs:
                if (attr == 'href' and value.find('enzyme_data1.pl?T1=') != -1):
                    self._current_id = value
                    self._current_id.replace('enzyme_data1.pl?T1=', '')
                    self._data[self._current_id] = []
                    self._current_textfield = ""
        if (self._current_id != None):
            if (tag == "td"):
                self._current_textfield = ""
            elif (tag == "sub"):
                self._current_textfield += ""
            elif (tag == "sup"):
                self._current_textfield += "-"
            elif (tag == "img" and attrs == ['SRC', './beta.gif']):
                self._current_textfield += "beta"

    def handle_endtag(self, tag):
        if (self._current_id != None):
            if (tag == "td"):
                self._data[self._current_id].append(self._current_textfield)
                self._current_textfield = ""
            elif (tag == "tr"):
                self._current_id = None
            elif (tag == "sub"):
                self._current_textfield += ""
            elif (tag == "sup"):
                self._current_textfield += ""

    def handle_data(self, data):
        if (self._current_id != None):
            self._current_textfield += data.replace('\n','').strip()

class EnzymeData(HTMLParser):
    def __init__(self, url):
        HTMLParser.__init__(self)
        self._data = []
        self._tables = []
        self._table_data = []
        self._num_cols = 0
        self._table_row = []
        self._table_cell = ""
        
        req = urlopen(url)
        self.feed(req.read())
        
    def handle_starttag(self, tag, attrs):
        if (tag == "table"):
            self._table_data = []
            self._table_row = []
            self._table_cell = ""
            self._num_cols = 0
        elif (tag == "tr"):
            self._table_row = []
            self._table_cell = ""
        elif (tag == "td"):
            self._table_cell = ""

    def handle_endtag(self, tag):
        if (tag == "html"):
            self._tables = []
        elif (tag == "table"):
            if (len(self._table_data) == 1):
                self._data += self._table_data[0]
            else:
                self._tables.append(self._table_data)
        elif (tag == "tr"):
            if (self._num_cols == 0):
                self._num_cols = len(self._table_row)
            elif (self._num_cols != len(self._table_row)):
                return # this row is invalid (number of cells doesn't match the first row)
            self._table_data.append(self._table_row)
        elif (tag == "td"):
            self._table_row.append(self._table_cell)
                
    def handle_data(self, data):
        self._table_cell += data.replace('\n','').strip()

def fetch_nist_data(db):
    domain = "http://xpdb.nist.gov/enzyme_thermodynamics"
    LOG = open('../res/nist.log', 'w')
    
    # First, gather all the links and reactions from the 'enzyme_show_all_data.pl' page
    url2reaction_map = {}
    col = 0
    logging.info("Parsing enzyme_show_all_data.pl")
    while True:
        col += 1
        print col,
        url = domain + "/enzyme_show_all_data.pl?col=%d.&R1=V2" % col
        spider = ShowAll(url)
        if len(spider._data) == 0:
            break
        for url, values in spider._data.iteritems():
            url2reaction_map[url] = values[1]
        
    # Then, for each URL, download the table of data and store in all_values
    all_field_names = set()
    all_values = []
    for url, reaction in url2reaction_map.iteritems():
        col = 0
        logging.info("Parsing URL: %s" % url)
        while True:
            col += 1
            url_curr = url.replace('enzyme_data1.pl?T1=','enzyme_data1.pl?col=%d.&T1=' % col)
            try:
                edata = EnzymeData(domain + "/" + url_curr)
                if (edata._tables == []):
                    break
                    
                data_map = {'Reaction' : reaction, 'URL' : url_curr}
                for val in edata._data:
                    i = val.find(':')
                    if (i == -1):
                        continue
                    data_map[val[0:i]] = val[(i+1):].strip()
                
                all_field_names = all_field_names.union(set(data_map.keys()))
                for table in edata._tables:
                    firstrow = table[0]
                    all_field_names = all_field_names.union(set(firstrow))
                    for row in table[1:]:
                        row_map = deepcopy(data_map)
                        for i in range(len(row)):
                            row_map[firstrow[i]] = row[i]
                        all_values.append(row_map)
                LOG.write(url_curr + ": OK\n")
            except HTMLParseError:
                LOG.write(url_curr + ": HTML Parse ERROR\n")
                logging.warning("HTML Parsing failed")
                break

    # Now, insert all the data into one big table in the DB (TABLE raw)
    n = len(all_field_names)
    db.CreateTable('nist_fields', ['col_number INT', 'name TEXT'])
    all_field_names = list(all_field_names)
    for i in xrange(n):
        db.Insert('nist_fields', [i, all_field_names[i]])
    db.CreateTable('nist_values', ['id INT'] + ["field%d TEXT" % i for i in xrange(n)])
    
    row_id = 0
    for row_map in all_values:
        db.Insert('nist_values', [row_id] + [row_map.get(key, None) for key in all_field_names])
        row_id += 1
    db.Commit()
    LOG.close()
    
################################################################################################################
#                                                 MAIN                                                         #        
################################################################################################################
if __name__ == "__main__":
    db = SqliteDatabase('../res/nist_raw.sqlite')
    fetch_nist_data(db)
