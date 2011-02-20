#!/usr/bin/python

import csv

class Elements:
    def __init__(self):
        csv_file = csv.reader(open('../data/thermodynamics/elements.csv', 'r'))
        csv_file.next()
        self.symbol_to_an = {}
        self.an_to_symbol = {}
        self.symbols = []
        for row in csv_file:
            # an = Atomic Number, mp = Melting Point, bp = Boiling Point, 
            # ec = Electron Configuration, ie = Ionization Energy
            an, unused_weight, unused_name, symbol, unused_mp, unused_bp, \
                unused_density, unused_earthcrust, unused_discovery, unused_ec, unused_ie = row
            self.symbol_to_an[symbol] = int(an)
            self.an_to_symbol[int(an)] = symbol
    
global ELEMENTS;
ELEMENTS = Elements()