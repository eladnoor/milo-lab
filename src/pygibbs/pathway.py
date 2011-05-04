#!/usr/bin/python

import re

from pygibbs.thermodynamic_constants import default_T, default_pH, default_I, default_pMg
from pygibbs.thermodynamic_constants import default_c0, default_c_mid, default_c_range


class PathwayConditions(object):
    
    def __init__(self,
                 temperature=default_T,
                 pH=default_pH,
                 pMg=default_pMg,
                 ionic_strength=default_I,
                 c0=default_c0,
                 media=None):
        self.T = temperature
        self.pH = pH
        self.pMg = pMg
        self.ionic_strength = ionic_strength
        self.c0 = c0
        self.media = media
    
    def _GetTemp(self):
        """Return the temperature associated with these conditions."""
        return self.T
    
    def _GetIonicStrength(self):
        """Return the ionic strength associated with these conditions."""
        return self.ionic_strength
    
    temp = property(_GetTemp)
    temperature = property(_GetTemp)
    I = property(_GetIonicStrength)

    @staticmethod
    def get_float_parameter(s, name, default_value):
        tokens = re.findall(name + "=([0-9\.e\+]+)", s)
        if len(tokens) == 0:
            return default_value
        if len(tokens) > 1:
            raise KeyError("The parameter %s appears more than once in %s" % (name, s))
        return float(tokens[0])

    @staticmethod
    def FromString(conditions_str):
        """Returns a list of conditions from a string in the KEGG file."""
        if not conditions_str or not conditions_str.strip():
            return [PathwayConditions()]
        
        conds = []
        for condition in conditions_str.split('\t'):
            media, pH, I, T, c0 = None, default_pH, default_I, default_T, default_c0
            media = re.findall("media=([a-zA-Z_]+)", condition)
            if media:
                media = media[0]
                if media == 'None':
                    media = None
                
            pH = PathwayConditions.get_float_parameter(condition, "pH", default_pH)
            I = PathwayConditions.get_float_parameter(condition, "I", default_I)
            T = PathwayConditions.get_float_parameter(condition, "T", default_T)
            c0 = PathwayConditions.get_float_parameter(condition, "c0", default_c0)
            conds.append(PathwayConditions(temperature=T,
                                           pH=pH, ionic_strength=I,
                                           c0=c0, media=media))
            
        return conds


class PathwayData(object):
    
    def __init__(self):
        self.field_map = {}
        self.name = None
        self.skip = False
        self.analysis_type = None
        self.redox_values = []
        self.conditions = []
        self.pH = default_pH
        self.pH_values = []
        self.I = default_I
        self.I_values = []
        self.T = default_T
        self.pMg = default_pMg
        self.c_mid = None
        self.c_range = None
        self.cid_mapping = {}
        self.kegg_module_id = None
        self.dG_methods = []
    
    @staticmethod
    def FromFieldMap(field_map):
        p = PathwayData()
        p.field_map = field_map
        p.analysis_type = field_map.GetStringField('TYPE')
        p.skip = field_map.GetBoolField('SKIP', default_value=False)
        p.name = field_map.get('NAME')
        
        
        p.conditions = PathwayConditions.FromString(field_map.get('CONDITIONS'))
        
        try:
            p.pH = field_map.GetFloatField('PH', default_pH)
        except ValueError:
            p.pH_values = field_map.GetVFloatField(
                'PH', default_value=[5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
        
        try:
            p.I = field_map.GetFloatField('I', default_I)
        except ValueError:
            p.I_values = field_map.GetVFloatField(
                'I', default_value=[0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4])
            
        p.T = field_map.GetFloatField('T', default_T)
        p.pMg = field_map.GetFloatField('PMG', default_pMg)
        p.c_mid = field_map.GetFloatField('C_MID', default_c_mid)
        p.redox_values = field_map.GetVFloatField("REDOX", [])
        
        c_range = field_map.GetVFloatField('C_RANGE', default_c_range)
        if c_range:
            p.c_range = tuple(c_range)
        
        if "MAP_CID" in field_map:
            for line in field_map["MAP_CID"].strip().split('\t'):
                cid_before, cid_after = [int(cid[1:]) for cid in line.split(None, 1)]
                p.cid_mapping[cid_before] = (cid_after, 1.0)
        
        if 'MODULE' in field_map:
            mid_str = field_map["MODULE"]
            if (mid_str[0] == 'M'):
                p.kegg_module_id = int(mid_str[1:])
            else:
                p.kegg_module_id = int(mid_str)
        
        if field_map.GetBoolField('MILO', default_value=True):
            p.dG_methods.append('MILO')
        if field_map.GetBoolField('HATZI', default_value=False):
            p.dG_methods.append('HATZI')
        
        return p
        
        
        
        

        
        
        
        