#!/usr/bin/python

import re

from pygibbs import kegg_parser
from pygibbs.thermodynamic_constants import default_T, default_pH, default_I
from pygibbs.thermodynamic_constants import default_pMg, default_c0
from pygibbs.thermodynamic_constants import default_c_range, default_c_mid

from pygibbs.metabolic_modelling import stoich_model
from pygibbs.metabolic_modelling import bounds
from pygibbs.kegg import Kegg



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
    def _RemapCompounds(S, cids, cid_mapping):
        """Remap the compounds in S.
        
        Args:
            S: the stoichiometric matrix.
            cids: the listing of cids in the matrix.
            cid_mapping: the mapping to replace them.
            
        Returns:
            Tuple of updated values (S, cids).
        """
        for i, cid in enumerate(list(cids)):
            if cid in cid_mapping:
                new_cid, coeff = cid_mapping[cid]
                cids[i] = new_cid
                S[:, i] *= coeff
        
        return S, cids
    
    def GetBounds(self):
        """Get a bounds.Bounds object."""
        lower_bounds = {1: 1}
        upper_bounds = {1: 1} # the default for H2O is 1
        
        bound_str = self.field_map.get("BOUND", "")
        
        for line in bound_str.strip().split('\t'):
            if not line:
                continue
            
            tokens = line.split(None)
            cid = int(tokens[0][1:])
            
            try:
                b_lower = float(tokens[1])
            except ValueError:
                b_lower = None    

            if len(tokens) == 2:
                b_upper = b_lower
            elif len(tokens) == 3:
                try:
                    b_upper = float(tokens[2])
                except ValueError:
                    b_upper = None    
            else:
                raise ValueError("Parsing error in BOUND definition for %s: %s" %
                                 (self.name, line))
                
            lower_bounds[cid] = b_lower
            upper_bounds[cid] = b_upper
        
        default_lb, default_ub = self.c_range
        return bounds.Bounds(lower_bounds, upper_bounds,
                             default_lb=default_lb, default_ub=default_ub)

    def GetStoichiometricModel(self, kegg):
        """Make a stoichiometric model from this pathway.
        
        Args:
            kegg: a Kegg instance.
        
        Returns:
            A StoichiometricModel instance.
        """
        # Explicitly map some of the CIDs to new ones.
        # This is useful, for example, when a KEGG module uses unspecific co-factor pairs,
        # like NTP => NDP, and we replace them with ATP => ADP 
        cid_mapping = self.cid_mapping
        field_map = self.field_map
        
        # Case 1: it's a module
        if self.kegg_module_id is not None:
            S, rids, fluxes, cids = kegg.get_module(self.kegg_module_id)
            S, cids = self._RemapCompounds(S, cids, cid_mapping)
        # Case 2: it's an explicitly defined pathway.
        else:
            S, rids, fluxes, cids = kegg.parse_explicit_module(field_map, cid_mapping) 
        
        model = stoich_model.StoichiometricModel(S, rids, cids,
                                                 fluxes=fluxes,
                                                 name=self.name)
        return model
        
    
    @staticmethod
    def FromFieldMap(field_map):
        p = PathwayData()
        p.field_map = field_map
        p.analysis_type = field_map.GetStringField('TYPE', default_value='None')
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
    
    def get_explicit_reactions(self, balance_water=True):
        kegg = Kegg.getInstance()
        return kegg.parse_explicit_module(self.field_map, self.cid_mapping,
                                          balance_water=balance_water) 

class KeggPathwayIterator(object):
    
    def __init__(self, parsed_kegg_file):
        self.parsed_kegg_file = parsed_kegg_file
        
    @staticmethod
    def FromFilename(fname):
        """Initialize a KeggPathwayIterator from a filename.
        
        Args:
            fname: a valid path to the file containing pathway definitions.
        """
        parsed = kegg_parser.ParsedKeggFile.FromKeggFile(fname)
        return KeggPathwayIterator(parsed)
    
    def __iter__(self):
        """Iterate over pathways."""
        for key in sorted(self.parsed_kegg_file.keys()):
            field_map = self.parsed_kegg_file[key]
            yield PathwayData.FromFieldMap(field_map)
