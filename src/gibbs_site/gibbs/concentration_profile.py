class _BaseConcentrationProfile(object):
    
    def Concentration(self, kegg_id):
        raise NotImplementedError


class MolarProfile(_BaseConcentrationProfile):
    
    def Concentration(self, kegg_id):
        return 1.0


class MilliMolarProfile(_BaseConcentrationProfile):
    
    def Concentration(self, kegg_id):
        if kegg_id == 'C00001':
            return 1.0
        return 0.001

PROFILES_BY_NAME = {'1M': MolarProfile(),
                    '1mM': MilliMolarProfile()}

def GetProfile(name):
    return PROFILES_BY_NAME.get(name, None)