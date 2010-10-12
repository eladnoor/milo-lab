import re

from django.db import models


class CommonName(models.Model):
    name = models.CharField(max_length=500)
    enabled = models.BooleanField()
    
    def __unicode__(self):
        return self.name
    

class Compound(models.Model):
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # The chemical formula.
    formula = models.CharField(max_length=500)
    
    # The molecular mass.
    mass = models.FloatField()  # In Daltons.
    
    
    def GetKeggLink(self):
        """Returns a link to the KEGG page for this compound."""
        if not self.kegg_id:
            return None
        
        return 'http://kegg.jp/dbget-bin/www_bget?cpd:%s' % self.kegg_id

    def GetHtmlFormattedFormula(self):
        """Returns the chemical formula with HTML formatted subscripts."""
        if not self.formula:
            return None
        
        return re.sub(r'(\d+)', r'<sub>\1</sub>', self.formula)
    
    def __unicode__(self):
        names = self.common_names.all()
        if names:
            return unicode(names[0])
        return unicode(self.formula)
    
    @staticmethod
    def AsLibrary():
        """Get's a map of common names to Compound objects from the DB."""
        all_compounds = Compound.objects.all()
        names_to_compounds = {}
        for compound in all_compounds:
            for name in compound.common_names.all():
                names_to_compounds[name] = compound
        
        return names_to_compounds
    
    @staticmethod
    def GetCompoundsByKeggId(kegg_ids):
        """Fetch compounds from a list of KEGG IDs.
        
        Args:
            kegg_ids: a list of KEGG IDs.
        
        Returns:
            A dictionary mapping KEGG ID to Compounds.
        """
        compounds = Compound.objects.filter(kegg_id__in=kegg_ids)
        return dict((c.kegg_id, c) for c in compounds if c != None)
        
        
    
