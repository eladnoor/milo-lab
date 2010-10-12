from django import forms
from gibbs import form_utils
from group_contribution import constants


class ReactionForm(forms.Form):
    reactantsId = form_utils.ListFormField()
    productsId = form_utils.ListFormField()
    reactantsCoeff = form_utils.ListFormField()
    productsCoeff = form_utils.ListFormField()
    
    temp = forms.FloatField(required=False)
    ph = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_reactantIds = property(lambda self: self.cleaned_data['reactantsId'])
    cleaned_productIds = property(lambda self: self.cleaned_data['productsId'])
    cleaned_reactantCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['reactantsCoeff']])
    cleaned_productCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['productsCoeff']])
    cleaned_temp = property(lambda self: float(self.cleaned_data['temp'] or
                                               constants.DEFAULT_TEMP))
    cleaned_ph = property(lambda self: float(self.cleaned_data['ph'] or
                                             constants.DEFAULT_PH))
    cleaned_ionic_strength = property(lambda self: float(self.cleaned_data['ionic_strength'] or
                                                         constants.DEFAULT_IONIC_STRENGTH))