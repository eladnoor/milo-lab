from django import forms
from gibbs import form_utils
from gibbs import constants


class ReactionForm(form_utils.BaseForm):
    reactantsId = form_utils.ListFormField()
    productsId = form_utils.ListFormField()
    reactantsCoeff = form_utils.ListFormField()
    productsCoeff = form_utils.ListFormField()
    
    temp = forms.FloatField(required=False)
    ph = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    concentrationProfile = forms.ChoiceField(required=False, choices=["1M", "1mM"])
    
    # Convenience accessors for clean data with defaults.
    cleaned_reactantIds = property(lambda self: self.cleaned_data['reactantsId'])
    cleaned_productIds = property(lambda self: self.cleaned_data['productsId'])
    cleaned_reactantCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['reactantsCoeff']])
    cleaned_productCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['productsCoeff']])
    cleaned_ph = property(lambda self: self._GetWithDefault('ph', constants.DEFAULT_PH))
    cleaned_ionic_strength = property(lambda self: self._GetWithDefault('ionic_strength',
                                                                        constants.DEFAULT_IONIC_STRENGTH))