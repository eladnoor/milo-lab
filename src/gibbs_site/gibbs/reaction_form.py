from django import forms
from gibbs import form_utils
from gibbs import constants


class ReactionForm(form_utils.BaseForm):
    reactantsId = form_utils.ListFormField()
    productsId = form_utils.ListFormField()
    reactantsCoeff = form_utils.ListFormField()
    productsCoeff = form_utils.ListFormField()
    reactantsName = form_utils.ListFormField()
    productsName = form_utils.ListFormField()
    reactantsConcentration = form_utils.ListFormField(required=False)
    productsConcentration = form_utils.ListFormField(required=False)

    ph = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    concentration_profile = forms.ChoiceField(required=False,
                                              choices=[('1M', '1M'),
                                                       ('1mM', '1mM'),
                                                       ('custom', 'custom')])
    
    query = forms.CharField(max_length=2048, required=False)
    balance_w_water = forms.BooleanField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_reactantIds = property(lambda self: self.cleaned_data['reactantsId'])
    cleaned_productIds = property(lambda self: self.cleaned_data['productsId'])
    cleaned_reactantCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['reactantsCoeff']])
    cleaned_productCoeffs = property(lambda self: [int(c) for c in self.cleaned_data['productsCoeff']])
    cleaned_reactantNames = property(lambda self: self.cleaned_data['reactantsName'])
    cleaned_productNames = property(lambda self: self.cleaned_data['productsName'])
    cleaned_reactantConcentrations = property(lambda self: [float(c) for c in self.cleaned_data['reactantsConcentration']])
    cleaned_productConcentrations = property(lambda self: [float(c) for c in self.cleaned_data['productsConcentration']])
    cleaned_ph = property(lambda self: self._GetWithDefault('ph', constants.DEFAULT_PH))
    cleaned_ionic_strength = property(lambda self: self._GetWithDefault('ionic_strength',
                                                                        constants.DEFAULT_IONIC_STRENGTH))
    cleaned_concentration_profile = property(lambda self: self.cleaned_data['concentration_profile'])
    cleaned_query = property(lambda self: self.cleaned_data['query'])
    cleaned_balance_w_water = property(lambda self: self._GetWithDefault('balance_w_water', False))
