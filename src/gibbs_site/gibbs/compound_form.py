from django import forms
from group_contribution import constants


class CompoundForm(forms.Form):
    compoundId = forms.CharField(max_length=50)
    temp = forms.FloatField(required=False)
    ph = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_compoundId = property(lambda self: self.cleaned_data['compoundId'])
    cleaned_temp = property(lambda self: float(self.cleaned_data['temp'] or
                                               constants.DEFAULT_TEMP))
    cleaned_ph = property(lambda self: float(self.cleaned_data['ph'] or
                                             constants.DEFAULT_PH))
    cleaned_ionic_strength = property(lambda self: float(self.cleaned_data['ionic_strength'] or
                                                         constants.DEFAULT_IONIC_STRENGTH))