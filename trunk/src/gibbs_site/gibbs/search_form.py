from django import forms
from gibbs import constants


class SearchForm(forms.Form):
    query = forms.CharField(max_length=2048)
    temp = forms.FloatField(required=False)
    ph = forms.FloatField(required=False)
    ionic_strength = forms.FloatField(required=False)
    
    # Convenience accessors for clean data with defaults.
    cleaned_query = property(lambda self: self.cleaned_data['query'])
    cleaned_temp = property(lambda self: (self.cleaned_data['temp'] or
                                          constants.DEFAULT_TEMP))
    cleaned_ph = property(lambda self: (self.cleaned_data['ph'] or
                                        constants.DEFAULT_PH))
    cleaned_ionic_strength = property(lambda self: (self.cleaned_data['ionic_strength'] or
                                                    constants.DEFAULT_IONIC_STRENGTH))