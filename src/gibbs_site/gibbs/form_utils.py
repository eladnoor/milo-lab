from django import forms
import logging

class ListFormField(forms.MultipleChoiceField):
    """A form field for a list of values that are unchecked.
    
    The Django MultipleChoiceField does *almost* what we want, except
    it validates that each choice is in a supplied list of choices, 
    even when that list is empty. We simply override the validation.
    """
    
    def valid_value(self, value):
        return True