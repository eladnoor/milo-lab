from django import forms

class SearchForm(forms.Form):
    query = forms.CharField(max_length=1024)
    
class FetchMailForm(forms.Form):
    email = forms.CharField(max_length=30)
    password = forms.CharField(max_length=30)
    
class UploadSigmaPdfForm(forms.Form):
    #title = forms.CharField(max_length=50)
    file  = forms.FileField()
    
class UpdateLocationForm(forms.Form):
    box = forms.IntegerField()
    row = forms.CharField(max_length=1)
    col = forms.IntegerField()
