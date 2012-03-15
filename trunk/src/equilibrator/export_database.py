#!/usr/bin/python

from util import django_utils

django_utils.SetupDjango()

from gibbs import models


def main():
    all_compounds = models.Compound.objects.all()
    for c in all_compounds:
        # Do something
        pass
    
    
if __name__ == '__main__':
    main()
 