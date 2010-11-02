
import unittest
from django.core.management import setup_environ
import settings

def SetupDjango():
    setup_environ(settings)
    