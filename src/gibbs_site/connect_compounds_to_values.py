#!/usr/bin/env python

import os
import sys
import StringIO
from django.core.management import setup_environ, call_command

print "Setting Up Environment... ",
try:
	import settings # Assumed to be in the same directory.
except ImportError:
	sys.stderr.write("Error: Can't find the file 'settings.py' in the directory containing %r.\n"
									 "(If the file settings.py does indeed exist, it's causing an ImportError somehow.)\n" % __file__)
	sys.exit(1)

setup_environ(settings)


from gibbs import models

print 'Connecting'
for e_formation in models.SpeciesFormationEnergy.objects.all():
	kegg_id = e_formation.specie.kegg_id
	try:
		compound = models.Compound.objects.get(kegg_id=kegg_id)
	except:
		print 'Missing kegg_id:', kegg_id
		continue

	pks = set(ef.pk for ef in compound.species_formation_energies.all())

	if e_formation.pk not in pks:
		compound.species_formation_energies.add(e_formation)
		compound.save()

print 'Done.'
