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
print "Done"

# Now we can import Django stuff!
from django.db import connection
from django.db.models import get_apps
from django.contrib.auth.management import create_superuser
from django.contrib.auth import models as auth_app
from django.db.models import signals
from django.dispatch import dispatcher

#print "Dumping Data... ",
#sys.stdout = open('data/dumped_data.json', 'w')
#call_command('dumpdata', format='json', indent=4)
#sys.stdout.close()
#sys.stdout = sys.__stdout__
#print "Done"

print "Deleting Tables... ",
app_labels = [app.__name__.split('.')[-2] for app in get_apps()]
sys.stdout = buffer = StringIO.StringIO()
call_command('sqlclear', *app_labels)
sys.stdout = sys.__stdout__

queries = buffer.getvalue().split(';')[1:-2]

cursor = connection.cursor()
for query in queries:
	cursor.execute(query.strip())
print "Done"

print "Sync-ing Database... "
#disable the "create a super user" question
#dispatcher.disconnect(create_superuser, sender=auth_app, signal=signals.post_syncdb)

call_command('syncdb')
print "Done"

#print "Loading Back Data... "
#call_command('loaddata', 'data/dumped_data.json')
#os.remove('dumped_data.json')
print "Done"
