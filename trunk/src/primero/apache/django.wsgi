import os
import sys
 
paths = ['/srv/www', '/srv/www/primero', '.']
for path in paths:
    if path not in sys.path:
        sys.path.insert(0, path)
 
os.environ['DJANGO_SETTINGS_MODULE'] = 'primero.settings'
 
import django.core.handlers.wsgi
application = django.core.handlers.wsgi.WSGIHandler()

