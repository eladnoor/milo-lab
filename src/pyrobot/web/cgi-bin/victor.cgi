#!/usr/bin/python

import cgi, os
import cgitb; cgitb.enable()
from victor_uploader import ImportFileToDatabase, CreateDummyDB

form = cgi.FieldStorage()

# Get filename here.
fileitem = form['filename']

# Test if the file was uploaded
if fileitem.filename:
    # strip leading path from file name to avoid 
    # directory traversal attacks
    fn = os.path.basename(fileitem.filename)
    db = CreateDummyDB()
    exp_id = ImportFileToDatabase(fileitem.file, db)
    message = 'The file "' + fn + '" was uploaded successfully'
   
else:
    message = 'No file was uploaded'
   
print """\
Content-Type: text/html\n
<html>
<body>
   <p>%s</p>
</body>
</html>
""" % (message,)
