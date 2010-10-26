import csv

def write_csv_query(cursor, filename, query, titles=None):
    writer = csv.writer(open(filename, 'w'))
    if (titles != None):
        writer.writerow(titles)
    for row in cursor.execute(query):
        writer.writerow(list(row))

def write_csv_table(cursor, filename, tablename, write_titles=True):
    column_names = []
    for row in cursor.execute("PRAGMA table_info(%s)" % tablename):
        #(index, name, type, can_be_null, default_value, stam) = row
        column_names.append(row[1])
    
    writer = csv.writer(open(filename, 'w'))
    if (write_titles):
        writer.writerow(column_names)

    cursor.execute("SELECT * FROM %s" % tablename)
    for row in cursor:
        writer.writerow(list(row))

def write_table_to_html(cursor, html_writer, tablename):
    column_names = []
    for row in cursor.execute("PRAGMA table_info(%s)" % tablename):
        #(index, name, type, can_be_null, default_value, stam) = row
        column_names.append(row[1])
    
    html_writer.write('<table border="1">\n')
    html_writer.write('  <tr>' + "".join([('<td><b>' + str(s) + '</b></td>') for s in column_names]) + '</tr>\n')
    for row in cursor.execute("SELECT * FROM %s" % tablename):
        html_writer.write('  <tr>' + "".join([('<td>' + str(s) + '</td>') for s in row]) + '</tr>\n')
    html_writer.write('</table>\n')

def write_query_to_html(cursor, html_writer, query, column_names=None):
    html_writer.write('<table border="1">\n')
    if (column_names != None):
        html_writer.write('  <tr>' + "".join([('<td><b>' + str(s) + '</b></td>') for s in column_names]) + '</tr>\n')
    for row in cursor.execute(query):
        html_writer.write('  <tr>' + "".join([('<td>' + str(s) + '</td>') for s in row]) + '</tr>\n')
    html_writer.write('</table>\n')
    html_writer.write('</p>\n')
