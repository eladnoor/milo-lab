import csv

def write_csv_query(cursor, filename, query, titles=None):
    writer = csv.writer(open(filename, 'w'))
    if (titles != None):
        writer.writerow(titles)
    for row in cursor.execute(query):
        writer.writerow(list(row))

def write_csv_table(cursor, filename, tablename, titles=True):
    writer = csv.writer(open(filename, 'w'))
    cursor.execute("SELECT * FROM " + tablename)
    if (titles):
        col_name_list = [tuple[0] for tuple in cursor.description]
        writer.writerow(col_name_list)
    for row in cursor:
        writer.writerow(list(row))
