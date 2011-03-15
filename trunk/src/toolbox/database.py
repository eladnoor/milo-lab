import sqlite3, csv
import os
import types
import logging

class Database(object):
    
    def __init__(self):
        raise Exception("class Database is an abstract class and cannot be instantiated")
    
    def CreateTable(self, table_name, columns, drop_if_exists=True):
        raise Exception("function not implemented")
            
    def Insert(self, table_name, list):
        raise Exception("function not implemented")
    
    def DoesTableExist(self, table_name):
        raise Exception("function not implemented")
    
class SQLDatabase(Database):
    
    def __init__(self, filename):
        raise Exception("class SQLDatabase is an abstract class and cannot be instantiated")
    
    def Execute(self, command, arguments=None):
        raise Exception("function not implemented")
    
    def Commit(self):
        raise Exception("function not implemented")
    
    def CreateTable(self, table_name, columns, drop_if_exists=True):
        if drop_if_exists:
            self.Execute("DROP TABLE IF EXISTS %s" % table_name)
        elif self.DoesTableExist(table_name):
            return
        
        if type(columns) == types.ListType:
            columns = ','.join(['%s TEXT' % col for col in columns])    
        
        self.Execute("CREATE TABLE %s (%s)" % (table_name, columns))
    
    def CreateIndex(self, index_name, table_name, columns, unique=True, drop_if_exists=True):
        if drop_if_exists:
            self.Execute("DROP INDEX IF EXISTS %s" % index_name)
        if unique:
            self.Execute("CREATE UNIQUE INDEX %s ON %s (%s)" % (index_name, table_name, columns))
        else:
            self.Execute("CREATE INDEX %s ON %s (%s)" % (index_name, table_name, columns))
            
    def Insert(self, table_name, list):
        sql_command = "INSERT INTO %s VALUES(%s)" % (table_name, ','.join(["?"]*len(list)))
        return self.Execute(sql_command, list)
    
    def DoesTableExist(self, table_name):
        for unused_row in self.Execute("SELECT name FROM sqlite_master WHERE name='%s'" % table_name):
            return True
        return False
    
    def Query2String(self, query, column_names=None):
        s = ""
        if column_names:
            s += ','.join(column_names) + "\n"
        for row in self.Execute(query):
            s += ','.join([str(x) for x in row]) + "\n"
        return s

    def Query2CSV(self, filename, query, column_names=None):
        writer = csv.writer(open(filename, 'w'))
        if column_names:
            writer.writerow(column_names)
        for row in self.Execute(query):
            writer.writerow(list(row))
    
    def Query2HTML(self, html_writer, query, column_names=None):
        html_writer.write('<table border="1">\n')
        if column_names:
            html_writer.write('  <tr>' + "".join([('<td><b>' + str(s) + '</b></td>') for s in column_names]) + '</tr>\n')
        for row in self.Execute(query):
            html_writer.write('  <tr>' + "".join([('<td>' + str(s) + '</td>') for s in row]) + '</tr>\n')
        html_writer.write('</table>\n')
        html_writer.write('</p>\n')

    def GetColumnNames(self, table_name):        
        column_names = []
        for unused_index, name, unused_type, unused_can_be_null, unused_default_value, unused_stam \
            in self.Execute("PRAGMA table_info(%s)" % table_name):
            
            column_names.append(name)
        
        return column_names
    
    def Table2String(self, table_name, write_titles=True):
        if write_titles:
            column_names = self.GetColumnNames(table_name)
        else:
            column_names = None
        return self.Query2String("SELECT * FROM %s" % table_name, column_names)
    
    def Table2CSV(self, filename, table_name, write_titles=True):
        if write_titles:
            column_names = self.GetColumnNames(table_name)
        else:
            column_names = None
        self.Query2CSV(filename, "SELECT * FROM %s" % table_name, column_names)
    
    def Table2HTML(self, html_writer, table_name, write_titles=True):
        if write_titles:
            column_names = self.GetColumnNames(table_name)
        else:
            column_names = None
        self.Query2HTML(html_writer, "SELECT * FROM %s" % table_name, column_names)

    def DictReader(self, table_name):
        titles = self.GetColumnNames(table_name)
        table_data = []
        for row in self.Execute("SELECT * FROM %s" % table_name):
            row_dict = {}
            for i, title in enumerate(titles):
                row_dict[title] = row[i]
            table_data.append(row_dict)
        return table_data    
        
class SqliteDatabase(SQLDatabase):
    
    def __init__(self, filename, flag='w'):
        """
            Connects to a Sqlite database.
            filename - the name of the file containing the database
            flag - can only be 'w' or 'r'.
                   use 'r' to verify that the file already exists (throws IOError otherwise)
                   use 'w' to create the file if it doesn't exist
        """
        if flag == 'r':
            if not os.path.exists(filename):
                raise IOError('No such file or directory: %s' % filename)
            self.comm = sqlite3.connect(filename)
        elif flag == 'w':
            self.comm = sqlite3.connect(filename)
        else:
            raise ValueError('Illegal flag: %s' % str(flag))
    
    def Execute(self, command, arguments=None):
        try:
            if arguments:
                return self.comm.execute(command, arguments)
            else:
                return self.comm.execute(command)
        except sqlite3.InterfaceError, e:
            logging.error('Failed to execute database command: ' + command)
            raise e
        
    def Commit(self):
        self.comm.commit()

    def __del__(self):
        self.comm.commit()
        self.comm.close()
