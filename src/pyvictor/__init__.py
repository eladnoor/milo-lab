import os

if not os.getcwd()[-4:] == '/src':
    os.chdir('..')