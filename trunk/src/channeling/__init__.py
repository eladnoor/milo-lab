import logging
import sys
import os

formatter = logging.Formatter('%(asctime)s %(filename)s %(lineno)s %(levelname)s  %(message)s')

stdout_handler = logging.StreamHandler(sys.stdout)
stderr_handler = logging.StreamHandler(sys.stderr)
stdout_handler.setFormatter(formatter)
stderr_handler.setFormatter(formatter)

logger = logging.getLogger('')
logger.addHandler(stdout_handler)
logger.setLevel(logging.INFO)

if not os.getcwd()[-4:] == '/src':
    os.chdir('..')