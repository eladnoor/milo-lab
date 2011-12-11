#!/bin/bash

mysqldump --databases djtest -u djangouser -p > data/sqldump.txt
