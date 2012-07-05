from toolbox.database import SqliteDatabase
import numpy as np

cond = "(cast(kgp.dGc1 as real) > %d AND cast(kgp.dGc2 as real) < %d)" % (40, -40)
query1 = """SELECT  p.*, pfi.score
            FROM (
                  SELECT  kgp.gene1 gene1, kgp.gene2 gene2, %s, count(*) ntot
                  FROM    kegg_gene_pairs kgp
                  GROUP BY kgp.gene1, kgp.gene2
                 ) p
            LEFT OUTER JOIN parkinson_functional_interactions pfi
            ON      (pfi.gene1 = p.gene1 AND pfi.gene2 = p.gene2
                     OR
                     pfi.gene1 = p.gene2 AND pfi.gene2 = p.gene1)""" % cond

db = SqliteDatabase('channeling/channeling.sqlite', 'w')

counter = 0
for row in db.Execute(query1):
    gene1, gene2, nqual, ntot, score = row
    if nqual > 0:
        counter += 1
        print counter, gene1, gene2

for row in db.Execute("select * from kegg_gene_pairs where gene1 = 'eco:b0074' and gene2 = 'eco:b0720'"):
    print row