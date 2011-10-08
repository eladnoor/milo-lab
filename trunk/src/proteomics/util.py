
import csv

def ReadProteinIDs(filename):
    f = open(filename)
    reader = csv.DictReader(f)
    gene_ids = dict((r.get("Identifier"), r.get("Name")) for r in reader)
    f.close()
    return gene_ids

def ReadProteinCounts(filename):
    f = open(filename)
    reader = csv.reader(f, dialect=csv.excel_tab)
    gene_counts = {}
    for row in reader:
        if not row or row[0].startswith('#'):
            continue
        
        _, id = row[1].split('.')
        gene_counts[id] = float(row[2])
    f.close()
    return gene_counts

def ExtractCounts(gene_counts, gene_ids):
    """Yields id, name, count for each gene in both sets."""
    
    for id, name in gene_ids.iteritems():
        count = gene_counts.get(id, None)
        if count is None:
            continue
        
        yield id, name, count