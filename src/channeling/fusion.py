######################################################################
#
# The class FusedChannel is designed to search (using blast) over
# a given range of molecules which pairs of them are suspected to
# undergo a fusion process in different organisms.
#
######################################################################


from Bio.Blast import NCBIWWW, NCBIXML
import re, logging


def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def getXMLInfo (xmlLine):
    regexp = re.match(r'<.+?>(.*?)<\/.+?>',xmlLine)
    if regexp:
        return regexp.group(1)
    return ''

class Gene (object):
    """
    Gene object holds further information about this gene.
    """

    _GENE_INFO_REGEXP1 = re.compile(r'gi\|(\d+).*')
    _GENE_INFO_REGEXP2 = re.compile(r'gi(\d+)')
    
    try:
        _GENE_INFO = eval(open('gene_info.txt').read())
    except:
        _GENE_INFO = {}

    def __init__ (self, idString, geneInfo=""):
        """
        Initiates a Gene object.
        idString - in the format of gi### or gi|<giNumber>|xxx|<accNumber>
        geneInfo - Information about this gene. Any string
        """
        m = self._GENE_INFO_REGEXP1.match(idString)
        if not m:
            m = self._GENE_INFO_REGEXP2.match(idString)
        if m:
            self.gi = long(m.group(1))
            self._idString = idString
            self._geneInfo = geneInfo
            if not geneInfo and self in Gene._GENE_INFO:
                self._geneInfo = Gene._GENE_INFO[self]
        else:
            raise Exception('Incorrect input: idString is not match')

    def __hash__ (self):
        return self.gi.__hash__()

    def __eq__ (self, other):
        return self.__hash__() == other.__hash__()

    def __repr__ (self):
        return "Gene('gi|%ld')"%self.gi

    __str__ = __repr__
    
    def getFullId (self):
        return self._idString

    def getInfo (self):
        return self._geneInfo
        
class FusedChannel (dict):
    """
    FusedChannel is a class that seeks for a substrate channeling
    mediated using fusion of two genes.
    This tool can get a list of whole genome from a specific organism
    and retrieve genes of external organism which have two or more
    genes from the current organism in each one.
    """
    
    # For other parameteres see http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml
    BLAST_TYPE = "blastp" # See http://www.ncbi.nlm.nih.gov/BLAST/blast_program.shtml
    DATABASE = "nr" # See http://www.ncbi.nlm.nih.gov/BLAST/blast_databases.shtml
    EVAL_CUTOFF = 1e-8
    ORGANISMS = ["Yeast"]
    MAX_OVERLAP = 0 # Maximum allowed overlaps between two genes over the fused
    
    _SEARCH_MODE = enum('NEW_ITER','ITER_ID','END_ITER','HIT_ID')

    def __init__ (self):
        """
        Note! After initiation, you may want to edit the constants, specially DATABASE and ORGANISMS
        """
        dict.__init__(self)
        self._usedGenes = set()
        self._pairs = {} # (original gene, fused gene) => (start_hit, end_hit) 
        
    def __getitem__ (self, item):
        if item not in self:
            dict.__setitem__(self,item,set())
        return dict.__getitem__(self,item)
        
    def match (self, giList, trace=False):
        """
        Gets list of whole genome of a specific organism to search over.
        The format of those genes are GI numbers.
        If trace is set, then a progress of the program will be returned.
        """
        logging.info('Starting a new match with %d genes.'%len(giList))
        if trace:
            counter = 0
            last_percentage = -1
        if self.ORGANISMS:
            orglist = " AND ".join([organism+"[Organism]" for organism in self.ORGANISMS])
        else:
            orglist = "(none)"
        for gi in giList:
            try:
                blastStream = NCBIWWW.qblast(self.BLAST_TYPE, self.DATABASE,gi,entrez_query=orglist,expect=self.EVAL_CUTOFF)
            except Exception, e:
                logging.error('BLAST exception at %s: %s'%(gi,e.message))
                continue
            self.parseBlast (blastStream)
            if trace:
                counter += 1
                percent = int(float(counter)/len(giList)*100)
                if percent > last_percentage:
                    print "%d%%"%(percent),
                    last_percentage = percent

    def parseBlast (self, stream):
        " Inner method, parses blast stream output "
        xml = NCBIXML.read(stream)
        currHitter = Gene(xml.query_id, xml.query)
        # Check whether it's already queried
        if currHitter in self._usedGenes:
            logging.warning('Duplication in gene %s. Given twice. Abort.'%currHitter)
        else:
            self._usedGenes.add(currHitter)
            for hit in xml.alignments:
                currGene = Gene(hit.hit_id, hit.hit_def)
                self[currGene].add(currHitter)
                hsp = hit.hsps[0]
                self._pairs[(currHitter,currGene)] = (hsp.sbjct_start,hsp.sbjct_end)
            logging.info('Gene: %s ended with %d results'%(currHitter,len(xml.alignments)))
    
    def getIntervals (self, fusedGene):
        " Returns the positions of all the hits of the blast "
        intervals = [self._pairs[(g, fusedGene)] for g in self[fusedGene]]
        intervals.sort()
        return intervals
            
    def getCandidatesForChanneling (self):
        " Return genes that are candidates to be fused genes"
        geneList = [gene for gene in self if len(self[gene])>1]
        candidates = []
        for gene in geneList:
            intervals = self.getIntervals(gene)
            overlap = min([max(intervals[i][1]-intervals[-1][0],0) for i in range(len(intervals)-1)])
            if overlap <= self.MAX_OVERLAP:
                candidates.append(gene)
        logging.info("Found %d candidates for channeling"%len(candidates))
        return candidates
    
    def clusterCandidates (self, candidateList):
        """
        Clusters similar candidates using the list of genes that hits different candidadtes,
        i.e. if two candidates where hitted using the same genes, they will cluster together.
        This function return a list of groups, where each group is a cluster of several genes.
        """
        cluster = {}
        for gene in candidateList:
            group = tuple(self[gene])
            if group in cluster:
                cluster[group].append(gene)
            else:
                cluster[group] = [gene]
        logging.info("Found %d clusters of candidates"%len(cluster))
        return cluster.values()
    
    def hasPair (self, gene, list_of_pairs):
        """
        Given list of pairs, the method checks whether exists a pair of genes that where fused
        together on a given gene.
        """
        for pair in list_of_pairs:
            if pair[0] in self[gene] and pair[1] in self[gene]:
                return True
        return False
    
    def countAllPairs (self, gene, list_of_pairs):
        """
        Given list of pairs, the method counts all possible pairs of genes that where
        found fused on a given gene are in the list of pairs.
        """
        gene_list = list(self[gene])
        counter = 0
        for i in range(len(gene_list)):
            for j in range(i+1,len(gene_list)):
                g1 = gene_list[i]
                g2 = gene_list[j]
                if (g1,g2) in list_of_pairs or (g2,g1) in list_of_pairs:
                    counter += 1
        return counter
    
    def findAllFusedPairs (self, candidates):
        """
        Return a list of all the pairs that where found fused on genes.
        The method gets list of candidates which is checked upon it.
        """
        pairs = set()
        counter = 0
        for gene in candidates:
            fused_genes = list(self[gene])
            intervals = [self._pairs[(g, gene)] for g in fused_genes]
            for i,first in enumerate(intervals):
                for j,second in enumerate(intervals[i+1:]):
                    if (first[1]-second[0]<=self.MAX_OVERLAP) or (second[1]-first[0]<=self.MAX_OVERLAP):
                        pairs.add((fused_genes[i],fused_genes[j+i+1]))
                        counter += 1
        # Remove pairs that apear twice
        new_pairs = set()
        for p1,p2 in pairs:
            if (p2,p1) not in new_pairs and p1!=p2:
                new_pairs.add((p1,p2))
        logging.info('FindAllFusedPairs: Over %d candidates, found %d pairs in which %d are different'
                     %(len(candidates), counter, len(new_pairs)))
        return new_pairs
        


def examples (number):
    f = FusedChannel()
    if number == 1:
        # Example one - Using list of gis from purine metabolism
        f.match(open('channeling/examples/eco_gis.txt').read().split('\n'),1)
    elif number == 2:
        # Example two - Using already found stream
        stream = open('examples/blastRes.xml')
        f.parseBlast(stream)
    elif number == 3:
        # Example three - GyrA and GyrB in e.coli fused into topoisomerase2 in yeast
        f.match(['GI|172073090','GI|91093601'])
            
        
if __name__ == '__main__':
    examples(-1) #None
    