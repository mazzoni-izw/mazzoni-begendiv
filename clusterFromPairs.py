import sys
import argparse

import networkx as nx
from Bio import SeqIO

def parseAllPairsTabOut(filePath):
    """
    Parse vsearch allpairs output and create a dictoanry that maps each 
    sequence to a list of tuples of matching sequences and match identities.
    
    Returns: {queryId: (matchingTargetId, matchIdentity), ...}
    """
    matches = {}
    for line in open(filePath, "rt"):
        query, target, ident = line.strip("\n").split("\t")
        try:
            matches[query].append((target, float(ident)))
        except KeyError:
            matches[query] = [(target, float(ident))]
        try:
            matches[target].append((query, float(ident)))
        except KeyError:
            matches[target] = [(query, float(ident))]
    return matches

def cluster(seqIds, pairs, prefix="clu", thr=99.0): # CHANGE IDENTITY THRESHOLD
    """
    Single linkage of match pairs implemented with the networkx library by 
    creating a adjacency graph and using the connetcetd_components function 
    to get clusters.
    
    Sequences are considered adjacent if the identity is equal or higher than
    the threshold (default: 99%).
    
    Returns a tuple with three dictonary 
        1. mapping cluster names to sequence Ids of their members
        2. mapping cluster names to minimum identity between clusters sequences
        3. mapping cluster names to maximum identity between clusters sequences
    ({cluId: [seqId, ...], ...},
     {cluId: [minIdent, ...], ...},
     {cluId: [maxIdent, ...], ...}
    )
    """
    graph = nx.Graph()
    #create nodes
    graph.add_nodes_from(seqIds)
    #create edges
    for query, matches in pairs.items():
        for target, ident in matches:
            if ident >= thr:
                graph.add_edge(query, target, weight=ident)
    rv = ({}, {}, {})
    for i,c in enumerate(sorted(nx.connected_components(graph), key=len, reverse=True)):
        #TODO start numbering from 1?
        clu = graph.subgraph(c)
        rv[0]["%s_%i" % (prefix, i)] = c
        if len(c) > 1:
            rv[1]["%s_%i" % (prefix, i)] = min([e[2] for e in clu.edges.data("weight")])
            rv[2]["%s_%i" % (prefix, i)] = max([e[2] for e in clu.edges.data("weight")])
        else:
            rv[1]["%s_%i" % (prefix, i)] = 0
            rv[2]["%s_%i" % (prefix, i)] = 0
        
    return rv

def main(vsearchPairsFile, fastqFile):

    pairs = parseAllPairsTabOut(vsearchPairsFile)

#    pairsFile = vsearchPairsFile.rsplit(".", 1)[0] + "_pairs.tsv"
#    with open(pairsFile, "wt") as out:
#        for query, matches in pairs.items():
#            for target, ident in matches:
#                out.write("%s\t%s\t%s\n" % (query, target, ident))

    seqSize = {}
    seqLen = {}
    for rec in SeqIO.parse(fastqFile, "fastq"):
        seqSize[rec.id] = int(rec.id.split(";")[-1].split("=")[1])
        seqLen[rec.id] = len(rec)

    sample = vsearchPairsFile.split("/")[-1].split(".")[0]

    clu99, clu99MinIdent, clu99MaxIdent = cluster(seqSize.keys(), pairs, sample, 90.0) # CHANGE IDENTITY THRESHOLD


    cluFile = vsearchPairsFile.rsplit(".", 1)[0] + "_clusters.tsv"
    with open(cluFile, "wt") as out:
        for cluId, cluMem in clu99.items():
            out.write("%s\t%i\t%s\n" % (cluId, len(cluMem), " ".join(cluMem)))

    cluSizeDist = {} #distribution of cluster sizes
    for cId, cluMem in clu99.items():
        try:
            cluSizeDist[len(cluMem)] += 1
        except KeyError:
            cluSizeDist[len(cluMem)] = 1

    for size, num in sorted(cluSizeDist.items(), key=lambda x: x[0]):
        sys.stderr.write("%i\t%i\n" % (size, num))


parser = argparse.ArgumentParser(description="builds clusters based on pairwise identites from vsearch")
parser.add_argument("-vs", "--vsearch", help="vsearch pairwise identity tabout file ")
parser.add_argument("-fq", "--fastq", help="fastq file used as an input for vsearch")
args = parser.parse_args()

if __name__ == "__main__":
    main(args.vsearch, args.fastq)
