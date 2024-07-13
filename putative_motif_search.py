import sys
import argparse

#-------------------------parameters start---------------------------------------------------------------------------
version = "v1.0"
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''
==========================================================
This script is designed to quickly screen multiple DNA sequences to identify 
shared kmers (subsequences of length k) among them. By detecting potential 
similar motif sequences, this script aids in narrowing down the research 
scope for subsequent detailed analysis.

Author: Youliang Pan, panyouliang@genomics.cn
Version: v1.0
Date: 2024-07-13, yyyy-mm-dd
==========================================================''')
parser.add_argument('-v', '--version',action='version',version=version)
parser.add_argument('-fa', metavar='fasta',type=str,required=True,help='Accepts input from a FASTA format file')
parser.add_argument('-k',metavar='k-mer length',type=int, default=10, help='Users can specify the length of the kmer, default [10]')
parser.add_argument('-max_error',metavar='tolerance',type=int, default=0, help='Users can set the tolerance rate for kmer matching,\
        default [0]')
parser.add_argument('-share_gene',metavar='share',type=int, default=2, help='Users can specify the minimum number of sequences a motif\
        must be shared in to be output, default [2].')
args = parser.parse_args()
#=========================parameters end=============================================================================



def find_shared_motifs(dicts, min_length, max_error):
    motif_set = {}

    motif_lis = []
    for seq in dicts.values():
        motifs = extract_kmers(seq,min_length,1)
        motif_lis += motifs

    motif_lst = list(set(motif_lis))

    for motif in motif_lst:
        for geneID,seqs in dicts.items():
            matches = find_matches(seqs,motif,max_error,'+')
            matches_rev = find_matches(seqs,rev(motif),max_error,'-')
            merge = matches+matches_rev
            if motif not in motif_set:
                motif_set[motif] = {geneID:merge}
            else:
                if geneID not in motif_set:
                    motif_set[motif][geneID] = merge
                else:
                    motif_set[motif][geneID].append(merge)

    return motif_set


def readfa(F):
    seqd = {}
    name,seq = '',[]
    with open(F,'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == '>':
                if name != '':
                    seqd[name] = ''.join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line.upper())
    seqd[name] = ''.join(seq)

    return seqd

def find_matches(long_read, query_sequence, max_errors, strand):
    matches = []
    for i in range(len(long_read) - len(query_sequence) + 1):
        errors = 0
        for j in range(len(query_sequence)):
            if long_read[i + j] != query_sequence[j]:
                errors += 1
                if errors > max_errors:
                    break
        if errors <= max_errors:
            matches.append((i, i + len(query_sequence) - 1, strand))
    return matches


def extract_kmers(sequence, k=10, step=10):
    kmers = []
    for i in range(0, len(sequence) - k + 1, step):
        kmer = sequence[i:i + k]
        kmers.append(kmer)
    return kmers


def rev(seq):
    base = {'A':'T','T':'A','G':'C','C':'G','N':'N','n':'n','a':'t','t':'a','c':'g','g':'c'}
    seq_list = list(reversed(seq))
    seq_rev = [base[k] for k in seq_list]
    seq_list = ''.join(seq_rev)
    return seq_list

def count_empty_values(d):
    empty_count = 0
    for key, value in d.items():
        if not value or all(not sublist for sublist in value):
            empty_count += 1
    return empty_count


if __name__ == '__main__':
    seqd = readfa(args.fa)
    min_length = args.k
    max_error = args.max_error
    min_share = args.share_gene

    shared_motifs = find_shared_motifs(seqd, min_length, max_error)
    out = open('putative_motif_result.txt','w')

    for k,v in shared_motifs.items():
        empty =  count_empty_values(v)
        if (len(v) - empty) >= min_share:
            out.write(f"motif:{k}, reverse:{rev(k)}, share_gene:{len(v)-empty}\n")
            for geneID,position in v.items():
                for loc in position:
                    if loc != None:
                        out.write(f"{geneID}\t{len(seqd[geneID])}\t{loc[0]}\t{loc[1]}\t{loc[2]}\n")
                    else:
                        out.write(f"{geneID}\n")

    out.close()




