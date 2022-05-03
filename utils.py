import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import bisect
import numpy as np 

def de_bruijn_ize(st,k):
    """create a de Bruijn graph. Currently unused."""
    edges=[]
    nodes=set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

def visualize_de_bruijn(st,k):
    """visualize a de Bruijn graph as a dot plot"""
    nodes, edges=de_bruijn_ize(st,k)
    dot_str='digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str+=' %s [label="%s"] ;\n' % (node, node)
    for src, dst in edges:
        dot_str+=' %s -> %s ;\n' %(src, dst)
    return dot_str+'}\n'

def translate(seq):
    """translate DNA sequences to amino acids"""
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i + 3]
        protein+= table[codon]
    return protein

def complement(seq):
    """gets complement sequence of DNA"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    for element in bases:
        if element not in complement:
            print(element) 
        letters = [complement[base] for base in element] 
        return ''.join(letters)

def reverse_complement(seq):
    """gets reverse complement of sequence"""
    return complement(seq[::-1])


def create_kmer_dictionary(tl, k):
    """ccreates dictionary from all possible kmers"""
    dict={}
    for t in tl:
        for m in range(2):
            for j in range(0,3):
                st=translate(t[j:])
                for i in range(len(st)-k):
                    if st[i:i+k] in dict:
                        dict[st[i:i+k]]+=1
                    else:
                        dict[st[i:i+k]]=1

            t=reverse_complement(t)
    dict={k: v for k, v in sorted(dict.items(), key=lambda item: item[1], reverse=True)}
    return dict

def read_fastq(filenames):
    flattened_filtered=[]
    """reads fastq file, splits at all non-DNA letters and returns only sequence"""
    lines=[]
    for file in filenames:
        with open(file) as f:
            lines.extend(f.readlines())
    split_lines=[re.split('[^GTCA]', l) for l in lines[1::4]]
    flattened = [val for sublist in split_lines for val in sublist]
    flattened_filtered.extend(list(filter(None, flattened)))
    return flattened_filtered

def flatten(A):
    """general util function to flatten nested lists"""
    rt = []
    for i in A:
        if isinstance(i,list): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

def find_path(dict, str,k,score):
    """recursive function to transverse a sequence, does not actually work"""
    # if last letter is stop sequence, return substring and remove from all kmers from dict
    if str[-1]=='_':
        # TODO: remove instances from dict
        return str,score/(len(str)-k+1)
    # else, find all possible paths
    else:
        next={key[-1]:val for key,val in dict.items() if key.startswith(str[-k+1:])}
        return [find_path(dict, str+n,k,score+s) for n,s in next.items()]

class contig:
    """represents the strand being transversed """
    def __init__(self, parent, kmer, abundance, min_length=10):
        self.min_length=min_length
        self.parent=parent
        self.abundance=abundance
        self.string=kmer
        self.length=len(kmer)
        self.can_end=False
        if bool(parent):
            self.length=parent.length+1-self.length
            self.can_end=parent.can_end

    def __str__(self):
        return "Str:"+str(self.string)+"|Ab:"+str(self.abundance)+"|Length:"+str(self.length)

    def is_stop(self):
        return self.string[-1]=='_'

    def is_start(self):
        return self.string[0]=='M'

    def extend(self, str,ab):
        self.string+=str
        self.length+=1

    def extend_back(self, ext):
        self.string=ext.string[0]+self.string
        self.abundance=+ext.abundance
        self.length+=1
        self.can_end= (self.can_end or (self.is_start() and self.length>self.min_length))

def clean_dict(dict, rate_top_kmers=0.001, rate_to_sig=0.05):
    key=list(dict.keys())
    val=list(dict.values())

    kmer_count_threshold=np.mean(val[:round(len(val)*rate_top_kmers)])*rate_to_sig
    index_to_keep=bisect.bisect_left(val[::-1], kmer_count_threshold)

    dict_clean={k:v for k,v in zip(key[0:-index_to_keep], val[0:-index_to_keep])}

    return dict_clean


def assemble_proteins(dict, k=8, min_length=30):
    start_sites=[contig(None,k, dict[k], min_length) for k in dict.keys() if k.endswith('_') and not '_' in k[:-1]]
    start_sites.sort(key=lambda x: x.abundance)
    protein_dict={}

    for i in range(len(start_sites)):
        queue=[]
        strand=start_sites.pop()
        count=0

        while bool(strand):
            # get all possible extensions backwards (must end with the previous sequence and not start with a stop codon)
            exts=[contig(strand, key, val, min_length) for key,val in dict.items() if (key.endswith(strand.string[:k-1]) and key[0]!='_')]
            # if there is one extension, add it on
            if len(exts)==1:
                strand.extend_back(exts[0])
                count+=1
                continue
            # otherwise if there is more than one possible extension add to queue, sort, and then extend with the most abundant kmer
            elif bool(exts):
                exts.sort(key=lambda x: x.abundance)
                strand=exts.pop()
                queue.extend(exts)
                queue.sort(key=lambda x: x.abundance)
            # if there is nothing, and there is still not stop point, try a different branch
            elif bool(queue) and not strand.can_end:
                strand=queue.pop()
            # then if there are no options whatsoever break
            else:
                break
        if strand.can_end:
            # remove amino acids from sequence
            protein=""
            parstrand=strand
            while parstrand.parent!=None:
                protein+=parstrand.string[:-k+1]
                parstrand=parstrand.parent
            protein+=parstrand.string

            protein='M'+protein.split('M',1)[1]
            protein_dict[protein]=sum([dict.pop(protein[j:j+k]) for j in range(0,len(protein)-k+1)])/(len(protein)-k+1)
            protein_dict={k: v for k, v in sorted(protein_dict.items(), key=lambda item: item[1], reverse=True)}


    return protein_dict

def write_seq(name,proteins):
    res = list(map(lambda i: i[ : -1], list(proteins.keys())))
    records=[]
    for (index, seq) in enumerate(res):
        records.append(SeqRecord(Seq(seq), str(index)))
    SeqIO.write(records, name+'.faa', "fasta")

    return res