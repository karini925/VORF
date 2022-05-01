import re

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

                    #if st[i:i+k-1] in dict:
                    #    if st[i+k] in dict[st[i:i+k-1]]:
                    #        dict[st[i:i+k-1]][st[i+k]]+=1
                    #    else:
                    #        dict[st[i:i+k-1]][st[i+k]]=1
                    #else:
                    #    dict[st[i:i+k-1]]={st[i+k]:1}
            t=reverse_complement(t)
    dict={k: v for k, v in sorted(dict.items(), key=lambda item: item[1], reverse=True)}
    return dict

def read_fastq(filename):
    """reads fastq file, splits at all non-DNA letters and returns only sequence"""
    with open(filename) as f:
        lines = f.readlines()
    split_lines=[re.split('[^GTCA]', l) for l in lines[1::4]]
    flattened = [val for sublist in split_lines for val in sublist]
    flattened_filtered=list(filter(None, flattened))
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
    def __init__(self, parent, kmer, abundance):
        self.parent=parent
        self.abundance=abundance
        self.string=kmer
        self.length=1
        if not parent==None:
            self.length=parent.length+1
    def __str__(self):
        return "Str:"+str(self.string)+"|Ab:"+str(self.abundance)+"|Length:"+str(self.length)

    def is_stop(self):
        return self.string[-1]=='_'

    def extend(self, str):
        self.string+=str
        self.length+=1