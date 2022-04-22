def de_bruijn_ize(st,k):
    edges=[]
    nodes=set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

def visualize_de_bruijn(st,k):
    nodes, edges=de_bruijn_ize(st,k)
    dot_str='digraph "DeBruijn graph" {\n'
    for node in nodes:
        dot_str+=' %s [label="%s"] ;\n' % (node, node)
    for src, dst in edges:
        dot_str+=' %s -> %s ;\n' %(src, dst)
    return dot_str+'}\n'

# translate into all possible frames (6)
# can jump into a different frame only between reads
# 

nodes, edges=de_bruijn_ize('ACGCGTCG',3)
print(nodes)
print(edges)


