#!/usr/bin/env python
# from https://gist.github.com/jhcepas/fd76041149a92ac9ba40
# Author Jaime Huerta-Cepas 
import sys
from collections import defaultdict
from ete3 import Tree

min_support=.95

try:
    t = Tree(sys.argv[1])
except IndexError:
    print('you need to provide a newick tree file as first argument\n\n', file=sys.stderr)
    print('Usage: Newick2FastTreeConstraints.py tree.nw > constraints.fa', file=sys.stderr)
    sys.exit(1)

n2content = t.get_cached_content()
all_leaves = n2content[t]
leaves = []
#for leaf in t:
#    leaves.append(leaf.name)
for node in t.traverse():
    if not node.is_leaf() and not node.is_root():
        if node.support < min_support:
            node.delete()
alg = defaultdict(list)
for n in t.traverse("postorder"):
    if len(n.children) > 1:
        ones = n2content[n]
        for leaf in all_leaves:
            if leaf in ones:
                alg[leaf].append("1")
            else:
                alg[leaf].append("0")

print("Number of leafs:", len(alg), file=sys.stderr)
print("Number of constraints:", set([len(v) for v in list(alg.values())]), file=sys.stderr)
for name, values in alg.items():
    print(">%s\n%s" %(name.name, ''.join(values)))