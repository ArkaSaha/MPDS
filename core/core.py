from os import system, path
from collections import deque
import sys

dataset = sys.argv[1]
output = sys.argv[2]
eta = sys.argv[3]
decomp = "output/(k-eta)-coreNumbers_eta={}_prec=100_{}".format(eta, dataset)
if not path.exists(decomp):
        system("rm -f core.class && javac core.java")
        system("java core {} {}".format(dataset, eta))
c = {}
mc = 0
with open(decomp, 'r') as f:
        for line in f:
                l = line.strip().split()
                v = int(l[0])
                d = int(l[1])
                c[v] = d
                if d > mc:
                        mc = d
core = set()
for v in c:
        if c[v] >= mc:
                core.add(v)
adj = {}
with open(dataset, 'r') as f:
        ll = f.readline()
        for line in f:
                l = line.strip().split()
                u = int(l[0])
                v = int(l[1])
                if u in core and v in core:
                        if u not in adj:
                                adj[u] = set()
                        adj[u].add(v)
                        if v not in adj:
                                adj[v] = set()
                        adj[v].add(u)
comp = []
visited = set()
for v in core:
        if v not in visited:
                c = set()
                q = deque()
                q.append(v)
                while len(q):
                        u = q.popleft()
                        if u not in visited:
                                visited.add(u)
                                c.add(u)
                                if u in adj:
                                        for w in adj[u]:
                                                if w not in visited:
                                                        q.append(w)
                comp.append(c)
with open(output, 'w') as f:
        for s in comp:
                for v in s:
                        f.write("{} ".format(v))
                f.write("{}\n".format(mc))