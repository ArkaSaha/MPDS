from os import system, path
import sys
from collections import deque

dataset = sys.argv[1][:-4]
output = sys.argv[2]
gamma = sys.argv[3]
decomp = "output/trussNumbers_gamma={}_prec=100_{}.txt".format(gamma, dataset)
if not path.exists(decomp):
        with open(dataset + ".txt", 'r') as fin, open(dataset + "_tmp.txt", 'w') as fout:
                ll = fin.readline()
                for line in fin:
                        l = line.strip().split()
                        fout.write("{},{},{}\n".format(l[0], l[1], l[2]))
        system("./truss {}_tmp.txt {} {}".format(dataset, gamma, decomp))
        system("rm -f {}_tmp.txt".format(dataset))
truss = set()
with open(decomp, 'r') as f:
        mt = int(f.readline().strip().split('=')[-1])
        for line in f:
                l = line.strip().split('\t')
                u = l[0]
                for e in l[1:-1]:
                        v, d = e.split(' ')
                        if int(d) >= mt:
                                truss.update([u, v])
adj = {}
with open(dataset + ".txt", 'r') as f:
        ll = f.readline()
        for line in f:
                l = line.strip().split()
                if l[0] in truss and l[1] in truss:
                        if l[0] not in adj:
                                adj[l[0]] = set()
                        adj[l[0]].add(l[1])
                        if l[1] not in adj:
                                adj[l[1]] = set()
                        adj[l[1]].add(l[0])
comp = []
visited = set()
for v in truss:
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
                f.write("{}\n".format(mt))
