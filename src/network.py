"""Simple undirected PPI graphs and observed-quantile degree matching."""
from collections import Counter
import numpy as np

def lcc_metrics(adj):
 """Bitset traversal, preserving first component in a largest-size tie."""
 remaining=(1<<len(adj))-1; best=0; edges=0
 while remaining:
  frontier=remaining & -remaining; component=0
  while frontier:
   bit=frontier & -frontier; frontier^=bit
   component|=bit; remaining&=~bit
   frontier|=adj[bit.bit_length()-1] & remaining & ~component
  size=component.bit_count()
  if size>best:
   best=size; degree_sum=0; bits=component
   while bits:
    bit=bits & -bits; bits^=bit
    degree_sum+=(adj[bit.bit_length()-1] & component).bit_count()
   edges=degree_sum//2
 return best,edges

def networks(genes, interaction_map, shared_masks):
 n=len(genes); bits=[1<<i for i in range(n)]
 second=[0]*n; combined=[0]*n
 full=[interaction_map[g] for g in genes]; shared=[shared_masks[g] for g in genes]
 for i,g in enumerate(genes):
  for j in range(i+1,n):
   s=bool(shared[i] & shared[j])
   if s:
    second[i]|=bits[j]; second[j]|=bits[i]
   if s or genes[j] in full[i] or g in full[j]:
    combined[i]|=bits[j]; combined[j]|=bits[i]
 return {'second_only':lcc_metrics(second),'first_and_second':lcc_metrics(combined)}

def sampling_setup(observed,background,interaction_map):
 degrees=[len(interaction_map[g]) for g in observed]
 edges=np.quantile(degrees,np.linspace(0,1,21))
 def which(degree): return min(int(np.searchsorted(edges[1:],degree,side='left')),19)
 counts=Counter(which(d) for d in degrees)
 pools={b:[g for g in background if which(len(interaction_map[g]))==b] for b in counts}
 for b,c in counts.items():
  if len(pools[b])<c: raise ValueError(f'Bin {b}: {len(pools[b])} unique genes for {c} slots')
 return [(pools[b],counts[b]) for b in sorted(counts)],edges.tolist(),dict(counts)

