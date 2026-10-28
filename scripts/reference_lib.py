"""Coordinate and sequence primitives. Genomic intervals are 1-based inclusive; transcript slices 0-based half-open."""
import hashlib
from Bio.Seq import Seq

def sha(sequence):return hashlib.sha256(sequence.encode()).hexdigest()
def positions(exons,strand='+'):
    out=[]
    for a,b in sorted(exons,reverse=strand=='-'):
        out.extend(range(a,b+1) if strand=='+' else range(b,a-1,-1))
    return out

def splice(genome,origin,exons,strand='+'):
    out=''.join(genome[a-origin:b-origin+1] for a,b in sorted(exons))
    return str(Seq(out).reverse_complement()) if strand=='-' else out

def blocks_for_slice(pos,start,end):
    selected=pos[start:end];blocks=[]
    for g in selected:
        if blocks and abs(g-blocks[-1][-1])==1:blocks[-1].append(g)
        else:blocks.append([g])
    return sorted([[min(b),max(b)] for b in blocks])

def intervals_from_positions(values):
    out=[]
    for g in sorted(set(values)):
        if out and g==out[-1][1]+1:out[-1][1]=g
        else:out.append([g,g])
    return out

def translate(nt):return str(Seq(nt).translate(table=1))

def find_orfs(seq,min_aa=8):
    """Every AUG start to its first in-frame standard-code stop; includes nested starts. No non-AUG or no-stop calls."""
    for start in range(len(seq)-2):
        if seq[start:start+3]!='ATG':continue
        for stop in range(start+3,len(seq)-2,3):
            if seq[stop:stop+3] in ('TAA','TAG','TGA'):
                if (stop-start)//3>=min_aa:yield start,stop+3
                break
