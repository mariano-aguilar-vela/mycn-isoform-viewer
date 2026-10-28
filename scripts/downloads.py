"""Derive public coordinate downloads from the viewer's canonical metadata."""
from pathlib import Path
import argparse,json
from urllib.parse import quote
ROOT=Path(__file__).resolve().parents[1]
def outputs(meta):
    txs=meta['transcripts']; starts=[e[0] for t in txs.values() for e in t['exons']];ends=[e[1] for t in txs.values() for e in t['exons']]
    gff=['##gff-version 3','##sequence-region chr2 15938000 15947200','chr2\tMYCN_viewer\tgene\t%d\t%d\t.\t+\t.\tID=MYCN;Name=MYCN'%(min(starts),max(ends))]
    def attrs(**values):return ';'.join(k+'='+quote(str(v),safe=':-_.') for k,v in values.items())
    def feature(kind,a,b,phase,att):gff.append('\t'.join(map(str,['chr2','MYCN_viewer',kind,a,b,'.','+',phase,att])))
    for key,t in txs.items():
        name='MYCN-'+key
        feature('mRNA',t['exons'][0][0],t['exons'][-1][1],'.',attrs(ID=name,Parent='MYCN',Name=name,transcript_id=t['transcript_id'],source_provenance=t['source']))
        for a,b in t['exons']:feature('exon',a,b,'.',attrs(Parent=name))
        n=0
        for a,b in t['cds']:
            feature('CDS',a,b,(3-n%3)%3,attrs(Parent=name));n+=b-a+1
    bed=['track name="MYCN ORF landscape" itemRgb="On"']
    for key,o in meta['orfs'].items():
        blocks=o['blocks'];start=blocks[0][0]-1;end=blocks[-1][1];color=o['color'].lstrip('#');rgb=','.join(str(int(color[i:i+2],16)) for i in [0,2,4])
        name=f"ORF{key}|{o['cls']}|{o['aa_len']}aa|{o['conserved']}"
        sizes=','.join(str(b-a+1) for a,b in blocks)+',';offsets=','.join(str(a-1-start) for a,b in blocks)+','
        bed.append('\t'.join(map(str,['chr2',start,end,name,0,'+',start,end,rgb,len(blocks),sizes,offsets])))
    return {'mycn_transcripts.gff3':'\n'.join(gff)+'\n','mycn_orfs.bed':'\n'.join(bed)+'\n'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--check',action='store_true');args=p.parse_args()
    for filename,text in outputs(json.loads((ROOT/'mycn_orfs.meta.json').read_text())).items():
        path=ROOT/filename
        if args.write:path.write_text(text)
        elif path.read_text()!=text:raise SystemExit('Stale coordinate download: '+filename)
    print('Coordinate downloads agree with metadata: eight transcripts, phased CDS, all ORFs.')
