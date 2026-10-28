"""Generate public sequence/summary/coordinate exports from the current canonical data."""
from pathlib import Path
import argparse,csv,json
from downloads import outputs
ROOT=Path(__file__).resolve().parents[1]
def generate(data):
    dest=ROOT/'data/current';dest.mkdir(parents=True,exist_ok=True)
    def fasta(kind):
        rows=[]
        for key,t in data['transcripts'].items():
            seq=t['sequence'] if kind=='transcripts' else t['protein_sequence']
            if not seq:continue
            name=t['transcript_id'] if kind=='transcripts' else t['protein_accession']
            header=f">{name} display={t['display_name']} lab_alias={t.get('lab_alias','')} source={t['source'].replace(' ','_')} release={t['gencode_release'].replace(' ','_')} assembly={data['assembly']} evidence={t['evidence_class']} app={data['application_version']} build={data['reference_build']}"
            rows.append(header+'\n'+'\n'.join(seq[i:i+60] for i in range(0,len(seq),60))+'\n')
        return ''.join(rows)
    for kind in ['transcripts','proteins']:(dest/(kind+'.fa')).write_text(fasta(kind))
    for kind,key in [('orfs-nt','nt'),('orfs-protein','aa')]:
        (dest/(kind+'.fa')).write_text(''.join(f">{o['build_orf_id']} evidence={o['evidence_class']} build={data['reference_build']}\n{o[key]}\n" for o in data['orfs'].values()))
    with (dest/'transcripts.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['display_name','lab_alias','accession','source','release','assembly','length_nt','exons','cds','protein_accession','protein_length','evidence','app_version','reference_build'])
        for t in data['transcripts'].values():w.writerow([t['display_name'],t.get('lab_alias',''),t['transcript_id'],t['source'],t['gencode_release'],data['assembly'],t['mlen'],json.dumps(t['exons']),json.dumps(t['cds']),t['protein_accession'],t['protein_length'],t['evidence_class'],data['application_version'],data['reference_build']])
    for name,text in outputs(data).items():(ROOT/name).write_text(text)
    from compare_reference_releases import compare
    previous=json.loads((ROOT/'data/releases/legacy-2026-07/reference.json').read_text());(dest/'release-diff.json').write_text(json.dumps(compare(previous,data),indent=2)+'\n')
if __name__=='__main__':generate(json.loads((ROOT/'data/current/reference.json').read_text()));print('Public sequence, summary, coordinate and history assets regenerated')
