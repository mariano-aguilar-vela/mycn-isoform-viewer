"""Fail on any current structure/sequence/annotation identity mismatch."""
import argparse,json
from pathlib import Path
from reference_lib import positions,splice,blocks_for_slice,translate,sha,find_orfs
CLASSES={'DATABASE_ANNOTATED','PRIMARY_LITERATURE_SUPPORTED','LAB_ANALYSIS_REPRODUCED','LAB_ANALYSIS_INHERITED_NOT_REPRODUCED','COMPUTATIONAL_PREDICTION','OWNER_SUPPLIED_NOMENCLATURE','UNSUPPORTED'}
def validate(data,loc):
    import jsonschema
    schema=json.loads((Path(__file__).resolve().parents[1]/'data/schemas/reference.schema.json').read_text())
    jsonschema.validate(data,schema)
    ids=set();assert data['status']=='CURRENT'
    for key,t in data['transcripts'].items():
        assert t['id']==t['transcript_id'] and t['id'] not in ids;ids.add(t['id'])
        assert t['mlen']==len(t['sequence'])==sum(b-a+1 for a,b in t['exons'])
        assert splice(loc['seq'],loc['start'],t['exons'],t['strand'])==t['sequence']
        assert sha(t['sequence'])==t['sequence_sha256']
        pos=positions(t['exons'],t['strand']);assert len(set(pos))==len(pos)
        cds=set(g for a,b in t['cds'] for g in range(a,b+1));utr=set(g for a,b in t['utr'] for g in range(a,b+1));assert not cds.intersection(utr)
        if t['protein_accession']:
            assert translate(splice(loc['seq'],loc['start'],t['cds'],t['strand']))==t['protein_sequence']
            assert len(t['protein_sequence'])==t['protein_length']
            a=t['cds_transcript_start'];b=t['cds_transcript_end'];assert t['sequence'][a:a+3]=='ATG' and t['sequence'][b:b+3] in ['TAA','TAG','TGA']
            assert any(o['transcript_records'].get(key,{}).get('relationship')=='ANNOTATED_CDS' for o in data['orfs'].values()),key
        assert t['evidence_class']=='DATABASE_ANNOTATED'
        if key=='XM':assert t['source'].startswith('RefSeq') and t['gencode_release']=='not a GENCODE transcript'
        else:assert t['source']==f"GENCODE v{data['gencode_release']}" and t['gencode_release']==f"v{data['gencode_release']}"
    observed={key:set() for key in data['transcripts']}
    for o in data['orfs'].values():
        assert len(o['nt'])==3*(o['aa_len']+1) and translate(o['nt'])==o['aa']+'*'
        assert o['nt'].startswith('ATG') and o['nt'][-3:] in ['TAA','TAG','TGA']
        assert sha(o['aa'])==o['property_provenance']['input_sequence_sha256']==o['protein_sequence_sha256']
        assert o['evidence_class'] in CLASSES and o['evidence_class']!='UNSUPPORTED'
        for t in o['carriers']:
            tx=data['transcripts'][t];r=o['transcript_records'][t];a=r['start']-1;b=r['end'];assert tx['sequence'][a:b]==o['nt'];assert blocks_for_slice(positions(tx['exons'],tx['strand']),a,b)==o['blocks'];observed[t].add((a,b))
        for d in o['domains']:assert d['input_sequence_sha256']==sha(o['aa']) and 1<=d['start']<=d['end']<=o['aa_len']
    for key,t in data['transcripts'].items():assert observed[key]==set(find_orfs(t['sequence'],8)),key
    assert set(data['products'])=={'nmyc','dmycn','mycnot','musep'}
    for product in data['products'].values():
        o=data['orfs'][product['orf_id']];assert o['aa']==product['protein_sequence'];assert sha(o['aa'])==product['protein_sequence_hash'];assert product['primary_transcript'] in o['carriers'];assert o['color']==product['color']
    features=data['relationships']['features']
    for unit,count in data['relationships']['counts'].items():
        chosen=[f for f in features if f['unit']==unit];assert count['all']==len(chosen);assert count['unique']==sum(f['relationship']=='unique' for f in chosen);assert count['shared']==sum(f['relationship']=='shared' for f in chosen)
    assert data['relationships']['counts']['protein_sequence']['unique']>0 and data['relationships']['counts']['protein_sequence']['shared']>0
    return {'transcripts':len(ids),'orfs':len(data['orfs']),'status':'PASS','coordinate_conventions':{'genomic':'1-based inclusive','transcript_public':'1-based inclusive','python_slices':'0-based half-open'},'checks':['complete enumerated AUG ORFs','exon sum','genomic splice equals official transcript','CDS translation equals official protein','start/stop','domain sequence hash','source release identity','RefSeq separation','carrier blocks','UTR/CDS nonoverlap']}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',default='data/current/reference.json');p.add_argument('--locus',default='mycn_locus_seq.json');a=p.parse_args();print(json.dumps(validate(json.loads(Path(a.reference).read_text()),json.loads(Path(a.locus).read_text())),indent=2))
