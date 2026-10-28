"""Curated product identity and explicitly scoped sequence/structural relationship layers."""
from pathlib import Path
from collections import defaultdict,Counter
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def digest(s):return hashlib.sha256(s.encode()).hexdigest()
def augment(data):
    registry=json.loads((ROOT/'data/curation/products.v1.json').read_text());products={}
    for product in registry['products']:
        p=dict(product);matches=[(key,o) for key,o in data['orfs'].items() if digest(o['aa'])==p['protein_sequence_hash'] and p['primary_transcript'] in o['carriers']]
        if len(matches)!=1:raise ValueError('Curated product has no unambiguous current ORF: '+p['product_id'])
        ident,o=matches[0];assert len(o['aa'])==p['protein_length'];p.update({'orf_id':ident,'build_orf_id':o['build_orf_id'],'protein_sequence':o['aa'],'nucleotide_sequence':o['nt'],'genomic_blocks':o['blocks'],'source_transcripts':[{'key':t,'canonical_name':data['transcripts'][t]['display_name'],'accession':data['transcripts'][t]['transcript_id'],'lab_alias':data['transcripts'][t]['lab_alias'],'translation_start':o['transcript_records'][t]['start'],'translation_end':o['transcript_records'][t]['end']} for t in o['carriers']],'active_status':'ACTIVE','release_last_reviewed':data['review_date'],'orf_evidence_class':o['evidence_class'],'domains':o['domains']})
        products[p['product_id']]=p
    for o in data['orfs'].values():
        mapped=[p for p in products.values() if p['orf_id']==str(o['id'])];o['product_ids']=[p['product_id'] for p in mapped];o['product_names']=[p['display_name'] for p in mapped]
        if mapped:o['color']=mapped[0]['color']
        o['relationship_class']='ORF_UNIQUE_TO_TRANSCRIPT' if len(o['carriers'])==1 else 'ORF_SHARED_ACROSS_TRANSCRIPTS';o['relationship_scope']='Current eight MYCN transcript models; exact genomic blocks and nucleotide sequence'
    # Project a canonical feature only when its complete amino-acid substring has one exact match.
    canonical=products['nmyc'];delta=products['dmycn'];mapped_domains=[]
    for d in canonical['domains']:
        if d['type'] not in ['Domain','Region'] or d['description'] not in ['bHLH','Leucine-zipper']:continue
        seq=canonical['protein_sequence'][d['start']-1:d['end']];offset=delta['protein_sequence'].find(seq)
        if offset>=0 and delta['protein_sequence'].count(seq)==1:mapped_domains.append({**d,'start':offset+1,'end':offset+len(seq),'source_coordinates':[d['start'],d['end']],'source_product':'N-Myc','evidence_class':'SEQUENCE_RELATIONSHIP_ONLY','input_sequence_sha256':delta['protein_sequence_hash'],'projection':'Exact complete-feature sequence match, not direct experimental domain validation in ΔMYCN'})
    delta['domains']=mapped_domains
    features=[];txs=data['transcripts'];tx_universe=[t['transcript_id'] for t in txs.values()];prod_universe=list(products)
    def add(unit,classification,carriers,**fields):
        identifier=f'{unit}-{sum(f["unit"]==unit for f in features)+1:04d}';features.append({'feature_id':identifier,'unit':unit,'classification':classification,'relationship':'unique' if '_UNIQUE_' in classification else 'shared','carriers':carriers,'universe':prod_universe if unit=='protein_sequence' else tx_universe,'evidence_status':'SEQUENCE_RELATIONSHIP_ONLY',**fields})
    for ident,o in data['orfs'].items():add('orf',o['relationship_class'],o['carriers'],label='ORF '+ident+(' · '+', '.join(o['product_names']) if o['product_names'] else ''),orf_id=ident,sequence=o['aa'],coordinates=o['blocks'],scope='Exact ORF model (genomic blocks + nucleotide sequence) among the eight current displayed transcripts; not protein-sequence or proteome uniqueness',minimum_length=8,length_unit='aa',matching='exact sequence and genomic blocks; transcript reading frame retained')
    def genomic_sequence(a,b,carrier):
        from reference_lib import positions
        tx=txs[carrier];mapping=dict(zip(positions(tx['exons'],tx['strand']),tx['sequence']))
        return ''.join(mapping[i] for i in range(a,b+1))
    for unit,field in [('exon','exons'),('junction','junctions')]:
        by=defaultdict(list)
        for key,t in txs.items():
            for interval in t[field]:by[tuple(interval)].append(key)
        for interval,carriers in sorted(by.items()):
            cat='TRANSCRIPT_'+('UNIQUE_' if len(carriers)==1 else 'SHARED_')+('EXON' if unit=='exon' else 'JUNCTION');add(unit,cat,carriers,label='chr2:'+str(interval[0])+'–'+str(interval[1]),coordinates=list(interval),sequence=genomic_sequence(*interval,carriers[0]) if unit=='exon' else '',scope='Exact '+('inclusive exon boundary pair' if unit=='exon' else 'splice donor/acceptor coordinate pair')+' among the eight current displayed transcripts; not substring uniqueness',minimum_length=None,length_unit='bp',matching='exact genomic coordinate pairs on the annotated + strand')
    boundaries=sorted({x for t in txs.values() for a,b in t['exons'] for x in [a,b+1]});regions=[]
    for a,b in zip(boundaries,boundaries[1:]):
        carriers=[key for key,t in txs.items() if any(lo<=a and b-1<=hi for lo,hi in t['exons'])]
        if not carriers:continue
        if regions and regions[-1]['carriers']==carriers and regions[-1]['end']==a-1:regions[-1]['end']=b-1
        else:regions.append({'start':a,'end':b-1,'carriers':carriers})
    for region in regions:add('nucleotide_region','TRANSCRIPT_'+('UNIQUE_' if len(region['carriers'])==1 else 'SHARED_')+'NUCLEOTIDE_REGION',region['carriers'],label='chr2:'+str(region['start'])+'–'+str(region['end']),coordinates=[region['start'],region['end']],sequence=genomic_sequence(region['start'],region['end'],region['carriers'][0]),scope='Exonic genomic-position coverage among eight current transcripts; not a genome-wide nucleotide substring comparison',minimum_length=1,length_unit='bp',matching='exact position coverage on the annotated strand; reverse-complement searching is not used')
    k=8;windows=defaultdict(set)
    for key,p in products.items():
        for i in range(len(p['protein_sequence'])-k+1):windows[p['protein_sequence'][i:i+k]].add(key)
    for key,p in products.items():
        seq=p['protein_sequence'];runs=[]
        for i in range(len(seq)-k+1):
            carriers=sorted(windows[seq[i:i+k]])
            if runs and runs[-1]['carriers']==carriers:runs[-1]['last_start']=i
            else:runs.append({'first_start':i,'last_start':i,'carriers':carriers})
        for run in runs:
            start=run['first_start']+1;end=run['last_start']+k;carriers=run['carriers'];add('protein_sequence','PRODUCT_'+('UNIQUE_' if len(carriers)==1 else 'SHARED_')+'SEQUENCE_REGION',carriers,label=p['display_name']+' aa '+str(start)+'–'+str(end),product_id=key,coordinates=[start,end],window_start_range=[run['first_start']+1,run['last_start']+1],window_count=run['last_start']-run['first_start']+1,sequence=seq[start-1:end],scope='Exact 8-aa windows among the four displayed MYCN products only; not the human proteome and not experimental peptide evidence',minimum_length=k,length_unit='aa',matching='exact, case-sensitive standard amino-acid sequence; no substitutions, reverse complement or nucleotide-frame comparison',display_note='Consecutive window starts with identical product carrier sets are merged; neighbouring displayed regions may overlap by 7 aa')
        p['unique_sequence_features']=[f['feature_id'] for f in features if f['unit']=='protein_sequence' and f['product_id']==key and f['relationship']=='unique'];p['shared_sequence_features']=[f['feature_id'] for f in features if f['unit']=='protein_sequence' and f['product_id']==key and f['relationship']=='shared']
    counts={}
    for unit in ['orf','exon','junction','nucleotide_region','protein_sequence']:counts[unit]={'all':sum(f['unit']==unit for f in features),'unique':sum(f['unit']==unit and f['relationship']=='unique' for f in features),'shared':sum(f['unit']==unit and f['relationship']=='shared' for f in features)}
    data.update({'schema_version':'2.1.0','application_version':json.loads((ROOT/'package.json').read_text())['version'],'structure_reference_build':'gencode-50_refseq-NM_001446166.1_mycn-1','reference_build':'gencode-50_refseq-NM_001446166.1_mycn-2','product_registry_version':registry['version'],'products':products,'relationships':{'version':'1.0.0','features':features,'counts':counts,'protein_window_length':k,'peptide_evidence_categories':'Not implemented: no retained experimental peptide-discrimination evidence is inferred from these sequence comparisons'},'evidence':{'curated_product_detection':'See each product; database annotation, literature detection and unresolved function are separate','relationship_calculations':'Exact current sequences/coordinates; comparison universes and length thresholds are stored per feature'}})
    return data
if __name__=='__main__':
    from build_mycn_reference import write
    from reference_lib import sha
    from build_public_assets import generate
    p=ROOT/'data/current/reference.json';data=augment(json.loads(p.read_text()));write(p,data);write(ROOT/'data/releases/gencode_v50/reference.json',data);manifest=json.loads((ROOT/'data/current/manifest.json').read_text());manifest.update({k:data[k] for k in ['application_version','schema_version','reference_build','structure_reference_build','product_registry_version']});manifest['reference_sha256']=sha(p.read_text());write(ROOT/'data/current/manifest.json',manifest);write(ROOT/f'data/manifests/{data["reference_build"]}.json',manifest);generate(data);print('Products mapped:',[(p['display_name'],p['orf_id']) for p in data['products'].values()]);print('Relationship counts:',data['relationships']['counts'])
