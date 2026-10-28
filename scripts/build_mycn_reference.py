"""Rebuild the complete current MYCN model from official GENCODE and NCBI inputs."""
import argparse,csv,gzip,json,re,sys
from pathlib import Path
from collections import defaultdict
import Bio
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from reference_lib import sha,positions,splice,blocks_for_slice,intervals_from_positions,translate,find_orfs
ROOT=Path(__file__).resolve().parents[1]

def attributes(text):
    out=defaultdict(list)
    for key,val in re.findall(r'(\w+) "([^"]*)"',text):out[key].append(val)
    return dict(out)
def one(a,k,default=''):return a.get(k,[default])[0]
def write(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def fasta(path,ids):
    result={}
    with gzip.open(path,'rt') as f:
        for r in SeqIO.parse(f,'fasta'):
            tokens=r.id.split('|');match=next((x for x in tokens if x in ids),None)
            if match:result[match]=str(r.seq).upper()
    return result

def build(cache,review):
    sm=json.loads((cache/'source-manifest.json').read_text());release=sm['gencode_release'];assembly=sm['assembly']
    if assembly!='GRCh38.p14':raise ValueError('New assembly needs reviewed locus remapping')
    # Verify immutable downloaded bytes before using any source.
    import hashlib
    for src in sm['sources']:
        h=hashlib.sha256()
        with (cache/src['file']).open('rb') as f:
            while b:=f.read(1024*1024):h.update(b)
        assert h.hexdigest()==src['sha256'],src['file']
    rows=[]
    with gzip.open(cache/f'gencode.v{release}.annotation.gtf.gz','rt') as f:
        for line in f:
            if line.startswith('#') or 'gene_id "ENSG00000134323.' not in line:continue
            fields=line.rstrip().split('\t');a=attributes(fields[8])
            if one(a,'gene_id').split('.')[0]=='ENSG00000134323':rows.append((fields,a))
    txrows={one(a,'transcript_id'):(f,a) for f,a in rows if f[2]=='transcript'}
    assert txrows,'No MYCN transcripts in official annotation'
    txseq=fasta(cache/f'gencode.v{release}.transcripts.fa.gz',set(txrows))
    proteinids={one(a,'protein_id') for f,a in rows if one(a,'protein_id')};proteins=fasta(cache/f'gencode.v{release}.pc_translations.fa.gz',proteinids)
    genome=SeqIO.read(cache/'refseq-genomic.gb','genbank');assert genome.id=='NC_000002.12' and len(genome)==9201
    loc={'chrom':'chr2','start':15938000,'end':15947200,'seq':str(genome.seq).upper(),'source':'NCBI NC_000002.12 GRCh38.p14 genomic slice','coordinate_convention':'1-based inclusive'}
    metadata={}
    for f in sorted(cache.glob(f'gencode.v{release}.metadata.*.gz')):
        selected=[]
        with gzip.open(f,'rt') as handle:
            for line in handle:
                key=line.split('\t')[0].split('.')[0]
                if key in {x.split('.')[0] for x in txrows}:selected.append(line.rstrip().split('\t'))
        metadata[f.name]=selected
    uniprot=json.loads((cache/'uniprot-P04198.json').read_text());uniseq=uniprot['sequence']['value'];assert uniprot['entryType'].startswith('UniProtKB reviewed')
    def domains(aa):
        if aa!=uniseq:return []
        return [{'type':f['type'],'description':f.get('description',''),'start':f['location']['start']['value'],'end':f['location']['end']['value'],'source':'UniProt P04198','source_version':uniprot['entryAudit'],'evidence_class':'DATABASE_ANNOTATED','evidences':f.get('evidences',[]),'input_sequence_sha256':sha(aa)} for f in uniprot['features'] if f['type'] in ['Domain','Region','Motif','Coiled coil','Compositional bias'] and f['location']['start'].get('modifier')=='EXACT' and f['location']['end'].get('modifier')=='EXACT']
    txs={}
    for accession,(f,a) in sorted(txrows.items(),key=lambda item:one(item[1][1],'transcript_name')):
        related=[(x,z) for x,z in rows if one(z,'transcript_id')==accession];strand=f[6]
        ex=sorted([(int(x[3]),int(x[4]),one(z,'exon_id'),one(z,'exon_number')) for x,z in related if x[2]=='exon'],reverse=strand=='-')
        exons=[[x[0],x[1]] for x in ex];pos=positions(exons,strand);seq=splice(loc['seq'],loc['start'],exons,strand)
        assert seq==txseq[accession],f'Transcript FASTA/genome mismatch {accession}'
        cds=sorted([[int(x[3]),int(x[4])] for x,z in related if x[2]=='CDS'],reverse=strand=='-')
        stop=sorted([[int(x[3]),int(x[4])] for x,z in related if x[2]=='stop_codon'],reverse=strand=='-')
        utr=sorted([[int(x[3]),int(x[4])] for x,z in related if x[2]=='UTR'],reverse=strand=='-')
        protein=next((one(z,'protein_id') for x,z in related if one(z,'protein_id')),'')
        aa=splice(loc['seq'],loc['start'],cds,strand) if cds else ''
        aa=translate(aa) if aa else ''
        if protein:assert aa==proteins[protein],f'Protein FASTA/CDS mismatch {accession}'
        name=one(a,'transcript_name');key=name[5:] if name.startswith('MYCN-') else accession
        tags=a.get('tag',[]);pindex={g:i for i,g in enumerate(pos)};codingpos=[pindex[g] for lo,hi in cds for g in range(lo,hi+1)]
        codingstart=min(codingpos) if codingpos else None;codingend=max(codingpos)+1 if codingpos else None
        stoppos=[pindex[g] for lo,hi in stop for g in range(lo,hi+1)]
        if cds:
            assert seq[codingstart:codingstart+3]=='ATG'
            assert stoppos and len(stoppos)==3 and min(stoppos)==codingend
            assert seq[codingend:codingend+3] in ['TAA','TAG','TGA']
        txs[key]={'id':accession,'transcript_id':accession,'display_name':name,'gene_id':one(a,'gene_id'),'chrom':f[0],'strand':strand,'exons':exons,'exon_ids':[e[2] for e in ex],'exon_order':[e[3] for e in ex],'cds':cds,'stop_codon_blocks':stop,'utr':utr,'utr5':blocks_for_slice(pos,0,codingstart) if cds else [],'utr3':blocks_for_slice(pos,codingend+3,len(pos)) if cds else [],'junctions':[[exons[i][1],exons[i+1][0]] for i in range(len(exons)-1)],'mlen':len(seq),'sequence':seq,'sequence_sha256':sha(seq),'protein_accession':protein,'protein_sequence':aa,'protein_length':len(aa),'protein_sha256':sha(aa) if aa else None,'cds_transcript_start':codingstart,'cds_transcript_end':codingend,'source':f'GENCODE v{release}','src':f'GENCODE v{release}','gencode_release':f'v{release}','biotype':one(a,'transcript_type'),'description':name+': '+one(a,'transcript_type'),'reference_release':f'GENCODE v{release}','tags':tags,'annotation_level':one(a,'level'),'transcript_support_level':one(a,'transcript_support_level','not supplied'),'appris':[t for t in tags if t.startswith('appris')],'mane':[t for t in tags if t.startswith('MANE')],'basic':'basic' in tags,'primary':'GENCODE_Primary' in tags,'annotation_attributes':a,'evidence_class':'DATABASE_ANNOTATED','domains':domains(aa),'metadata':{k:[r for r in v if r[0].split('.')[0]==accession.split('.')[0]] for k,v in metadata.items()}}
    # Current curated NCBI transcript; place it using its current or explicitly replaced genomic annotation, requiring exact exon/sequence identity.
    ref=SeqIO.read(cache/'refseq-transcript.gb','genbank');assert ref.id.startswith('NM_'),'RefSeq record is not curated mRNA'
    assert 'REVIEWED REFSEQ' in ref.annotations.get('comment',''),'Unexpected RefSeq review status'
    previous=re.findall(r'XM_\d+\.\d+',ref.annotations.get('comment',''))
    placements=[x for x in genome.features if x.type=='mRNA' and x.qualifiers.get('transcript_id',[''])[0] in [ref.id]+previous and str(x.extract(genome.seq))==str(ref.seq)]
    assert len(placements)==1,'Ambiguous or absent exact RefSeq genomic placement'
    exons=[[15938000+int(p.start),15938000+int(p.end)-1] for p in placements[0].location.parts]
    pos=positions(exons,'+');seq=str(ref.seq);assert splice(loc['seq'],loc['start'],exons)==seq
    cdsfeat=next(x for x in ref.features if x.type=='CDS');cs=int(cdsfeat.location.start);ce=int(cdsfeat.location.end);aa=cdsfeat.qualifiers['translation'][0];protein=cdsfeat.qualifiers['protein_id'][0]
    assert translate(seq[cs:ce])==aa+'*';assert seq[cs:cs+3]=='ATG'
    if (cache/'refseq-protein.gb').exists():assert str(SeqIO.read(cache/'refseq-protein.gb','genbank').seq)==aa
    txs['XM']={'id':ref.id,'transcript_id':ref.id,'display_name':'MYCN-XM','gene_id':'GeneID:4613','chrom':'chr2','strand':'+','exons':exons,'exon_ids':[],'exon_order':['1','2'],'cds':blocks_for_slice(pos,cs,ce-3),'stop_codon_blocks':blocks_for_slice(pos,ce-3,ce),'utr':blocks_for_slice(pos,0,cs)+blocks_for_slice(pos,ce,len(seq)),'utr5':blocks_for_slice(pos,0,cs),'utr3':blocks_for_slice(pos,ce,len(seq)),'junctions':[[exons[i][1],exons[i+1][0]] for i in range(len(exons)-1)],'mlen':len(seq),'sequence':seq,'sequence_sha256':sha(seq),'protein_accession':protein,'protein_sequence':aa,'protein_length':len(aa),'protein_sha256':sha(aa),'cds_transcript_start':cs,'cds_transcript_end':ce-3,'source':'RefSeq reviewed transcript variant 4','src':'RefSeq','gencode_release':'not a GENCODE transcript','refseq_id':ref.id,'refseq_obsolete_id':previous[0] if previous else None,'refseq_history':previous,'placement_annotation_accession':placements[0].qualifiers['transcript_id'][0],'placement_note':'Current reviewed mRNA is byte-identical to the explicitly replaced genomic mRNA annotation; every exon boundary and complete spliced sequence were verified.','biotype':'protein_coding','tags':['REVIEWED REFSEQ'],'annotation_attributes':{},'transcript_support_level':'NCBI reviewed','appris':[],'mane':[],'basic':False,'primary':False,'evidence_class':'DATABASE_ANNOTATED','domains':domains(aa),'description':ref.description,'reference_release':'RefSeq '+ref.id+'; record '+ref.annotations.get('date','not supplied'),'metadata':{}}
    # Recompute ORFs independently on every current transcript. Merge only identical genomic blocks AND nucleotide strings.
    merged={}
    for key,tx in txs.items():
        pos=positions(tx['exons'],tx['strand']);seq=tx['sequence'];cstart=tx['cds_transcript_start'];cend=tx['cds_transcript_end']
        for start,end in find_orfs(seq,8):
            nt=seq[start:end];aa=translate(nt);assert aa.endswith('*') and '*' not in aa[:-1]
            blocks=blocks_for_slice(pos,start,end);identity=json.dumps([blocks,nt,tx['strand']]);o=merged.setdefault(identity,{'blocks':blocks,'nt':nt,'aa':aa[:-1],'aa_len':len(aa)-1,'carriers':[],'transcript_records':{},'annotated_on':[]})
            annotated=cstart==start and cend is not None and end==cend+3 and tx['protein_sequence']==aa[:-1]
            if annotated:relation='ANNOTATED_CDS';o['annotated_on'].append(key)
            elif cstart is None:relation='ALTERNATIVE_ORF'
            elif end<=cstart:relation='UPSTREAM_ORF'
            elif start>=cend+3:relation='DOWNSTREAM_ORF'
            else:relation='OVERLAPPING_ORF'
            overlap=sum(1 for g in pos[start:end] if any(a<=g<=b for a,b in tx['cds']))
            o['carriers'].append(key);o['transcript_records'][key]={'start':start+1,'end':end,'coordinate_convention':'1-based inclusive transcript','frame':start%3,'relationship':relation,'cds_overlap_nt':overlap,'start_codon':nt[:3],'stop_codon':nt[-3:],'prediction_method':'AUG_FIRST_STOP_MIN8_STANDARD_CODE_V1','evidence_class':'DATABASE_ANNOTATED' if annotated else 'COMPUTATIONAL_PREDICTION'}
    # IDs are local to this scientific build; do not imply continuity with legacy ORF numbers.
    values=sorted(merged.values(),key=lambda o:(not bool(o['annotated_on']),-o['aa_len'],o['blocks'],o['nt']))
    orfs={}
    for n,o in enumerate(values,1):
        key=str(n);host=o['carriers'][0];start=o['transcript_records'][host]['start'];end=o['transcript_records'][host]['end'];an=bool(o['annotated_on']);up=all(x['relationship']=='UPSTREAM_ORF' for x in o['transcript_records'].values());cls='canonical-CDS' if an else 'novel-uORF' if up else 'novel-other';pa=ProteinAnalysis(o['aa'])
        props={'mw':round(pa.molecular_weight(),3),'pI':round(pa.isoelectric_point(),4),'gravy':round(pa.gravy(),5),'instab':round(pa.instability_index(),4),'netq':round(pa.charge_at_pH(7),5),'composition':pa.count_amino_acids()}
        accessions=sorted({txs[t]['protein_accession'] for t in o['annotated_on']});ds=domains(o['aa']) if an else []
        o.update({'id':n,'build_orf_id':f'gencode-{release}-orf-{n:04d}','cls':cls,'annotation_class':'ANNOTATED_CDS' if an else 'PREDICTED_ORF','evidence_class':'DATABASE_ANNOTATED' if an else 'COMPUTATIONAL_PREDICTION','acc':'; '.join(accessions),'protein_accessions':accessions,'host':host,'m1':start,'m2':end,'span':[o['blocks'][0][0],o['blocks'][-1][1]],'location':','.join(sorted({v['relationship'] for v in o['transcript_records'].values()})),'source':f'GENCODE v{release} and current reviewed RefSeq sequence' if an else 'Sequence-based prediction, no expression or functional inference','color':'#4A6FA5' if an else '#5B7D71' if up else '#657D96','start_codon':'ATG','stop_codon':o['nt'][-3:],'protein_sequence_sha256':sha(o['aa']),'nt_sha256':sha(o['nt']),'property_provenance':{'software':'Biopython','version':Bio.__version__,'algorithm':'Bio.SeqUtils.ProtParam.ProteinAnalysis','parameters':{'genetic_code':1,'net_charge_pH':7,'average_molecular_mass':True},'input_sequence_sha256':sha(o['aa']),'calculation_date':review,'evidence_class':'COMPUTATIONAL_PREDICTION'},'domains':ds,'conserved':'NOT_ASSESSED',**props})
        orfs[key]=o
    buildid=f'gencode-{release}_refseq-{ref.id}_mycn-1'
    data={'schema_version':'2.0.0','application_version':'2.0.0','reference_build':buildid,'status':'CURRENT','review_date':review,'assembly':assembly,'ensembl_release':sm['ensembl_release'],'gencode_release':release,'refseq_accession':ref.id,'transcript_order':list(txs),'transcripts':txs,'orfs':orfs,'regions':{},'palette':{'nmyc':'#4A6FA5','dmycn':'#D98A3D','mycnot':'#C0504D','musep':'#B0568C','ri_line':'#2F7A57','ri_bg':'#C3E0D4','pred_ext':'#8FA9C9','utr':'#C7CFCC','term':'#8A9791','gold':'#C9A227'},'prediction_policy':{'method':'AUG_FIRST_STOP_MIN8_STANDARD_CODE_V1','start_codons':['ATG'],'minimum_amino_acids':8,'genetic_code':1,'include_nested_starts':True,'include_no_stop':False,'non_AUG_calls':'not enumerated','experimental_claims':'none; annotation and prediction do not establish expression or function'},'official_source_manifest':sm,'official_metadata':metadata,'claim_policy':{'structures':'DATABASE_ANNOTATED','ORF_enumeration':'COMPUTATIONAL_PREDICTION','physicochemical_properties':'COMPUTATIONAL_PREDICTION','domains':'DATABASE_ANNOTATED only for exact sequence match to UniProt P04198','laboratory_aliases':'OWNER_SUPPLIED_NOMENCLATURE','legacy_MS_RiboSeq_conservation':'removed from current factual claims; not independently reproduced'}}
    return data,loc

def prepare_layers(data):
    aliases=json.loads(re.search(r'var aliases = Object.freeze\((\{.*?\})\);',(ROOT/'lab-aliases.js').read_text(),re.S).group(1))
    previous=json.loads((ROOT/'data/releases/legacy-2026-07/reference.json').read_text())['transcripts']
    for k,t in data['transcripts'].items():t['lab_alias']=aliases.get(t['display_name'],'');t['history']=previous.get(k)
    from build_product_layers import augment
    return augment(data)

def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',required=True);p.add_argument('--review-date',required=True);p.add_argument('--application-commit');args=p.parse_args()
    data,loc=build(Path(args.cache),args.review_date)
    import subprocess,platform,importlib.metadata
    data['pipeline_versions']={'python':platform.python_version(),'biopython':Bio.__version__,'jsonschema':importlib.metadata.version('jsonschema')}
    data['application_commit']=args.application_commit or subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    prepare_layers(data)
    from validate_reference import validate
    validate(data,loc)
    write(ROOT/'data/current/reference.json',data);write(ROOT/f'data/releases/gencode_v{data["gencode_release"]}/reference.json',data);write(ROOT/'mycn_locus_seq.json',loc)
    manifest={k:data[k] for k in ['schema_version','application_version','reference_build','status','review_date','assembly','ensembl_release','gencode_release','refseq_accession','prediction_policy','application_commit','pipeline_versions','structure_reference_build','product_registry_version']};manifest['reference_sha256']=sha((ROOT/'data/current/reference.json').read_text());manifest['sources']=data['official_source_manifest'];write(ROOT/'data/current/manifest.json',manifest);write(ROOT/f'data/manifests/{data["reference_build"]}.json',manifest)
    # Retain exact official locus subsets with complete-file source hashes in the manifest.
    cache=Path(args.cache);release=data['gencode_release'];sub=[]
    with gzip.open(cache/f'gencode.v{release}.annotation.gtf.gz','rt') as f:
        for line in f:
            if 'ENSG00000134323.' in line:sub.append(line)
    dest=ROOT/f'data/releases/gencode_v{release}';(dest/'MYCN.annotation.gtf').write_text(''.join(sub))
    for filename,key in [('MYCN.transcripts.fa','sequence'),('MYCN.proteins.fa','protein_sequence')]:
        (dest/filename).write_text(''.join('>'+tx['transcript_id']+(' protein='+tx['protein_accession'] if key=='protein_sequence' else '')+'\n'+tx[key]+'\n' for tx in data['transcripts'].values() if tx[key]))
    from build_public_assets import generate
    generate(data)
    print('Rebuilt',len(data['transcripts']),'transcripts and',len(data['orfs']),'unique ORF predictions/annotations; build',data['reference_build'])
if __name__=='__main__':main()
