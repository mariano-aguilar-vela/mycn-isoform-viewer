"""Compare structural releases without recalculating historical study estimates."""
import argparse,json
from pathlib import Path
from reference_lib import sha
ROOT=Path(__file__).resolve().parents[1]
def compare(old,new):
    rows=[]
    for key in sorted(set(old['transcripts'])|set(new['transcripts'])):
        a=old['transcripts'].get(key);b=new['transcripts'].get(key);changes=[]
        if not a or not b:changes=['ANNOTATION_CHANGE']
        else:
            if a['exons']!=b['exons']:changes+=['BOUNDARY_CHANGE']
            if len(a['exons'])<len(b['exons']):changes+=['EXON_ADDED']
            if len(a['exons'])>len(b['exons']):changes+=['EXON_REMOVED']
            if a.get('junctions')!=b.get('junctions'):changes+=['JUNCTION_CHANGE']
            if a.get('sequence_sha256')!=b.get('sequence_sha256'):changes+=['SEQUENCE_CHANGE']
            if a['cds']!=b['cds']:changes+=['CDS_CHANGE']
            if a.get('protein_sequence')!=b.get('protein_sequence'):changes+=['PROTEIN_CHANGE']
            if a.get('transcript_id')!=b.get('transcript_id') and not changes:changes+=['VERSION_ONLY']
            if not changes:changes=['UNCHANGED']
        rows.append({'canonical_name':b['display_name'] if b else a['display_name'],'classifications':changes,'previous':a,'current':b,'previous_orf_count':sum(key in o.get('carriers',[]) for o in old.get('orfs',{}).values()),'current_orf_count':sum(key in o['carriers'] for o in new['orfs'].values()),'domain_change':'Current domains require exact sequence identity to reviewed UniProt; prior claims not carried forward','evidence_change':['EVIDENCE_CHANGE','REMOVED_UNSUPPORTED_CLAIM'],'note':'ORF counts use different enumerations and are not directly a gain/loss of expressed proteins.'})
    return {'previous_status':'HISTORICAL','current_reference_build':new['reference_build'],'transcripts':rows,'claim_changes':[{'class':'REMOVED_UNSUPPORTED_CLAIM','scope':'Legacy record-level MS, ribosome-profiling, conservation and functional assertions','reason':'Underlying analyses not recovered/reproduced for current ORFs; no implicit evidence transfer'},{'class':'NEW_SUPPORTED_CLAIM','scope':'Current GENCODE/RefSeq structures, exact sequence-matched UniProt domains and explicit computational properties','reason':'Direct primary-database retrieval and sequence validation'}]}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',default='data/releases/legacy-2026-07/reference.json');p.add_argument('--new',default='data/current/reference.json');p.add_argument('--output',default='data/current/release-diff.json');a=p.parse_args();Path(a.output).write_text(json.dumps(compare(json.loads(Path(a.old).read_text()),json.loads(Path(a.new).read_text())),indent=2)+'\n');print('Release diff generated')
