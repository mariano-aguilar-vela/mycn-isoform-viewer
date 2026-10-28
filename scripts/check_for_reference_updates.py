"""DETECT and report; never change data, push commits or deploy."""
import argparse,json,re,urllib.request
from pathlib import Path
from datetime import datetime,timezone
from Bio import SeqIO
from io import StringIO
from reference_lib import sha

def check(manifest,page,record):
    release=int(re.search(r'Release\s+(\d+)\s*\(',page).group(1));ref=SeqIO.read(StringIO(record),'genbank');change=release!=manifest['gencode_release'] or ref.id!=manifest['refseq_accession']
    return {'checked_utc':datetime.now(timezone.utc).isoformat(),'status':'UPDATE_DETECTED' if change else 'NO_UPDATE','deployed_gencode':manifest['gencode_release'],'official_gencode':release,'deployed_refseq':manifest['refseq_accession'],'official_refseq':ref.id,'refseq_sequence_sha256':sha(str(ref.seq)),'record_reviewed':'REVIEWED REFSEQ' in ref.annotations.get('comment',''),'next_step':'FETCH_DIFF_REVIEW_TEST_DEPLOY' if change else 'REVIEW_EVIDENCE_AT_NEXT_SCHEDULE','automatic_publication':False}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',default='data/current/manifest.json');p.add_argument('--output',default='reference-update-report.json');p.add_argument('--fixtures');p.add_argument('--manifest-url');p.add_argument('--report-only',action='store_true');p.add_argument('--candidate-output',default='reference-update-candidate');p.add_argument('--candidate-cache',default='.reference-update-cache');a=p.parse_args()
    if a.fixtures:
        d=Path(a.fixtures);page=(d/'gencode-current.html').read_text();record=(d/'refseq-transcript.gb').read_text()
    else:
        page=urllib.request.urlopen('https://www.gencodegenes.org/human/',timeout=40).read().decode();record=urllib.request.urlopen('https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NM_001446166&db=nuccore&report=genbank&retmode=text',timeout=40).read().decode()
    manifest=json.load(urllib.request.urlopen(a.manifest_url,timeout=40)) if a.manifest_url else json.loads(Path(a.manifest).read_text())
    report=check(manifest,page,record);report['manifest_source']=a.manifest_url or 'local reviewed manifest';
    if report['status']=='UPDATE_DETECTED' and not a.report_only and not a.fixtures:
        import subprocess,sys
        subprocess.run([sys.executable,str(Path(__file__).with_name('fetch_reference_sources.py')),'--cache',a.candidate_cache],check=True)
        from build_mycn_reference import build
        from validate_reference import validate
        from compare_reference_releases import compare
        candidate,loc=build(Path(a.candidate_cache),datetime.now(timezone.utc).date().isoformat());validate(candidate,loc)
        dest=Path(a.candidate_output);dest.mkdir(parents=True,exist_ok=True)
        (dest/'candidate-reference.json').write_text(json.dumps(candidate,indent=2)+'\n')
        current=json.load(urllib.request.urlopen(a.manifest_url.rsplit('/',1)[0]+'/reference.json',timeout=60)) if a.manifest_url else json.loads(Path(a.manifest).with_name('reference.json').read_text())
        (dest/'structural-diff.json').write_text(json.dumps(compare(current,candidate),indent=2)+'\n');report['candidate_status']='GENERATED_FOR_REVIEW_ONLY'
    Path(a.output).write_text(json.dumps(report,indent=2)+'\n');print(report['status'])
