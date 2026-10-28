"""Fetch official inputs to a private/external cache, retaining full-file SHA-256 provenance."""
import argparse,concurrent.futures,hashlib,json,re,urllib.parse,urllib.request,html
from pathlib import Path
from datetime import datetime,timezone
ALLOWED={'www.gencodegenes.org','ftp.ebi.ac.uk','www.ncbi.nlm.nih.gov','eutils.ncbi.nlm.nih.gov','api.ncbi.nlm.nih.gov','rest.uniprot.org','www.ebi.ac.uk','ftp.ncbi.nlm.nih.gov'}
def fetch(url,dest,refresh=False):
    u=urllib.parse.urlparse(url)
    if u.scheme!='https' or u.hostname not in ALLOWED:raise ValueError('Unapproved reference URL')
    dest=Path(dest);dest.parent.mkdir(parents=True,exist_ok=True)
    if refresh or not dest.exists():
        tmp=dest.with_suffix(dest.suffix+'.part')
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'MYCN-reference-builder/2.0'}),timeout=120) as r,tmp.open('wb') as out:
            if urllib.parse.urlparse(r.geturl()).hostname not in ALLOWED:raise ValueError('Unapproved redirect')
            while b:=r.read(1024*1024):out.write(b)
        tmp.replace(dest)
    h=hashlib.sha256()
    with dest.open('rb') as f:
        while b:=f.read(1024*1024):h.update(b)
    return {'file':dest.name,'url':url,'sha256':h.hexdigest(),'bytes':dest.stat().st_size,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'source_identity':u.hostname}
def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',required=True);args=p.parse_args();d=Path(args.cache);d.mkdir(parents=True,exist_ok=True)
    current=fetch('https://www.gencodegenes.org/human/',d/'gencode-current.html',refresh=True);page=(d/'gencode-current.html').read_text();release=int(re.search(r'Release\s+(\d+)\s*\(',page).group(1));assembly=re.search(r'Release\s+\d+\s*\(([^)]+)\)',page).group(1)
    base=f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{release}/'
    names=[f'gencode.v{release}.annotation.gtf.gz',f'gencode.v{release}.transcripts.fa.gz',f'gencode.v{release}.pc_translations.fa.gz']+[f'gencode.v{release}.metadata.{name}.gz' for name in ['RefSeq','SwissProt','Transcript_source','Annotation_remark','Transcript_supporting_feature','Exon_supporting_feature']]
    urls={n:base+n for n in names}
    urls.update({'MD5SUMS':base+'MD5SUMS','refseq-protein.gb':'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NP_001433095.1&db=protein&report=genbank&retmode=text','gencode-releases.html':'https://www.gencodegenes.org/human/releases.html','refseq-transcript.gb':'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NM_001446166&db=nuccore&report=genbank&retmode=text','refseq-genomic.gb':'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=NC_000002.12&db=nuccore&report=genbank&retmode=text&from=15938000&to=15947200','refseq-gene.json':'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/4613','uniprot-P04198.json':'https://rest.uniprot.org/uniprotkb/P04198.json'})
    sources=[current];errors=[]
    def run(item):
        name,url=item
        try:
            r=fetch(url,d/name,refresh=not name.endswith(".gz"));print('Fetched',name,r['bytes'],flush=True);return r
        except Exception as e:print('FAILED',name,str(e),flush=True);return {'file':name,'url':url,'error':str(e)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for r in pool.map(run,urls.items()):
            (errors if 'error' in r else sources).append(r)
    ensembl=None
    for row in re.findall(r'<tr[^>]*>.*?</tr>',(d/'gencode-releases.html').read_text(),re.S):
        cells=[html.unescape(re.sub('<[^>]*>','',c)).strip() for c in re.findall(r'<td[^>]*>(.*?)</td>',row,re.S)]
        if len(cells)>5 and cells[1]==str(release):ensembl=int(cells[5])
    assert ensembl,'Missing official Ensembl release'
    expected={line.split()[-1].lstrip('*./'):line.split()[0] for line in (d/'MD5SUMS').read_text().splitlines() if len(line.split())>=2}
    for src in sources:
        if src['file'] in expected:
            h=hashlib.md5()
            with (d/src['file']).open('rb') as f:
                while b:=f.read(1024*1024):h.update(b)
            assert h.hexdigest()==expected[src['file']],src['file']
            src['publisher_md5']=expected[src['file']]
    manifest={'gencode_release':release,'assembly':assembly,'ensembl_release':ensembl,'sources':sources,'unavailable_sources':errors}
    (d/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Source manifest complete; unavailable',len(errors),flush=True)
if __name__=='__main__':main()
