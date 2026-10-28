"""Create an allowlisted static artifact; never publish archives or audit material."""
from pathlib import Path
import argparse, hashlib, json, shutil, re
parser=argparse.ArgumentParser();parser.add_argument('--output',default='dist');args=parser.parse_args()
root=Path(__file__).resolve().parents[1];dest=Path(args.output).resolve()
if dest==root or root in dest.parents and dest.name!='dist':raise SystemExit('Choose dist or an external output directory')
files=['index.html','app.js','styles.css','lab-aliases.js','favicon.svg','404.html','mycn_orfs.meta.json','mycn_locus_seq.json','chr2_cytoband.json','mycn_transcripts.gff3','mycn_orfs.bed','vendor/d3.v7.9.0.min.js','vendor/D3-LICENSE']
from downloads import outputs
for filename,content in outputs(json.loads((root/'mycn_orfs.meta.json').read_text())).items():
 if (root/filename).read_text()!=content:raise SystemExit('Stale coordinate download: '+filename)
patterns=[r'/(?:home|mnt|scratch|lustre)/[\w/.-]+',r'gh[pousr]_[A-Za-z0-9]{20,}',r'github_pat_[A-Za-z0-9_]+',r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY',r'(?i)(?:patient_id|sample_id)\s*["\x27]?\s*[:=]',r'https?://(?:localhost|127\.0\.0\.1)']
for f in files:
 text=(root/f).read_text()
 for pat in patterns:
  if re.search(pat,text):raise SystemExit('Disclosure check failed: '+f+' / '+pat)
if dest.exists():
 extras=[str(p.relative_to(dest)) for p in dest.rglob('*') if p.is_file() and str(p.relative_to(dest)) not in files+['.nojekyll','asset-manifest.json']]
 if extras:raise SystemExit('Unexpected existing output files; refusing mixed artifact: '+', '.join(extras))
dest.mkdir(parents=True,exist_ok=True)
manifest={}
for f in files:
 p=dest/f;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(root/f,p);manifest[f]=hashlib.sha256(p.read_bytes()).hexdigest()
(dest/'.nojekyll').write_text('');(dest/'asset-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Built',len(files),'allowlisted assets; disclosure checks passed')
