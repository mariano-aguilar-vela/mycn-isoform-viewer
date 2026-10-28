# Living reference policy

The default is CURRENT: the complete MYCN set extracted from the current official comprehensive GENCODE chromosome annotation, plus the reviewed NCBI MYCN transcript variant displayed as MYCN-XM. Laboratory aliases are owner-supplied nomenclature. Additional official models receive their official name and no invented laboratory alias.

The monthly `Review reference updates` workflow and `npm run reference:check` follow DETECT → DIFF → REVIEW → TEST → DEPLOY. They retrieve the current official release and RefSeq record and report changes. When a change is detected, the checker builds a candidate in an isolated directory and produces a structural diff. It never modifies current assets, pushes changes or publishes a scientific release. Read-only permissions and pinned GitHub Actions are used. Monthly execution starts only once the workflow reaches the repository default branch.

Rebuild with Python 3.13 and `pip install -r requirements-reference.txt`:

```
python scripts/fetch_reference_sources.py --cache /path/to/external-reference-cache
python scripts/build_mycn_reference.py --cache /path/to/external-reference-cache --review-date YYYY-MM-DD
python scripts/compare_reference_releases.py
python scripts/downloads.py --write
npm test
npm run build
npm run test:browser
```

The cache stores complete official source files and publisher MD5/SHA-256 provenance. It must remain outside the public source tree. The generated release manifest contains public source identities, URLs and hashes, never machine-specific cache paths. Genomic intervals are 1-based inclusive; Python slices are 0-based half-open; exported transcript coordinates are 1-based inclusive; BED is zero-based half-open. CDS intervals exclude stop codons; each complete ORF nucleotide sequence includes its terminal stop.

Every AUG start reaching its first in-frame standard-code stop is enumerated, including nested starts, at a minimum length of eight amino acids. This explicit prediction policy excludes non-AUG and no-stop ORFs. IDs are local to a scientific build and are not inherited from the former viewer. Annotated CDSs must match official protein FASTA. Computational ORFs do not establish expression, translation or function. Protein properties are Biopython predictions. UniProt feature coordinates are assigned only to an exact matching amino-acid sequence and retain source evidence and sequence hashes.

Application releases use semantic versions; scientific releases use explicit GENCODE/RefSeq build identifiers. CURRENT becomes SUPERSEDED after a reviewed successor. Historical structural snapshots and Git commits remain HISTORICAL or ARCHIVED. Historical study estimates are not recalculated by a structural update; the viewer contains no Kelly abundance payload. A release with changed assembly, ambiguous placement, incomplete official source identity or failed sequence validation must stop for review.

Publication uses the existing branch-based GitHub Pages mechanism. After all scientific, accessibility, Chromium/Firefox/WebKit, export, disclosure and clean-build checks pass, verify authenticated Pages settings, record a rollback identity, publish the reviewed commit through the established mechanism, and verify the public artifact and site. Never push the separate private audit branch.
