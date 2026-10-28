# MYCN living transcript and ORF viewer

Application **2.0.0** displays the complete **GENCODE v50 / GRCh38.p14 / Ensembl 116** MYCN transcript set and the reviewed **RefSeq NM_001446166.1** addition, shown as **MYCN-XM (TV4)**. Data were scientifically reviewed on **2026-10-28**. The current build is `gencode-50_refseq-NM_001446166.1_mycn-1`.

[Public site](https://mariano-aguilar-vela.github.io/mycn-isoform-viewer/) · [Reference-update policy](REFERENCE_UPDATE_POLICY.md) · [Current source manifest](data/current/manifest.json) · [Scientific release changes](data/current/release-diff.json)

The structures are regenerated from official GENCODE annotation, transcript FASTA and protein FASTA, and NCBI transcript/protein/genomic records. Every transcript sequence matches concatenated genomic exons; annotated CDS translations match the official proteins. MYCN-202, MYCN-204 and MYCN-205 use the current `.2` accessions with their corresponding rebuilt exon boundaries, not relabelled v49 models.

TV1–TV6, TV2+, and NC are owner-supplied laboratory aliases, maintained centrally in `lab-aliases.js`. Current database names and accessions remain authoritative. MYCN-XM is a study display label for the reviewed RefSeq transcript variant 4, not a GENCODE name. Its replaced XM accessions and exact genomic placement verification are retained in provenance.

## ORFs and evidence

The current sequences yield **42 distinct ORFs** under a declared computational policy: every AUG start to its first in-frame standard-code stop, at least 8 amino acids, including nested starts. Non-AUG and no-stop ORFs are not enumerated. These build-specific IDs do not correspond to the old 56-record ORF numbering. Counts do not represent expressed proteins.

Annotated CDSs carry `DATABASE_ANNOTATED`; other ORFs and Biopython physicochemical values carry `COMPUTATIONAL_PREDICTION`. Domain/region coordinates are transferred only when the protein sequence exactly matches reviewed UniProt P04198, with sequence hashes, source versions and database evidence retained. No domain assigned does not mean no domain exists. Protein instability index and other properties are predictions, not experimental measurements.

The old record-level mass-spectrometry, ribosome-profiling, conservation and functional assertions are removed from the current factual view because their underlying analyses have not been reproduced against these new ORFs. Historical structures and the previous public application commit remain reproducible; historical abundance estimates are neither included nor recalculated.

## Build and test

Use Node 20+, Python 3.13, and the pinned reference requirements:

```
pip install -r requirements-reference.txt
npm ci
npm test
npm run build
```

The static build uses generated assets in `data/current/`, checks them against the JSON schema and scientific sequence rules, rejects stale coordinate exports and scans the allowlisted artifact for defined disclosure patterns. Serve `dist/` at `/mycn-isoform-viewer/`. D3 7.9.0 is vendored locally with its upstream license. There are no runtime CDN, analytics, cookie or persistent-storage dependencies.

`npm run test:browser` uses Playwright and axe with `MYCN_TEST_URL` and `MYCN_TEST_OUTPUT`. All Chromium, Firefox and WebKit engines must pass. A compatible Playwright container is suitable when host WebKit libraries are unavailable. The existing GitHub Pages deployment method is unchanged; the separate monthly workflow only detects and prepares reference updates for review.

See the reference-update policy for full-source download, checksum verification, regeneration and release-diff commands. Reference caches and private audit material belong outside the public repository. Never push private audit branches.

## Interaction, exports and accessibility

Select a transcript to expand its ORFs; the ORF lens shows annotations and predictions separately. Search matches canonical names, exact laboratory aliases, versioned/stable accessions, protein accessions and transcript descriptions. TV2 and TV2+ stay distinct. Use the locus field, zoom buttons, Reset, pointer navigation, keyboard transcript/exon controls and structured data table.

Reference and evidence shows current provenance, build versions, research-use limitations and sequence/coordinate downloads. History links the previous structural snapshot and scientific changes. Data JSON contains complete transcript/ORF structures, sequences, translations, laboratory aliases, sources, evidence statuses and build versions. Standard FASTA/GFF3/BED/TSV downloads provide format-specific records; the JSON and manifest are their full provenance sidecars. SVG and PNG exports show the current reference and build; SVG also carries detailed machine-readable metadata.

The viewer is for research and education, not clinical diagnosis. Annotation and laboratory aliases do not establish sample-specific expression or function. Automated accessibility scores and keyboard/browser checks do not constitute formal WCAG certification.
