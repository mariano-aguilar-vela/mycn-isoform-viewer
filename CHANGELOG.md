# Changes — 2026-10-28

- Added the eight owner-supplied laboratory aliases through one shared mapping, exact alias search, accessible labels and figure metadata.
- Added a keyboard-accessible data dialog and transcript-structure table, clearer frozen-reference provenance, research-use information, visible focus, and keyboard exon/shared-ORF controls.
- Fixed transcript filtering being overridden by CSS, toolbar overflow at intermediate widths, narrow shared-ORF labels, and mobile drawer focus handling.
- Regenerated coordinate downloads from the canonical metadata: corrected spliced CDS phases in GFF3 and included all 56 ORFs in BED.
- Reconciled the terminal stop codons in eight nucleotide strings with their existing coordinate blocks and the retained locus sequence; amino-acid strings are unchanged.
- Clarified that the 593-aa extension is predicted, that missing CDS annotation does not establish noncoding function, and that the below-30-aa instability-index threshold is a retained reporting policy rather than a mathematical limit.
- Bundled D3 7.9.0 locally with its license, added reproducible static packaging, data/browser regression checks, metadata and a 404 page.

The frozen GENCODE v49/v50 plus RefSeq models remain pinned. A uniform v50 update needs separate model and sequence reconciliation for MYCN-202, -204 and -205. No new biological validation is claimed.
