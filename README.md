# MYCN isoform & ORF viewer

An interactive viewer of the **MYCN** locus showing its **8 transcript isoforms** and the **48 encoded ORFs**,
with per-transcript ORF specificity (unique vs shared), base-level sequence, and the encoded protein / uORF
landscape. Built as a custom d3/SVG genome-browser-style viewer — transcript-first, with an ideogram, a
genomic ruler, a reference-sequence strip, zoom/pan, and figure export.

**Live site:** https://mariano-aguilar-vela.github.io/mycn-isoform-viewer/

## Data

- Reference: **GRCh38**; transcript annotation from **GENCODE v49** (transcript 207 from v50; XM from RefSeq).
- Transcript models and the ORF set were derived in the thesis project (MYCN transcript-isoform
  characterisation). The staged locus reference sequence is the GRCh38 slice `chr2:15,938,000–15,947,200`.
- Per-ORF metrics shown: amino-acid length, Kozak context, PhyloCSF, phyloP, MW / pI / GRAVY, net charge,
  instability index, conservation, the **per-resource translation-evidence block**, **MS / protein evidence**,
  accession, and the carrier transcripts.
- **Instability index is blank below 30 aa**, with its reason shown in place. The dipeptide-weighted index is
  undefined at that length, so for those ORFs the viewer renders *not computed* rather than a number. A blank
  is not a zero, and GRAVY (uncalibrated, but defined) is retained for every ORF.
- **Conservation caveat:** conservation is a *weak discriminator for short uORFs* — the known functional uORFs
  **MYCNOT** (ORF9) and **MUSEP** (ORF10) are both conservation-negative yet functional. The viewer therefore
  presents conservation as evidence, not a verdict.
- **MS / protein-evidence caveat:** the MS axis is **lysate mass spectrometry** (OpenProt), the lower-confidence
  MS assay. Detection is shown at three states — *detected* (≥2 unique peptides), *present-not-detected*
  (in the catalogue, 0 peptides), and *absent-from-MS-catalogue* (checked, not found) — kept distinct from
  *not-assessed*. **MYCNOT is MS-detected at 7 unique peptides**, the control proving the assay is live at this
  locus. No HLA-immunopeptidomics was tested, so *no immunopeptidomic support* is never an earned claim.
- **Ribo-seq is not one axis, and its absences are not all negatives.** Evidence is reported per resource, in
  five states that are never conflated: *DETECTED*, *SUB-THRESHOLD*, *PRESENT-NOT-DETECTED*, *ABSENT*, and
  *NOT-ASSESSED* (always with its reason). Only **nuORFdb v1.2** is both powered and in scope here, so only its
  absences are real negatives. **sORFs.org** is non-exhaustive and scoped to ≤100 aa; **Ribo-uORF** covers the
  uORF/uTIS class only; and **five consensus catalogues return zero records anywhere at the MYCN locus** — they
  are **UNPOWERED and NOT-ASSESSED, not negative**. A *no Ribo-seq support* claim across all of them would be
  false, and is not made.
- **Status of the ORFs:** apart from the annotated proteins (N-Myc, ΔMYCN), the known functional uORFs
  (MYCNOT, MUSEP), and **ORF 24 — OpenProt IP_083082, which now carries direct protein-level evidence (2 unique
  peptides by lysate MS, two-executor confirmed, frame-disjoint from N-Myc)** — the ORFs shown are **predictions**.
  ORF 24 is a bounded protein-level positive at the detection floor, not a discovery.

## How to use

- **Click a transcript** (left gutter or its row) to expand and see the ORFs it encodes, ranked by specificity:
  *unique to that transcript* first, then a collapsible *shared* group.
- **By-ORF lens** (top-left toggle): browse the 48 ORFs grouped by class; selecting an ORF **highlights it on
  every carrier transcript** and dims non-carriers, with full protein detail in the right rail.
- **Navigate:** *scroll to zoom* (cursor-centred), *drag to move*, *drag the ruler to zoom to a region*; the
  `chr2:start-end` box and the +/− / Reset buttons also work.
- **Zoom in** far enough to reveal the **base-level reference sequence** (A/C/G/T) and, for a selected ORF, its
  codons and amino acids (start green / stop red), with intron-spanning codons handled via the spliced sequence.
- **Export** the current view as **SVG or PNG** (ideogram + ruler + reference + stack).

## Built with

- [d3 v7](https://d3js.org/) (loaded from CDN). No backend — the site is fully static (HTML/CSS/JS + JSON data).

## Provenance

Part of **Mariano Aguilar Vela's** PhD (QUT / APCRC-Q), on MYCN transcript-isoform characterisation.
