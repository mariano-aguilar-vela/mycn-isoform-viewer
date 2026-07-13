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
- **MS / protein evidence — the axis has TWO limits, and both bound every cell.**
  The MS axis is **lysate mass spectrometry** (OpenProt 2.2), the lower-confidence MS assay.
  - **COVERAGE.** **MYCNOT is MS-detected at 7 unique peptides** — the control proving the **assay fires** at
    this locus. **It does not extend the search space.** OpenProt catalogues only **9 of the 763** ORFs; an ORF
    with no OpenProt accession was **never in the search space and was therefore never searched** →
    **NOT-ASSESSED, not a negative.** *Absent from the catalogue is not absent from the lysate.* The only real
    MS negatives are the **3** ORFs that **are** catalogued and returned 0 peptides (*present-not-detected*).
  - **ATTRIBUTION.** **No ORF at this locus is start-discriminated by lysate MS.** The locus is a dense nest of
    same-stop ORFs, and **every MS-detected ORF has a same-stop host that contains it in its entirety** (N-Myc
    18 hosts, MUSEP 5, ORF 24 two, ΔMYCN 1, **MYCNOT 1 — the control included**), so no peptide of it can ever
    exclude that host. All four detections are therefore **`DETECTED (GROUP-LEVEL)`**, and the viewer draws the
    evidence mark **once per stop group, never once per ORF** — a mark per ORF would let a reader count *N*
    detections where the evidence supports **one**. **The named proteins are identified by ANNOTATION, not by
    MS.** The honest ceiling of the axis is: ***"a protein from this stop group is translated."***
  - **Peptide uniqueness is bounded.** Proteome-uniqueness **excludes a canonical-fragment explanation**; it
    **cannot attribute the start**. The uniqueness search compared against four named proteins plus UniProt
    (UP000005640, 147,506 seqs) — **and five of the six Region-1 candidates are absent from every reference
    proteome**, so it ran against a database that **did not contain the competitors**. It was never capable of
    discriminating them.
  - No HLA-immunopeptidomics was tested, so *no immunopeptidomic support* is never an earned claim.
  - Accessions follow **Table B's namespace** (`XP_047300390.1`, `ENSP00000491476.1`). The RefSeq ids
    (`NP_001280157.1`, `NP_001280160.1`) are **aliases of the same proteins**, listed on one OpenProt record —
    a database-namespace difference, **not a conflict**. Recorded so it is not "corrected" back.
- **THE RIBO-SEQ AXIS HAS NO DEFENSIBLE NEGATIVES AT THIS LOCUS.** Evidence is reported per resource, in states
  that are never conflated: *DETECTED*, *SUB-THRESHOLD*, *PRESENT-NOT-DETECTED* and *NOT-ASSESSED* (always with
  its reason). **`ABSENT` is now UNUSED on every axis — no ORF qualifies for it.**
  - **A catalogue's coverage is not an assay's power.** Every Ribo-seq resource here is a **catalogue of
    positives**: it publishes the ORFs it **called**, not the candidate set it **tested**. An absence therefore
    cannot distinguish *"tested, no evidence"* from *"never a candidate"* — **and neither is provable from the
    artifact**, so an absence is unfalsifiable.
  - **nuORFdb v1.2** — its **single** MYCN record (**NC-NEST-090**) is **not one of the ORFs in this table**.
    Emitting one record proves the pipeline *can* call an ORF here; it says nothing about the search space.
  - **sORFs.org** — **unpowered for negatives, and demonstrably so**: **MUSEP** (86 aa, a known functional uORF,
    squarely in scope at ≤100 aa) is **missed entirely**, and all 17 hits fall in one ~575 bp island. A catalogue
    proven to miss an in-scope true positive at this locus cannot supply a reliable negative here.
  - **Ribo-uORF** covers the **uORF/uTIS class only**, and its records are per-sample **calls**, not candidates.
  - **Five consensus catalogues return zero records anywhere at the MYCN locus** — **UNPOWERED, NOT-ASSESSED.**
  - **What the axis has is POSITIVES:** 17 sORFs.org records (7 terminating ORFs across 6 regions, **one study**,
    one ~575 bp island, **out-of-frame phasing where tested**), 1 nuORFdb record, 1 Ribo-uORF record.
    **It can say nothing about anything it did not find.**
- **The two axes agree on the ceiling.** MS: *"a protein from this stop group is translated"* — no start
  discrimination. Ribo-seq: *"reads are present at these coordinates"* — no defensible absence. **The ORF
  landscape is a search space with a handful of positives and no reliable negatives on either axis.**
- **Status of the ORFs:** apart from the annotated proteins (N-Myc, ΔMYCN) and the known functional uORFs
  (MYCNOT, MUSEP), the ORFs shown are **predictions**. **ORF 24** (OpenProt IP_083082) carries protein-level
  evidence — 2 unique peptides by lysate MS, two-executor confirmed, frame-disjoint from N-Myc — **but at
  GROUP LEVEL, not as a start claim**: those peptides are contained in 6 of the 8 Region-1 members and in 16
  ORFs across the full 763-ORF catalogue, so **a protein IS translated from that stop group, but WHICH START
  produced it is NOT established**. Detection is not function; and here, detection is not attribution.
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
