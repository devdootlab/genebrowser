# genebrowser

A browser for **non-coding variants** across the coagulation cascade — every clotting factor and
natural anticoagulant — ranked by predicted splice disruption.

**Factors** F2 · F3 · F5 · F7 · F8 · F9 · F10 · F11 · F12 · F13A1 · F13B · FGA · FGB · FGG · VWF
**Contact** KLKB1 · KNG1
**Anticoagulants** SERPINC1 · PROC · PROS1 · THBD · PROCR · TFPI · SERPIND1 · SERPINA5 · PROZ
**Fibrinolysis** PLG · PLAT · SERPINE1 · SERPINF2
**Related** ADAMTS13 · VKORC1 · GGCX · LMAN1 · MCFD2
**Immune** OAS1 · GSDMB

```bash
npx http-server . -p 8089 -c-1     # then open /
node scripts/test.mjs              # 481 acceptance checks
```

## What it does

Pulls every gnomAD variant in each gene, attaches SpliceAI / Pangolin / CADD / phyloP, and draws a
4-panel AlphaGenome assessment for the ones predicted to disrupt splicing.

The premise is that a language model can only return variants that have been *written about*, so
anything it offers is by construction already studied. gnomAD returns what nobody has looked at.

## Numbers

| | |
|---|---:|
| genes | 37 |
| gnomAD variants | 614,941 |
| non-coding | 563,259 |
| with a SpliceAI score | 536,806 |
| SpliceAI ≥ 0.5 | 2,699 |
| AlphaGenome panels | 555 |

## The interface

Table first, the way gnomAD and ClinVar are — nobody finds a variant by looking at a lollipop.

- sortable on SpliceAI, phyloP, AlphaGenome effect, frequency, position, **variant type**
- a **density strip** with a real y-axis; click any bin to filter the table to it
- a **zoom** showing GRCh38 bases, polypyrimidine tracts, splice dinucleotides and exon edges —
  **all four computed on the transcribed strand** (see below)
- **indels drawn as indels**: ● substitution, ▼ deletion with a band over the bases it removes,
  ▲ insertion, ■ multi-base substitution
- **density outliers classified** as calling artefact (homopolymer, tandem repeat, segmental
  duplication) or real GC/CpG-driven mutability — with the evidence in the tooltip
- coding variants hidden by default, so every gene shows the same population

## Splice motifs are read on the transcribed strand

Until 2026-08-18 this page searched the **plus strand** for `GT` and `AG` on every gene. On a
minus-strand gene a donor `GT` reads `AC` on the plus strand and an acceptor `AG` reads `CT`, so it
boxed dinucleotides that are not splice sites and drew nothing at the ones that are. The
polypyrimidine tract detector had the same bug — it counted C/T on the plus strand, which on a
minus-strand gene is the purines.

`scripts/check_strand_motifs.mjs` is the known-answer run. ~99% of human introns are GT-AG **in
transcript orientation**, which is a fact about biology rather than about this repo:

```
                       introns  strand-aware(fixed)  plus-strand-only(old)
18 minus-strand genes      236              236/236                  0/236
18 plus-strand genes       173              169/173                171/173

strand-aware    : 405/409 canonical GT-AG  (99.0%)
plus-strand-only: 171/409 canonical GT-AG  (41.8%)
```

**236 of 409 introns — 57.7% of every intron here — were annotated wrong**, split perfectly clean
along strand. The 4 residual misses are GC-AG minor introns, ~1% of human introns.

For a minus-strand gene the zoom now draws a second sequence row: the complement, in place, read
right to left.

## What the numbers mean, and what they do not

**Frequency is measured.** SpliceAI, Pangolin, CADD and phyloP come from gnomAD.

**The AlphaGenome effect is a prediction**, reported as `splice_maxabsdiff` — the largest change in
the splice-site track within 256 bp of the variant, in the single most-affected track of 667.

That metric was chosen by experiment, not taste. Against 30 known-positives (canonical ±1/±2 splice
variants) and 30 allele-frequency-matched known-negatives (>500 bp from any boundary, SpliceAI 0):

| metric | positives | negatives | overlap |
|---|---:|---:|---:|
| `splice_log2fc` | 0.427 | 0.203 | **100%** |
| `delta_log2fc` (RNA) | 0.287 | 0.046 | 80% |
| **`splice_maxabsdiff`** | **0.992** | **0.004** | **0%** |

The log-ratio metrics cannot distinguish a destroyed splice site from a random intronic base. The
control set is in `data/control_set.json` and is a regression test: any change must keep that
overlap at zero.

SpliceAI and AlphaGenome agree at **r = 0.72** across 503 variants on `splice_maxabsdiff`.

**There is no null distribution.** Nothing here tells you whether a given effect size is unusual for
a random intronic base at that locus. Treat the ranking as ordinal.

**A null gnomAD genome AF means the variant was called in the exomes, not the genomes.** Deep in an
intron that means off-target capture — low depth, error-prone. It is 14–56% of rows gene-wide and up
to 82% inside a density peak. Those rows are drawn hollow and there is a *genome-called only* filter.

**gnomAD is ~800,000 largely healthy people.** This finds variants that exist in the population and
are predicted to disrupt splicing — not variants that cause disease in a patient. It is a discovery
tool, not a diagnostic filter.

## Adding a gene

```bash
node scripts/add_genes.mjs SERPINC1 PROC PROS1     # structure + variants + predictors, one pass
PY=/path/to/python scripts/run_peak_panels.sh      # one panel per gene, deep-intronic peaks
```

The Ensembl ID is looked up **from the symbol**, never hand-entered, and every variant must fall
inside its own gene's span or the file is not written. Both guards exist because a hand-entered ID
once put 62 SSTR2 variants into F2 with figures drawn at coordinates belonging to neither gene.

## Refusing to draw

`make_panels.py` runs preflight checks before any API call:

- chromosome and position inside the named gene
- the stated reference base **is** the base at that coordinate
- alleles are single ACGTN (AlphaGenome scores substitutions only)

A failure is refused with its reason printed, and `drawn` / `refused` / `failed` are recorded
separately — a missing figure is never indistinguishable from one nobody asked for. The manifest
**merges** across runs, because it used to be overwritten and a 16-gene loop drew 16 figures and
recorded 1, with no error.

The canonical ±1/±2 bucket is the known-answer set: those destroy a splice site, so the run fails
loudly if their median effect is near zero. It once reported 0.005 and caught two bugs that had
already shipped 40 figures.

## Known-wrong

1. **The deep-intronic peak selection does not do what its docstring says.** `--bucket peak` claims
   to pick "strongest SpliceAI, then commonest" per peak; none of the 16 picks has a SpliceAI score
   and most have no AF, so every sort key was 0 and it silently fell through to distance order. The
   16 panels in `data/panels_peak.json` are real and preflighted, but **the selection is not what is
   documented and the results should not be read as a finding** until it is fixed and run against a
   matched null.
2. **Exons render empty in the density strip.** Coding is hidden by default, but the strip uses the
   filtered list for the backdrop as well as the overlay.
3. **`--queue` is not implemented.** The right-click queue copies a command `make_panels.py` cannot
   parse.
4. **5 of 11 curated rows have a reference base that is not at that coordinate.** Two are strand
   complements, three are not — different fixes. See `data/ref_mismatch.json`.
5. **9 of 16 density outliers are unexplained** — no repeat, no GC/CpG elevation. There is no
   significance test behind the 4× threshold.
6. **F8 is on chrX.** Males are hemizygous, so gnomAD AF is over mixed ploidy and the frequency bands
   are calibrated for autosomes.
7. **21 of 37 genes have no AlphaGenome panel at all.**

## Retracted

"SpliceAI and AlphaGenome do not agree, r = −0.05" was measured with `splice_log2fc`, which has 100%
positive/negative overlap on the control set above. The correct metric gives **r = 0.72**.

## Figures

555 panel PNGs, **not committed** — regenerable from `data/panels*.json` via `scripts/make_panels.py`.
Needs an AlphaGenome API key in `ALPHAGENOME_API_KEY`.

## Data sources

gnomAD v4.1 (GraphQL) · Ensembl REST (GRCh38) · UCSC `genomicSuperDups` · AlphaGenome
