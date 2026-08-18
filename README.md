# genebrowser

A browser for **non-coding variants** in ten coagulation and immune genes, ranked by predicted
splice disruption.

**F2 · F5 · F8 · VWF · FGA · FGB · FGG · SERPINC1 · OAS1 · GSDMB**

```bash
npx http-server . -p 8089 -c-1     # then open /browse.html
node scripts/test.mjs              # acceptance checks
```

## What it does

Pulls every gnomAD variant in each gene, attaches SpliceAI / Pangolin / CADD / phyloP, and draws a
4-panel AlphaGenome assessment for the ones predicted to disrupt splicing.

The premise is that a language model can only return variants that have been *written about*, so
anything it offers is by construction already studied. gnomAD returns what nobody has looked at.

## Numbers

| | |
|---|---:|
| genes | 10 |
| gnomAD variants | ~190,000 |
| with a SpliceAI score | ~180,000 |
| SpliceAI ≥ 0.5 | ~780 |
| AlphaGenome panels | 539 |

## The interface

Table first, the way gnomAD and ClinVar are — nobody finds a variant by looking at a lollipop.

- sortable on SpliceAI, phyloP, AlphaGenome effect, frequency, position
- a **density strip** with a real y-axis; click any bin to filter the table to it
- a **zoom** showing GRCh38 bases, polypyrimidine tracts, GT/AG dinucleotides and exon edges
- **density outliers classified** as calling artefact (homopolymer, tandem repeat, segmental
  duplication) or real GC/CpG-driven mutability — with the evidence in the tooltip
- coding variants hidden by default, so every gene shows the same population

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

**gnomAD is ~800,000 largely healthy people.** This finds variants that exist in the population and
are predicted to disrupt splicing — not variants that cause disease in a patient. It is a discovery
tool, not a diagnostic filter.

## Adding a gene

```bash
node scripts/add_genes.mjs SERPINC1        # structure + variants + predictors, one pass
python scripts/make_panels.py --all --bucket spliceai --spliceai 0.5 --limit 400
node scripts/build-readme.mjs
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
separately — a missing figure is never indistinguishable from one nobody asked for.

The canonical ±1/±2 bucket is the known-answer set: those destroy a splice site, so the run fails
loudly if their median effect is near zero. It once reported 0.005 and caught two bugs that had
already shipped 40 figures.

## Known-wrong

1. **Exons render empty in the density strip.** Coding is hidden by default, but the strip uses the
   filtered list for the backdrop as well as the overlay.
2. **`--queue` is not implemented.** The right-click queue copies a command `make_panels.py` cannot
   parse.
3. **5 of 11 curated rows have a reference base that is not at that coordinate.** Two are strand
   complements, three are not — different fixes. See `data/ref_mismatch.json`.
4. **9 of 16 density outliers are unexplained** — no repeat, no GC/CpG elevation. There is no
   significance test behind the 4× threshold.
5. **F8 is on chrX.** Males are hemizygous, so gnomAD AF is over mixed ploidy and the frequency bands
   are calibrated for autosomes.

## Figures

539 panel PNGs, ~74 MB, **not committed** — regenerable from `data/panels.json` via
`scripts/make_panels.py`. Needs an AlphaGenome API key in `ALPHAGENOME_API_KEY`.

## Data sources

gnomAD v4.1 (GraphQL) · Ensembl REST (GRCh38) · UCSC `genomicSuperDups` · AlphaGenome
