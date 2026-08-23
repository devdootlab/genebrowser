// cb++ — variants built to survive a reviewer whose null hypothesis is "you picked something dumb".
//
// The four ways a deep-intronic pick is dumb, and what each becomes here:
//
//   1. IT IS NOT REAL          -- a calling artefact                     -> HARD GATE
//   2. IT HAS NO FUNCTION      -- would come back flat in an MPRA        -> ranked term
//   3. IT IS TOO MINOR         -- no human relevance                     -> ranked term
//   4. IT IS OBVIOUS           -- a donor/acceptor hit, we already knew   -> HARD GATE
//
// 1 and 4 are DISQUALIFYING. No allele frequency should rescue a homopolymer stutter, and no
// conservation score should make a canonical +1 interesting again. So they gate rather than score.
// 2 and 3 are matters of degree, so they rank.
//
// WHAT THIS IS NOT: a measurement. It is a shortlisting heuristic with named parts, and every part
// is stored per row so the table can show WHY something ranked, not just that it did. A single
// opaque number would let a high-frequency-low-conservation variant and a low-frequency-high-
// conservation one land on the same score with no way to tell them apart.
//
// THE TENSION WORTH SAYING OUT LOUD: abundance fights function. A common deep-intronic variant in a
// coagulation gene is, by purifying selection, probably not severely deleterious. cb++ hits will
// therefore skew towards QUANTITATIVE TRAIT MODIFIERS -- factor levels, VWF antigen, PT/aPTT --
// rather than disease alleles. That is still coagulopathic physiology and still human-relevant, but
// it is a different claim from "this causes haemophilia" and must be framed as such.
//
// WHAT THIS SCRIPT ESTABLISHED, AND WHY ITS OUTPUT IS NOT USED
// ------------------------------------------------------------
// Run with abundance weighted at 0.40 as originally specified, the top 2,000 came back with median
// phyloP -0.02 -- literally unconserved -- and median AF 16%, plus 12 rows at AF >= 0.999 where the
// REFERENCE is the minor allele. Ranking on abundance produces a list that fails the very objection
// it was built to answer ("no obvious function, would fail an MPRA").
//
// The trade-off underneath it is real and is the reason: of 108,755 variants that pass the reality
// and non-obviousness gates, only 1,436 reach phyloP >= 2, and their median AF is 2.6e-5. Every
// route to "this has function" lands on variants at about 1 in 40,000. That is purifying selection,
// observed. It is why 2,000 cb++ variants do not exist as a static list.
//
// So this script is kept as the DIAGNOSTIC that established the trade-off, and its ranked output is
// deliberately not written to disk any more -- scripts/build_cbpp_screen.mjs supersedes it, turning
// 2,000 from an answer into an AlphaGenome screening budget with a matched null arm.
//
//   node scripts/score_cbpp.mjs            # prints the gate census and the trade-off
//   node scripts/score_cbpp.mjs --top 500
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const WEIGHT = JSON.parse(fs.readFileSync(R + 'data/gene_weights.json', 'utf8'));
const TOP = (() => { const i = process.argv.indexOf('--top'); return i > 0 ? +process.argv[i + 1] : 2000; })();

// Every AlphaGenome result we already have, keyed the way the page keys them.
const PANEL = {};
for (const f of fs.readdirSync(R + 'data').filter(f => /^panels.*\.json$/.test(f))) {
  for (const d of (JSON.parse(fs.readFileSync(R + 'data/' + f, 'utf8')).drawn || []))
    PANEL[d.gene + ':' + d.pos + ':' + d.ref + d.alt] = d.splice_maxabsdiff;
}

/* ---- gate helpers -------------------------------------------------------------------------- */
const runAt = (s, i) => {            // length of the homopolymer run containing index i
  const c = s[i]; if (!c) return 0;
  let a = i, b = i;
  while (a > 0 && s[a - 1] === c) a--;
  while (b < s.length - 1 && s[b + 1] === c) b++;
  return b - a + 1;
};
function strAround(s, i, half = 24) { // longest 2-6mer tandem repeat overlapping index i
  const lo = Math.max(0, i - half), w = s.slice(lo, i + half);
  let best = 0;
  for (let u = 2; u <= 6; u++) for (let k = 0; k + u * 4 <= w.length; k++) {
    const m = w.substr(k, u); if (/^(.)\1*$/.test(m)) continue;
    let n = 1; while (w.substr(k + n * u, u) === m) n++;
    if (n >= 4) { const st = lo + k, en = st + n * u; if (i >= st && i < en && n * u > best) best = n * u; }
  }
  return best;
}

/* ---- ranked terms, each 0..1 ---------------------------------------------------------------- */
// Abundance. The user's headline term. AF spans six orders of magnitude, so rank on log10 and put
// the knee where a variant starts being findable in a real cohort: 1e-4 scores ~0, 1e-2 scores ~0.67.
const fAbund = af => Math.min(1, Math.max(0, (Math.log10(af) + 4) / 3));
// Conservation. phyloP 100-way. The best proxy available offline for "there is something here to
// break". Anything at or below 0 is unconstrained; 4+ is strongly conserved.
const fCons = p => p == null ? 0 : Math.min(1, Math.max(0, p / 4));
// CADD, on its usual PHRED-like scale. 20 = top 1% of the genome.
const fCadd = c => c == null ? 0 : Math.min(1, Math.max(0, c / 25));
// Subtlety. Zero at the splice-region boundary, full marks by 200 bp in and flat thereafter. Flat,
// not falling: deep regulatory elements exist at every distance, and a decay term would be an
// invented prior about where introns keep their function.
const fSubtle = d => Math.min(1, Math.max(0, (d - 100) / 100));
// The ML term. Scaled against the control run: 0.1313 is the largest value 30 AF-matched
// deep-intronic negatives produced, 0.9338 the smallest a canonical splice-killer produced.
const fML = v => v == null ? null : Math.min(1, Math.max(0, (v - 0.1313) / (0.9338 - 0.1313)));

// Relative weights. Abundance dominates because that is what was asked for. ML is scored only when
// present; the remaining weights are renormalised so a variant is never punished for not having had
// an API call spent on it yet.
const W = { abund: 0.40, cons: 0.20, cadd: 0.15, subtle: 0.10, ml: 0.15 };

const rows = [];
const reject = { noAF: 0, notSNV: 0, coding: 0, exonic: 0, tooCloseToSplice: 0, homopolymer: 0, str: 0, multiallelic: 0, noWeight: 0, afFloor: 0 };
let considered = 0;

for (const [sym, g] of Object.entries(GENES)) {
  const gw = WEIGHT[sym];
  if (!gw) { continue; }                                    // weight 0 or absent = out of domain
  const pack = JSON.parse(fs.readFileSync(R + `data/variants_${sym}.json`, 'utf8'));
  const ix = n => pack.cols.indexOf(n);
  const [Irs, Ipos, Iref, Ialt, Icons, Iaf, , Isa, Ipg, Icadd, Iphy, Icod] =
    ['rsid','pos','ref','alt','cons','af','curated','spliceai','pangolin','cadd','phylop','coding'].map(ix);

  const seq = fs.existsSync(R + `data/seq_${sym}.txt`)
    ? fs.readFileSync(R + `data/seq_${sym}.txt`, 'utf8').trim() : null;

  const ex = g.exons.map(e => [Math.min(...e), Math.max(...e)]).sort((a, b) => a[0] - b[0]);
  const bounds = []; for (let i = 0; i < ex.length - 1; i++) bounds.push(ex[i][1], ex[i + 1][0]);

  // multi-allelic pile-ups: one mutable site throwing a ladder of alleles is the rs753240633
  // signature, and it is the single clearest tell that a position is a repeat, not a variant.
  const perPos = {};
  for (const r of pack.rows) perPos[r[Ipos]] = (perPos[r[Ipos]] || 0) + 1;

  for (const r of pack.rows) {
    considered++;
    const pos = r[Ipos], ref = String(r[Iref]), alt = String(r[Ialt]), af = r[Iaf];

    /* ---- GATE 1: is it real? ---- */
    if (af == null)                          { reject.noAF++; continue; }          // exome-only call
    if (af < 1e-5)                           { reject.afFloor++; continue; }       // singleton territory
    if (ref.length !== 1 || alt.length !== 1 || !'ACGT'.includes(ref) || !'ACGT'.includes(alt))
                                             { reject.notSNV++; continue; }
    if (perPos[pos] > 3)                     { reject.multiallelic++; continue; }
    if (seq) {
      const i = pos - g.start;
      if (i >= 0 && i < seq.length) {
        if (runAt(seq, i) >= 8)              { reject.homopolymer++; continue; }
        if (strAround(seq, i) >= 12)         { reject.str++; continue; }
      }
    }

    /* ---- GATE 2: is it non-obvious? ---- */
    if (r[Icod])                             { reject.coding++; continue; }
    let intron = -1;
    for (let i = 0; i < ex.length - 1; i++) if (ex[i][1] < pos && pos < ex[i + 1][0]) { intron = i; break; }
    if (intron < 0)                          { reject.exonic++; continue; }
    const dist = Math.min(...bounds.map(b => Math.abs(pos - b)));
    if (dist < 100)                          { reject.tooCloseToSplice++; continue; }

    /* ---- ranked terms ---- */
    const ml = fML(PANEL[sym + ':' + pos + ':' + ref + alt]);
    const parts = {
      abund:  fAbund(af),
      cons:   fCons(r[Iphy]),
      cadd:   fCadd(r[Icadd]),
      subtle: fSubtle(dist),
      ml:     ml
    };
    // renormalise over the terms actually available
    let num = 0, den = 0;
    for (const k of Object.keys(W)) { if (parts[k] == null) continue; num += W[k] * parts[k]; den += W[k]; }
    const score = gw * (num / den);

    rows.push({
      gene: sym, rsid: r[Irs], pos, ref, alt, af,
      phylop: r[Iphy], cadd: r[Icadd], spliceai: r[Isa], pangolin: r[Ipg],
      dist_splice: dist, intron: (intron + 1) + '/' + (ex.length - 1),
      gene_weight: gw, ml_raw: PANEL[sym + ':' + pos + ':' + ref + alt] ?? null,
      parts, score: +score.toFixed(5)
    });
  }
}

rows.sort((a, b) => b.score - a.score);
const top = rows.slice(0, TOP);

/* ---- report -------------------------------------------------------------------------------- */
const pc = n => ((100 * n) / considered).toFixed(1) + '%';
console.log(`considered ${considered.toLocaleString()} rows across ${Object.keys(GENES).filter(g => WEIGHT[g]).length} in-domain genes\n`);
console.log('rejected by gate:');
for (const [k, v] of Object.entries(reject).sort((a, b) => b[1] - a[1]))
  if (v) console.log(`  ${k.padEnd(20)} ${String(v).padStart(8)}  ${pc(v).padStart(6)}`);
console.log(`\n${rows.length.toLocaleString()} passed every gate  (${pc(rows.length)})`);
if (rows.length < TOP)
  console.log(`\n  *** only ${rows.length} candidates exist, fewer than the ${TOP} asked for.\n` +
              `  *** Reporting all of them rather than padding the list with rows that failed a gate.`);

const withML = rows.filter(r => r.parts.ml != null).length;
console.log(`${withML} of them already have an AlphaGenome score; the rest are ranked on the other four terms.`);

const q = (arr, p) => arr.length ? arr[Math.floor((arr.length - 1) * p)] : null;
const afs = top.map(r => r.af).sort((a, b) => a - b);
const phys = top.map(r => r.phylop).filter(x => x != null).sort((a, b) => a - b);
console.log(`\ntop ${top.length}:`);
console.log(`  allele frequency   median ${q(afs, .5).toExponential(2)}   range ${q(afs, 0).toExponential(2)} – ${q(afs, 1).toExponential(2)}`);
console.log(`  phyloP             median ${q(phys, .5)?.toFixed(2)}   90th pct ${q(phys, .9)?.toFixed(2)}`);
console.log(`  distance to splice median ${q(top.map(r => r.dist_splice).sort((a, b) => a - b), .5).toLocaleString()} bp`);
const byGene = {}; for (const r of top) byGene[r.gene] = (byGene[r.gene] || 0) + 1;
console.log('  genes: ' + Object.entries(byGene).sort((a, b) => b[1] - a[1]).slice(0, 8).map(([g, n]) => g + ' ' + n).join(', '));

// The ranked list is deliberately NOT written. See the header: this weighting produces a top-2,000
// whose median phyloP is -0.02, which fails the objection the score exists to answer. Leaving the
// file on disk would let something downstream read a list that has been argued against, and a
// retracted result that is still loadable is a retracted result that will be used.
console.log(`\nno file written. This script is the diagnostic that established the abundance/function`);
console.log(`trade-off; the screening set that supersedes it comes from scripts/build_cbpp_screen.mjs.`);
