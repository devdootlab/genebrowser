// Pick a known-positive and a known-negative set, so AlphaGenome can be asked a question whose
// answer we already know.
//
// WHY
// r = -0.05 between SpliceAI and AlphaGenome across 170 variants says one of three things: they
// measure different quantities, one is wrong, or our extraction is wrong. None of those can be told
// apart from a correlation. What separates them is a set where the truth is not in dispute.
//
//   POSITIVE  canonical +-1/+-2, SpliceAI >= 0.9. These destroy a splice site. Any usable metric
//             must light up here.
//   NEGATIVE  deep intronic, >500 bp from any boundary, SpliceAI <= 0.01, AND matched on allele
//             frequency to the positives so the two sets differ in position, not in rarity.
//
// If the two sets do not separate, the metric is noise and every ranking built on it is arbitrary.
// If they separate on splice_maxabsdiff but not splice_log2fc, the pseudocount is the problem.
//
//   node scripts/synthetic_control.mjs        -> data/control_set.json
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
// repo root, derived from this file's location -- no absolute path baked in
const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));

const POS = [], NEG = [];
for (const gene of Object.keys(GENES)) {
  const g = GENES[gene];
  const p = JSON.parse(fs.readFileSync(R + `data/variants_${gene}.json`, 'utf8'));
  const i = Object.fromEntries(p.cols.map((c, k) => [c, k]));
  const ex = g.exons.map(e => [Math.min(...e), Math.max(...e)]).sort((a, b) => a[0] - b[0]);
  const bounds = [];
  for (let k = 0; k < ex.length - 1; k++) bounds.push(ex[k][1], ex[k + 1][0]);

  for (const r of p.rows) {
    const pos = r[i.pos], ref = r[i.ref], alt = r[i.alt];
    if (ref.length !== 1 || alt.length !== 1 || !'ACGT'.includes(ref) || !'ACGT'.includes(alt)) continue;
    let intron = -1;
    for (let k = 0; k < ex.length - 1; k++) if (pos > ex[k][1] && pos < ex[k + 1][0]) { intron = k; break; }
    if (intron < 0) continue;
    const d = Math.min(...bounds.map(b => Math.abs(pos - b)));
    const sa = r[i.spliceai], af = r[i.af];
    const rec = { gene, chrom: p.chrom, pos, ref, alt, af, spliceai: sa, dist: d,
                  rsid: r[i.rsid], intron: intron + 1, nintrons: ex.length - 1 };
    if (d <= 2 && sa != null && sa >= 0.9) POS.push(rec);
    else if (d > 500 && sa != null && sa <= 0.01) NEG.push(rec);
  }
}

// Match on allele frequency. Without this the negatives are systematically rarer than the
// positives and any difference could be a frequency effect rather than a splice one.
const band = v => v == null ? 'na' : v >= 1e-3 ? 'a' : v >= 1e-4 ? 'b' : v >= 1e-5 ? 'c' : 'd';
const byBand = {};
NEG.forEach(n => (byBand[band(n.af)] ||= []).push(n));
// deterministic pick: sort within band, take evenly spaced, so re-runs give the same control set
for (const k of Object.keys(byBand)) byBand[k].sort((a, b) => a.pos - b.pos);
const want = {};
POS.forEach(p => { want[band(p.af)] = (want[band(p.af)] || 0) + 1; });

const picked = [];
for (const [b, n] of Object.entries(want)) {
  const pool = byBand[b] || [];
  if (!pool.length) { console.log(`  no negative available in AF band ${b} (needed ${n})`); continue; }
  const step = Math.max(1, Math.floor(pool.length / n));
  for (let k = 0; k < n && k * step < pool.length; k++) picked.push(pool[k * step]);
}

const CAP = 30;
const pos = POS.slice(0, CAP), neg = picked.slice(0, CAP);
console.log(`positives available ${POS.length}, negatives available ${NEG.length}`);
console.log(`taking ${pos.length} positive, ${neg.length} negative (AF-matched)`);
console.log('\nAF band distribution:');
for (const b of ['a', 'b', 'c', 'd', 'na']) {
  const P = pos.filter(x => band(x.af) === b).length, N = neg.filter(x => band(x.af) === b).length;
  if (P || N) console.log(`  ${b}: positive ${P}  negative ${N}`);
}
console.log('\npositive: median dist to splice site', pos.map(x => x.dist).sort((a,b)=>a-b)[pos.length>>1],
            ' median SpliceAI', pos.map(x => x.spliceai).sort((a,b)=>a-b)[pos.length>>1]);
console.log('negative: median dist to splice site', neg.map(x => x.dist).sort((a,b)=>a-b)[neg.length>>1],
            ' median SpliceAI', neg.map(x => x.spliceai).sort((a,b)=>a-b)[neg.length>>1]);

fs.writeFileSync(R + 'data/control_set.json', JSON.stringify({ positive: pos, negative: neg }, null, 1));
console.log('\nwrote data/control_set.json');
console.log('next: python scripts/make_panels.py --control');
