// Read the cb++ screen and answer the only question that matters about it:
// does the model flag hint-enriched variants more often than AF-matched variants with no hint?
//
// Without the null arm, "AlphaGenome scored 40 of our 2,000 candidates highly" is unfalsifiable --
// there is no way to tell it from a model that scores 2% of anything highly. Arm C is the answer,
// and it is AF-matched so the comparison is about the HINT rather than about frequency.
//
// The threshold is not chosen here. It comes from data/panels_control.json: 0.1313 is the largest
// splice_maxabsdiff that any of 30 AF-matched deep-intronic negatives produced, so it is the
// ceiling of "nothing is happening" as measured by this pipeline on this kind of variant.
//
//   node scripts/analyse_cbpp.mjs
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const P = JSON.parse(fs.readFileSync(R + 'data/panels_cbpp.json', 'utf8'));
const C = JSON.parse(fs.readFileSync(R + 'data/panels_control.json', 'utf8'));
const negs = C.drawn.filter(d => d.label === 'negative').map(d => d.splice_maxabsdiff).sort((a, b) => a - b);
const THRESH = negs[negs.length - 1];
const posMin = Math.min(...C.drawn.filter(d => d.label === 'positive').map(d => d.splice_maxabsdiff));

const rows = P.drawn;
const arms = ['abundance', 'hint', 'null'];
const by = a => rows.filter(r => r.arm === a);
const med = xs => { if (!xs.length) return null; const s = xs.slice().sort((x, y) => x - y); return s[Math.floor(s.length / 2)]; };

console.log(`cb++ screen: ${rows.length} scored of ${P.drawn.length + P.refused.length + P.failed.length} attempted` +
            `  (${P.refused.length} refused, ${P.failed.length} failed)`);
console.log(`threshold ${THRESH.toFixed(4)} = the highest value 30 AF-matched deep-intronic negatives produced.`);
console.log(`for scale, the weakest canonical splice-killer in that same run scored ${posMin.toFixed(4)}.\n`);

console.log('arm          n     median AG    n > threshold    hit rate');
const hits = {};
for (const a of arms) {
  const s = by(a); if (!s.length) { console.log(`  ${a.padEnd(12)} 0`); continue; }
  const v = s.map(r => r.splice_maxabsdiff);
  hits[a] = s.filter(r => r.splice_maxabsdiff > THRESH);
  console.log(`  ${a.padEnd(11)} ${String(s.length).padStart(4)}   ${med(v).toFixed(4).padStart(9)}   ` +
    `${String(hits[a].length).padStart(13)}    ${(100 * hits[a].length / s.length).toFixed(1)}%`);
}

/* ---- the comparison the null arm exists for ------------------------------------------------- */
function fisher2x2(a, b, c, d) {          // one-sided p, hypergeometric tail
  const lgamma = z => { // Lanczos
    const g = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059,
               12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    z -= 1; let x = 0.99999999999980993;
    for (let i = 0; i < g.length; i++) x += g[i] / (z + i + 1);
    const t = z + g.length - 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  };
  const lfact = n => lgamma(n + 1);
  const p = (a, b, c, d) => Math.exp(lfact(a+b)+lfact(c+d)+lfact(a+c)+lfact(b+d)-lfact(a)-lfact(b)-lfact(c)-lfact(d)-lfact(a+b+c+d));
  let s = 0; const n = a + b + c + d;
  for (let i = a; i <= Math.min(a + b, a + c); i++) s += p(i, a + b - i, a + c - i, d - (i - a));
  return Math.min(1, s);
}

console.log('');
for (const arm of ['hint', 'abundance']) {
  const A = by(arm), N = by('null');
  if (!A.length || !N.length) continue;
  const ah = hits[arm].length, an = A.length - ah, nh = hits['null'].length, nn = N.length - nh;
  const rr = (ah / A.length) / (nh / N.length || 1e-9);
  const p = fisher2x2(ah, an, nh, nn);
  console.log(`${arm} vs null:  ${ah}/${A.length} (${(100*ah/A.length).toFixed(1)}%) vs ` +
    `${nh}/${N.length} (${(100*nh/N.length).toFixed(1)}%)   ` +
    `rate ratio ${rr.toFixed(2)}x   Fisher one-sided p = ${p < 1e-4 ? p.toExponential(1) : p.toFixed(4)}`);
}

/* ---- AF-stratified, because that is the confound ---------------------------------------------- */
console.log('\nhit rate within allele-frequency bands (the null arm is AF-matched, so these compare like with like):');
console.log('AF band            abundance        hint          null');
const bands = [[1e-5,1e-3,'1e-5 – 0.1%'],[1e-3,1e-2,'0.1% – 1%'],[1e-2,5e-2,'1% – 5%'],[5e-2,1,'>= 5%']];
for (const [lo, hi, lab] of bands) {
  const cell = a => { const s = by(a).filter(r => r.af >= lo && r.af < hi);
    return s.length ? `${s.filter(r => r.splice_maxabsdiff > THRESH).length}/${s.length}` : '—'; };
  console.log(`  ${lab.padEnd(16)} ${cell('abundance').padStart(10)} ${cell('hint').padStart(13)} ${cell('null').padStart(13)}`);
}

/* ---- the actual cb++ list ---------------------------------------------------------------------- */
const cb = rows.filter(r => r.splice_maxabsdiff > THRESH)
               .sort((a, b) => b.splice_maxabsdiff - a.splice_maxabsdiff);
console.log(`\n=== cb++ : ${cb.length} variants the model flags above the negative-control ceiling ===`);
console.log('gene       rsid / pos            AF        dist    phyloP   SpliceAI   AlphaGenome  arm');
for (const r of cb.slice(0, 30))
  console.log(`${r.gene.padEnd(10)} ${String(r.rsid || r.pos).padEnd(20)} ${r.af.toExponential(1).padStart(8)}  ` +
    `${String(r.dist_splice).padStart(6)}  ${String(r.phylop ?? '—').padStart(7)}  ${String(r.spliceai_gnomad ?? '—').padStart(8)}   ` +
    `${r.splice_maxabsdiff.toFixed(4).padStart(10)}  ${r.arm}`);
if (cb.length > 30) console.log(`  ... and ${cb.length - 30} more`);

fs.writeFileSync(R + 'data/cbpp.json', JSON.stringify({
  generated_from: 'scripts/analyse_cbpp.mjs',
  threshold: THRESH, threshold_source: 'max splice_maxabsdiff over 30 AF-matched deep-intronic negatives in data/panels_control.json',
  positive_min: posMin, screened: rows.length,
  arm_counts: Object.fromEntries(arms.map(a => [a, by(a).length])),
  arm_hits: Object.fromEntries(arms.map(a => [a, (hits[a] || []).length])),
  rows: cb
}));
console.log(`\nwrote data/cbpp.json  (${cb.length} rows, ${(fs.statSync(R + 'data/cbpp.json').size / 1024).toFixed(0)} KB)`);
