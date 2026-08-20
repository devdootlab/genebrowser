// What is actually IN each deep-intronic density peak, before anything is claimed about it.
//
// A peak is a window holding several times the gene's typical variant count. Two very different
// things produce that shape and they look identical in the strip:
//   real   -- a mutable sequence (CpG), or a region under no constraint
//   artefact -- a homopolymer or tandem repeat, where alignment and length-calling produce a ladder
//               of alleles at one anchor (rs753240633 is 16 rows at ONE position for this reason)
// So this reports the composition and the sequence context and does not decide.
//
// The AlphaGenome column is scaled against data/panels_control.json, which is the known-answer run:
// 30 canonical +-1/+-2 positives (median 0.9915, min 0.9338) and 30 AF-matched deep-intronic
// negatives with SpliceAI 0 (median 0.0040, MAX 0.1313). A score above 0.1313 is outside everything
// 30 matched controls produced.
//
//   node scripts/peak_composition.mjs
import fs from 'node:fs';

import path from 'node:path';
import { fileURLToPath } from 'node:url';
// repo root, derived from this file's location -- no absolute path baked in
const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const P = JSON.parse(fs.readFileSync(R + 'data/panels_peak.json', 'utf8'));
const C = JSON.parse(fs.readFileSync(R + 'data/panels_control.json', 'utf8'));
const negMax = Math.max(...C.drawn.filter(d => d.label === 'negative').map(d => d.splice_maxabsdiff));
const posMin = Math.min(...C.drawn.filter(d => d.label === 'positive').map(d => d.splice_maxabsdiff));
const sleep = ms => new Promise(r => setTimeout(r, ms));

async function seq(chrom, a, b) {
  for (let t = 0; t < 4; t++) {
    try {
      const r = await fetch(`https://rest.ensembl.org/sequence/region/human/${String(chrom).replace('chr','')}:${a}..${b}?content-type=application/json`,
        { headers: { accept: 'application/json' } });
      if (r.ok) return (await r.json()).seq.toUpperCase();
    } catch (e) { /* retry */ }
    await sleep(900 * (t + 1));
  }
  return null;
}
const longestRun = s => { let b = 0, c = 1; for (let i = 1; i < s.length; i++) { c = s[i] === s[i-1] ? c+1 : 1; if (c > b) b = c; } return b; };
// a 2-6mer repeated >=4 times back to back
function longestSTR(s) {
  let best = 0, unit = '';
  for (let u = 2; u <= 6; u++) for (let i = 0; i + u * 4 <= s.length; i++) {
    const m = s.substr(i, u); if (/^(.)\1*$/.test(m)) continue;
    let n = 1; while (s.substr(i + n * u, u) === m) n++;
    if (n >= 4 && n * u > best) { best = n * u; unit = `(${m})x${n}`; }
  }
  return { len: best, unit };
}

const rows = [];
for (const d of P.drawn) {
  const g = GENES[d.gene];
  const pk = JSON.parse(fs.readFileSync(R + `data/variants_${d.gene}.json`, 'utf8'));
  const ix = n => pk.cols.indexOf(n);
  // the peak window that produced this pick: recover it from make_panels' own bin arithmetic
  const binw = 200, half = 600;
  const lo = d.pos - half, hi = d.pos + half;
  const inWin = pk.rows.filter(r => r[ix('pos')] >= lo && r[ix('pos')] <= hi && !r[ix('coding')]);
  const nullAF = inWin.filter(r => r[ix('af')] == null).length;
  const indel = inWin.filter(r => String(r[ix('ref')]).length !== 1 || String(r[ix('alt')]).length !== 1).length;
  const posns = new Set(inWin.map(r => r[ix('pos')]));
  const s = await seq(g.chrom, lo, hi);
  const str = s ? longestSTR(s) : { len: 0, unit: '' };
  rows.push({
    id: d.id, gene: d.gene, pos: d.pos, dist: d.dist_splice, ag: d.splice_maxabsdiff,
    n: inWin.length, uniq: posns.size, nullAF: nullAF, indel,
    hp: s ? longestRun(s) : 0, str: str.len, unit: str.unit,
    cpg: s ? (s.match(/CG/g) || []).length : 0
  });
  process.stdout.write(`\r  ${rows.length}/${P.drawn.length}`);
  await sleep(250);
}
console.log('\n');
console.log(`AlphaGenome scale from the control run: negatives max ${negMax.toFixed(4)}, positives min ${posMin.toFixed(4)}\n`);
console.log('peak (±600bp)          AGmax   vs ctrl   rows  uniq pos  alleles/pos  %noAF  %indel  homopol  STR');
rows.sort((a, b) => b.ag - a.ag);
for (const r of rows) {
  const verdict = r.ag >= posMin ? 'SPLICE-LIKE' : r.ag > negMax ? 'above all' : 'within';
  console.log(
    (r.gene + ':' + r.pos).padEnd(22) +
    r.ag.toFixed(4).padStart(7) + '  ' + verdict.padEnd(11) +
    String(r.n).padStart(5) + String(r.uniq).padStart(9) +
    (r.n / r.uniq).toFixed(2).padStart(12) +
    (100 * r.nullAF / r.n).toFixed(0).padStart(7) +
    (100 * r.indel / r.n).toFixed(0).padStart(8) +
    String(r.hp).padStart(9) + '  ' + (r.str >= 12 ? r.str + 'bp ' + r.unit : '—'));
}
console.log('\nalleles/pos > 1.5 or a long homopolymer/STR means the count is inflated by one mutable');
console.log('site, not by many independent variants. That does not make the AlphaGenome score wrong —');
console.log('the score is about the single substitution tested — but it does mean the PEAK is not');
console.log('evidence of anything on its own.');
fs.writeFileSync(R + 'data/peak_composition.json', JSON.stringify({ negMax, posMin, rows }, null, 1));
console.log('\nwrote data/peak_composition.json');
