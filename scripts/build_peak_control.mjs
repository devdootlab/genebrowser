// A DENSITY-MATCHED control for the peak run.
//
// The peak representatives were compared against the cb++ screen's null arm, which is AF-matched
// deep-intronic variants -- matched on frequency but NOT on the thing under test. "Peaks beat that
// null" could therefore be any property of dense regions rather than density itself, which is the
// same mistake the cb++ hint arm made and it should not be made twice.
//
// This draws, for each peak, one variant from the SAME GENE, the SAME allele-frequency band, and
// the same deep-intronic distance rule, but from a bin that did NOT clear the density threshold.
// Everything is held constant except the one variable the filter selects on.
//
//   node scripts/build_peak_control.mjs      # writes data/peak_control_set.json
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const PEAK = JSON.parse(fs.readFileSync(R + 'data/panels_peak.json', 'utf8'));
const MINDIST = 300, BINW = 200, SEED = 20260819;

// same seeded generator as the cb++ screen -- an unreproducible control is not a control
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const rnd = mulberry32(SEED);
const band = af => af == null ? 'none' : String(Math.round(Math.log10(af) * 2) / 2);

const byGene = {};
for (const d of PEAK.drawn) (byGene[d.gene] ||= []).push(d);

const out = [];
const shortfall = [];
for (const [sym, peaks] of Object.entries(byGene)) {
  const g = GENES[sym];
  const p = JSON.parse(fs.readFileSync(R + `data/variants_${sym}.json`, 'utf8'));
  const i = n => p.cols.indexOf(n);
  const ex = g.exons.map(e => [Math.min(...e), Math.max(...e)]).sort((a, b) => a[0] - b[0]);
  if (ex.length < 2) continue;
  const bounds = []; for (let k = 0; k < ex.length - 1; k++) bounds.push(ex[k][1], ex[k + 1][0]);

  // rebuild the deep set and the bin counts exactly as peaks_for does
  const deep = [];
  for (const r of p.rows) {
    const pos = r[i('pos')];
    if (r[i('coding')]) continue;
    if (String(r[i('ref')]).length !== 1 || String(r[i('alt')]).length !== 1) continue;
    if (!'ACGT'.includes(r[i('ref')]) || !'ACGT'.includes(r[i('alt')])) continue;
    let ok = false; for (let k = 0; k < ex.length - 1; k++) if (ex[k][1] < pos && pos < ex[k + 1][0]) { ok = true; break; }
    if (!ok) continue;
    if (Math.min(...bounds.map(b => Math.abs(pos - b))) < MINDIST) continue;
    deep.push({ pos, ref: r[i('ref')], alt: r[i('alt')], rsid: r[i('rsid')], af: r[i('af')],
                spliceai: r[i('spliceai')], phylop: r[i('phylop')] });
  }
  const bins = {};
  for (const v of deep) (bins[Math.floor(v.pos / BINW)] ||= new Set()).add(v.pos);
  const npos = {}; for (const k in bins) npos[k] = bins[k].size;
  const counts = Object.values(npos).sort((a, b) => a - b);
  const med = counts[counts.length >> 1];
  const mad = counts.map(c => Math.abs(c - med)).sort((a, b) => a - b)[counts.length >> 1] * 1.4826;
  const thresh = Math.max(med + 4 * mad, 5);
  // every bin the peak run used, INCLUDING the gene-floor ones, is off limits as a control
  const used = new Set();
  for (const d of peaks) for (let b = Math.floor(d.peak_start / BINW); b < Math.floor(d.peak_end / BINW); b++) used.add(b);
  const coldPool = deep.filter(v => {
    const b = Math.floor(v.pos / BINW);
    return !used.has(b) && npos[b] < thresh;
  });

  for (const d of peaks) {
    const want = band(d.af);
    let cand = coldPool.filter(v => band(v.af) === want && !out.some(o => o.pos === v.pos && o.gene === sym));
    let matched = 'af-band';
    if (!cand.length) { cand = coldPool.filter(v => !out.some(o => o.pos === v.pos && o.gene === sym)); matched = 'gene-only'; }
    if (!cand.length) { shortfall.push(sym + ':' + d.pos); continue; }
    const pick = cand[Math.floor(rnd() * cand.length)];
    out.push({ gene: sym, ...pick, matched_to: d.id, matched_on: matched,
               peak_via: d.picked_via, peak_af: d.af });
  }
}

console.log(`${PEAK.drawn.length} peaks -> ${out.length} density-matched controls`);
const m = {}; for (const o of out) m[o.matched_on] = (m[o.matched_on] || 0) + 1;
console.log('  matching quality: ' + JSON.stringify(m));
if (shortfall.length) console.log(`  ${shortfall.length} peaks had no usable cold-bin control: ${shortfall.slice(0, 6).join(', ')}`);
const withSA = out.filter(o => o.spliceai != null && o.spliceai >= 0.05).length;
console.log(`  controls with SpliceAI >= 0.05: ${withSA} (${(100 * withSA / out.length).toFixed(0)}%)`);
fs.writeFileSync(R + 'data/peak_control_set.json', JSON.stringify({
  generated_from: 'scripts/build_peak_control.mjs', seed: SEED,
  design: 'one variant per peak: same gene, same AF band, deep intronic >=300bp, from a bin BELOW the density threshold and not used by any peak',
  n: out.length, rows: out }, null, 1));
console.log(`wrote data/peak_control_set.json`);
console.log('');
console.log('queue string for scoring them:');
console.log(out.map(o => (o.rsid || o.pos) + ':' + o.ref + o.alt).join(','));
