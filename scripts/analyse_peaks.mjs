// Does variant DENSITY predict an AlphaGenome splice effect?
//
// Three scorings of the same 172 peaks, plus a matched control, because each one alone is
// uninterpretable:
//
//   peak-auto     representative = highest SpliceAI in the peak.  CIRCULAR: SpliceAI and
//                 AlphaGenome agree at r = 0.72, so picking on one and scoring the other mostly
//                 re-measures that. 44% of these carry a SpliceAI call.
//   peak-centre   representative = nearest the peak midpoint, no score consulted. UNBIASED.
//   peak-control  one variant per peak from the SAME gene and AF band, deep intronic, but from a
//                 bin BELOW the density threshold. DENSITY-MATCHED: everything held constant
//                 except the thing the filter selects on. 2% carry a SpliceAI call.
//
// The only test that means anything is peak-centre vs peak-control. Comparing peak-auto to anything
// tests the SpliceAI correlation, and comparing either to the cb++ null arm tests AF-matching only,
// which is the mistake this exists to avoid repeating.
//
//   node scripts/analyse_peaks.mjs
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const J = f => JSON.parse(fs.readFileSync(R + f, 'utf8'));
const has = f => fs.existsSync(R + f);

const CB = J('data/cbpp.json');
const T = CB.threshold;
const auto = J('data/panels_peak.json').drawn;
const centre = has('data/panels_peak_centre.json') ? J('data/panels_peak_centre.json').drawn : [];
const ctl = has('data/panels_peakctl.json') ? J('data/panels_peakctl.json').drawn : [];
const CTLSET = has('data/peak_control_set.json') ? J('data/peak_control_set.json').rows : [];

const hit = r => r.splice_maxabsdiff > T;
const pct = (a, b) => b ? (100 * a / b).toFixed(1) + '%' : '—';
const med = a => { if (!a.length) return null; const s = a.slice().sort((x, y) => x - y); return s[s.length >> 1]; };

function fisher(a, b, c, d) {
  const lg = z => { const g = [676.5203681218851,-1259.1392167224028,771.32342877765313,-176.61502916214059,
    12.507343278686905,-0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lg(1 - z);
    z -= 1; let x = 0.99999999999980993; for (let i = 0; i < g.length; i++) x += g[i] / (z + i + 1);
    const t = z + g.length - 0.5; return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x); };
  const lf = n => lg(n + 1);
  const pr = (a,b,c,d) => Math.exp(lf(a+b)+lf(c+d)+lf(a+c)+lf(b+d)-lf(a)-lf(b)-lf(c)-lf(d)-lf(a+b+c+d));
  let s = 0; for (let i = a; i <= Math.min(a+b, a+c); i++) s += pr(i, a+b-i, a+c-i, d-(i-a));
  return Math.min(1, s);
}

console.log(`threshold ${T.toFixed(4)} — the 99th percentile of the cb++ screen's 600-variant matched null.\n`);
console.log('set             n    median AG   above threshold');
for (const [lab, rows] of [['peak-auto', auto], ['peak-centre', centre], ['peak-control', ctl]]) {
  if (!rows.length) { console.log(`  ${lab.padEnd(14)} — not run yet`); continue; }
  const v = rows.map(r => r.splice_maxabsdiff);
  console.log(`  ${lab.padEnd(14)}${String(rows.length).padStart(4)}   ${med(v).toFixed(4).padStart(9)}   ` +
    `${String(rows.filter(hit).length).padStart(3)}   ${pct(rows.filter(hit).length, rows.length)}`);
}

if (centre.length && ctl.length) {
  const ch = centre.filter(hit).length, kh = ctl.filter(hit).length;
  const p = fisher(ch, centre.length - ch, kh, ctl.length - kh);
  console.log('\nTHE TEST — unbiased pick in a dense bin vs unbiased pick in a matched sparse bin:');
  console.log(`  peak-centre  ${ch}/${centre.length} (${pct(ch, centre.length)})`);
  console.log(`  peak-control ${kh}/${ctl.length} (${pct(kh, ctl.length)})`);
  console.log(`  rate ratio ${kh ? (ch / centre.length / (kh / ctl.length)).toFixed(2) + 'x' : 'n/a'}   Fisher one-sided p = ${p < 1e-4 ? p.toExponential(1) : p.toFixed(4)}`);
  console.log(`  -> ${p < 0.01 ? 'DENSITY PREDICTS AN EFFECT' : 'no evidence that density predicts an effect'}`);
}

// internal check: peaks the statistic found vs peaks the per-gene floor supplied. If the floor
// group matches or beats the threshold group, the threshold is not selecting for anything.
if (centre.length) {
  console.log('\ninternal check — did clearing the density threshold matter at all?');
  for (const [lab, rows] of [['peak-auto', auto], ['peak-centre', centre]]) {
    const t = rows.filter(r => r.picked_via === 'threshold'), f = rows.filter(r => r.picked_via === 'gene-floor');
    if (!t.length || !f.length) continue;
    const th = t.filter(hit).length, fh = f.filter(hit).length;
    console.log(`  ${lab.padEnd(12)} threshold-found ${th}/${t.length} (${pct(th, t.length)})   ` +
      `gene-floor ${fh}/${f.length} (${pct(fh, f.length)})   p = ${fisher(th, t.length-th, fh, f.length-fh).toFixed(3)}`);
  }
  console.log('  A gene-floor peak did NOT clear the threshold. If it scores as well, the threshold is inert.');
}

if (ctl.length && CTLSET.length) {
  const withSA = r => (CTLSET.find(c => c.pos === r.pos && c.gene === r.gene) || {}).spliceai;
  const sa = ctl.filter(r => { const s = withSA(r); return s != null && s >= 0.05; });
  console.log(`\ncontrols carrying a SpliceAI call: ${sa.length}/${ctl.length} (${pct(sa.length, ctl.length)})` +
    `  — against 44% of peak-auto, which is why peak-auto cannot be the comparison.`);
}
