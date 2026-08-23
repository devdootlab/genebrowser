// Read the cb++ screen and answer the only question that matters about it:
// does AlphaGenome flag deep-intronic variants for any reason OTHER than SpliceAI already flagging them?
//
// TWO THINGS THIS SCRIPT LEARNED THE HARD WAY, both of which inverted the headline result.
//
// 1. THE THRESHOLD WAS UNDER-POWERED. It was originally 0.1313, the maximum over the 30 AF-matched
//    negatives in data/panels_control.json. Against the 600-variant null arm here, 0.1313 turns out
//    to be the 99.0th percentile -- so "about 1% of every arm exceeds it" is the DEFINITION of that
//    threshold, not a discovery. A max-of-30 cannot set a cut-off for 2,000 tests. The null arm is
//    now the control distribution: 20x the size, AF-matched by construction, same pipeline.
//
// 2. THE HINT ARM'S EFFECT WAS ENTIRELY SpliceAI. Arm B admits a variant on ANY of phyloP>=1,
//    CADD>=10, or SpliceAI>=0.05. Reported whole, it looked spectacular -- 7.4% vs 1.0%, a 7.4x
//    rate ratio at p = 2.4e-9. Split by which hint let each variant in:
//        SpliceAI-entered   44/114  38.6%
//        no SpliceAI hint    8/586   1.4%
//        null arm            6/600   1.0%
//    All of it was AlphaGenome agreeing with SpliceAI, which is the r = 0.72 correlation this repo
//    already measured, re-expressed through a selection criterion. phyloP and CADD contribute
//    nothing detectable. Any arm whose entry criterion overlaps the readout must be split this way
//    before it is believed.
//
//   node scripts/analyse_cbpp.mjs
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const P = JSON.parse(fs.readFileSync(R + 'data/panels_cbpp.json', 'utf8'));
const C = JSON.parse(fs.readFileSync(R + 'data/panels_control.json', 'utf8'));
const posMin = Math.min(...C.drawn.filter(d => d.label === 'positive').map(d => d.splice_maxabsdiff));

const rows = P.drawn;
const by = a => rows.filter(r => r.arm === a);
const hasSA = r => (r.spliceai_gnomad ?? -1) >= 0.05;
const nulv = by('null').map(r => r.splice_maxabsdiff).sort((a, b) => a - b);
const q = p => nulv[Math.min(nulv.length - 1, Math.floor(nulv.length * p))];
const THRESH = q(0.99);
const hit = r => r.splice_maxabsdiff > THRESH;
const pct = (a, b) => b ? (100 * a / b).toFixed(1) + '%' : '—';

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

console.log(`cb++ screen: ${rows.length} scored (${P.refused.length} refused, ${P.failed.length} failed)\n`);
console.log(`null arm, n=${nulv.length}, AF-matched, no functional hint:`);
console.log(`  median ${q(.5).toFixed(5)}  90th ${q(.90).toFixed(4)}  95th ${q(.95).toFixed(4)}  99th ${q(.99).toFixed(4)}  max ${nulv[nulv.length-1].toFixed(4)}`);
console.log(`  THRESHOLD = ${THRESH.toFixed(4)} (99th percentile of this null).`);
console.log(`  a canonical splice-killer scores ${posMin.toFixed(2)}+, so the whole scale here is far below "destroys a splice site".\n`);

console.log('arm          n   median AG   flagged    rate');
for (const a of ['abundance', 'hint', 'null']) {
  const s = by(a), h = s.filter(hit).length;
  const v = s.map(r => r.splice_maxabsdiff).sort((x, y) => x - y);
  console.log(`  ${a.padEnd(11)}${String(s.length).padStart(4)}   ${v[v.length>>1].toFixed(4).padStart(9)}   ${String(h).padStart(7)}   ${pct(h, s.length).padStart(6)}`);
}

/* ---- the split that decides whether any of this is new -------------------------------------- */
const saHint = by('hint').filter(hasSA), freeHint = by('hint').filter(r => !hasSA(r)), nul = by('null');
const a1 = saHint.filter(hit).length, b1 = freeHint.filter(hit).length, c1 = nul.filter(hit).length;
console.log('\nWHICH HINT DID THE WORK — arm B split by entry criterion');
console.log(`  SpliceAI >= 0.05 entered   ${String(a1).padStart(4)}/${String(saHint.length).padEnd(5)} ${pct(a1, saHint.length).padStart(7)}`);
console.log(`  phyloP / CADD only         ${String(b1).padStart(4)}/${String(freeHint.length).padEnd(5)} ${pct(b1, freeHint.length).padStart(7)}`);
console.log(`  null arm                   ${String(c1).padStart(4)}/${String(nul.length).padEnd(5)} ${pct(c1, nul.length).padStart(7)}`);
const pFree = fisher(b1, freeHint.length - b1, c1, nul.length - c1);
console.log(`\n  SpliceAI-free hint vs null: ${pct(b1, freeHint.length)} vs ${pct(c1, nul.length)}   ` +
            `Fisher one-sided p = ${pFree.toFixed(3)}  ->  ${pFree < 0.01 ? 'INDEPENDENT SIGNAL' : 'NO INDEPENDENT SIGNAL'}`);
console.log('  Read: conservation and CADD do not predict an AlphaGenome flag once SpliceAI is removed.');
console.log('  The model is, on this data, largely restating SpliceAI.');

/* ---- what survives as a candidate list --------------------------------------------------------- */
// cb++ is redefined here as the CONCORDANT set: two independent splice models agree, the variant is
// deep intronic, real, and in a coagulation gene. That satisfies "not artefact", "not the obvious
// donor/acceptor", "human-relevant gene", and gives the functional criterion the best evidence
// available -- but it is CONFIRMATION by ML, not DISCOVERY by ML, and must not be sold as the latter.
const cb = rows.filter(r => hit(r) && hasSA(r)).sort((a, b) => (b.af || 0) - (a.af || 0));
const novel = rows.filter(r => hit(r) && !hasSA(r)).sort((a, b) => b.splice_maxabsdiff - a.splice_maxabsdiff);
console.log(`\n=== cb++ (concordant): ${cb.length} variants — SpliceAI >= 0.05 AND AlphaGenome > ${THRESH.toFixed(4)} ===`);
console.log('gene       variant                AF      dist   phyloP  SpliceAI  AlphaGenome');
for (const r of cb.slice(0, 25))
  console.log(`${r.gene.padEnd(10)} ${String(r.rsid || r.pos).padEnd(20)} ${(r.af*100).toFixed(2).padStart(6)}% ${String(r.dist_splice).padStart(6)}  ` +
    `${String(r.phylop ?? '—').padStart(6)}  ${String(r.spliceai_gnomad ?? '—').padStart(7)}  ${r.splice_maxabsdiff.toFixed(4).padStart(11)}`);
if (cb.length > 25) console.log(`  ... and ${cb.length - 25} more`);
console.log(`\n${novel.length} further variants clear the threshold WITHOUT a SpliceAI call, but that count is`);
console.log(`consistent with the null rate (${pct(c1, nul.length)}), so they are not distinguishable from background.`);

fs.writeFileSync(R + 'data/cbpp.json', JSON.stringify({
  generated_from: 'scripts/analyse_cbpp.mjs',
  threshold: THRESH,
  threshold_source: '99th percentile of the 600-variant AF-matched null arm of this same screen',
  threshold_note: 'the previous 0.1313 was the max of only 30 controls, which is the 99.0th percentile here — a ~1% hit rate against it was arithmetic, not signal',
  positive_min: posMin, screened: rows.length,
  arm_counts: Object.fromEntries(['abundance','hint','null'].map(a => [a, by(a).length])),
  arm_hits: Object.fromEntries(['abundance','hint','null'].map(a => [a, by(a).filter(hit).length])),
  spliceai_split: { entered: saHint.length, entered_hits: a1, free: freeHint.length, free_hits: b1,
                    null_n: nul.length, null_hits: c1, free_vs_null_p: pFree },
  independent_signal: pFree < 0.01,
  rows: cb, unconfirmed: novel
}));
console.log(`\nwrote data/cbpp.json  (${cb.length} concordant, ${novel.length} unconfirmed)`);
