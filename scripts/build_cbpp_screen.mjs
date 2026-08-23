// Build the 2,000-variant cb++ AlphaGenome screening set, in three arms.
//
// The point of the arms is that the result has to survive "your ML is just finding noise". One
// screening set with no internal comparison cannot answer that. Three can:
//
//   A  abundance   700   AF >= 1%, gate-passers, ranked by AF x gene weight.
//                        The variants a reviewer would call human-relevant. Low prior of function.
//   B  hint        700   phyloP>=1 OR CADD>=10 OR SpliceAI>=0.05, ranked by AF x gene weight.
//                        Higher prior, but rarer -- the abundance/function trade-off, made explicit.
//   C  null        600   NO functional hint, drawn at random, AF-MATCHED to the pooled A+B.
//
// Why C is AF-matched rather than uniform: a uniform draw from the gated pool has median AF ~2.6e-5
// while arm A is >=1%, so "the model flags more of A than C" would be an allele-frequency effect and
// would be read as a functional one. Matching on AF makes the comparison about the HINT, which is
// the thing actually under test. This is the same reasoning behind data/control_set.json, whose
// AF-matched negatives are what gave the panels pipeline a real 0% overlap instead of a flattering one.
//
// Arms are disjoint, assignment order A then B then C, so no variant is its own control.
//
// A per-gene cap keeps this from becoming a VWF screen: an unconstrained abundance ranking put 737
// of 2,000 in VWF purely because VWF is 176 kb. A finding that only replicates in the biggest gene
// is not a finding about coagulation.
//
//   node scripts/build_cbpp_screen.mjs
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const WEIGHT = JSON.parse(fs.readFileSync(R + 'data/gene_weights.json', 'utf8'));

const N_A = 700, N_B = 700, N_C = 600;
const GENE_CAP = 0.15;                 // no gene may exceed this share of an arm
const SEED = 20260818;                 // fixed: the draw must be reproducible

// Math.random() cannot be reproduced, and an unreproducible control arm is not a control.
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const rnd = mulberry32(SEED);

const runAt = (s, i) => { const c = s[i]; if (!c) return 0; let a = i, b = i;
  while (a > 0 && s[a-1] === c) a--; while (b < s.length-1 && s[b+1] === c) b++; return b-a+1; };
function strAround(s, i, half = 24) {
  const lo = Math.max(0, i - half), w = s.slice(lo, i + half); let best = 0;
  for (let u = 2; u <= 6; u++) for (let k = 0; k + u*4 <= w.length; k++) {
    const m = w.substr(k, u); if (/^(.)\1*$/.test(m)) continue;
    let n = 1; while (w.substr(k + n*u, u) === m) n++;
    if (n >= 4) { const st = lo+k, en = st + n*u; if (i >= st && i < en && n*u > best) best = n*u; }
  }
  return best;
}

/* ---- gate every variant, keeping the survivors ---------------------------------------------- */
const pool = [];
for (const [sym, g] of Object.entries(GENES)) {
  const gw = WEIGHT[sym]; if (!gw) continue;
  const p = JSON.parse(fs.readFileSync(R + `data/variants_${sym}.json`, 'utf8'));
  const ix = n => p.cols.indexOf(n);
  const [Irs,Ipos,Iref,Ialt,Iaf,Isa,Ipg,Icadd,Iphy,Icod] =
    ['rsid','pos','ref','alt','af','spliceai','pangolin','cadd','phylop','coding'].map(ix);
  const seq = fs.existsSync(R + `data/seq_${sym}.txt`) ? fs.readFileSync(R + `data/seq_${sym}.txt`,'utf8').trim() : null;
  const ex = g.exons.map(e => [Math.min(...e), Math.max(...e)]).sort((a,b)=>a[0]-b[0]);
  const bounds = []; for (let i=0;i<ex.length-1;i++) bounds.push(ex[i][1], ex[i+1][0]);
  const perPos = {}; for (const r of p.rows) perPos[r[Ipos]] = (perPos[r[Ipos]]||0)+1;

  for (const r of p.rows) {
    const af = r[Iaf], pos = r[Ipos], ref = String(r[Iref]), alt = String(r[Ialt]);
    if (af == null || af < 1e-5) continue;
    // AF >= 0.5 means the REFERENCE is the minor allele. Calling that "a variant whose effect we
    // predict" is a statement about the reference assembly, not about a person.
    if (af >= 0.5) continue;
    if (ref.length !== 1 || alt.length !== 1 || !'ACGT'.includes(ref) || !'ACGT'.includes(alt)) continue;
    if (perPos[pos] > 3 || r[Icod]) continue;
    let intron = -1;
    for (let i=0;i<ex.length-1;i++) if (ex[i][1] < pos && pos < ex[i+1][0]) { intron = i; break; }
    if (intron < 0) continue;
    const dist = Math.min(...bounds.map(b => Math.abs(pos - b)));
    if (dist < 100) continue;
    if (seq) { const i = pos - g.start;
      if (i >= 0 && i < seq.length) { if (runAt(seq,i) >= 8) continue; if (strAround(seq,i) >= 12) continue; } }

    const phylop = r[Iphy], cadd = r[Icadd], spliceai = r[Isa];
    const hint = (phylop != null && phylop >= 1) || (cadd != null && cadd >= 10) || (spliceai != null && spliceai >= 0.05);
    pool.push({ gene: sym, rsid: r[Irs], pos, ref, alt, af, phylop, cadd, spliceai,
                pangolin: r[Ipg], dist_splice: dist, intron: (intron+1)+'/'+(ex.length-1),
                gene_weight: gw, hint });
  }
}
console.log(`gated pool: ${pool.length.toLocaleString()}   with a functional hint: ${pool.filter(v=>v.hint).length.toLocaleString()}`);

/* ---- arm selection --------------------------------------------------------------------------- */
const taken = new Set();
const key = v => v.gene + ':' + v.pos + ':' + v.ref + v.alt;

function pick(cands, n, label) {
  const cap = Math.ceil(n * GENE_CAP), perGene = {}, out = [];
  for (const v of cands) {
    if (out.length >= n) break;
    if (taken.has(key(v))) continue;
    if ((perGene[v.gene] || 0) >= cap) continue;
    perGene[v.gene] = (perGene[v.gene] || 0) + 1;
    taken.add(key(v)); out.push({ ...v, arm: label });
  }
  return out;
}

// A — abundance first
const armA = pick(pool.filter(v => v.af >= 0.01).sort((a,b) => (b.af*b.gene_weight) - (a.af*a.gene_weight)), N_A, 'abundance');
// B — functional hint, then abundance
const armB = pick(pool.filter(v => v.hint).sort((a,b) => (b.af*b.gene_weight) - (a.af*a.gene_weight)), N_B, 'hint');

// C — AF-matched null, no hint.
//
// MATCH TO EACH ARM SEPARATELY, NOT TO THEIR POOL. The first version matched the pooled A+B, and
// that pool is BIMODAL: arm A has median AF 32%, arm B has 0.6%. A distribution matched to the
// mixture sat at 18% and therefore matched neither -- KS distance 0.333 against the pool it was
// supposed to mirror, 30x too common for arm B. The comparison it was built to support was
// confounded by the very variable it was meant to control.
//
// Matching each arm and concatenating gives a null whose AF distribution is the same MIXTURE, so
// every band contains controls drawn at that band's frequency. Band-stratified analysis then
// compares like with like, which pooled matching cannot do when the pool is multi-modal.
const bin = af => Math.round(Math.log10(af) * 2) / 2;          // half-decade bins
// how many controls each arm needs, per band, kept SEPARATE
const wantPerBin = {};
for (const [arm, n] of [[armA, N_C * armA.length / (armA.length + armB.length)],
                        [armB, N_C * armB.length / (armA.length + armB.length)]]) {
  const per = {};
  for (const v of arm) per[bin(v.af)] = (per[bin(v.af)] || 0) + 1;
  for (const k of Object.keys(per)) wantPerBin[k] = (wantPerBin[k] || 0) + per[k] * (n / arm.length);
}
const noHintByBin = {};
for (const v of pool) if (!v.hint && !taken.has(key(v))) (noHintByBin[bin(v.af)] ||= []).push(v);

const armC = [];
const shortfall = {};
for (const b of Object.keys(wantPerBin)) {
  const want = Math.round(wantPerBin[b]);
  const avail = noHintByBin[b] || [];
  // Fisher-Yates with the seeded generator, so the same seed gives the same control arm
  for (let i = avail.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [avail[i], avail[j]] = [avail[j], avail[i]]; }
  const got = avail.slice(0, want);
  if (got.length < want) shortfall[b] = want - got.length;
  for (const v of got) { taken.add(key(v)); armC.push({ ...v, arm: 'null' }); }
}
// A shortfall in a bin means the pool has no unhinted variants at that frequency -- itself a
// finding, and it must be reported rather than quietly back-filled from a different bin.
if (Object.keys(shortfall).length)
  console.log('  AF bins where the null arm could not be filled: ' +
    Object.entries(shortfall).map(([b,n]) => `10^${b} short ${n}`).join(', '));

const rows = armA.concat(armB, armC);
console.log(`\narm A abundance : ${armA.length}   median AF ${median(armA.map(v=>v.af)).toExponential(2)}   ${new Set(armA.map(v=>v.gene)).size} genes`);
console.log(`arm B hint      : ${armB.length}   median AF ${median(armB.map(v=>v.af)).toExponential(2)}   ${new Set(armB.map(v=>v.gene)).size} genes   ` +
  `phyloP>=1 ${armB.filter(v=>(v.phylop??-9)>=1).length}, CADD>=10 ${armB.filter(v=>(v.cadd??-9)>=10).length}, SpliceAI>=0.05 ${armB.filter(v=>(v.spliceai??-1)>=0.05).length}`);
console.log(`arm C null      : ${armC.length}   median AF ${median(armC.map(v=>v.af)).toExponential(2)}   ${new Set(armC.map(v=>v.gene)).size} genes`);
console.log(`\ntotal ${rows.length}   all disjoint: ${new Set(rows.map(key)).size === rows.length}`);

function median(a){ const s=a.slice().sort((x,y)=>x-y); return s[Math.floor(s.length/2)]; }

fs.writeFileSync(R + 'data/cbpp_screen.json', JSON.stringify({
  generated_from: 'scripts/build_cbpp_screen.mjs', seed: SEED,
  design: { A: 'abundance: AF>=1%, ranked AF x gene weight',
            B: 'hint: phyloP>=1 OR CADD>=10 OR SpliceAI>=0.05, ranked AF x gene weight',
            C: 'null: no hint, random, AF-matched to pooled A+B' },
  gene_cap: GENE_CAP, pool_size: pool.length, rows
}, null, 1));
console.log(`wrote data/cbpp_screen.json  (${(fs.statSync(R+'data/cbpp_screen.json').size/1024).toFixed(0)} KB)`);
