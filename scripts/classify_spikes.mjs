// Classify every density outlier as calling noise or real context-driven mutability.
//
// A tall bin means one of two very different things and the chart cannot tell them apart:
//
//   ARTIFACT   the bin sits on a homopolymer, a short tandem repeat, or a segmental duplication.
//              Sequencers slip in repeats and aligners misplace reads in duplicated sequence, so
//              spurious calls pile up. OAS1's tallest bin is 18 consecutive T.
//   REAL       the bin is GC- or CpG-elevated with no repeat structure. CpG dinucleotides mutate
//              roughly 10x faster than the genome average because methylated cytosine deaminates
//              to thymine, so a genuine density spike over a CpG island is expected biology.
//
// Repeat structure is computed from the sequence itself; segmental duplications come from UCSC's
// genomicSuperDups track, which cannot be derived from sequence alone.
//
//   node scripts/classify_spikes.mjs
import fs from 'node:fs';
const R = 'C:/DevLab/GitFolder/daily/13b-genebrowser/';
const NB = 220, MULT = 4;                       // same binning and threshold the strip uses
const sleep = ms => new Promise(r => setTimeout(r, ms));

async function get(u, tries = 4) {
  for (let i = 0; i < tries; i++) {
    try { const r = await fetch(u, { headers: { accept: 'application/json' } }); if (r.ok) return r.json(); }
    catch (e) { /* retry */ }
    await sleep(900 * (i + 1));
  }
  return null;
}

// longest run of one base
const homopolymer = s => {
  let best = 0, cur = 1;
  for (let i = 1; i < s.length; i++) { if (s[i] === s[i - 1]) cur++; else cur = 1; if (cur > best) best = cur; }
  return best;
};
// short tandem repeat: a unit of 1-6 bp repeated at least 5 times consecutively
const str = s => {
  let best = { unit: '', copies: 0 };
  for (let k = 1; k <= 6; k++) {
    for (let i = 0; i + k * 2 <= s.length; i++) {
      const u = s.substr(i, k);
      let c = 1, j = i + k;
      while (s.substr(j, k) === u) { c++; j += k; }
      if (c >= 5 && c * k > best.copies * best.unit.length) best = { unit: u, copies: c };
    }
  }
  return best;
};
const gc = s => (s.match(/[GC]/g) || []).length / s.length;
// CpG observed/expected. Above ~0.6 with GC>50% is the standard CpG-island signature, and CpG
// sites mutate ~10x faster than average, so density there is biology not noise.
const cpgOE = s => {
  const c = (s.match(/C/g) || []).length, g = (s.match(/G/g) || []).length;
  const cg = (s.match(/CG/g) || []).length;
  return (c && g) ? (cg * s.length) / (c * g) : 0;
};

const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const out = {};
for (const gene of Object.keys(GENES)) {
  const G = GENES[gene];
  const p = JSON.parse(fs.readFileSync(R + `data/variants_${gene}.json`, 'utf8'));
  const i = Object.fromEntries(p.cols.map((c, k) => [c, k]));
  const span = G.end - G.start;
  const bins = new Array(NB).fill(0);
  for (const r of p.rows) { const k = Math.floor((r[i.pos] - G.start) / span * NB); if (k >= 0 && k < NB) bins[k]++; }
  const mean = p.rows.length / NB;
  const hot = bins.map((n, k) => ({ k, n })).filter(x => x.n > MULT * mean);
  if (!hot.length) { console.log(`${gene}: no outlier bins`); continue; }

  // segmental duplications for the whole gene, one request
  const segs = await get(`https://api.genome.ucsc.edu/getData/track?genome=hg38;track=genomicSuperDups;` +
    `chrom=${G.chrom};start=${G.start};end=${G.end}`);
  const segList = (segs && (segs.genomicSuperDups || segs[G.chrom] || [])) || [];
  await sleep(400);

  out[gene] = [];
  for (const h of hot) {
    const b0 = Math.round(G.start + h.k / NB * span), b1 = Math.round(G.start + (h.k + 1) / NB * span);
    const sq = await get(`https://rest.ensembl.org/sequence/region/human/` +
      `${String(G.chrom).replace('chr', '')}:${b0}..${b1}?content-type=application/json`);
    await sleep(300);
    const seq = (sq && sq.seq || '').toUpperCase();
    if (!seq) { out[gene].push({ bin: h.k, start: b0, end: b1, n: h.n, verdict: 'unknown', why: 'no sequence' }); continue; }

    const hp = homopolymer(seq), st = str(seq), G_ = gc(seq), oe = cpgOE(seq);
    const inSeg = segList.some(s => s.chromStart < b1 && s.chromEnd > b0);

    let verdict, why;
    if (inSeg) { verdict = 'artifact'; why = 'inside a segmental duplication — reads misplace between copies'; }
    else if (hp >= 8) { verdict = 'artifact'; why = `${hp}-base homopolymer — sequencers slip in these`; }
    else if (st.copies >= 5 && st.unit.length >= 2) { verdict = 'artifact'; why = `short tandem repeat (${st.unit})x${st.copies}`; }
    else if (oe >= 0.6 && G_ >= 0.5) { verdict = 'real'; why = `CpG island signature (GC ${(G_*100).toFixed(0)}%, CpG o/e ${oe.toFixed(2)}) — CpG mutates ~10x faster`; }
    else if (G_ >= 0.6) { verdict = 'real'; why = `GC-rich (${(G_*100).toFixed(0)}%) with no repeat structure`; }
    else { verdict = 'unexplained'; why = `no repeat, GC ${(G_*100).toFixed(0)}%, CpG o/e ${oe.toFixed(2)}`; }

    out[gene].push({ bin: h.k, start: b0, end: b1, n: h.n, fold: +(h.n / mean).toFixed(1),
                     homopolymer: hp, str: st.copies >= 5 ? `${st.unit}x${st.copies}` : null,
                     gc: +(G_ * 100).toFixed(0), cpg_oe: +oe.toFixed(2), segdup: inSeg, verdict, why });
    console.log(`  ${gene} bin ${h.k}  ${b0}-${b1}  ${h.n} variants (${(h.n/mean).toFixed(1)}x)  ${verdict.toUpperCase()}  ${why}`);
  }
}
fs.writeFileSync(R + 'data/spike_classification.json', JSON.stringify(out, null, 1));
const all = Object.values(out).flat();
const c = v => all.filter(x => x.verdict === v).length;
console.log(`\n${all.length} outlier bins: ${c('artifact')} artifact, ${c('real')} real, ${c('unexplained')} unexplained, ${c('unknown')} unknown`);
console.log('wrote data/spike_classification.json');
