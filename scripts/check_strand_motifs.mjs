// Known-answer test for the strand-aware splice annotation in browse.html.
//
// ~99% of human introns begin GT and end AG *in transcript orientation*. That is a fact about
// biology, not about this repo, which makes it a real known-answer set: run the page's own
// dinucleotide logic over every internal exon boundary of every gene and it must hit ~99%.
//
// The bug this catches: browse.html searched the PLUS strand for GT and AG regardless of strand.
// On a minus-strand gene the donor GT reads AC on the plus strand and the acceptor AG reads CT, so
// the page marked non-sites and drew nothing at the real ones. Run both ways below -- the old
// plus-strand-only rule should score near chance on minus-strand genes, the fix near 100%.
//
//   node scripts/check_strand_motifs.mjs
import fs from 'node:fs';

import path from 'node:path';
import { fileURLToPath } from 'node:url';
// repo root, derived from this file's location -- no absolute path baked in
const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const sleep = ms => new Promise(r => setTimeout(r, ms));
const CO = { A: 'T', C: 'G', G: 'C', T: 'A' };
const txBase = c => CO[c] || 'N';
// exactly the rule now in browse.html
const txDi = (s, i, strand) => strand === -1 ? txBase(s[i + 1]) + txBase(s[i]) : s[i] + s[i + 1];

async function seq(chrom, a, b) {
  for (let t = 0; t < 4; t++) {
    try {
      const r = await fetch(`https://rest.ensembl.org/sequence/region/human/${String(chrom).replace('chr','')}:${a}..${b}?content-type=application/json`,
        { headers: { accept: 'application/json' } });
      if (r.ok) return (await r.json()).seq.toUpperCase();
    } catch (e) { /* retry */ }
    await sleep(1000 * (t + 1));
  }
  return null;
}

let fixOK = 0, fixBad = 0, oldOK = 0, oldBad = 0, missed = 0;
const perGene = [];

for (const [sym, g] of Object.entries(GENES)) {
  const ex = (g.exons || []).map(e => [Math.min(...e), Math.max(...e)]).sort((a, b) => a[0] - b[0]);
  if (ex.length < 2) continue;
  // introns are the gaps between consecutive exons, in GENOMIC order
  const introns = [];
  for (let k = 0; k < ex.length - 1; k++) introns.push([ex[k][1] + 1, ex[k + 1][0] - 1]);
  let gf = 0, go = 0, n = 0, gm = 0;
  for (const [a, b] of introns) {
    if (b - a < 4) continue;
    // only the 2 bp at each END of the intron matters. Fetching whole introns pulled 30 kb
    // per F8 intron for 4 informative bases.
    const head = await seq(g.chrom, a, a + 1), tail = await seq(g.chrom, b - 1, b);
    // A dropped fetch used to `continue` silently, so the intron vanished from the DENOMINATOR and
    // the run still printed a confident percentage. Two runs of this script reported FGB as 7 and
    // then 5 introns when it has 7 -- same code, same gene, different silent losses. Count them.
    if (!head || !tail || head.length < 2 || tail.length < 2) { missed++; gm++; continue; }
    const s = head + tail;
    // In TRANSCRIPT order the intron starts with a donor and ends with an acceptor. On a minus-
    // strand gene the transcript start of the intron is its HIGH genomic end.
    const minus = g.strand === -1;
    const donorI    = minus ? 2 : 0;
    const acceptorI = minus ? 0 : 2;
    const dFix = txDi(s, donorI, g.strand), aFix = txDi(s, acceptorI, g.strand);
    // the rule the page used before: plus strand only, GT at the low end, AG at the high end
    const dOld = head, aOld = tail;
    const okFix = dFix === 'GT' && aFix === 'AG';
    const okOld = dOld === 'GT' && aOld === 'AG';
    okFix ? fixOK++ : fixBad++; okOld ? oldOK++ : oldBad++;
    if (okFix) gf++; if (okOld) go++;
    n++;
    await sleep(200);
  }
  perGene.push({ sym, strand: g.strand === -1 ? '-' : '+', n, fixed: gf, old: go, missed: gm, want: introns.filter(([x,y])=>y-x>=4).length });
  process.stdout.write(`\r  ${perGene.length}/${Object.keys(GENES).length} genes`);
}
console.log('\n');
console.log('gene        strand  introns  strand-aware(fixed)  plus-strand-only(old)  unfetched');
for (const p of perGene) {
  console.log(`${p.sym.padEnd(11)} ${p.strand}      ${String(p.n).padStart(6)}   ` +
    `${String(p.fixed + '/' + p.n).padStart(18)}   ${String(p.old + '/' + p.n).padStart(20)}   ` +
    (p.missed ? `${p.missed} of ${p.want} DROPPED` : ''));
}
const pc = (a, b) => (100 * a / (a + b)).toFixed(1) + '%';
console.log(`\nstrand-aware   : ${fixOK}/${fixOK + fixBad} canonical GT-AG  (${pc(fixOK, fixBad)})`);
console.log(`plus-strand-only: ${oldOK}/${oldOK + oldBad} canonical GT-AG  (${pc(oldOK, oldBad)})`);
if (missed) {
  console.log(`\n${missed} intron(s) could not be fetched from Ensembl and are NOT in either ` +
    `denominator. The percentages above describe ${fixOK + fixBad} introns, not every intron in ` +
    `these genes. Re-run to pick them up -- the failures are transient 500s, not missing data.`);
}
const passed = fixOK / (fixOK + fixBad) > 0.95;
// A run that dropped introns still reports, but does not pass. A silently incomplete green tick is
// the exact failure mode this whole session has been about.
console.log(passed && !missed ? '\nPASS — the strand-aware rule recovers the canonical motif.'
          : passed ? `\nINCOMPLETE — rule looks right (${pc(fixOK, fixBad)}) but ${missed} intron(s) went unchecked.`
          : '\nFAIL — the strand-aware rule is still wrong.');
process.exit(passed && !missed ? 0 : 1);
