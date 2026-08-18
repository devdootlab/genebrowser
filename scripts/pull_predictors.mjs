// Pull SpliceAI, Pangolin, CADD and phyloP from gnomAD for every variant, and merge into the
// per-gene files.
//
// These were available in the pull we already ran and nobody asked for them. gnomAD's variant
// record carries in_silico_predictors: spliceai_ds_max, pangolin_largest_ds, cadd, phylop. That
// makes the whole "how do we get SpliceAI" question a batching problem, not an infrastructure one.
//
//   node scripts/pull_predictors.mjs [GENE]
import fs from 'node:fs';

const R = 'C:/DevLab/GitFolder/daily/13b-genebrowser/';
const ONLY = process.argv[2] || null;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const sleep = ms => new Promise(r => setTimeout(r, ms));

// Region query returns every variant in one request with its predictors attached. Per-variant
// queries would be 37,655 requests; this is a few dozen.
const Q = (chrom, start, stop) => JSON.stringify({
  query: `{ region(chrom: "${chrom}", start: ${start}, stop: ${stop}, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) { variant_id pos in_silico_predictors { id value } } } }`
});

async function post(body, tries = 5) {
  for (let i = 0; i < tries; i++) {
    try {
      const r = await fetch('https://gnomad.broadinstitute.org/api',
        { method: 'POST', headers: { 'content-type': 'application/json' }, body });
      if (r.status === 429) { await sleep(20000); continue; }   // rate limited, not failed
      const j = await r.json();
      if (j.errors && !j.data) { await sleep(6000); continue; }
      return j;
    } catch (e) { await sleep(4000 * (i + 1)); }
  }
  return null;
}

for (const gene of Object.keys(GENES)) {
  if (ONLY && gene !== ONLY) continue;
  const g = GENES[gene];
  const chrom = String(g.chrom).replace('chr', '');
  const f = R + `data/variants_${gene}.json`;
  if (!fs.existsSync(f)) { console.log(`${gene}: no per-gene file, skipping`); continue; }
  const p = JSON.parse(fs.readFileSync(f, 'utf8'));
  const idx = Object.fromEntries(p.cols.map((c, i) => [c, i]));

  // chunk the gene so no single response is enormous
  const CH = 10000;
  const found = new Map();
  let reqs = 0;
  for (let s0 = g.start; s0 <= g.end; s0 += CH) {
    const s1 = Math.min(g.end, s0 + CH - 1);
    const j = await post(Q(chrom, s0, s1));
    reqs++;
    const vs = (j && j.data && j.data.region && j.data.region.variants) || [];
    for (const v of vs) {
      const pr = {};
      for (const x of (v.in_silico_predictors || [])) pr[x.id] = x.value;
      found.set(v.variant_id, pr);
    }
    process.stdout.write(`\r  ${gene}: ${reqs} requests, ${found.size} variants with predictors`);
    await sleep(6500);                       // ~9 requests/min, inside gnomAD's limit
  }
  console.log('');

  // merge. variant_id is CHROM-POS-REF-ALT, which is exactly what we can rebuild.
  const add = ['spliceai', 'pangolin', 'cadd', 'phylop'];
  for (const a of add) if (!p.cols.includes(a)) p.cols.push(a);
  let hit = 0;
  p.rows = p.rows.map(r => {
    const key = `${chrom}-${r[idx.pos]}-${r[idx.ref]}-${r[idx.alt]}`;
    const pr = found.get(key) || {};
    const vals = [
      pr.spliceai_ds_max == null ? null : +pr.spliceai_ds_max,
      pr.pangolin_largest_ds == null ? null : +pr.pangolin_largest_ds,
      pr.cadd == null ? null : +pr.cadd,
      pr.phylop == null ? null : +pr.phylop
    ];
    if (vals[0] != null) hit++;
    const base = r.slice(0, p.cols.length - add.length);
    return base.concat(vals);
  });
  fs.writeFileSync(f, JSON.stringify(p));
  const sp = p.rows.map(r => r[p.cols.indexOf('spliceai')]).filter(x => x != null);
  console.log(`  ${gene}: ${hit} of ${p.rows.length} matched a SpliceAI score` +
    (sp.length ? `  |  >=0.2: ${sp.filter(x => x >= 0.2).length}  >=0.5: ${sp.filter(x => x >= 0.5).length}  >=0.8: ${sp.filter(x => x >= 0.8).length}` : ''));
}
console.log('\ndone. browse.html can now rank on spliceai instead of distance-to-splice-site.');
