// Add genes end-to-end: Ensembl structure -> gnomAD variants + predictors -> per-gene file.
//
// One script for the whole chain, because the previous state of this repo was a pile of scripts
// whose order lived only in someone's head, and that is how F2 ended up queried with SSTR2's
// Ensembl ID. Here the ID is looked up FROM THE SYMBOL and asserted against the coordinates that
// come back, so a wrong ID cannot survive step one.
//
//   node scripts/add_genes.mjs F8 VWF FGA FGB FGG
import fs from 'node:fs';

const R = 'C:/DevLab/GitFolder/daily/13b-genebrowser/';
const SYMBOLS = process.argv.slice(2);
if (!SYMBOLS.length) { console.error('give gene symbols'); process.exit(1); }
const sleep = ms => new Promise(r => setTimeout(r, ms));

async function ens(path, tries = 4) {
  for (let i = 0; i < tries; i++) {
    try {
      const r = await fetch('https://rest.ensembl.org' + path, { headers: { accept: 'application/json' } });
      if (r.ok) return r.json();
    } catch (e) { /* retry */ }
    await sleep(1200 * (i + 1));
  }
  return null;
}
async function gql(body, tries = 5) {
  for (let i = 0; i < tries; i++) {
    try {
      const r = await fetch('https://gnomad.broadinstitute.org/api',
        { method: 'POST', headers: { 'content-type': 'application/json' }, body });
      if (r.status === 429) { await sleep(20000); continue; }
      const j = await r.json();
      if (j.data) return j;
    } catch (e) { /* retry */ }
    await sleep(5000 * (i + 1));
  }
  return null;
}

const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));

for (const sym of SYMBOLS) {
  console.log(`\n=== ${sym} ===`);

  // 1. structure, from the SYMBOL. Never hand-enter an Ensembl ID.
  const g = await ens(`/lookup/symbol/homo_sapiens/${sym}?expand=1;content-type=application/json`);
  if (!g || !g.id) { console.log(`  could not resolve ${sym} in Ensembl — skipping`); continue; }
  const canon = (g.Transcript || []).find(t => t.is_canonical) || (g.Transcript || [])[0];
  if (!canon || !(canon.Exon || []).length) { console.log(`  no canonical transcript with exons — skipping`); continue; }
  const exons = canon.Exon.map(e => [e.start, e.end]).sort((a, b) => a[0] - b[0]);

  // the assertion that would have caught the SSTR2/F2 substitution immediately
  if (!exons.every(([a, b]) => a >= g.start && b <= g.end)) { console.log('  exons fall outside the gene span — refusing'); continue; }

  const rec = {
    symbol: sym, ensembl_id: g.id, canonical_transcript: canon.id,
    chrom: 'chr' + g.seq_region_name, strand: g.strand, start: g.start, end: g.end,
    length: g.end - g.start + 1, exons, exon_count: exons.length,
    description: (g.description || '').replace(/\s*\[Source.*$/, '')
  };
  console.log(`  ${g.id}  ${rec.chrom}:${g.start.toLocaleString()}-${g.end.toLocaleString()}  ` +
    `${(rec.length / 1000).toFixed(0)} kb  strand ${g.strand === 1 ? '+' : '-'}  ${exons.length} exons`);

  // 2. variants + predictors, chunked so no single response is enormous
  const chrom = String(g.seq_region_name);
  const CH = 10000;
  const rows = [];
  let reqs = 0;
  for (let s0 = g.start; s0 <= g.end; s0 += CH) {
    const s1 = Math.min(g.end, s0 + CH - 1);
    const j = await gql(JSON.stringify({
      query: `{ region(chrom: "${chrom}", start: ${s0}, stop: ${s1}, reference_genome: GRCh38) {
        variants(dataset: gnomad_r4) { variant_id pos ref alt rsids consequence
          genome { af } in_silico_predictors { id value } } } }` }));
    reqs++;
    for (const v of ((j && j.data && j.data.region && j.data.region.variants) || [])) {
      const pr = {};
      for (const x of (v.in_silico_predictors || [])) pr[x.id] = x.value;
      rows.push([
        (v.rsids && v.rsids[0]) || String(v.pos),
        v.pos, v.ref, v.alt,
        (v.consequence || '').replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase()),
        (v.genome && typeof v.genome.af === 'number') ? v.genome.af : null,
        0,
        pr.spliceai_ds_max == null ? null : +pr.spliceai_ds_max,
        pr.pangolin_largest_ds == null ? null : +pr.pangolin_largest_ds,
        pr.cadd == null ? null : +pr.cadd,
        pr.phylop == null ? null : +pr.phylop,
        // scope flag, so a whole-gene pull can be shown alongside non-coding-only ones
        /missense|synonymous|stop|frameshift|inframe|start lost|coding sequence|protein altering/i
          .test(v.consequence || '') ? 1 : 0
      ]);
    }
    process.stdout.write(`\r  ${reqs} requests, ${rows.length} variants`);
    await sleep(6500);
  }
  console.log('');

  // every row must sit inside the gene it is filed under
  const outside = rows.filter(r => r[1] < g.start || r[1] > g.end).length;
  if (outside) { console.log(`  ${outside} variants outside the gene span — refusing to write`); continue; }

  const pack = { gene: sym, chrom: rec.chrom,
    cols: ['rsid', 'pos', 'ref', 'alt', 'cons', 'af', 'curated', 'spliceai', 'pangolin', 'cadd', 'phylop', 'coding'], rows };
  fs.writeFileSync(R + `data/variants_${sym}.json`, JSON.stringify(pack));
  rec.n_variants = rows.length;
  GENES[sym] = rec;
  fs.writeFileSync(R + 'data/genes.json', JSON.stringify(GENES, null, 1));

  const sa = rows.map(r => r[7]).filter(x => x != null);
  console.log(`  wrote ${rows.length.toLocaleString()} variants  (${(fs.statSync(R + `data/variants_${sym}.json`).size / 1048576).toFixed(1)} MB)`);
  console.log(`  SpliceAI scored ${sa.length.toLocaleString()}  >=0.2: ${sa.filter(x => x >= 0.2).length}` +
    `  >=0.5: ${sa.filter(x => x >= 0.5).length}  >=0.8: ${sa.filter(x => x >= 0.8).length}`);
}
console.log('\ndone. reload browse.html — new genes appear as buttons automatically.');
