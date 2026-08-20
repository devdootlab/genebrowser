// Pull each gene's full reference span from Ensembl once, to data/seq_<GENE>.txt.
//
// The zoom fetched sequence from Ensembl live, per window. Even with block caching that is a
// cross-origin round-trip the first time you look anywhere new, so the panel sat on "fetching
// reference sequence…" for a second or more on every fresh region -- on a 40 bp window, waiting on
// the network for 40 characters.
//
// The whole reference for all 37 genes is 1.6 Mb, and only one gene is ever on screen: the median
// gene is 27 kb (~7 KB gzipped), the worst is F8 at 191 kb (~47 KB), which is smaller than any of
// the variant files sitting next to it. Shipping it makes zoom instant and removes a runtime
// dependency on Ensembl being up and fast.
//
//   node scripts/fetch_sequences.mjs            # only genes with no file yet
//   node scripts/fetch_sequences.mjs --force    # refetch everything
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
const FORCE = process.argv.includes('--force');
const sleep = ms => new Promise(r => setTimeout(r, ms));
const CHUNK = 50000;   // well inside Ensembl's limit, and a failure costs one chunk not a gene

async function region(chrom, a, b) {
  const c = String(chrom).replace('chr', '');
  for (let t = 0; t < 5; t++) {
    try {
      const r = await fetch(`https://rest.ensembl.org/sequence/region/human/${c}:${a}..${b}?content-type=application/json`,
        { headers: { accept: 'application/json' } });
      if (r.ok) {
        const j = await r.json();
        if (j && j.seq) return j.seq.toUpperCase();
      }
    } catch (e) { /* retry */ }
    await sleep(1200 * (t + 1));
  }
  return null;
}

let wrote = 0, skipped = 0, failed = [];
for (const [sym, g] of Object.entries(GENES)) {
  const f = R + `data/seq_${sym}.txt`;
  const want = g.end - g.start + 1;
  if (!FORCE && fs.existsSync(f) && fs.readFileSync(f, 'utf8').trim().length === want) { skipped++; continue; }

  let seq = '';
  for (let s = g.start; s <= g.end; s += CHUNK) {
    const e = Math.min(g.end, s + CHUNK - 1);
    const part = await region(g.chrom, s, e);
    if (part == null) { seq = null; break; }
    if (part.length !== e - s + 1) { console.log(`  ${sym}: chunk ${s}-${e} came back ${part.length} bp, expected ${e - s + 1}`); seq = null; break; }
    seq += part;
    process.stdout.write(`\r  ${sym}: ${(seq.length / 1000).toFixed(0)}/${(want / 1000).toFixed(0)} kb`);
    await sleep(250);
  }
  // Length is the whole check that matters: a short file silently shifts every base in the zoom,
  // which would look like a real sequence and be wrong by an offset.
  if (seq == null || seq.length !== want) { failed.push(sym); console.log(`\r  ${sym}: FAILED`); continue; }
  fs.writeFileSync(f, seq);
  wrote++;
  console.log(`\r  ${sym.padEnd(10)} ${(want / 1000).toFixed(0).padStart(4)} kb  ${(fs.statSync(f).size / 1024).toFixed(0)} KB on disk`);
}
console.log(`\n${wrote} written, ${skipped} already present, ${failed.length} failed${failed.length ? ': ' + failed.join(' ') : ''}`);
if (failed.length) process.exit(1);
