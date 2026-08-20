// Backfill the three fields scripts/test.mjs requires onto genes pulled before add_genes.mjs wrote
// them: pull_scope, n_coding/n_noncoding, and a non-blank consequence label.
//
// gnomAD returns consequence: null on a minority of rows. add_genes.mjs turned that into '', which
// renders as a blank filter chip that looks like a category and matches nothing.
//
//   node scripts/backfill_scope.mjs
import fs from 'node:fs';

import path from 'node:path';
import { fileURLToPath } from 'node:url';
// repo root, derived from this file's location -- no absolute path baked in
const R = path.dirname(path.dirname(fileURLToPath(import.meta.url))) + path.sep;
const GENES = JSON.parse(fs.readFileSync(R + 'data/genes.json', 'utf8'));
let touched = 0, blanks = 0;

for (const sym of Object.keys(GENES)) {
  const f = R + `data/variants_${sym}.json`;
  if (!fs.existsSync(f)) { console.log(`  ${sym}: no variants file, skipping`); continue; }
  const p = JSON.parse(fs.readFileSync(f, 'utf8'));
  const ix = n => p.cols.indexOf(n);
  const ci = ix('coding'), xi = ix('cons');
  let fixed = 0;
  for (const r of p.rows) {
    if (!String(r[xi] || '').trim()) { r[xi] = 'Unannotated'; fixed++; }
    if (r[ci] !== 0 && r[ci] !== 1) r[ci] = 0;
  }
  if (fixed) { fs.writeFileSync(f, JSON.stringify(p)); blanks += fixed; }

  const g = GENES[sym];
  const nc = p.rows.filter(r => r[ci] === 1).length;
  const before = JSON.stringify([g.pull_scope, g.n_coding, g.n_noncoding]);
  if (typeof g.pull_scope !== 'string')
    g.pull_scope = 'whole gene (coding and non-coding), gnomAD v4.1 genomes, region query';
  g.n_coding = nc;
  g.n_noncoding = p.rows.length - nc;
  g.n_variants = p.rows.length;
  if (JSON.stringify([g.pull_scope, g.n_coding, g.n_noncoding]) !== before) touched++;
  console.log(`  ${sym.padEnd(10)} ${String(p.rows.length).padStart(7)} rows  ` +
    `${String(nc).padStart(5)} coding  ${String(p.rows.length - nc).padStart(7)} non-coding` +
    (fixed ? `  (${fixed} blank consequences labelled)` : ''));
}
fs.writeFileSync(R + 'data/genes.json', JSON.stringify(GENES, null, 1));
console.log(`\n${touched} gene records updated, ${blanks} blank consequence labels filled`);
