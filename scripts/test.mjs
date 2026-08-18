// Acceptance checks for 13b-genebrowser. Run: node scripts/test.mjs
//
// 13b had no tests at all, which is how five chr17 variants shipped as F2, a 4-track prediction
// shipped as a 667-track one, and a caption claimed pin height was phyloP conservation when every
// pin sat at one of two heights. Each check below exists because that specific thing was wrong.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const R = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const J = f => JSON.parse(fs.readFileSync(path.join(R, f), 'utf8'));
const has = f => fs.existsSync(path.join(R, f));

let pass = 0; const fails = [];
const ok = (c, m) => { if (c) pass++; else fails.push(m); };

const GENES = J('data/genes.json');
const GN = Object.keys(GENES);

/* ---- gene records ---------------------------------------------------------------------------- */
ok(GN.length >= 4, `genes: at least 4 genes, got ${GN.length}`);
ok(GN.every(g => GENES[g].exons && GENES[g].exons.length >= 2), 'genes: every gene has >=2 exons');
ok(GN.every(g => [1, -1].includes(GENES[g].strand)), 'genes: strand is +1 or -1');
// exons are [start,end] PAIRS. Reading .start off an array gave NaN and rendered an invisible
// gene model on all four genes, silently.
ok(GN.every(g => GENES[g].exons.every(e => Array.isArray(e) && e.length === 2 && e.every(Number.isFinite))),
   'genes: exons are [start,end] pairs of finite numbers');
ok(GN.every(g => GENES[g].exons.every(([a, b]) => a >= GENES[g].start && b <= GENES[g].end)),
   'genes: every exon lies inside its own gene span');
ok(GN.every(g => typeof GENES[g].n_variants === 'number' && GENES[g].n_variants > 0),
   'genes: every gene carries a variant count for the button label');

/* ---- per-gene variant files ------------------------------------------------------------------ */
for (const g of GN) {
  const f = `data/variants_${g}.json`;
  ok(has(f), `data: ${f} exists (per-gene split; the 23.6 MB monolith is withdrawn)`);
  if (!has(f)) continue;
  const p = J(f);
  const i = n => p.cols.indexOf(n);
  ok(p.rows.length === GENES[g].n_variants, `data: ${g} row count matches genes.json`);
  ok(p.chrom === GENES[g].chrom, `data: ${g} chromosome matches its gene record`);
  // five F2 entries were chr17 coordinates filed as chr11 F2, with figures drawn for all five
  const outside = p.rows.filter(r => r[i('pos')] < GENES[g].start || r[i('pos')] > GENES[g].end);
  ok(outside.length === 0, `data: ${g} has no variant outside its gene span (${outside.length} found)`);
  const afs = p.rows.map(r => r[i('af')]).filter(x => typeof x === 'number');
  ok(afs.every(x => x >= 0 && x <= 1), `data: ${g} allele frequencies are proportions, not percentages`);
  ok(new Set(afs.map(x => x.toFixed(8))).size > 20, `data: ${g} AF is not one repeated constant`);
  // Budget scales with the gene. VWF is 178 kb with 52 exons and legitimately carries 66,351
  // variants; a flat cap set when 746 KB was the worst case is simply wrong about it.
  const budget = Math.max(1.5e6, GENES[g].length * 45);
  const sz = fs.statSync(path.join(R, f)).size;
  ok(sz < budget, `data: ${g} is ${(sz/1e6).toFixed(1)} MB, budget ${(budget/1e6).toFixed(1)} MB`);
}

/* ---- the panels ------------------------------------------------------------------------------- */
if (has('data/panels.json')) {
  const P = J('data/panels.json');
  ok(P.drawn.length > 0, 'panels: something was drawn');
  ok(P.failed.length === 0, `panels: nothing failed (${P.failed.length} did)`);
  // ontology_terms=[] returns the full track set. Passing one cell type returned 4, and picking
  // "the most-affected track" from 4 is not the same operation as picking from 667.
  ok(P.drawn.every(r => r.n_tracks > 600), 'panels: every panel used the full track set (>600)');
  ok(P.drawn.every(r => typeof r.splice_log2fc === 'number' && Number.isFinite(r.splice_log2fc)),
     'panels: every panel carries a finite splice-site score');
  // canonical +-1/+-2 destroys a splice site. A near-zero median means the pipeline is wrong.
  if ((P.buckets||[]).includes('canonical')) {
    const med = P.drawn.map(r => Math.abs(r.splice_log2fc)).sort((a, b) => a - b)[Math.floor(P.drawn.length / 2)];
    ok(med > 0.1, `panels: canonical median |splice log2FC| is ${med.toFixed(3)}, must exceed 0.1 ` +
                  `(the RNA-only, 4-track build reported 0.005)`);
  }
  // one rsID can carry two alts; a gene_rsid stem overwrote the first and 40 drawn became 39 files
  const stems = P.drawn.map(r => r.id);
  ok(new Set(stems).size === stems.length, 'panels: figure filenames are unique across alleles');
  // Figures are gitignored — regenerable from the manifest. In a fresh clone the directory is
  // empty and that is correct, so only check the files when some are present.
  const anyFigs = has('images/panels') && fs.readdirSync(path.join(R,'images/panels')).some(f=>f.endsWith('.png'));
  if (anyFigs) ok(P.drawn.every(r => has(r.rich_img)), 'panels: every manifest entry has its file on disk');
  const onDisk = has('images/panels') ? fs.readdirSync(path.join(R, 'images/panels')).filter(f => f.endsWith('.png')) : [];
  if (anyFigs) ok(onDisk.length === P.drawn.length,
     `panels: ${onDisk.length} files on disk vs ${P.drawn.length} in the manifest`);
}

// Four genes were pulled non-coding-only and five whole-gene. Same screen, different populations.
// Every row now carries a coding flag and every gene records its pull scope.
for (const g of GN) {
  const p = J(`data/variants_${g}.json`);
  const ci = p.cols.indexOf('coding'), xi = p.cols.indexOf('cons');
  ok(typeof GENES[g].pull_scope === 'string', `scope: ${g} records how it was pulled`);
  ok(ci >= 0, `scope: ${g} rows carry a coding flag`);
  ok(p.rows.every(r => r[ci] === 0 || r[ci] === 1), `scope: ${g} coding flag is 0 or 1`);
  ok(p.rows.every(r => String(r[xi] || '').trim()), `scope: ${g} has no blank consequence label`);
  ok(GENES[g].n_coding + GENES[g].n_noncoding === p.rows.length, `scope: ${g} coding counts sum to its rows`);
}
ok(/SHOWCODING/.test(fs.readFileSync(path.join(R,'browse.html'),'utf8')),
   'scope: the page filters coding out by default so all genes show one population');

/* ---- honesty of the page ---------------------------------------------------------------------- */
const page = fs.readFileSync(path.join(R, 'browse.html'), 'utf8');
ok(!/gnomad_noncoding_variants\.json/.test(page),
   'page: does not load the withdrawn 23.6 MB monolith');
ok(/consOf/.test(page),
   'page: consequences go through one canonical mapping, not raw strings from two vocabularies');
ok(/variants \/ bin/.test(page), 'page: the density strip has a labelled y-axis');
ok(/rect class="bin"/.test(page) || /rect\.bin/.test(page), 'page: density bins are clickable');
ok(/homopolymer/.test(page), 'page: density outliers are flagged as possible artifacts');
// The old build printed 'PIN HEIGHT REPRESENTS PHYLOP 100-WAY CONSERVATION AT THAT BASE' while every
// pin sat at one of two heights. Mentioning phyloP is fine -- claiming to SHOW it is not.
ok(!/(height|pin)[^.]{0,40}(represents|is|=)[^.]{0,20}phyloP/i.test(page),
   'page: never claims a drawn height represents phyloP');
if (has('data/ref_mismatch.json')) {
  ok(Object.keys(J('data/ref_mismatch.json')).length > 0 && /REFBAD/.test(page),
     'page: variants with a wrong reference base carry a warning on their panel');
}

/* ---- secrets ------------------------------------------------------------------------------------ */
ok(has('.gitignore'), 'repo: a .gitignore exists');
if (has('.gitignore')) {
  const gi = fs.readFileSync(path.join(R, '.gitignore'), 'utf8');
  ok(/^\.env$/m.test(gi), 'repo: .env is ignored (a live AlphaGenome key lives there)');
}
const tracked = ['browse.html', 'index.html', 'scripts/make_panels.py', 'scripts/test.mjs'];
ok(tracked.filter(has).every(f => !/AIza[0-9A-Za-z_-]{20,}|sk-[A-Za-z0-9]{20,}/.test(
     fs.readFileSync(path.join(R, f), 'utf8'))), 'repo: no key-shaped string in tracked source');

console.log(fails.length ? `\n${pass} passed, ${fails.length} FAILED\n\n` + fails.map(f => '  x ' + f).join('\n')
                         : `\n${pass} passed, 0 failed — all acceptance checks pass`);
process.exit(fails.length ? 1 : 0);
