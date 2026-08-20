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
// THBD (thrombomodulin) is genuinely intronless -- one exon, no splice sites, nothing for an
// intronic browser to show. That is biology, not a bad pull, so the check is >=1 and the page has
// to SAY the gene has no introns rather than silently showing an empty splice track.
ok(GN.every(g => GENES[g].exons && GENES[g].exons.length >= 1), 'genes: every gene has >=1 exon');
{
  const single = GN.filter(g => GENES[g].exons.length < 2);
  ok(single.every(g => GENES[g].exon_count === 1), `genes: single-exon genes record it (${single.join(',') || 'none'})`);
  ok(!single.length || /no introns/.test(fs.readFileSync(path.join(R, 'browse.html'), 'utf8')),
     'genes: the page says so when a gene has no introns');
}
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
  // Budget scales with the VARIANT COUNT, not the gene length. The file is a list of variants, and
  // 45 bytes/bp mis-sized both ends: ADAMTS13 is dense for its length (2.4 MB against a 2.0 MB
  // length-based budget) while a long sparse gene got a budget it could never use.
  // A row is ~90 bytes: ["rs1371497467",139553812,"C","A","Intron Variant",null,0,null,null,1.2,0.4,0]
  // 120 leaves headroom for long rsIDs and long indel alleles while still catching the regression
  // this exists for -- the 23.6 MB all-genes monolith, and any accidental duplication of rows.
  const budget = Math.max(1.5e6, p.rows.length * 120);
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
  // Count against EVERY manifest, not just panels.json. There are several buckets now (canonical,
  // control, peak) and comparing the whole directory to one of them reported orphans that were
  // simply another bucket's figures. A genuine orphan -- a PNG no manifest points at -- is still a
  // failure: it is a figure nobody can trace back to a variant, which is how a stale pre-fix figure
  // survives a re-run and gets read as current.
  const onDisk = has('images/panels') ? fs.readdirSync(path.join(R, 'images/panels')).filter(f => f.endsWith('.png')) : [];
  const claimed = new Set();
  for (const f of fs.readdirSync(path.join(R, 'data')).filter(f => /^panels.*\.json$/.test(f)))
    for (const r of (J('data/' + f).drawn || [])) claimed.add(path.basename(r.rich_img || ''));
  const orphans = onDisk.filter(f => !claimed.has(f));
  if (anyFigs) ok(orphans.length === 0,
     `panels: ${orphans.length} figures on disk that no manifest claims (${orphans.slice(0, 4).join(', ')}${orphans.length > 4 ? ', …' : ''})`);
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

/* ---- strand ------------------------------------------------------------------------------------ */
// Splicing reads the TRANSCRIBED strand. browse.html searched the plus strand for GT/AG on every
// gene, which scored 0/10, 0/24, 0/4, 0/8, 0/50, 0/25 and 0/6 canonical motifs on the seven
// minus-strand genes -- every mark wrong, every real site unmarked. scripts/check_strand_motifs.mjs
// is the known-answer run (150/152, 98.7%); these checks stop the rule silently reverting.
ok(/txBase|txDi/.test(page), 'strand: the zoom has a transcript-orientation helper');
ok(/strand === -1 \? txBase/.test(page), 'strand: dinucleotides are read in transcript orientation');
ok(/annotateSeq\([^)]*strand\)/.test(page), 'strand: annotateSeq is passed the gene strand');
// the polypyrimidine tract is pyrimidines ON THE TRANSCRIPT -- same bug, second place
ok(/var tx = seq\.split/.test(page), 'strand: polypyrimidine runs on the transcribed strand');
ok(/minus-strand gene/.test(page), 'strand: the page says so when the gene is on the minus strand');

/* ---- indels ------------------------------------------------------------------------------------ */
// Every row rendered as one circle at r.pos, so a 13 bp deletion, a 6 bp insertion and a point
// substitution were the same mark -- and rs753240633's 16 length alleles were one dot.
ok(/function vtype/.test(page), 'indels: rows are classified by type');
ok(/deletion/.test(page) && /insertion/.test(page), 'indels: the zoom legend names both');
ok(/k:'vtype'/.test(page), 'indels: the table has a sortable type column');
ok(/ONLY === 'indel'/.test(page), 'indels: there is an indels-only filter');
ok(/AlphaGenome scores substitutions only/.test(page),
   'indels: an indel says WHY it has no panel rather than a bare "no panel"');
// the classifier itself, on the real SERPINC1 poly-A ladder
{
  const vt = new Function('r', page.match(/function vtype\(r\)\{[^]*?\n\}/)[0] + '; return vtype(r);');
  const cases = [['C','T','SNV'],['CA','C','del'],['C','CA','ins'],['CAAAAAAAAAAAAA','C','del'],['AT','GC','MNV']];
  ok(cases.every(c => vt({ref:c[0],alt:c[1]}).kind === c[2]), 'indels: vtype classifies the poly-A ladder correctly');
  ok(vt({ref:'CAAAAAAAAAAAAA',alt:'C'}).d === -13, 'indels: a 13-base deletion reports -13, not -14 or -1');
}

/* ---- favicon ----------------------------------------------------------------------------------- */
// A favicon that 404s looks identical to no favicon at all: a blank page glyph, which is exactly
// what makes a tab unfindable. Both files must exist and both must be referenced.
ok(has('favicon.svg'), 'icon: favicon.svg exists');
ok(has('favicon.png'), 'icon: favicon.png exists (fallback for browsers without SVG icon support)');
ok(/rel="icon"[^>]*favicon\.svg/.test(page), 'icon: the page links favicon.svg');
ok(/favicon\.png/.test(page), 'icon: the page links the PNG fallback');
{
  // SVG is XML, and this repo writes `--` in comments everywhere. A double hyphen inside an XML
  // comment is a hard parse error, so the first version of favicon.svg served HTTP 200 with the
  // right content-type and would not decode. Nothing looked wrong: a favicon that fails to parse
  // renders exactly like no favicon.
  const svg = fs.readFileSync(path.join(R, 'favicon.svg'), 'utf8');
  const comments = svg.match(/<!--[\s\S]*?-->/g) || [];
  ok(comments.every(c => !/--/.test(c.slice(4, -3))),
     'icon: no double hyphen inside an SVG comment (hard XML parse error)');
  ok(/^\s*<svg[\s>]/.test(svg) && /<\/svg>\s*$/.test(svg), 'icon: favicon.svg opens and closes an <svg> root');
  // crude well-formedness: every non-void tag that opens must close, and vice versa
  const opens = (svg.match(/<(svg|g|path|circle|title|defs|use)\b(?![^>]*\/>)/g) || []).length;
  const closes = (svg.match(/<\/(svg|g|path|circle|title|defs|use)>/g) || []).length;
  ok(opens === closes, `icon: favicon.svg tags balance (${opens} open, ${closes} close)`);

  // a real PNG, not a pasted blob nobody rendered
  const b = fs.readFileSync(path.join(R, 'favicon.png'));
  ok(b.slice(0, 8).toString('hex') === '89504e470d0a1a0a', 'icon: favicon.png has a valid PNG signature');
  ok(b.readUInt32BE(16) === 32 && b.readUInt32BE(20) === 32, 'icon: favicon.png is 32x32');
  ok(b[25] === 6, 'icon: favicon.png is RGBA, so the corners are transparent not black');
}
ok(!/href="data:image\/png;base64/.test(page),
   'icon: no inline base64 image (an icon nobody rendered is indistinguishable from a corrupt one)');

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
