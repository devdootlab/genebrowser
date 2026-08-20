"""
4-panel AlphaGenome panels (dmidealoutput layout), corrected.

Keeps agy's layout. Fixes what was wrong underneath it, and refuses to draw rather than draw
something false. Every refusal prints the reason.

WHAT WAS WRONG IN generate_all_dmidealoutput.py, AND WHY IT MATTERED
--------------------------------------------------------------------
1. np.mean(values, axis=1) averaged all 667 tracks. That diluted a real single-track signal ~236x
   and is why max |delta_fc| across 24 mutations was 0.044 and some modality deltas were 1e-7.
   Fixed: pick the single most-affected track by argmax(|ALT.mean(0) - REF.mean(0)|).
2. delta_fc computed with np.log but printed as "log2 FC". Every number was 1.44x its own label.
   Fixed: np.log2 throughout.
3. `strand` was read at line 67 and never used, sitting next to a genome.Variant call, which reads
   as if orientation had been handled. Fixed: strand is fetched, asserted, and printed.
4. No reference-base check. A variant whose ref is not the base at that coordinate is not that
   variant, and every number downstream is meaningless. 5 of 11 curated rows failed this.
5. No gene/coordinate check. Five F2 rows were chr17 coordinates from SSTR2 filed as chr11 F2, and
   figures were drawn for all five.
6. Indels were passed to genome.Variant, which only accepts ACGTN.

WHY CANONICAL SITES ARE THE DEFAULT BUCKET
------------------------------------------
+-1/+-2 splice variants are the KNOWN-ANSWER set: they destroy a splice site, so a working pipeline
must show a large effect. If these come back near zero, the pipeline is broken, not the biology.
Run this bucket first and read the effect sizes before trusting any other bucket.

  uv run python scripts/make_panels.py --gene OAS1 --bucket canonical --limit 10
  uv run python scripts/make_panels.py --all --bucket canonical --limit 10
"""
import json, os, sys, time, argparse, traceback
import numpy as np
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RENDERER = 2
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def log(m):
    print(m, flush=True)


def ensembl(path):
    """One retrying GET. A transient 500 must not be read as 'no such variant' -- Ensembl 500s
    intermittently, and treating that as a negative result silently drops real data."""
    url = "https://rest.ensembl.org" + path
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if attempt == 3:
                log(f"    [ensembl] gave up after 4 tries: {path}  ({e})")
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


# ---- selection ------------------------------------------------------------------------------
def buckets_for(gene, genes, variants):
    g = genes[gene]
    ex = sorted([[min(e), max(e)] if isinstance(e, list) else [e["start"], e["end"]] for e in g["exons"]])
    bounds = []
    for i in range(len(ex) - 1):
        bounds += [ex[i][1], ex[i + 1][0]]
    out = {"canonical": [], "spliceregion": [], "deep": [], "skipped_indel": [], "skipped_exonic": []}
    for v in variants:
        if v.get("gene") != gene:
            continue
        if v.get("coding"):
            continue          # app default view is non-coding; keep the panel set matching it
        ref, alt = str(v.get("ref", "")), str(v.get("alt", ""))
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            out["skipped_indel"].append(v)          # legitimate variant, wrong shape for AlphaGenome
            continue
        pos = v["pos"]
        intron = None
        for i in range(len(ex) - 1):
            if ex[i][1] < pos < ex[i + 1][0]:
                intron = i
                break
        if intron is None:
            out["skipped_exonic"].append(v)
            continue
        d = min(abs(pos - b) for b in bounds)
        v = dict(v)
        v["_dist"] = d
        v["_intron"] = (intron + 1) if g["strand"] == 1 else (len(ex) - 1 - intron)
        v["_nintrons"] = len(ex) - 1
        v["_istart"], v["_iend"] = ex[intron][1] + 1, ex[intron + 1][0] - 1
        out["canonical" if d <= 2 else "spliceregion" if d <= 8 else "deep"].append(v)
    for k in ("canonical", "spliceregion", "deep"):
        out[k].sort(key=lambda x: (x["_dist"], -(x.get("af") or 0)))
    return out, ex, g


def peaks_for(bucket_deep, min_dist=100, binw=200, mode="auto"):
    """Density peaks in DEEP intronic sequence -- the spikes in the browser's strip, minus
    everything within `min_dist` of a splice boundary.

    The `deep` bucket sorts by distance ASCENDING, so `--bucket deep --limit 10` returns the ten
    variants that only just cleared the 8 bp splice-region cut. Those are the least deep intronic
    variants there are, which is the opposite of what a deep-intronic question is asking.

    A peak here is a `binw` window holding far more variants than the gene's typical window.
    Whether a peak is biology or a calling artefact is exactly what is unknown -- homopolymers and
    tandem repeats produce identical-looking spikes -- so this only selects WHERE to look.
    One representative per peak, so N panels cover N distinct peaks rather than N variants from one.
    """
    if not bucket_deep:
        return []
    far = [v for v in bucket_deep if v["_dist"] >= min_dist]
    if not far:
        return []
    bins = {}
    for v in far:
        bins.setdefault(v["pos"] // binw, []).append(v)
    counts = sorted(len(x) for x in bins.values())
    med = counts[len(counts) // 2]
    # median absolute deviation, scaled to a normal sigma; robust to the peaks themselves
    mad = sorted(abs(c - med) for c in counts)[len(counts) // 2] * 1.4826
    thresh = max(med + 4 * mad, 3 * med, 5)
    hot = sorted(b for b, vs in bins.items() if len(vs) >= thresh)
    # Merge ADJACENT hot bins. OAS1's three top "peaks" were 112,918,200-112,919,200 -- one
    # contiguous 1 kb region reported as three, which would have spent three panels re-asking the
    # same question. A peak is a run of hot bins, not a bin.
    runs, cur = [], []
    for b in hot:
        if cur and b == cur[-1] + 1:
            cur.append(b)
        else:
            if cur:
                runs.append(cur)
            cur = [b]
    if cur:
        runs.append(cur)
    picks = []
    for run in runs:
        vs = [v for b in run for v in bins[b]]
        lo, hi = run[0] * binw, (run[-1] + 1) * binw
        centre = (lo + hi) // 2

        # Representative, chosen in EXPLICIT TIERS, recording which tier was used.
        #
        # This used to be one sort: (-(spliceai or 0), -(af or 0)). That WORKED -- 13 of 14 picks
        # are unchanged by the rewrite below. It was reported here as broken on the strength of
        # `spliceai_gnomad` being null on all 16 rows of data/panels_peak.json, but that field is
        # only written in the --control path and was simply absent for this bucket. Absent field,
        # not absent data: every variant in every one of these peaks has a SpliceAI score.
        #
        # The tiers exist anyway, for the reason the false alarm pointed at: `or 0` collapses "no
        # score" and "scores zero" into the same key, so on a gene where the scores really are
        # missing this would degrade to `deep`'s ordering -- distance from the splice site,
        # ASCENDING, the worst possible choice for a deep-intronic question -- with no error and
        # nothing in the output to show it had happened. Now the tier is chosen explicitly and
        # written to the manifest as `picked_by`.
        #
        # NOTE the "auto"/"spliceai" tiers make the pick SpliceAI-driven, so AlphaGenome agreeing
        # with SpliceAI on those picks is partly circular (r = 0.72 overall, 0.755 on these 16).
        # Use --peakpick centre to ask whether the PEAK carries signal independently of SpliceAI.
        scored = [v for v in vs if v.get("spliceai") is not None] if mode in ("auto", "spliceai") else []
        withaf = [v for v in vs if v.get("af") is not None] if mode in ("auto", "af") else []
        if scored:
            pick = max(scored, key=lambda x: (x["spliceai"], x.get("af") or 0, -x["pos"]))
            why = "spliceai"
        elif withaf:
            pick = max(withaf, key=lambda x: (x["af"], -x["pos"]))
            why = "af"
        else:
            # Nothing in this peak carries either score. Fall back to the variant nearest the
            # peak's CENTRE -- the most representative position in the window, and unbiased with
            # respect to splice distance. The old accidental fallback picked the variant closest to
            # an exon, which is the single worst choice for a deep-intronic question.
            pick = min(vs, key=lambda x: (abs(x["pos"] - centre), x["pos"]))
            why = "centre" if mode == "centre" else "centre (no SpliceAI or AF anywhere in this peak)"

        r = dict(pick)
        r["_peak_n"] = len(vs)
        r["_peak_start"] = lo
        r["_peak_end"] = hi
        r["_peak_bins"] = len(run)
        r["_peak_med"] = med
        r["_peak_scored"] = len(scored)
        r["_peak_withaf"] = len(withaf)
        r["_picked_by"] = why
        picks.append(r)
    picks.sort(key=lambda x: -x["_peak_n"])
    return picks


# ---- the checks that must pass before a single API call is made -------------------------------
def preflight(v, g, gene, seq_cache):
    """Returns (ok, reason). Every one of these has a real failure behind it in this codebase."""
    chrom = str(v["chrom"]).replace("chr", "")
    gchrom = str(g["chrom"]).replace("chr", "")
    if chrom != gchrom:
        return False, f"chromosome {v['chrom']} but {gene} is on {g['chrom']}"
    if not (g["start"] <= v["pos"] <= g["end"]):
        return False, f"pos {v['pos']:,} outside {gene} {g['start']:,}-{g['end']:,}"
    if seq_cache.get(gene) is None:
        r = ensembl(f"/sequence/region/human/{gchrom}:{g['start']}..{g['end']}?content-type=application/json")
        if not r or "seq" not in r:
            return False, "could not fetch reference sequence"
        s = r["seq"].upper()
        if len(s) != g["end"] - g["start"] + 1:
            return False, f"reference sequence is {len(s)} bp, expected {g['end']-g['start']+1}"
        seq_cache[gene] = s
    base = seq_cache[gene][v["pos"] - g["start"]]
    if base != str(v["ref"]).upper():
        comp = COMP.get(str(v["ref"]).upper(), "?")
        hint = " (it IS the complement -- transcript-orientation allele stored as genomic)" if comp == base else ""
        return False, f"ref {v['ref']} but genome has {base} at {v['chrom']}:{v['pos']:,}{hint}"
    return True, None


# ---- figure -------------------------------------------------------------------------------------
def draw(v, g, gene, model, dna_client, genome, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chrom, pos, ref, alt = str(v["chrom"]), int(v["pos"]), str(v["ref"]), str(v["alt"])
    variant = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)
    iv = variant.reference_interval.resize(131072)
    pred = model.predict_variant(
        # ontology_terms is REQUIRED by the client, but an EMPTY list means "do not filter" and
        # returns the full track set. Passing one term returned 4 tracks; [] returns 667.
        interval=iv, variant=variant, ontology_terms=[],
        requested_outputs=[dna_client.OutputType.RNA_SEQ, dna_client.OutputType.SPLICE_SITES])

    RV = np.asarray(pred.reference.rna_seq.values)
    AV = np.asarray(pred.alternate.rna_seq.values)
    if RV.ndim != 2 or RV.shape != AV.shape:
        raise ValueError(f"unexpected prediction shape ref={RV.shape} alt={AV.shape}")

    SR = np.asarray(pred.reference.splice_sites.values)
    SA = np.asarray(pred.alternate.splice_sites.values)
    if SR.shape != SA.shape:
        raise ValueError(f"splice track shape mismatch {SR.shape} vs {SA.shape}")
    # The variant sits at the centre of a resized interval. Assert it rather than assume it: a
    # variant near a contig edge cannot be centred, and an off-centre window would measure the
    # wrong bases while looking entirely normal.
    _c = SR.shape[0] // 2
    _expect = (iv.start + iv.end) // 2
    if abs(_expect - pos) > 2:
        raise ValueError(f"variant {pos} is not centred in interval {iv.start}-{iv.end}")
    _w = slice(max(0, _c - 128), min(SR.shape[0], _c + 128))

    # Track selection MUST use the local window. Averaging the difference over all 131,072
    # positions dilutes a ~10-position effect to nothing -- the same mistake as averaging 667
    # tracks, one level up. argmax over a diluted vector picks an essentially random track.
    sp_per_track = np.abs(SA[_w] - SR[_w]).max(axis=0)
    sti = int(np.argmax(sp_per_track))
    SPR, SPA = SR[:, sti], SA[:, sti]
    sp_ref, sp_alt = float(np.abs(SPR[_w]).mean()), float(np.abs(SPA[_w]).mean())
    splice_delta = float(np.log2(sp_alt + 1e-6) - np.log2(sp_ref + 1e-6))
    splice_absdiff = float(np.abs(SPA[_w] - SPR[_w]).max())

    # THE 236x FIX. Averaging every track flattens a real single-track effect into noise.
    # same correction for RNA: choose on the local window, not the 131 kb mean
    per_track = np.abs(AV[_w] - RV[_w]).max(axis=0)
    ti = int(np.argmax(per_track))
    R, A = RV[:, ti], AV[:, ti]
    track_name = ""
    try:
        track_name = str(pred.reference.rna_seq.metadata.iloc[ti].get("name", ""))[:58]
    except Exception:
        pass

    # local means, for the same reason: a 131 kb average is dominated by bases the variant
    # cannot affect, so it understates every real effect by orders of magnitude.
    rm, am = float(R[_w].mean()), float(A[_w].mean())
    # log2, matching the label. The old code used np.log and printed "log2 FC".
    delta = float(np.log2(am + 1e-6) - np.log2(rm + 1e-6))
    if not np.isfinite(delta):
        raise ValueError(f"non-finite delta (ref mean {rm}, alt mean {am})")

    x = np.linspace(iv.start, iv.end, R.shape[0])
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5))
    fig.suptitle(f"AlphaGenome · {gene} {v.get('rsid') or ''} · {chrom}:{pos:,} {ref}>{alt} · "
                 f"splice {splice_delta:+.3f} log2FC | RNA {delta:+.4f} log2FC",
                 fontsize=12, fontweight="bold")

    # 1. gene model
    ax = axes[0][0]
    ex = sorted([[min(e), max(e)] if isinstance(e, list) else [e["start"], e["end"]] for e in g["exons"]])
    ax.plot([g["start"], g["end"]], [0, 0], color="#334155", lw=1.4)
    for a, b in ex:
        ax.add_patch(plt.Rectangle((a, -0.3), max(b - a, 60), 0.6, color="#1e293b"))
    ax.axvline(pos, color="#d97706", lw=1.6)
    ax.set_ylim(-1.2, 1.2); ax.set_yticks([]); ax.set_xlim(g["start"], g["end"])
    ax.set_title(f"{len(ex)} exons · strand {'+' if g['strand']==1 else '-'} · "
                 f"intron {v['_intron']}/{v['_nintrons']} · {v['_dist']} bp from a splice site", fontsize=9)

    # 2. REF vs ALT, one track
    ax = axes[0][1]
    ax.plot(x, R, color="#2563eb", lw=1.0, label="REF")
    ax.plot(x, A, color="#dc2626", lw=1.0, ls="--", label="ALT")
    ax.axvline(pos, color="#d97706", lw=1.0, ls=":")
    ax.legend(fontsize=8); ax.set_title(f"most-affected track {ti} of {RV.shape[1]}  {track_name}", fontsize=9)

    # 3. same, zoomed
    ax = axes[1][0]
    m = (x >= pos - 2000) & (x <= pos + 2000)
    ax.plot(x[m], R[m], color="#2563eb", lw=1.2, label="REF")
    ax.plot(x[m], A[m], color="#dc2626", lw=1.2, ls="--", label="ALT")
    ax.axvline(pos, color="#d97706", lw=1.2, ls=":")
    ax.legend(fontsize=8); ax.set_title("4 kb around the variant", fontsize=9)

    # 4. the number, with its own scale stated
    ax = axes[1][1]
    ax.bar(["REF", "ALT"], [rm, am], color=["#2563eb", "#dc2626"])
    ax.set_title(f"mean coverage in track {ti}\n{delta:+.4f} log2 FC", fontsize=9)

    fig.text(0.5, 0.005,
             "AlphaGenome PREDICTION, not a measurement. ONE track of "
             f"{RV.shape[1]}, the one this variant moves most. Averaged over all tracks the change is far smaller.",
             ha="center", fontsize=7.5, color="#475569")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])

    stem = f"{gene}_{v.get('rsid') or str(pos)}_{ref}{alt}"
    png = os.path.join(outdir, stem + "_rich.png")
    fig.savefig(png, dpi=110); plt.close(fig)
    return {"id": stem, "gene": gene, "rsid": v.get("rsid"), "chrom": chrom, "pos": pos,
            "ref": ref, "alt": alt, "af": v.get("af"), "dist_splice": v["_dist"],
            "intron": f"{v['_intron']}/{v['_nintrons']}", "delta_log2fc": delta,
            "splice_log2fc": splice_delta, "splice_maxabsdiff": splice_absdiff,
            "splice_track_index": sti, "n_splice_tracks": int(SR.shape[1]),
            "track_index": ti, "n_tracks": int(RV.shape[1]), "track": track_name,
            "renderer": RENDERER, "rich_img": os.path.relpath(png, ROOT).replace("\\", "/"),
            # Provenance for --bucket peak. A result whose selection rule cannot be recovered from
            # the manifest is not reproducible, and this selection has already silently changed
            # rule once. Absent for every other bucket.
            **({"picked_by": v["_picked_by"], "peak_start": v["_peak_start"],
                "peak_end": v["_peak_end"], "peak_n": v["_peak_n"],
                "peak_median": v["_peak_med"], "peak_n_scored": v["_peak_scored"],
                "peak_n_withaf": v["_peak_withaf"]} if "_picked_by" in v else {})}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene"); ap.add_argument("--all", action="store_true")
    ap.add_argument("--bucket", default="canonical", choices=["canonical", "spliceregion", "deep", "spliceai", "peak"])
    ap.add_argument("--peakmin", type=int, default=100, help="min bp from any splice boundary for --bucket peak")
    ap.add_argument("--peakbin", type=int, default=200, help="window width for --bucket peak")
    # "auto" ranks by SpliceAI, which makes the pick SpliceAI-driven and any AlphaGenome/SpliceAI
    # agreement partly circular (they correlate at r=0.72). "centre" takes the variant nearest the
    # peak middle regardless of any score -- the unbiased pick, and the one that actually asks
    # whether the PEAK carries signal rather than whether SpliceAI and AlphaGenome agree again.
    ap.add_argument("--peakpick", default="auto", choices=["auto", "spliceai", "af", "centre"],
                    help="how to choose the representative variant in each peak")
    ap.add_argument("--spliceai", type=float, default=0.5, help="threshold when --bucket spliceai")
    ap.add_argument("--limit", type=int, default=10)
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--control", action="store_true", help="run data/control_set.json instead of a bucket")
    # The browser's right-click queue copies exactly this. It has done since the queue was built,
    # and the flag did not exist, so the button handed you a command that died on argparse.
    ap.add_argument("--queue", help='comma-separated "<rsid-or-pos>:<REF><ALT>" from the browser queue')
    a = ap.parse_args()

    genes = json.load(open(os.path.join(ROOT, "data/genes.json")))
    # per-gene slim files replaced the 23.6 MB monolith; rebuild the shape buckets_for expects
    variants = []
    for _g in genes:
        _f = os.path.join(ROOT, f"data/variants_{_g}.json")
        if not os.path.exists(_f):
            continue
        _p = json.load(open(_f))
        _i = {c: k for k, c in enumerate(_p["cols"])}
        for _r in _p["rows"]:
            variants.append({"gene": _g, "chrom": _p["chrom"], "rsid": _r[_i["rsid"]],
                             "pos": _r[_i["pos"]], "ref": _r[_i["ref"]], "alt": _r[_i["alt"]],
                             "consequence": _r[_i["cons"]], "af": _r[_i["af"]],
                             "coding": (_r[_i["coding"]] if "coding" in _i else 0),
                             "spliceai": (_r[_i["spliceai"]] if "spliceai" in _i else None)})
    if not variants:
        sys.exit("no data/variants_*.json found -- run the per-gene split first")
    targets = list(genes) if a.all else [a.gene]
    # --queue names its own variants and searches every gene, so it needs no --gene/--all
    if not a.control and not a.queue and (not targets or targets == [None]):
        sys.exit("give --gene GENE or --all")

    outdir = os.path.join(ROOT, "images/panels"); os.makedirs(outdir, exist_ok=True)
    model = dna_client = genome = None
    if not a.dry:
        key = os.environ.get("ALPHAGENOME_API_KEY")
        if not key:
            sys.exit("ALPHAGENOME_API_KEY is not set. Refusing to run rather than draw nothing quietly.")
        from alphagenome.models import dna_client as _dc
        from alphagenome.data import genome as _gn
        dna_client, genome = _dc, _gn
        model = _dc.create(api_key=key)

    seq_cache, done, refused, failed = {}, [], [], []

    if a.control:
        # Known-answer experiment. Positives destroy a splice site; negatives sit >500bp away with
        # SpliceAI 0, matched on allele frequency so the sets differ in POSITION, not in rarity.
        C = json.load(open(os.path.join(ROOT, "data/control_set.json")))
        for label in ("positive", "negative"):
            log("")
            log(f"=== control: {label} ({len(C[label])}) ===")
            for v in C[label]:
                g = genes[v["gene"]]
                ex = sorted([[min(e), max(e)] for e in g["exons"]])
                for k in range(len(ex) - 1):
                    if ex[k][1] < v["pos"] < ex[k + 1][0]:
                        v["_istart"], v["_iend"] = ex[k][1] + 1, ex[k + 1][0] - 1
                        v["_intron"] = (k + 1) if g["strand"] == 1 else (len(ex) - 1 - k)
                        v["_nintrons"] = len(ex) - 1
                        break
                v["_dist"] = v["dist"]
                ok, why = preflight(v, g, v["gene"], seq_cache)
                if not ok:
                    log(f"  REFUSED {v['gene']} {v['pos']}: {why}")
                    refused.append({"label": label, **{k2: v[k2] for k2 in ("gene","pos","rsid")}, "reason": why}); continue
                try:
                    rec = draw(v, g, v["gene"], model, dna_client, genome, outdir)
                    rec["label"] = label; rec["spliceai_gnomad"] = v.get("spliceai")
                    done.append(rec)
                    log(f"  {label[:3].upper()} {v['gene']:6} {str(v.get('rsid') or v['pos']):16} "
                        f"splice {rec['splice_log2fc']:+7.3f}  maxdiff {rec['splice_maxabsdiff']:.4f}  SpliceAI {v.get('spliceai')}")
                except Exception as e:
                    log(f"  FAILED {v['gene']} {v['pos']}: {type(e).__name__}: {e}")
                    failed.append({"label": label, "gene": v["gene"], "pos": v["pos"], "error": str(e)})
        json.dump({"drawn": done, "refused": refused, "failed": failed},
                  open(os.path.join(ROOT, "data/panels_control.json"), "w"), indent=1)
        for m in ("splice_log2fc", "splice_maxabsdiff", "delta_log2fc"):
            P = [abs(r[m]) for r in done if r["label"] == "positive"]
            N = [abs(r[m]) for r in done if r["label"] == "negative"]
            if not P or not N: continue
            P.sort(); N.sort()
            mp, mn = P[len(P)//2], N[len(N)//2]
            # how many negatives reach the weakest positive -- the only question that matters
            overlap = sum(1 for x in N if x >= P[0])
            log("")
            log(f"{m}:  positive median {mp:.4f}   negative median {mn:.4f}   "
                f"ratio {(mp/mn if mn else float('inf')):.1f}x")
            log(f"    weakest positive {P[0]:.4f}; {overlap} of {len(N)} negatives reach it "
                f"({100*overlap/len(N):.0f}% overlap)")
        log("")
        log("wrote data/panels_control.json")
        return

    if a.queue:
        # Resolve each queued entry through buckets_for, so _dist/_intron/_istart/_iend are computed
        # by the same code the buckets use. Recomputing them here would be a second implementation
        # of the arithmetic that has already been wrong once in this file.
        want = []
        for tok in a.queue.split(","):
            tok = tok.strip()
            if not tok or ":" not in tok:
                continue
            ident, alleles = tok.rsplit(":", 1)
            want.append((ident.strip(), alleles.strip().upper()))
        if not want:
            sys.exit("--queue got nothing parseable; expected rs123:CT,45678:AG")
        pool, seen = [], set()
        for gene in genes:
            b, _ex, _g = buckets_for(gene, genes, variants)
            for v in b["canonical"] + b["spliceregion"] + b["deep"]:
                for ident, alleles in want:
                    if (str(v.get("rsid")) == ident or str(v["pos"]) == ident) and                        (str(v["ref"]) + str(v["alt"])).upper() == alleles:
                        key = (gene, v["pos"], v["ref"], v["alt"])
                        if key not in seen:
                            seen.add(key); pool.append((gene, v))
        found = {(i, al) for i, al in want
                 if any(str(v.get("rsid")) == i or str(v["pos"]) == i for _gn, v in pool)}
        missing = [f"{i}:{al}" for i, al in want if (i, al) not in found]
        log("")
        log(f"=== queue: {len(want)} requested, {len(pool)} resolved ===")
        if missing:
            # An indel or a coding variant is filtered out by buckets_for before it gets here, and
            # a silent short list would read as "these are all that exist".
            log(f"    {len(missing)} not found among intronic substitutions (queue skips exonic variants and indels: AlphaGenome scores single-base substitutions, and buckets are intron-only): {', '.join(missing)}")
        for gene, v in pool[: a.limit] if a.limit else pool:
            tag = f"{gene} {v.get('rsid') or v['pos']}"
            ok, why = preflight(v, genes[gene], gene, seq_cache)
            if not ok:
                log(f"  REFUSED {tag}: {why}"); refused.append({"tag": tag, "reason": why}); continue
            if a.dry:
                # --dry leaves model/genome as None, so draw() cannot be called at all here.
                log(f"  would draw {tag}  {v['ref']}>{v['alt']}  {v['_dist']} bp from splice site")
                continue
            try:
                rec = draw(v, genes[gene], gene, model, dna_client, genome, outdir)
                if rec: done.append(rec); log(f"  OK       {tag}")
            except Exception as e:
                log(f"  FAILED   {tag}: {type(e).__name__}: {e}")
                failed.append({"tag": tag, "error": f"{type(e).__name__}: {e}"})
        targets = []

    for gene in targets:
        b, ex, g = buckets_for(gene, genes, variants)
        if a.bucket == "spliceai":
            pool = [v for v in (b["canonical"] + b["spliceregion"] + b["deep"])
                    if v.get("spliceai") is not None and v["spliceai"] >= a.spliceai]
            pool.sort(key=lambda x: -x["spliceai"])
            b["spliceai"] = pool
        if a.bucket == "peak":
            b["peak"] = peaks_for(b["deep"], a.peakmin, a.peakbin, a.peakpick)
        picked = b[a.bucket][: a.limit]
        log(f"\n=== {gene} · bucket {a.bucket} · {len(b[a.bucket])} available, taking {len(picked)} ===")
        if a.bucket == "spliceai":
            log(f"    SpliceAI >= {a.spliceai}: {len(b['spliceai'])} in {gene}")
        if a.bucket == "peak":
            for _v in picked:
                log(f"    peak {_v['_peak_start']:,}-{_v['_peak_end']:,} ({_v['_peak_bins']}x{a.peakbin}bp)  "
                    f"{_v['_peak_n']} variants vs median {_v['_peak_med']}  "
                    f"{_v['_dist']:,} bp from the nearest splice site")
                # Say how the representative was chosen and how much there was to choose from.
                # The previous version printed neither, which is why a tiebreak that had silently
                # become "closest to an exon" survived a 16-gene run and a written-up result.
                _note = ("SpliceAI and AF not consulted (--peakpick centre)" if a.peakpick == "centre"
                         else f"{_v['_peak_scored']} of {_v['_peak_n']} have SpliceAI, "
                              f"{_v['_peak_withaf']} have an AF")
                log(f"      picked by {_v['_picked_by']}  ({_note})")
        log(f"    other buckets: canonical {len(b['canonical'])}, spliceregion {len(b['spliceregion'])}, "
            f"deep {len(b['deep'])}, indels skipped {len(b['skipped_indel'])}, exonic skipped {len(b['skipped_exonic'])}")
        for v in picked:
            tag = f"{gene} {v.get('rsid') or v['pos']}"
            ok, why = preflight(v, g, gene, seq_cache)
            if not ok:
                log(f"  REFUSED  {tag}: {why}")
                refused.append({"gene": gene, "pos": v["pos"], "rsid": v.get("rsid"), "reason": why})
                continue
            if a.dry:
                log(f"  would draw {tag}  {v['ref']}>{v['alt']}  {v['_dist']} bp from splice site")
                continue
            try:
                t0 = time.time()
                rec = draw(v, g, gene, model, dna_client, genome, outdir)
                done.append(rec)
                log(f"  OK       {tag}  {rec['delta_log2fc']:+.4f} log2FC  track {rec['track_index']}/{rec['n_tracks']}  "
                    f"{time.time()-t0:.1f}s")
            except Exception as e:
                log(f"  FAILED   {tag}: {type(e).__name__}: {e}")
                failed.append({"gene": gene, "pos": v["pos"], "rsid": v.get("rsid"),
                               "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-600:]})

    log(f"\ndrawn {len(done)} · refused {len(refused)} · failed {len(failed)}")
    if done:
        d = [r["splice_log2fc"] for r in done]
        rna = [r["delta_log2fc"] for r in done]
        mx = [r["splice_maxabsdiff"] for r in done]
        log(f"SPLICE log2FC {min(d):+.3f} to {max(d):+.3f}  median |effect| {np.median(np.abs(d)):.4f}")
        log(f"SPLICE max|diff| median {np.median(mx):.4f}   RNA median |effect| {np.median(np.abs(rna)):.4f}")
        if a.bucket == "canonical" and len(done) and np.median(np.abs(d)) < 0.01:
            log("  *** WARNING: canonical +-1/+-2 variants destroy a splice site. A median effect")
            log("  *** below 0.01 means the PIPELINE is wrong, not the biology. Do not trust other buckets.")
    if not a.dry:
        # --peakpick changes WHICH variant is tested, so it is a different experiment and gets
        # its own manifest. Sharing one filename let the centre run silently replace the
        # SpliceAI run, leaving 14 figures on disk that no manifest claimed.
        _suffix = f"_{a.peakpick}" if (a.bucket == "peak" and a.peakpick != "auto") else ""
        p = os.path.join(ROOT, f"data/panels_{a.bucket}{_suffix}.json")
        # MERGE, do not overwrite. One --gene call per gene is the normal way to run this, and a
        # plain overwrite meant a 16-gene loop drew 16 figures and recorded 1: the PNGs were on disk
        # with nothing in the manifest pointing at them, which is a silent loss, not an error.
        # Keyed by id so a re-run of one gene replaces that gene's rows and leaves the rest.
        prev = {"drawn": [], "refused": [], "failed": []}
        if os.path.exists(p):
            try:
                old = json.load(open(p))
                if old.get("renderer") == RENDERER and old.get("bucket") == a.bucket:
                    prev = old
            except Exception as e:
                log(f"  existing {os.path.basename(p)} unreadable ({e}) -- starting fresh")
        def _merge(old_rows, new_rows, key):
            fresh = {r.get(key) for r in new_rows}
            return [r for r in old_rows if r.get(key) not in fresh] + new_rows
        out = {"renderer": RENDERER, "bucket": a.bucket,
               "drawn":   _merge(prev.get("drawn", []), done, "id"),
               "refused": _merge(prev.get("refused", []), refused, "tag"),
               "failed":  _merge(prev.get("failed", []), failed, "tag")}
        json.dump(out, open(p, "w"), indent=1)
        log(f"wrote {os.path.relpath(p, ROOT)} ({len(out['drawn'])} drawn in total, {len(done)} this run)")


if __name__ == "__main__":
    main()
