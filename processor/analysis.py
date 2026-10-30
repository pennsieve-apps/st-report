"""
What the report says: per section, the tables, the numbers, and the observations (plain
sentences drawn from the data) and caveats (what to keep in mind before trusting a number,
mostly the limitations each upstream app's README states). Figures are drawn in figures.py
from what this returns.

The text describes what the data show; it does not tell the reader what to do.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import PurePath
from typing import Any

import numpy as np
import pandas as pd

from processor.inputs import Inputs, Run

# ---------------------------------------------------------------- formatting


def pct(frac: Any, digits: int | None = None) -> str:
    """A 0-1 fraction as a percentage with sensible precision: 19%, 1.8%, 0.06%."""
    try:
        v = 100.0 * float(frac)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(v):
        return "–"
    if digits is not None:
        return f"{v:.{digits}f}%"
    a = abs(v)
    if a == 0:
        return "0%"
    if a >= 10:
        return f"{v:.0f}%"
    if a >= 1:
        return f"{v:.1f}%"
    if a >= 0.1:
        return f"{v:.2f}%"
    return f"{v:.1g}%" if a >= 0.001 else "<0.001%"


def pval(p: Any) -> str:
    """6e-93, 0.003, 0.04."""
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(p):
        return "–"
    if p == 0:
        return "<1e-300"
    if p < 0.001:
        return f"{p:.0e}".replace("e-0", "e-")
    return f"{p:.2g}"


def num(value: Any, digits: int = 2, signed: bool = False) -> str:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "–" if value is None else str(value)
    if math.isnan(f):
        return "–"
    return f"{f:{'+' if signed else ''}.{digits}f}"


def count(value: Any) -> str:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(f):
        return "–"
    return f"{int(round(f)):,}"


def sig(value: Any, digits: int = 2) -> str:
    """A number with `digits` significant figures, thousands separated: 0.18, 6,453, 13.1."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return "–"
    if math.isnan(f):
        return "–"
    if f == 0:
        return "0"
    if abs(f) >= 10 ** digits:
        return count(f)
    d = max(0, digits - 1 - int(math.floor(math.log10(abs(f)))))
    return f"{f:,.{d}f}"


def joined(items: list[Any], conj: str = "and") -> str:
    items = [str(i) for i in items]
    if len(items) <= 1:
        return "".join(items)
    return ", ".join(items[:-1]) + f" {conj} " + items[-1]


def plural(n: int, word: str, many: str | None = None) -> str:
    return f"{count(n)} {word if n == 1 else (many or word + 's')}"


def clean(obj: Any) -> Any:
    """JSON-safe: NaN -> None, numpy scalars -> Python, DataFrames -> records."""
    if isinstance(obj, pd.DataFrame):
        return clean(obj.to_dict(orient="records"))
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return None if math.isnan(f) or math.isinf(f) else f
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, PurePath):
        return str(obj)
    return obj


def truthy(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().lower() in ("true", "1", "yes", "t")


def unit_of(kind: str | None, n: int = 2) -> str:
    """'spots', 'bins', 'cells', 'observations'."""
    word = {"spot": "spot", "bin": "bin", "cell": "cell"}.get(str(kind or ""), "observation")
    return word if n == 1 else word + "s"


def section() -> dict[str, Any]:
    return {"observations": [], "caveats": []}


# ---------------------------------------------------------------- study design

REQUIRED = ("sample_id", "path", "platform", "kind", "n_obs", "n_genes", "n_obs_in_tissue")
PLATFORM_NAMES = {"visium": "Visium", "visium-hd": "Visium HD", "slide-seq": "Slide-seq", "xenium": "Xenium",
                  "merscope": "MERSCOPE", "cosmx": "CosMx", "h5ad": "AnnData", "spatialdata": "SpatialData"}
GROUPING_HINTS = ("donor", "patient", "subject", "individual", "participant", "batch", "slide", "block",
                  "replicate", "section")


def condition_column(inp: Inputs) -> str | None:
    """The grouping the comparison used (st-compare's CONDITION_COLUMN), else 'condition'."""
    for run in inp.many("compare"):
        col = (run.summary or {}).get("condition_column")
        if col and inp.samples is not None and col in inp.samples.columns:
            return str(col)
    if inp.samples is not None and "condition" in inp.samples.columns:
        return "condition"
    return None


def sample_meta(inp: Inputs) -> pd.DataFrame | None:
    if inp.samples is None:
        return None
    return inp.samples.set_index("sample_id")


def condition_of(inp: Inputs) -> dict[str, str]:
    col = condition_column(inp)
    if not col:
        return {}
    return {str(k): str(v) for k, v in zip(inp.samples["sample_id"], inp.samples[col])}


def kind_of(inp: Inputs) -> str | None:
    if inp.samples is not None and "kind" in inp.samples.columns:
        kinds = inp.samples["kind"].dropna().astype(str).unique().tolist()
        if len(kinds) == 1:
            return kinds[0]
    qc = inp.one("qc")
    if qc is not None:
        df = qc.csv("qc_summary.csv", "kind")
        if df is not None and df["kind"].nunique() == 1:
            return str(df["kind"].iloc[0])
    return None


def demo_datasets(inp: Inputs) -> list[dict[str, Any]]:
    from processor import citations as C

    paths: list[tuple[str, str, int]] = []
    if inp.samples is not None and "path" in inp.samples.columns:
        for _, r in inp.samples.iterrows():
            paths.append((str(r.get("path", "")), str(r.get("platform", "")), int(r.get("n_genes", 0) or 0)))
    elif inp.prep:
        for s in inp.prep.get("samples", {}).values():
            if isinstance(s, dict):
                paths.append((str(s.get("source", "")), str(s.get("platform", "")), int(s.get("n_genes", 0) or 0)))
    out = []
    for ds in C.DATASETS:
        hits = [p for p, plat, ng in paths if plat == ds["platform"] and any(m in p.split("/") or m in p for m in ds["match"])]
        if ds["key"] == "xenium":  # Rep1/Rep2 alone is a weak signal: the 313-gene breast panel too
            hits = [p for p, plat, ng in paths if plat == "xenium" and ng == 313 and p.split("/")[-1] in ds["match"]]
        if hits:
            d = dict(ds)
            d["n_samples"] = len(hits)
            d["reference"] = C.reference(ds["ref"]) if ds.get("ref") else None
            d["extra_reference"] = C.reference(ds["extra_ref"]) if ds.get("extra_ref") else None
            out.append(d)
    return out


def study(inp: Inputs) -> dict[str, Any] | None:
    s = inp.samples
    prep = inp.prep or {}
    if s is None and not prep:
        return None
    out = section()
    if s is None:  # prep.json only
        rows = [{"sample_id": k, **{c: v.get(c) for c in ("platform", "kind", "n_obs", "n_genes")}}
                for k, v in prep.get("samples", {}).items() if isinstance(v, dict)]
        s = pd.DataFrame(rows)
    n = len(s)
    platforms = [PLATFORM_NAMES.get(p, p) for p in s["platform"].astype(str).unique()] if "platform" in s else []
    kind = str(s["kind"].iloc[0]) if "kind" in s and s["kind"].nunique() == 1 else None
    word = unit_of(kind)
    extra = [c for c in s.columns if c not in REQUIRED]
    out.update(n_samples=n, platforms=platforms, kind=kind, unit=word, extra_columns=extra,
               samples=s.to_dict(orient="records"), columns=[c for c in s.columns if c != "path"])
    sentence = f"{plural(n, 'sample')}" + (f", {joined(platforms)}" if platforms else "")
    if "n_obs" in s:
        obs = pd.to_numeric(s["n_obs"], errors="coerce")
        out["n_obs"] = float(obs.sum())
        sentence += (f": {count(obs.sum())} {word} in all" +
                     (f" ({count(obs.min())}–{count(obs.max())} per sample)" if n > 1 else ""))
    if "n_genes" in s:
        g = pd.to_numeric(s["n_genes"], errors="coerce")
        sentence += (f", {count(g.iloc[0])} genes" if g.nunique() == 1 else f", {count(g.min())}–{count(g.max())} genes")
        out["n_genes"] = [float(g.min()), float(g.max())]
        if g.max() <= 1000:
            out["targeted_panel"] = True
            out["caveats"].append(f"A targeted panel of {count(g.max())} genes: anything the panel does not measure "
                                  "(cell types, pathways, ligand–receptor pairs) is invisible to every step.")
    out["observations"].append(sentence + ".")

    cond = condition_column(inp) if inp.samples is not None else None
    out["condition_column"] = cond
    if cond:
        counts = s[cond].astype(str).value_counts().sort_index()
        out["conditions"] = {k: int(v) for k, v in counts.items()}
        out["observations"].append(f"Conditions (column '{cond}'): " +
                                   joined([f"{k} ({v} sample{'s' if v != 1 else ''})" for k, v in counts.items()]) + ".")
        if len(counts) > 1 and counts.min() < 3:
            out["caveats"].append(
                f"{counts.idxmin()} has only {counts.min()} sample(s). Tests between conditions use the sample as "
                "the replicate, so with this few sections only large, consistent differences can be detected.")
        # Columns that carry the same grouping, or nest in it (donor, slide, batch ...).
        for c in extra:
            if c == cond or s[c].isna().all():
                continue
            pairs = s[[c, cond]].astype(str).drop_duplicates()
            if s[c].nunique() == counts.size and len(pairs) == counts.size:
                out["caveats"].append(f"'{c}' and '{cond}' group the samples identically: a difference between "
                                      f"the {cond} levels is also a difference in {c}.")
            elif (any(h in c.lower() for h in GROUPING_HINTS) and s[c].nunique() > 1
                  and pairs.groupby(c).size().max() == 1 and len(counts) > 1):
                out["caveats"].append(f"Every '{c}' value occurs in one {cond} level only, so {cond} and {c} are "
                                      f"confounded: a difference between conditions can be a difference between {c}s.")
    if kind == "spot":
        out["caveats"].append("Visium spots (55 µm) each hold several cells (typically 1–10): expression, domains "
                              "and spatially variable genes describe spot mixtures, not single cells.")
    # prep.json details
    samples_prep = {k: v for k, v in (prep.get("samples") or {}).items() if isinstance(v, dict)}
    if samples_prep:
        sub = [k for k, v in samples_prep.items() if (v.get("observations") or {}).get("n_subsampled_to")]
        if sub:
            out["caveats"].append(
                f"{joined(sub)} {'was' if len(sub) == 1 else 'were'} subsampled by st-prep (MAX_OBS_PER_SAMPLE): "
                "neighbourhood graphs are sparser and counts smaller than on the full data; use such a run for a "
                "quick look, not for final numbers.")
            out["subsampled"] = sub
        ann = [(k, (v.get("annotation") or {}).get("n_annotated")) for k, v in samples_prep.items()]
        ann = [(k, a) for k, a in ann if a]
        if ann and "n_obs" in s:
            total = sum(a for _, a in ann)
            out["observations"].append(f"A user annotation (annotations.csv) labels {count(total)} {word} in "
                                       f"{plural(len(ann), 'sample')}.")
        dropped: dict[str, int] = {}
        for v in samples_prep.values():
            for typ, d in ((v.get("features") or {}).get("dropped") or {}).items():
                if isinstance(d, dict):
                    dropped[typ] = max(dropped.get(typ, 0), int(d.get("n_features") or 0))
        if dropped:
            out["observations"].append("Control features summed into control_counts and dropped from the panel: " +
                                       joined([f"{n} {t}" for t, n in dropped.items()]) + ".")
        warns = [(k, w) for k, v in samples_prep.items() for w in (v.get("warnings") or [])]
        out["warnings"] = [{"sample_id": k, "warning": w} for k, w in warns]
        units = {v.get("coordinate_unit") for v in samples_prep.values() if v.get("coordinate_unit")}
        if units:
            out["coordinate_unit"] = joined(sorted(units))
        out["images"] = sum(1 for v in samples_prep.values() if (v.get("image") or {}).get("lowres"))
        skipped = prep.get("skipped") or []
        if skipped:
            out["caveats"].append(f"st-prep skipped {plural(len(skipped), 'entry', 'entries')} of the input "
                                  "(sheet rows matching nothing, duplicates or failed samples); see prep.json.")
    return out


# ---------------------------------------------------------------- QC


def qc(inp: Inputs) -> dict[str, Any] | None:
    run = inp.one("qc")
    if run is None:
        return None
    df = run.csv("qc_summary.csv", "sample_id")
    if df is None:
        return None
    out = section()
    df = df.copy()
    df["sample_id"] = df["sample_id"].astype(str)
    cond = condition_of(inp)
    df["condition"] = df["sample_id"].map(cond) if cond else None
    kind = str(df["kind"].iloc[0]) if "kind" in df and df["kind"].nunique() == 1 else kind_of(inp)
    word = unit_of(kind)
    out["unit"] = word
    ok = df[df["status"].astype(str) == "ok"] if "status" in df else df
    failed = df[df["status"].astype(str) != "ok"] if "status" in df else df.iloc[0:0]
    before, after = float(ok["n_obs_before"].sum()), float(ok["n_obs_after"].sum())
    removed = before - after
    out.update(rows=df.to_dict(orient="records"), n_before=before, n_after=after, n_removed=removed,
               filter_mode=str(df["filter_mode"].iloc[0]) if "filter_mode" in df else None)
    pr = pd.to_numeric(ok["pct_removed"], errors="coerce")
    out["observations"].append(
        f"QC removed {count(removed)} of {count(before)} {word} ({pct(removed / before if before else 0)})" +
        (f"; per sample {pct(pr.min() / 100)} ({ok.loc[pr.idxmin(), 'sample_id']}) to {pct(pr.max() / 100)} "
         f"({ok.loc[pr.idxmax(), 'sample_id']})" if len(ok) > 1 else "") + ".")
    reasons = {"low counts": "n_removed_low_counts", "few genes": "n_removed_low_genes",
               "high mitochondrial %": "n_removed_mt", "high control-probe fraction": "n_removed_control"}
    by = {k: int(pd.to_numeric(ok[c], errors="coerce").fillna(0).sum()) for k, c in reasons.items() if c in ok}
    out["by_reason"] = by
    nz = {k: v for k, v in by.items() if v}
    if nz:
        out["observations"].append("Failing each test: " + joined([f"{k} {count(v)}" for k, v in nz.items()]) +
                                   " (one observation can fail several).")
    if len(ok) > 2:
        med = float(pr.median())
        high = ok[(pr > max(2 * med, med + 2)) & (pr > 3)]
        for _, r in high.iterrows():
            out["observations"].append(f"{r['sample_id']} lost more than the others: {pct(r['pct_removed'] / 100)} "
                                       f"against a median of {pct(med / 100)}.")
    if "median_counts" in ok and len(ok) > 1:
        mc = pd.to_numeric(ok["median_counts"], errors="coerce")
        if mc.min() > 0:
            fold = mc.max() / mc.min()
            out["depth_fold"] = float(fold)
            text = (f"Median counts per {unit_of(kind, 1)} after QC range from {sig(mc.min(), 3)} "
                    f"({ok.loc[mc.idxmin(), 'sample_id']}) to {sig(mc.max(), 3)} ({ok.loc[mc.idxmax(), 'sample_id']})")
            if fold >= 3:
                out["observations"].append(text + f", a {fold:.0f}-fold difference in depth.")
                out["caveats"].append("Samples differ several-fold in sequencing depth or capture; normalization scales "
                                      "library size per observation but does not equalize detection of weakly "
                                      "expressed genes. Weigh depth before reading differences between samples.")
            else:
                out["observations"].append(text + ".")
    if "median_pct_mt" in ok:
        mt = pd.to_numeric(ok["median_pct_mt"], errors="coerce").fillna(0)
        if (mt == 0).all():
            out["observations"].append("Median mitochondrial % is 0 in every sample: the panel or probe set has no "
                                       "(or almost no) mitochondrial genes, so only the count, gene and control cuts "
                                       "could act.")
    if "n_genes_before" in ok and "n_genes_after" in ok:
        gd = pd.to_numeric(ok["n_genes_before"]) - pd.to_numeric(ok["n_genes_after"])
        if gd.max() > 0:
            out["observations"].append(f"Genes detected in too few observations were dropped per sample "
                                       f"({count(gd.min())}–{count(gd.max())}), so samples keep different gene sets; "
                                       "steps that combine samples use the genes shared by all.")
    for _, r in failed.iterrows():
        out["observations"].append(f"{r['sample_id']} failed QC and was skipped: {r['status']}.")
    if out["filter_mode"] == "mad":
        mads = (inp.one("qc").provenance or {}).get("parameters", {}).get("MADS", 5)
        out["caveats"].append(f"Thresholds adapt to each sample ({mads} MADs from its own median on the log scale): a "
                              "sample that is poor throughout keeps most of its observations. Its medians in the "
                              "table below are where to notice it.")
    if kind == "cell":
        out["caveats"].append("Segmentation errors (two cells merged, fragments) are not detected by QC; low-count "
                              "small cells such as lymphocytes are kept unless they fall below the thresholds.")
    th = run.csv("thresholds.csv", "sample_id", "metric", "threshold")
    out["thresholds"] = th.to_dict(orient="records") if th is not None else []
    out["plots"] = [p for p in run.pngs() if p.name.endswith(("_qc.png", "_spatial_counts.png"))]
    return out


# ---------------------------------------------------------------- domains


def domains(inp: Inputs, max_items: int) -> dict[str, Any] | None:
    run = inp.one("domains")
    if run is None:
        return None
    ds = run.csv("domain_summary.csv", "sample_id", "domain", "n_obs")
    sm = run.summary or {}
    if ds is None and not sm:
        return None
    out = section()
    word = unit_of(kind_of(inp))
    cond = condition_of(inp)
    min_obs = int(sm.get("min_obs_present") or 20)
    if ds is not None:
        ds["sample_id"] = ds["sample_id"].astype(str)
        ds["domain"] = ds["domain"].astype(str)
        if "fraction" not in ds:
            ds["fraction"] = ds["n_obs"] / ds.groupby("sample_id")["n_obs"].transform("sum")
        size = ds.groupby("domain")["n_obs"].sum().sort_values(ascending=False)
        total = float(size.sum())
        present = ds[ds["n_obs"] >= min_obs].groupby("domain")["sample_id"].apply(list)
        samples = list(dict.fromkeys(ds["sample_id"]))
        rows = []
        for d, n in size.items():
            ps = present.get(d, [])
            conds = sorted({cond.get(x, "") for x in ps} - {""})
            rows.append({"domain": d, "n_obs": int(n), "fraction": n / total if total else None,
                         "n_samples": len(ps), "samples": ps, "conditions": conds})
        out["coverage"] = ds.pivot_table(index="domain", columns="sample_id", values="fraction", fill_value=0) \
            .reindex(index=list(size.index), columns=samples)
        out["coverage_counts"] = ds.pivot_table(index="domain", columns="sample_id", values="n_obs", fill_value=0) \
            .reindex(index=list(size.index), columns=samples)
    else:
        sizes = sm.get("domain_sizes") or {}
        total = float(sum(sizes.values()))
        cov = sm.get("domain_sample_coverage") or {}
        rows = [{"domain": d, "n_obs": n, "fraction": n / total if total else None, "n_samples": cov.get(d),
                 "samples": [], "conditions": []} for d, n in sizes.items()]
        samples = list((sm.get("per_sample") or {}).keys())
    n_samples = len(samples)
    markers = run.csv("domain_markers.csv", "domain", "gene")
    top: dict[str, list[str]] = {}
    if markers is not None:
        markers["domain"] = markers["domain"].astype(str)
        sort_cols = ["domain", "rank"] if "rank" in markers else ["domain"]
        for d, g in markers.sort_values(sort_cols).groupby("domain", sort=False):
            top[d] = g["gene"].astype(str).head(5).tolist()
    elif sm.get("top_markers"):
        top = {k: list(v)[:5] for k, v in sm["top_markers"].items()}
    for r in rows:
        r["markers"] = top.get(r["domain"], [])
    out.update(rows=rows, n_domains=len(rows), samples=samples, markers=top, unit=word, min_obs=min_obs,
               lambda_=sm.get("lambda"), resolution=sm.get("resolution"), integration=sm.get("integration"),
               n_genes_banksy=sm.get("n_genes_banksy"), gene_selection=sm.get("gene_selection"))
    shared = [r for r in rows if (r["n_samples"] or 0) >= 2]
    in_all = [r for r in rows if n_samples > 1 and (r["n_samples"] or 0) == n_samples]
    single = [r["domain"] for r in rows if (r["n_samples"] or 0) == 1]
    none = [r["domain"] for r in rows if (r["n_samples"] or 0) == 0]
    text = f"BANKSY found {plural(len(rows), 'domain')}"
    if sm.get("lambda") is not None:
        text += f" (λ = {sm['lambda']}, " + (f"resolution {sm['resolution']}, " if sm.get("resolution") else "") + \
            f"{sm.get('integration') or 'no'} integration)"
    if n_samples > 1:
        text += (f"; {len(shared)} are present (≥ {min_obs} {word}) in at least two samples" +
                 (f", {len(in_all)} in all {n_samples}" if in_all else ""))
    out["observations"].append(text + ".")
    if single:
        out["observations"].append(f"Confined to a single sample: {joined(single)}" +
                                   (f"; {joined(none)} {'has' if len(none) == 1 else 'have'} fewer than {min_obs} "
                                    f"{word} in every sample" if none else "") + ".")
    elif none:
        out["observations"].append(f"{joined(none)} {'has' if len(none) == 1 else 'have'} fewer than {min_obs} "
                                   f"{word} in every sample.")
    for r in rows[:3]:
        out["observations"].append(f"{r['domain']} covers {pct(r['fraction'])} of all {word}" +
                                   (f" (markers {', '.join(r['markers'][:4])})" if r["markers"] else "") + ".")
    levels = sorted(set(cond.values())) if cond else []
    if len(levels) > 1:
        only = [(r["domain"], r["conditions"][0], r["samples"]) for r in rows
                if len(r["conditions"]) == 1 and r["n_samples"] >= 1]
        if only:
            out["observations"].append("Present in one condition only: " + joined(
                [f"{d} ({c}: {', '.join(ss)})" for d, c, ss in only[:8]]) + ".")
            out["condition_specific"] = [{"domain": d, "condition": c, "samples": ss} for d, c, ss in only]
    ann = sm.get("annotation") or {}
    ag = run.csv("agreement.csv", "sample_id", "ari", "nmi")
    if ann or ag is not None:
        per = []
        if ag is not None:
            ag["sample_id"] = ag["sample_id"].astype(str)
            pooled = ag[ag["sample_id"] == "all"]
            per = ag[ag["sample_id"] != "all"].to_dict(orient="records")
            ari = float(pooled["ari"].iloc[0]) if len(pooled) else ann.get("ari")
            nmi = float(pooled["nmi"].iloc[0]) if len(pooled) else ann.get("nmi")
        else:
            ari, nmi = ann.get("ari"), ann.get("nmi")
            per = [{"sample_id": k, **v} for k, v in (ann.get("per_sample") or {}).items()]
        out["agreement"] = {"ari": ari, "nmi": nmi, "per_sample": per, "ignored": ann.get("ignored") or [],
                            "n_obs_used": ann.get("n_obs_used")}
        text = f"Agreement with the user annotation: ARI {num(ari)}, NMI {num(nmi)} pooled"
        if per:
            p = pd.DataFrame(per)
            lo, hi = p.loc[p["ari"].idxmin()], p.loc[p["ari"].idxmax()]
            text += f"; per sample ARI {num(lo['ari'])} ({lo['sample_id']}) to {num(hi['ari'])} ({hi['sample_id']})"
        if ann.get("n_obs_used"):
            text += f", over {count(ann['n_obs_used'])} annotated {word}"
            if ann.get("ignored"):
                text += f" ({joined(ann['ignored'])} left out)"
        out["observations"].append(text + ".")
        bm = run.csv("annotation_match.csv", "annotation", "domain", "precision", "recall")
        best = bm.to_dict(orient="records") if bm is not None else (ann.get("best_match") or [])
        out["best_match"] = best
        good = [b for b in best if (b.get("f1") or 0) >= 0.3]
        weak = [b for b in best if (b.get("f1") or 0) < 0.3]
        if good:
            out["observations"].append("Best-matching domains: " + joined(
                [f"{b['annotation']} → {b['domain']} (precision {num(b['precision'])}, recall {num(b['recall'])})"
                 for b in sorted(good, key=lambda b: -(b.get("f1") or 0))]) + ".")
        if weak:
            out["observations"].append("No domain matches " + joined([b["annotation"] for b in weak], "or") +
                                       " well (best F1 below 0.3).")
        out["caveats"].append("ARI rewards getting the large classes right, not small structures; a modest ARI is "
                              "expected when the annotation has a few coarse classes and many unannotated "
                              f"{word}.")
    out["caveats"].append("Domains are clusters, not anatomy: BANKSY makes clusters spatially smooth but does not "
                          "know what a structure is. One region can split over several domains and a domain can mix "
                          "tissues that share a neighbourhood profile; look at the maps and markers before naming them.")
    if sm.get("integration") == "harmony" and n_samples > 1:
        out["caveats"].append("Harmony removed differences between samples before clustering, including real ones: "
                              "a region that exists in only some sections, or that a condition changes, can be pulled "
                              "into another domain.")
    if sm.get("lambda") is not None and float(sm["lambda"]) >= 0.5:
        out["caveats"].append(f"With λ = {sm['lambda']} the neighbourhood dominates: structures only a few "
                              f"{word} across are smoothed into their surroundings, and boundaries between regions "
                              "can form thin domains of their own.")
    fig_dir = run.dir / "figures"
    out["maps"] = [p for p in run.pngs("figures") if p.name.endswith("_domains.png")]
    out["upstream_figures"] = {n: fig_dir / f"{n}.png" for n in ("umap", "markers_dotplot", "confusion",
                                                                 "domain_composition") if (fig_dir / f"{n}.png").is_file()}
    out["markers_table"] = markers
    return out


# ---------------------------------------------------------------- spatially variable genes


def marker_index(inp: Inputs) -> dict[str, list[str]]:
    """gene -> the domains / clusters it marks (top 10 markers of each), for cross-references."""
    idx: dict[str, list[str]] = {}
    d = inp.one("domains")
    if d is not None:
        m = d.csv("domain_markers.csv", "domain", "gene")
        if m is not None:
            m = m[m["rank"] <= 10] if "rank" in m else m.groupby("domain").head(10)
            for g, dom in zip(m["gene"].astype(str), m["domain"].astype(str)):
                idx.setdefault(g, []).append(dom)
    c = inp.one("celltype")
    if c is not None:
        m = c.csv("cluster_markers.csv", "cluster", "gene")
        ann = c.csv("cluster_annotation.csv", "cluster", "cell_type")
        names = dict(zip(ann["cluster"].astype(str), ann["cell_type"].astype(str))) if ann is not None else {}
        if m is not None:
            m = m[m["rank"] <= 10] if "rank" in m else m.groupby("cluster").head(10)
            for g, cl in zip(m["gene"].astype(str), m["cluster"].astype(str)):
                idx.setdefault(g, []).append(f"{cl} {names.get(cl, '')}".strip())
    return idx


def svg(inp: Inputs, fdr: float, max_items: int) -> dict[str, Any] | None:
    run = inp.one("svg")
    if run is None:
        return None
    sm = run.summary or {}
    cons = run.csv("svg_consensus.csv", "gene")
    per_gene = run.csv("svg.csv", "sample_id", "gene", "morans_i")
    if cons is None and per_gene is None and not sm:
        return None
    out = section()
    word = unit_of(kind_of(inp))
    up_fdr = float(sm.get("fdr") or (run.provenance or {}).get("parameters", {}).get("FDR") or 0.05)
    out["fdr"], out["upstream_fdr"] = fdr, up_fdr
    rows = []
    if per_gene is not None:
        per_gene["sample_id"] = per_gene["sample_id"].astype(str)
        if "pval_adj" in per_gene:
            per_gene["sig"] = pd.to_numeric(per_gene["pval_adj"], errors="coerce") < fdr
        else:
            per_gene["sig"] = per_gene["significant"].map(truthy)
        for sid, g in per_gene.groupby("sample_id", sort=False):
            ps = (sm.get("per_sample") or {}).get(sid, {})
            top = g.sort_values("morans_i", ascending=False)
            rows.append({"sample_id": sid, "n_obs": ps.get("n_obs"), "n_tested": len(g), "n_significant": int(g["sig"].sum()),
                         "fraction_significant": g["sig"].mean(), "median_morans_i": float(g["morans_i"].median()),
                         "max_morans_i": float(g["morans_i"].max()),
                         "top_genes": top["gene"].astype(str).head(6).tolist(),
                         "graph_mean_degree": ps.get("graph_mean_degree"), "condition": ps.get("condition")})
    else:
        for sid, ps in (sm.get("per_sample") or {}).items():
            rows.append({"sample_id": sid, "n_obs": ps.get("n_obs"), "n_tested": ps.get("n_genes_tested"),
                         "n_significant": ps.get("n_significant"), "fraction_significant": ps.get("fraction_significant"),
                         "median_morans_i": ps.get("median_morans_i"), "max_morans_i": None,
                         "top_genes": (ps.get("top_genes") or [])[:6], "graph_mean_degree": ps.get("graph_mean_degree"),
                         "condition": ps.get("condition")})
    out["per_sample"] = rows
    idx = marker_index(inp)
    consensus = []
    if cons is not None:
        if per_gene is not None and fdr != up_fdr:
            cons = cons.drop(columns=[c for c in ("n_significant",) if c in cons]).merge(
                per_gene.groupby("gene")["sig"].sum().rename("n_significant").reset_index(), on="gene", how="left")
            cons = cons.sort_values(["n_significant", "mean_morans_i"], ascending=False).reset_index(drop=True)
            cons["rank"] = np.arange(1, len(cons) + 1)
        for _, r in cons.head(max_items).iterrows():
            consensus.append({"gene": r["gene"], "n_samples": r.get("n_samples"), "n_significant": r.get("n_significant"),
                              "mean_morans_i": r.get("mean_morans_i"), "median_rank": r.get("median_rank"),
                              "rank": r.get("rank"), "markers": idx.get(str(r["gene"]), [])})
    out["consensus"] = consensus
    n_samp = len(rows)
    if per_gene is not None:
        n_tested = per_gene["gene"].nunique()
        sig_any = per_gene[per_gene["sig"]]["gene"].nunique()
        sig_all = int((per_gene.groupby("gene")["sig"].sum() == n_samp).sum()) if n_samp > 1 else sig_any
        out.update(n_genes_tested=n_tested, n_significant_any=sig_any, n_significant_all=sig_all)
    else:
        out.update(n_genes_tested=sm.get("n_genes_tested"), n_significant_any=sm.get("n_genes_significant_any"),
                   n_significant_all=sm.get("n_genes_significant_all"))
    stat = "Moran's I" if str(sm.get("statistic", "moran")).startswith("moran") else str(sm.get("statistic"))
    if rows:
        r = pd.DataFrame(rows)
        if n_samp > 1:
            lo, hi = r.loc[r["n_significant"].idxmin()], r.loc[r["n_significant"].idxmax()]
            out["observations"].append(
                f"Spatially variable genes ({stat}, FDR {fdr:g}) per sample: {count(lo['n_significant'])} "
                f"({lo['sample_id']}) to {count(hi['n_significant'])} ({hi['sample_id']}) of the "
                f"{count(r['n_tested'].median())} genes tested in each.")
        else:
            out["observations"].append(f"{count(r['n_significant'].iloc[0])} of {count(r['n_tested'].iloc[0])} genes "
                                       f"tested are spatially variable ({stat}, FDR {fdr:g}).")
    if n_samp > 1 and out.get("n_genes_tested"):
        out["observations"].append(f"{count(out['n_significant_any'])} of {count(out['n_genes_tested'])} genes tested "
                                   f"are significant in at least one sample, {count(out['n_significant_all'])} in "
                                   + ("both." if n_samp == 2 else f"all {n_samp}."))
    if consensus:
        top = consensus[:8]
        out["observations"].append("Most consistently spatially variable: " + joined(
            [f"{c['gene']} ({count(c['n_significant'])}/{count(c['n_samples'])}, mean I {num(c['mean_morans_i'])})"
             for c in top]) + " (significant in / tested in samples).")
        marked = [c for c in consensus[:max_items] if c["markers"]]
        if marked:
            out["observations"].append(
                f"{len(marked)} of the top {len(consensus)} are also top-10 markers of a domain or cluster, e.g. " +
                joined([f"{c['gene']} ({', '.join(c['markers'][:2])})" for c in marked[:4]]) +
                ": their spatial pattern largely follows where that region or cell type lies.")
    if per_gene is not None and per_gene["sig"].any():
        sigs = per_gene[per_gene["sig"]]
        small = float((sigs["morans_i"] < 0.05).mean())
        out["small_effect_fraction"] = small
        if small > 0.1:
            out["observations"].append(f"{pct(small)} of the significant gene × sample results have I < 0.05: "
                                       "statistically detectable, but a weak pattern.")
    out["caveats"].append(f"Significance is not effect size: with thousands of {word} the analytic test calls genes "
                          "with I of 0.02–0.03 significant. Read I and the rank, not only the count of significant genes.")
    if kind_of(inp) in ("spot", "bin"):
        out["caveats"].append("In spot data a gene is spatially variable mostly because the cells expressing it are "
                              "arranged in space (immunoglobulins follow plasma cells); this does not separate "
                              "composition from regulation within a cell type.")
    out["caveats"].append(f"{stat} depends on the neighbour graph and density: compare values between samples only "
                          "when they were built the same way. It measures one global pattern; a gene patterned in a "
                          "small region of a large section scores low.")
    if per_gene is not None and per_gene.groupby("sample_id")["gene"].apply(frozenset).nunique() > 1:
        out["caveats"].append("Each sample tested its own highly variable genes, so gene sets differ; 'tested in' "
                              "says in how many samples a gene was tested.")
    out["upstream_figures"] = {n: run.dir / f"{n}.png" for n in ("consensus_top", "morans_distribution")
                               if (run.dir / f"{n}.png").is_file()}
    out["maps"] = [p for p in run.pngs() if p.name.endswith("_svg_top.png")]
    return out


# ---------------------------------------------------------------- condition comparison


def min_rank_p(n1: int, n2: int) -> float:
    """The smallest two-sided Mann-Whitney p-value possible with n1 vs n2 samples (exact)."""
    if n1 < 1 or n2 < 1:
        return 1.0
    return min(1.0, 2.0 / math.comb(n1 + n2, n1))


def compare_one(run: Run, inp: Inputs, fdr: float, max_items: int) -> dict[str, Any] | None:
    sm = run.summary or {}
    ct = run.csv("composition_tests.csv", "label")
    dsum = run.csv("de_summary.csv", "label", "status")
    if ct is None and dsum is None:
        return None
    out = section()
    up_fdr = float((run.provenance or {}).get("parameters", {}).get("FDR") or 0.05)
    ref = sm.get("reference")
    levels = [str(l) for l in (sm.get("levels") or [])]
    others = [l for l in levels if l != ref]
    per_level = {str(k): int(v) for k, v in (sm.get("n_samples_per_level") or {}).items()}
    label_col = sm.get("label_column") or (run.provenance or {}).get("parameters", {}).get("LABEL_COLUMN") or run.label
    cond_col = sm.get("condition_column")
    out.update(label=run.label, label_column=label_col, condition_column=cond_col, reference=ref, levels=levels,
               per_level=per_level, fdr=fdr, upstream_fdr=up_fdr, block=sm.get("block_column"))
    design = joined([f"{k} {v}" for k, v in per_level.items()])
    out["design"] = design
    word = unit_of(kind_of(inp))

    # composition
    if ct is not None:
        ct = ct.copy()
        ct["label"] = ct["label"].astype(str)
        tested = ct["tested"].map(truthy) if "tested" in ct else pd.Series(True, index=ct.index)
        padj = pd.to_numeric(ct.get("parametric_p_adj"), errors="coerce")
        ct["significant_report"] = tested & (padj < fdr)
        rows = []
        for _, r in ct.sort_values("parametric_p" if "parametric_p" in ct else "label").iterrows():
            row = {"label": r["label"], "comparison": r.get("comparison"), "clr_difference": r.get("clr_difference"),
                   "p": r.get("parametric_p"), "p_adj": r.get("parametric_p_adj"), "rank_p_adj": r.get("rank_p_adj"),
                   "tested": bool(tested.loc[_]), "significant": bool(r["significant_report"]),
                   "means": {l: r.get(f"mean_fraction_{l}") for l in levels}}
            rows.append(row)
        out["composition"] = rows
        n_sig = int(ct["significant_report"].sum())
        n_tested = int(tested.sum())
        out["composition_n_significant"], out["composition_n_tested"] = n_sig, n_tested
        comp = (f"{others[0]} vs {ref}" if len(others) == 1 else f"across {joined(levels)}") if ref else "between conditions"
        if n_sig:
            sig_rows = [r for r in rows if r["significant"]]
            out["observations"].append(
                f"Composition ({comp}): {plural(n_sig, 'label')} of {n_tested} differ at FDR {fdr:g}: " +
                joined([f"{r['label']} " + (f"{pct(r['means'].get(ref))} → {pct(r['means'].get(others[0]))}"
                                             if len(others) == 1 else "") + f" (adj. p {pval(r['p_adj'])})"
                        for r in sig_rows[:6]]) + ".")
        elif n_tested:
            best = min((r for r in rows if r["tested"]), key=lambda r: r["p_adj"] if r["p_adj"] == r["p_adj"] else 9)
            why = []
            if len(others) == 1 and ref in per_level and others[0] in per_level:
                n1, n2 = per_level[ref], per_level[others[0]]
                why.append(f"with {n1} vs {n2} samples the rank test cannot go below p = {pval(min_rank_p(n1, n2))} "
                           f"before correction over {n_tested} labels")
                if min(n1, n2) < 5:
                    why.append("and the t-test has few degrees of freedom")
            out["observations"].append(
                f"Composition ({comp}): no label's share differs significantly at FDR {fdr:g}. The closest is "
                f"{best['label']}" + (f" ({pct(best['means'].get(ref))} of {word} in {ref} vs "
                                       f"{pct(best['means'].get(others[0]))} in {others[0]}" if len(others) == 1 else " (")
                + f"; adj. p {pval(best['p_adj'])})." + (" Not significant here is not evidence of no change: " +
                                                       " ".join(why) + "." if why else ""))
            out["nothing_significant"] = True
        else:
            out["observations"].append("Composition: no label could be tested (every level needs at least two samples).")
    # DE
    if dsum is not None:
        de = run.csv("de.csv", "label", "gene", "p_adj")
        dsum = dsum.copy()
        dsum["label"] = dsum["label"].astype(str)
        if de is not None:
            de["label"] = de["label"].astype(str)
            de["p_adj"] = pd.to_numeric(de["p_adj"], errors="coerce")
            de["log2_fold_change"] = pd.to_numeric(de.get("log2_fold_change"), errors="coerce")
            de["sig"] = de["p_adj"] < fdr
            if fdr != up_fdr:  # recount at the report's FDR
                g = de[de["sig"]].groupby(["label", "comparison"])
                n = g.size().rename("n_significant")
                up = de[de["sig"] & (de["log2_fold_change"] > 0)].groupby(["label", "comparison"]).size().rename("n_up")
                dn = de[de["sig"] & (de["log2_fold_change"] < 0)].groupby(["label", "comparison"]).size().rename("n_down")
                cnt = pd.concat([n, up, dn], axis=1).fillna(0).reset_index()
                dsum = dsum.drop(columns=[c for c in ("n_significant", "n_up", "n_down") if c in dsum]).merge(
                    cnt, on=["label", "comparison"], how="left")
                tested_mask = dsum["status"] == "tested"
                for c in ("n_significant", "n_up", "n_down"):
                    dsum.loc[tested_mask, c] = dsum.loc[tested_mask, c].fillna(0)
        rows = dsum.to_dict(orient="records")
        out["de_summary"] = rows
        tested_rows = [r for r in rows if r["status"] == "tested"]
        skipped = [r for r in rows if r["status"] != "tested"]
        total = int(sum((r.get("n_significant") or 0) for r in tested_rows if r.get("n_significant") == r.get("n_significant")))
        out["de_n_significant"] = total
        out["de_n_tested"] = len(tested_rows)
        if tested_rows:
            most = max(tested_rows, key=lambda r: r.get("n_significant") or 0)
            comp = most.get("comparison") or ""
            out["observations"].append(
                f"Pseudobulk DE ({comp}): {count(total)} significant gene × label results at FDR {fdr:g} over "
                f"{plural(len(tested_rows), 'tested label')}; most in {most['label']} ({count(most.get('n_significant'))}: "
                f"{count(most.get('n_up'))} up, {count(most.get('n_down'))} down in {others[0] if len(others) == 1 else 'the level'}).")
            zero = [r["label"] for r in tested_rows if not (r.get("n_significant") or 0)]
            if zero:
                out["observations"].append(f"No significant gene in {joined(zero)}.")
        if skipped:
            reasons: dict[str, list[str]] = {}
            for r in skipped:
                reason = str(r.get("reason") or r["status"])
                reason = re.sub(r"has \d+ pseudobulks \(< MIN_SAMPLES_PER_GROUP=(\d+)\)",
                                r"has fewer than \1 pseudobulks", reason)
                reasons.setdefault(reason, []).append(r["label"])
            out["observations"].append(f"Not tested ({plural(len(skipped), 'label')}): " + "; ".join(
                f"{joined(ls)} ({reason})" for reason, ls in reasons.items()) + ".")
        if de is not None and de["sig"].any():
            s = de[de["sig"]].copy()
            s["abs"] = s["log2_fold_change"].abs()
            top = s.sort_values(["p_adj", "abs"], ascending=[True, False]).head(max_items)
            out["top_genes"] = top[[c for c in ("label", "comparison", "gene", "log2_fold_change", "lfc_se", "base_mean",
                                                "p_adj") if c in top]].to_dict(orient="records")
            # genes changing the same way in most tested labels: an organ- or condition-wide shift
            n_lab = len(tested_rows)
            if n_lab >= 3:
                s["dir"] = np.sign(s["log2_fold_change"])
                shared = s.groupby(["gene", "dir"])["label"].nunique().reset_index()
                shared = shared[shared["label"] >= max(3, math.ceil(n_lab * 0.6))].sort_values("label", ascending=False)
                if len(shared):
                    out["shared_genes"] = shared.head(max_items).rename(columns={"label": "n_labels"}).to_dict(orient="records")
                    ups = shared[shared["dir"] > 0]["gene"].head(6).tolist()
                    dns = shared[shared["dir"] < 0]["gene"].head(6).tolist()
                    out["observations"].append(
                        f"{plural(len(shared), 'gene')} change the same way in at least {max(3, math.ceil(n_lab * 0.6))} of "
                        f"the {n_lab} tested labels" + (f" (up: {', '.join(ups)}" if ups else " (") +
                        (f"; down: {', '.join(dns)})" if dns else ")") +
                        ": a shift of the whole tissue between conditions rather than of one region or cell type.")
    # caveats
    if per_level and min(per_level.values()) < 5:
        out["caveats"].append(f"The sample is the replicate ({design} samples). With this few, only large and "
                              "consistent changes are detectable, and the rank test cannot reach significance in small "
                              "designs.")
    out["caveats"].append(f"Each {cond_col or 'condition'} level is a different set of sections; anything else that "
                          "differs between them (organ, tumour type, donor, processing) is part of the comparison "
                          + ("unless it is the block." if out["block"] else "(no BLOCK_COLUMN was given)."))
    if kind_of(inp) == "spot":
        out["caveats"].append(f"Visium pseudobulks mix the cell types of each {label_col}: a DE gene can reflect a "
                              "change in the cell-type mix within the label rather than in expression.")
    out["caveats"].append("Log fold changes are not shrunken and the FDR is controlled within each label, not across "
                          "labels; composition tests are CLR + t-test per label (scCODA is the fuller model).")
    out["figures_upstream"] = {n: run.dir / f"{n}.png" for n in ("composition", "de_heatmap") if (run.dir / f"{n}.png").is_file()}
    out["volcanoes"] = [p for p in run.pngs() if p.name.startswith("volcano_")]
    return out


def compare(inp: Inputs, fdr: float, max_items: int) -> dict[str, Any] | None:
    runs = [r for r in (compare_one(run, inp, fdr, max_items) for run in inp.many("compare")) if r]
    # the workflow's own labels (domains, cell types) before a user annotation
    runs.sort(key=lambda r: (str(r["label_column"]) == "annotation", str(r["label"])))
    if not runs:
        return None
    return {"runs": runs, "observations": [], "caveats": []}


# ---------------------------------------------------------------- cell types (st-celltype)


def celltypes(inp: Inputs, max_items: int) -> dict[str, Any] | None:
    run = inp.one("celltype")
    figs = inp.one("celltype-figures")
    if run is None and figs is None:
        return None
    out = section()
    sm = (run.summary if run else None) or {}
    word = unit_of(kind_of(inp))
    out["unit"] = word
    if run is not None:
        ann = run.csv("cluster_annotation.csv", "cluster", "cell_type")
        summ = run.csv("cell_type_summary.csv", "sample_id", "cell_type", "n_obs")
        markers = run.csv("cluster_markers.csv", "cluster", "gene")
        top: dict[str, list[str]] = {}
        if markers is not None:
            sort_cols = ["cluster", "rank"] if "rank" in markers else ["cluster"]
            for c, g in markers.sort_values(sort_cols).groupby("cluster", sort=False):
                top[str(c)] = g["gene"].astype(str).head(5).tolist()
        clusters = []
        if ann is not None:
            for _, r in ann.iterrows():
                score, second = r.get("score"), r.get("second_score")
                margin = float(score) - float(second) if pd.notna(score) and pd.notna(second) else None
                clusters.append({"cluster": str(r["cluster"]), "cell_type": r["cell_type"], "n_obs": r.get("n_obs"),
                                 "score": score, "second_best": r.get("second_best"), "second_score": second,
                                 "margin": margin, "close": margin is not None and margin < 0.25,
                                 "markers": top.get(str(r["cluster"]), [])})
        out["clusters"] = clusters
        if summ is not None:
            summ["sample_id"] = summ["sample_id"].astype(str)
            if "fraction" not in summ:
                summ["fraction"] = summ["n_obs"] / summ.groupby("sample_id")["n_obs"].transform("sum")
            comp = summ.pivot_table(index="sample_id", columns="cell_type", values="fraction", fill_value=0)
            order = summ.groupby("cell_type")["n_obs"].sum().sort_values(ascending=False).index.tolist()
            comp = comp.reindex(columns=order)
            out["composition"] = comp
            tot = summ.groupby("cell_type")["n_obs"].sum().reindex(order)
            out["types"] = [{"cell_type": t, "n_obs": int(n), "fraction": n / tot.sum(),
                             "min": float(comp[t].min()), "max": float(comp[t].max())} for t, n in tot.items()]
        n_cl = sm.get("n_clusters") or len(clusters)
        n_ct = sm.get("n_cell_types") or (len(out.get("types", [])) or None)
        unas = sm.get("fraction_unassigned")
        out["observations"].append(
            f"{plural(int(n_cl or 0), 'cluster')} named as {plural(int(n_ct or 0), 'cell type')}" +
            (f" from {count(sm.get('n_obs'))} {word}" if sm.get("n_obs") else "") +
            (f"; {pct(unas)} Unassigned" if unas is not None else "") + ".")
        if out.get("types"):
            ts = out["types"]
            out["observations"].append("Largest: " + joined(
                [f"{t['cell_type']} {pct(t['fraction'])}" + (f" ({pct(t['min'])}–{pct(t['max'])} per sample)"
                                                             if len(out['composition']) > 1 else "") for t in ts[:4]]) + ".")
            if len(out["composition"]) > 1:
                diffs = sorted(ts, key=lambda t: -(t["max"] - t["min"]))
                d = diffs[0]
                if d["max"] - d["min"] >= 0.05:
                    out["observations"].append(f"The share that varies most between samples is {d['cell_type']} "
                                               f"({pct(d['min'])} to {pct(d['max'])}).")
        multi = {}
        for c in clusters:
            multi.setdefault(c["cell_type"], []).append(c["cluster"])
        multi = {k: v for k, v in multi.items() if len(v) > 1}
        if multi:
            out["observations"].append("Several clusters share a call: " + joined(
                [f"{k} ({', '.join(v)})" for k, v in multi.items()]) +
                "; they differ in expression in ways the marker table does not name.")
        close = [c for c in clusters if c["close"]]
        if close:
            out["observations"].append("Close calls (runner-up within 0.25 z): " + joined(
                [f"{c['cluster']} {c['cell_type']} ({num(c['score'])}) vs {c['second_best']} ({num(c['second_score'])})"
                 for c in close]) + ".")
        ms = sm.get("marker_sets") or {}
        if ms.get("n_given") is not None and ms.get("n_used") is not None:
            out["marker_sets"] = {"source": ms.get("source"), "n_given": ms.get("n_given"), "n_used": ms.get("n_used"),
                                  "skipped": sorted((ms.get("skipped") or {}).keys()) if isinstance(ms.get("skipped"), dict)
                                  else list(ms.get("skipped") or [])}
            if ms["n_used"] < ms["n_given"]:
                out["observations"].append(
                    f"{ms['n_used']} of {ms['n_given']} {ms.get('source') or ''} marker sets had enough genes on the "
                    f"panel to be scored; the other {ms['n_given'] - ms['n_used']} cell types cannot be called.")
        ag = run.csv("agreement.csv", "ari")
        if ag is not None:
            out["agreement"] = ag.to_dict(orient="records")
        out["caveats"].append("Marker-based labels are heuristic: a call says which marker set a cluster expresses most "
                              "relative to the other clusters of this dataset. Every cluster gets its best set even "
                              "when none fits well (check the score and the runner-up); it does not tell tumour from "
                              "normal epithelium, states or subtypes the table lacks.")
        out["caveats"].append("Rare populations (below ~1–2% of cells) usually merge into a larger cluster at the "
                              "default resolution; their absence here is not evidence that they are absent.")
        if kind_of(inp) == "cell":
            out["caveats"].append("Segmentation spill-over (transcripts of a neighbouring cell assigned to a cell) blends "
                                  "signatures, e.g. epithelial genes in tumour-infiltrating macrophages.")
        params = (run.provenance or {}).get("parameters", {})
        if params.get("INTEGRATION", sm.get("embedding", {}).get("integration") if isinstance(sm.get("embedding"), dict) else None) == "harmony":
            out["caveats"].append("Harmony removed slide-level differences before clustering, including real ones.")
    if figs is not None:
        out["upstream_figures"] = {n: figs.dir / f"{n}.png" for n in ("umap_cell_type", "umap_cluster",
                                                                     "cell_type_by_sample", "dotplot_marker_sets")
                                   if (figs.dir / f"{n}.png").is_file()}
        out["maps"] = [p for p in figs.pngs() if p.name.endswith("_cell_type.png") and not p.name.startswith("umap_")]
    return out


# ---------------------------------------------------------------- deconvolution (cell2location)


def deconvolution(inp: Inputs, max_items: int) -> dict[str, Any] | None:
    run = inp.one("deconvolution")
    if run is None:
        return None
    sm = run.summary or {}
    ts = run.csv("type_summary.csv", "sample_id", "cell_type")
    if ts is None and not sm:
        return None
    out = section()
    word = unit_of(kind_of(inp))
    out["unit"] = word
    if ts is not None:
        ts["sample_id"] = ts["sample_id"].astype(str)
        val = "fraction_of_cells" if "fraction_of_cells" in ts else ("total_abundance" if "total_abundance" in ts else None)
        if val:
            comp = ts.pivot_table(index="sample_id", columns="cell_type", values=val, fill_value=0)
            if val == "total_abundance":
                comp = comp.div(comp.sum(axis=1), axis=0)
            order = comp.mean().sort_values(ascending=False).index.tolist()
            out["composition"] = comp[order]
            tot = (ts.groupby("cell_type")["total_abundance"].sum() if "total_abundance" in ts
                   else comp.mean()).reindex(order)
            frac = tot / tot.sum()
            dom = ts.groupby("cell_type")["n_dominant"].sum().reindex(order) if "n_dominant" in ts else None
            out["types"] = [{"cell_type": t, "fraction": float(frac[t]), "min": float(comp[t].min()),
                             "max": float(comp[t].max()), "n_dominant": int(dom[t]) if dom is not None else None}
                            for t in order]
    ref = sm.get("reference") or {}
    mapping = sm.get("mapping") or {}
    out.update(n_obs=sm.get("n_obs"), n_cell_types=sm.get("n_cell_types") or len(out.get("types", [])),
               reference=ref, mapping=mapping, genes=sm.get("genes") or {}, n_mixed=sm.get("n_mixed"),
               min_dominant_fraction=sm.get("min_dominant_fraction"), dominant_counts=sm.get("dominant_type_counts") or {})
    text = f"cell2location estimated the abundance of {plural(int(out['n_cell_types'] or 0), 'cell type')}"
    if sm.get("n_obs"):
        text += f" in {count(sm['n_obs'])} {word}"
    if mapping.get("mean_cells_per_obs"):
        text += f" (on average {num(mapping['mean_cells_per_obs'], 1)} cells per {unit_of(kind_of(inp), 1)})"
    out["observations"].append(text + ".")
    if ref:
        out["observations"].append(
            f"Reference: {count(ref.get('n_cells_used'))} cells of {plural(int(ref.get('n_types') or 0), 'type')}"
            + (f"; dropped for too few cells: {joined(ref['dropped_types'])}" if ref.get("dropped_types") else "")
            + (f"; {count((sm.get('genes') or {}).get('n_matched'))} genes shared with the spatial data"
               if (sm.get("genes") or {}).get("n_matched") else "") + ".")
    if out.get("types"):
        t = out["types"]
        out["observations"].append("Largest estimated shares of cells: " + joined(
            [f"{x['cell_type']} {pct(x['fraction'])}" for x in t[:5]]) + ".")
    dc = out["dominant_counts"]
    if dc:
        n = sum(dc.values())
        top = sorted(((k, v) for k, v in dc.items() if k != "Mixed"), key=lambda kv: -kv[1])
        out["observations"].append(
            "Dominant type per spot: " + joined([f"{k} {count(v)} ({pct(v / n)})" for k, v in top[:4]]) +
            (f"; {count(dc.get('Mixed', 0))} {word} ({pct(dc.get('Mixed', 0) / n)}) are Mixed (top type below "
             f"{pct(out['min_dominant_fraction'])} of the spot)" if dc.get("Mixed") else "") + ".")
        if out.get("types"):
            hidden = [x["cell_type"] for x in out["types"] if x["fraction"] >= 0.05 and not dc.get(x["cell_type"])]
            if hidden:
                out["observations"].append(f"{joined(hidden)} make up a sizeable share of cells but dominate no "
                                           f"{unit_of(kind_of(inp), 1)}: they are invisible to anything built on "
                                           "dominant_type (st-niches, st-compare).")
    coloc = (sm.get("colocation") or {}).get("groups") or []
    cl = run.csv("colocation.csv", "group", "cell_type", "weight")
    groups = []
    if cl is not None:
        for g, d in cl.groupby("group"):
            d = d.sort_values("weight", ascending=False)
            groups.append({"group": g, "top": [(r.cell_type, float(r.weight)) for r in d.head(4).itertuples() if r.weight >= 0.2]})
    elif coloc:
        groups = [{"group": g["group"], "top": [(t.split(" (")[0], None) for t in g.get("top_types", [])[:3]]} for g in coloc]
    out["colocation"] = groups
    if groups:
        out["observations"].append(f"{plural(len(groups), 'co-location group')}: " + joined(
            [f"{g['group']} " + " + ".join(t for t, _ in g["top"][:3]) for g in groups]) + ".")
    if mapping.get("log1p_correlation") is not None:
        out["observations"].append(f"Reconstruction check: the model's expected counts correlate with the observed "
                                   f"ones at r = {num(mapping['log1p_correlation'])} (log1p scale).")
        if float(mapping["log1p_correlation"]) < 0.6:
            out["caveats"].append("The reconstruction correlation is low: the reference may miss cell types present in "
                                  "the tissue, or the training did not converge (see training_history.png).")
    out["caveats"].append("Abundances are model estimates, not counts: cell2location explains each spot's counts with "
                          "the reference signatures, so cell types absent from the reference cannot be found and "
                          "similar types (subsets of T cells) trade abundance between them. Benchmarks (Spotless, "
                          "Sang-aram et al. 2024) find it among the most accurate methods, with accuracy falling "
                          "for rare and transcriptionally similar types.")
    out["caveats"].append("The dominant type keeps only the most abundant type of each spot; minority types and mixed "
                          "spots are lost to steps that use it.")
    if mapping.get("epochs") and int(mapping["epochs"]) < 10000:
        out["caveats"].append(f"The spatial model was trained for {count(mapping['epochs'])} epochs (cell2location "
                              "recommends ~30,000): a quick run; abundances may not have converged.")
    fig_dir = run.dir / "figures"
    names = ("type_composition", "colocation_weights", "reconstruction_qc", "training_history", "reference_gene_filter")
    out["upstream_figures"] = {n: fig_dir / f"{n}.png" for n in names if (fig_dir / f"{n}.png").is_file()}
    out["maps"] = [p for p in run.pngs("figures") if p.name.endswith(("_abundance.png", "_dominant_type.png",
                                                                      "_colocation.png"))]
    return out


# ---------------------------------------------------------------- niches


def _pair(a: dict[str, Any]) -> str:
    return a["label_a"] if a.get("homotypic") or a["label_a"] == a["label_b"] else f"{a['label_a']}–{a['label_b']}"


def niches_one(run: Run, inp: Inputs, max_items: int) -> dict[str, Any] | None:
    sm = run.summary or {}
    nc = run.csv("niche_composition.csv", "niche", "label", "fraction")
    ne = run.csv("nhood_enrichment.csv", "sample_id", "label_a", "label_b", "zscore")
    if nc is None and ne is None and not sm:
        return None
    out = section()
    label_col = sm.get("label_column") or (run.provenance or {}).get("parameters", {}).get("LABEL_COLUMN") or run.label
    out.update(label=run.label, label_column=label_col, n_obs=sm.get("n_obs"), n_samples=sm.get("n_samples"))
    word = unit_of(kind_of(inp))
    het = sm.get("strongest_enriched_heterotypic") or []
    dep = sm.get("strongest_depleted") or []
    homo = [p for p in (sm.get("strongest_enriched") or []) if p.get("homotypic")]
    if ne is not None and not het:
        pooled = ne[ne["sample_id"].astype(str) == "pooled"]
        pooled = pooled[pooled["label_a"] < pooled["label_b"]]
        het = pooled.sort_values("zscore", ascending=False).head(5).to_dict(orient="records")
        dep = pooled.sort_values("zscore").head(5).to_dict(orient="records")
    out["enriched"], out["depleted"], out["homotypic"] = het[:max_items], dep[:max_items], homo[:max_items]
    if ne is not None:
        pooled = ne[ne["sample_id"].astype(str) == "pooled"].copy()
        if len(pooled):
            pooled = pooled[pooled["label_a"].astype(str) <= pooled["label_b"].astype(str)]
            out["pairs"] = pooled.sort_values("zscore", ascending=False)
    if homo:
        strongest = (sm.get("strongest_enriched") or [{}])[0]
        lead = ("The strongest enrichment is between neighbours of the same label" if strongest.get("homotypic")
                else "Same-label contacts")
        out["observations"].append(lead + " (z: " + joined(
            [f"{p['label_a']} {num(p['zscore'], 0)}" for p in homo[:4]]) + ").")
    if het:
        out["observations"].append("Strongest mixed contacts (pooled neighbourhood-enrichment z): " + joined(
            [f"{_pair(p)} {num(p['zscore'], 0)}" for p in het[:5] if float(p['zscore']) > 2]) + ".")
    if dep:
        neg = [p for p in dep[:5] if float(p["zscore"]) < -2]
        if neg:
            out["observations"].append("Most segregated (avoiding each other): " + joined(
                [f"{_pair(p)} {num(p['zscore'], 0)}" for p in neg]) + ".")
    if sm.get("n_pairs_enriched") is not None:
        out["n_pairs_enriched"], out["n_pairs_depleted"] = sm.get("n_pairs_enriched"), sm.get("n_pairs_depleted")
    niches = []
    if nc is not None:
        for n, g in nc.groupby("niche"):
            g = g.sort_values("fraction", ascending=False)
            niches.append({"niche": n, "n_obs": int(g["n_obs"].iloc[0]) if "n_obs" in g else None,
                           "top": [{"label": r.label, "fraction": float(r.fraction),
                                    "enrichment": float(getattr(r, "enrichment", float("nan")))} for r in g.head(3).itertuples()]})
        tot = sum(n["n_obs"] or 0 for n in niches)
        for n in niches:
            n["fraction"] = (n["n_obs"] / tot) if tot and n["n_obs"] is not None else None
        out["composition"] = nc.pivot_table(index="niche", columns="label", values="fraction", fill_value=0)
    elif sm.get("niches"):
        niches = [{"niche": n["niche"], "n_obs": n.get("n_obs"), "fraction": n.get("fraction"),
                   "top": n.get("top_labels", [])} for n in sm["niches"]]
    out["niches"] = niches
    sel = sm.get("niche_selection") or {}
    if niches:
        text = f"{plural(len(niches), 'niche')}"
        sil = sel.get("silhouette") or {}
        if sel.get("chosen_by") == "silhouette" and sil:
            k = str(len(niches))
            best = sil.get(k)
            alts = sorted(((kk, v) for kk, v in sil.items() if kk != k), key=lambda kv: -kv[1])
            text += f", chosen by silhouette ({num(best)}" + (f"; k = {alts[0][0]} scores {num(alts[0][1])}" if alts else "") + ")"
        elif sel.get("requested"):
            text += " (N_NICHES)"
        out["observations"].append(text + ": " + "; ".join(
            f"{n['niche']} ({pct(n['fraction'])} of {word}): " + ", ".join(
                f"{t['label']} {pct(t['fraction'])}" for t in n["top"][:3]) for n in niches[:8]) + ".")
    ns = run.csv("niche_summary.csv", "sample_id", "niche", "fraction")
    if ns is not None:
        ns["sample_id"] = ns["sample_id"].astype(str)
        piv = ns.pivot_table(index="niche", columns="sample_id", values="fraction", fill_value=0)
        out["per_sample"] = piv
        if piv.shape[1] > 1:
            spread = (piv.max(axis=1) - piv.min(axis=1)).sort_values(ascending=False)
            n0 = spread.index[0]
            if spread.iloc[0] >= 0.05:
                out["observations"].append(f"{n0} varies most between samples: {pct(piv.loc[n0].min())} "
                                           f"({piv.loc[n0].idxmin()}) to {pct(piv.loc[n0].max())} ({piv.loc[n0].idxmax()}).")
    per = sm.get("per_sample") or {}
    sparse = [k for k, v in per.items() if isinstance(v, dict) and v.get("mean_degree") is not None
              and float(v["mean_degree"]) < 4 and v.get("graph") != "grid"]
    if sparse:
        out["caveats"].append(f"The neighbour graph of {joined(sparse)} is sparse (mean degree below 4), typical of a "
                              "subsample of cells: many cells have few or no neighbours under the distance cap. Use "
                              "the full data for final niches and enrichment.")
    if (sm.get("n_obs") or 0) > 20000:
        out["caveats"].append(f"Enrichment z-scores grow with the number of {word}: with {count(sm['n_obs'])} almost "
                              "every pair is 'significant'. Compare z-scores between pairs, and read the interaction "
                              "fractions for effect size; they are not corrected for multiple tests.")
    out["caveats"].append("The pooled z-score treats the samples as one tissue (labels shuffled within each sample); "
                          "it is not a test across replicates or conditions.")
    out["caveats"].append(f"Everything here is about the {label_col} labels, not expression: it is only as good as the "
                          "labelling upstream; mislabelling at boundaries shows up as enrichment between neighbouring labels.")
    if str(label_col) == "dominant_type":
        out["caveats"].append("dominant_type keeps only the most abundant cell type per spot: these niches are "
                              "neighbourhoods of spot types, and minority types are invisible.")
    out["caveats"].append("Niches are KMeans clusters of neighbourhood compositions with a silhouette-chosen number: a "
                          "convenient partition, not a unique one.")
    out["upstream_figures"] = {n: run.dir / f"{n}.png" for n in ("nhood_enrichment", "niche_composition")
                               if (run.dir / f"{n}.png").is_file()}
    out["maps"] = [p for p in run.pngs() if p.name.endswith(("_niches.png", "_nhood_enrichment.png"))]
    out["cooccurrence"] = [p for p in run.pngs() if p.name.startswith("co_occurrence_")]
    return out


def niches(inp: Inputs, max_items: int) -> dict[str, Any] | None:
    runs = [r for r in (niches_one(run, inp, max_items) for run in inp.many("niches")) if r]
    return {"runs": runs, "observations": [], "caveats": []} if runs else None


# ---------------------------------------------------------------- cell-cell communication


def ccc_one(run: Run, inp: Inputs, max_items: int) -> dict[str, Any] | None:
    sm = run.summary or {}
    cons = run.csv("lr_consensus.csv", "source", "target", "ligand_complex", "receptor_complex")
    sp = run.csv("lr_spatial.csv", "sample_id", "ligand", "receptor", "global_morans")
    if cons is None and sp is None and not sm:
        return None
    out = section()
    label_col = sm.get("label_column") or (run.provenance or {}).get("parameters", {}).get("LABEL_COLUMN") or run.label
    out.update(label=run.label, label_column=label_col, resource=sm.get("resource"),
               n_interactions=sm.get("n_interactions"), n_consensus=sm.get("n_consensus_interactions"),
               n_every=sm.get("n_consensus_in_every_sample"), n_samples=sm.get("n_samples"))
    per = sm.get("samples") or sm.get("per_sample") or {}
    cover = {k: v.get("n_lr_pairs_in_data") for k, v in per.items() if isinstance(v, dict) and v.get("n_lr_pairs_in_data") is not None}
    out["coverage"] = cover
    spatial = None
    if sp is not None:
        sp = sp.copy()
        sp["pos_sig"] = (pd.to_numeric(sp["global_morans"], errors="coerce") > 0) & \
                        (pd.to_numeric(sp.get("global_pval"), errors="coerce") <= 0.05)
        g = sp.groupby(["ligand", "receptor"])
        spatial = pd.DataFrame({"mean_global_morans": g["global_morans"].mean(), "n_samples": g.size(),
                                "n_samples_significant": g["pos_sig"].sum(),
                                "n_obs_significant": g["n_obs_significant"].sum() if "n_obs_significant" in sp else np.nan}) \
            .reset_index().sort_values("mean_global_morans", ascending=False)
        out["spatial"] = spatial
    top = []
    if cons is not None:
        cons = cons.sort_values("rank") if "rank" in cons else cons
        sp_r = {(r.ligand, r.receptor): r.mean_global_morans for r in spatial.itertuples()} if spatial is not None else {}
        for _, r in cons.head(max_items).iterrows():
            top.append({"rank": r.get("rank"), "source": r["source"], "target": r["target"],
                        "ligand": r["ligand_complex"], "receptor": r["receptor_complex"],
                        "n_samples": r.get("n_samples"), "n_samples_top": r.get("n_samples_top"),
                        "mean_magnitude_rank": r.get("mean_magnitude_rank"),
                        "mean_specificity_rank": r.get("mean_specificity_rank"),
                        "spatial_r": sp_r.get((r["ligand_complex"], r["receptor_complex"]))})
    out["top"] = top
    n_s = int(sm.get("n_samples") or (len(per) if per else 0) or 0)
    if top:
        out["observations"].append("Top-ranked interactions across samples (LIANA+ consensus): " + joined(
            [f"{t['source']} → {t['target']} {t['ligand']} → {t['receptor']}" +
             (f" (top-{(sm.get('parameters') or {}).get('TOP_N', 20)} in {t['n_samples_top']}/{t['n_samples']} samples)"
              if t.get("n_samples_top") == t.get("n_samples_top") and t.get("n_samples") else "") for t in top[:4]]) + ".")
    if out["n_consensus"]:
        out["observations"].append(f"{count(out['n_consensus'])} interactions scored across samples" +
                                   (f", {count(out['n_every'])} of them in every sample" if out.get("n_every") is not None and n_s > 1 else "") + ".")
    if spatial is not None and len(spatial):
        pos = spatial[spatial["n_samples_significant"] > 0]
        out["observations"].append(
            f"Spatially co-expressed (bivariate Moran's R > 0, p ≤ 0.05) for {plural(len(pos), 'pair')} of the "
            f"{len(spatial)} scored: strongest " + joined(
                [f"{r.ligand} → {r.receptor} (R {num(r.mean_global_morans)}, {int(r.n_samples_significant)}/{int(r.n_samples)} samples)"
                 for r in spatial.head(4).itertuples()]) + ".")
        anti: dict[str, float] = {}
        for t in top[:max_items]:
            if t.get("spatial_r") is not None and t["spatial_r"] < -0.05:
                anti.setdefault(f"{t['ligand']} → {t['receptor']}", t["spatial_r"])
        if anti:
            out["anti_correlated"] = list(anti)
            out["observations"].append(
                "Ranked high by expression but anti-correlated in space: " +
                joined([f"{k} (R {num(r)})" for k, r in list(anti.items())[:4]]) +
                "; the cells expressing ligand and receptor sit apart, so the expression rank alone overstates these.")
    if cover:
        lo = min(cover.values())
        out["observations"].append(f"Ligand–receptor pairs of the resource with all genes measured: {count(lo)}" +
                                   (f"–{count(max(cover.values()))}" if max(cover.values()) != lo else "") + " per sample.")
        if lo < 200:
            out["caveats"].append(f"Only {count(lo)} ligand–receptor pairs of the {out.get('resource') or ''} resource "
                                  "are fully covered by the measured genes: most pathways are invisible to this panel, "
                                  "and a missing interaction says nothing about its absence.")
    out["caveats"].append("Inference, not measurement: ligand and receptor mRNA in the right cells is necessary but not "
                          "sufficient for signalling (protein, secretion, co-receptors and the response are not "
                          "observed). Treat the list as ranked hypotheses; ranks are not p-values.")
    if kind_of(inp) == "cell":
        out["caveats"].append("Segmentation spill-over between touching cells inflates the spatial scores of pairs whose "
                              "genes mark neighbouring cell types.")
    elif kind_of(inp) == "spot":
        out["caveats"].append("Spots mix several cells, so autocrine and paracrine signals cannot be told apart.")
    out["caveats"].append("The consensus counts samples; it does not test conditions against each other.")
    out["dotplots"] = [p for p in run.pngs() if p.name.endswith("_dotplot.png")]
    out["source_target"] = [p for p in run.pngs() if p.name.endswith("_source_target.png")]
    out["local"] = [p for p in run.pngs() if "_local_" in p.name]
    return out


def ccc(inp: Inputs, max_items: int) -> dict[str, Any] | None:
    runs = [r for r in (ccc_one(run, inp, max_items) for run in inp.many("ccc")) if r]
    return {"runs": runs, "observations": [], "caveats": []} if runs else None


# ---------------------------------------------------------------- methods & citations


def _param(v: Any) -> str:
    if isinstance(v, (list, dict)):
        return json.dumps(v, separators=(",", ":"))
    if v is None:
        return "null"
    return str(v)


def _method_text(m: Any) -> list[str]:
    if isinstance(m, str):
        return [m]
    if isinstance(m, list):
        return [str(x) for x in m]
    if isinstance(m, dict):
        return [f"{k}: {v if isinstance(v, str) else json.dumps(v)}" for k, v in m.items()]
    return []


def methods(inp: Inputs, own: dict[str, Any]) -> dict[str, Any]:
    from processor import citations as C

    groups: dict[str, dict[str, Any]] = {}
    for path, prov in inp.provenances:
        app = str(prov.get("app"))
        if app == own["app"]:
            continue
        version = prov.get("version") or prov.get("appVersion") or ""
        tools = prov.get("tools") or prov.get("tool_versions") or {}
        params = prov.get("parameters") or {}
        g = groups.setdefault(app, {"app": app, "runs": [], "step": C.APPS.get(app, {}).get("step", ""),
                                    "cite_keys": C.keys_for(app), "given_citations": []})
        if not g["cite_keys"] and isinstance(prov.get("citations"), list):
            g["given_citations"] = [c for c in prov["citations"] if isinstance(c, dict)]
        key = json.dumps([version, params], sort_keys=True, default=str)
        for run in g["runs"]:
            if run["key"] == key:
                run["files"].append(inp.rel(path))
                for k, v in tools.items():
                    run["tools"].setdefault(k, v)
                break
        else:
            g["runs"].append({"key": key, "version": version, "tools": dict(tools), "parameters": params,
                              "files": [inp.rel(path)], "method": _method_text(prov.get("method"))})
    groups[own["app"]] = {"app": own["app"], "step": C.APPS[own["app"]]["step"], "cite_keys": C.keys_for(own["app"]),
                          "given_citations": [],
                          "runs": [{"key": "", "version": own["version"], "tools": own["tools"],
                                    "parameters": own["parameters"], "files": [], "method": []}]}
    ordered = sorted(groups.values(), key=lambda g: (C.ORDER.index(g["app"]) if g["app"] in C.ORDER else 99, g["app"]))
    # one numbered, deduplicated reference list, in order of first citation
    keys: list[str] = []
    for g in ordered:
        for k in g["cite_keys"]:
            if k not in keys:
                keys.append(k)
    refs = [C.reference(k) for k in keys]
    number = {r["key"]: i + 1 for i, r in enumerate(refs)}
    for g in ordered:
        g["cite_numbers"] = [number[k] for k in g["cite_keys"]]
        for run in g["runs"]:
            run["parameters_text"] = ", ".join(f"{k}={_param(v)}" for k, v in run["parameters"].items())
            run["tools_text"] = ", ".join(f"{k} {v}" for k, v in run["tools"].items() if v)
    return {"apps": ordered, "references": refs}
