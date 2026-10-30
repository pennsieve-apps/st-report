"""
The HTML page (Jinja2, inline CSS, no scripts, no external resources: it works offline from
a dataset folder), the summary tables and summary.json.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from processor import analysis as A

TEMPLATES = Path(__file__).parent / "templates"


def environment() -> Environment:
    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=select_autoescape(["html", "j2"]),
                      trim_blocks=True, lstrip_blocks=True)
    env.filters["pct"] = A.pct
    env.filters["pval"] = A.pval
    env.filters["num"] = A.num
    env.filters["count"] = A.count
    env.filters["sig"] = A.sig
    env.filters["pctv"] = lambda v, d=1: "–" if v is None or pd.isna(v) else f"{float(v):.{d}f}%"
    env.filters["isnum"] = lambda v: v is not None and not (isinstance(v, float) and pd.isna(v))
    env.tests["present"] = lambda v: v is not None and not (isinstance(v, float) and pd.isna(v)) and v != ""
    return env


def write_index(report_dir: Path, context: dict[str, Any]) -> Path:
    page = report_dir / "index.html"
    page.write_text(environment().get_template("index.html.j2").render(**context), encoding="utf-8")
    return page


def _join(v: Any) -> Any:
    return ";".join(map(str, v)) if isinstance(v, (list, tuple)) else v


def write_tables(tables_dir: Path, ctx: dict[str, Any]) -> list[str]:
    written: list[str] = []

    def save(name: str, df: pd.DataFrame | None) -> None:
        if df is None or df.empty:
            return
        df = df.copy()
        for c in df.columns:
            if df[c].map(lambda v: isinstance(v, (list, tuple))).any():
                df[c] = df[c].map(_join)
        df.to_csv(tables_dir / name, index=False)
        written.append(name)

    st, qc, dom, sv, cmp_, ct, dc, ni, cc = (ctx.get(k) for k in ("st", "qc", "dom", "sv", "cmp", "ct", "dc", "ni", "cc"))
    if st:
        save("samples.csv", pd.DataFrame(st["samples"]))
    if qc:
        save("qc.csv", pd.DataFrame(qc["rows"]))
    if dom:
        save("domains.csv", pd.DataFrame(dom["rows"]))
        cov, cnt = dom.get("coverage"), dom.get("coverage_counts")
        if cov is not None:
            long = cov.stack().rename("fraction").reset_index().rename(columns={"level_1": "sample_id"})
            long = long.merge(cnt.stack().rename("n_obs").reset_index().rename(columns={"level_1": "sample_id"}),
                              on=["domain", "sample_id"])
            save("domain_coverage.csv", long[["sample_id", "domain", "n_obs", "fraction"]])
        if dom.get("agreement"):
            ag = dom["agreement"]
            rows = list(ag["per_sample"]) + [{"sample_id": "all", "ari": ag["ari"], "nmi": ag["nmi"]}]
            save("domain_agreement.csv", pd.DataFrame(rows))
    if sv:
        save("svg_samples.csv", pd.DataFrame(sv["per_sample"]))
        save("svg_consensus.csv", pd.DataFrame(sv["consensus"]).rename(columns={"markers": "top10_marker_of"}))
    if cmp_:
        comp, des, top = [], [], []
        for r in cmp_["runs"]:
            for row in r.get("composition") or []:
                row = dict(row)
                means = row.pop("means")
                comp.append({"label_set": r["label"], "label_column": r["label_column"], **row,
                             **{f"mean_fraction_{k}": v for k, v in means.items()}})
            des += [{"label_set": r["label"], **row} for row in r.get("de_summary") or []]
            top += [{"label_set": r["label"], **row} for row in r.get("top_genes") or []]
        save("composition_tests.csv", pd.DataFrame(comp))
        save("de_summary.csv", pd.DataFrame(des))
        save("de_top_genes.csv", pd.DataFrame(top))
    if ct:
        if ct.get("composition") is not None:
            long = ct["composition"].stack().rename("fraction").reset_index()
            long.columns = ["sample_id", "cell_type", "fraction"]
            save("cell_types.csv", long)
        save("clusters.csv", pd.DataFrame(ct.get("clusters") or []))
    if dc:
        if dc.get("composition") is not None:
            long = dc["composition"].stack().rename("fraction_of_cells").reset_index()
            long.columns = ["sample_id", "cell_type", "fraction_of_cells"]
            save("deconvolution_composition.csv", long)
        save("deconvolution_types.csv", pd.DataFrame(dc.get("types") or []))
    if ni:
        nrows, prow = [], []
        for r in ni["runs"]:
            for n in r.get("niches") or []:
                nrows.append({"label_set": r["label"], "niche": n["niche"], "n_obs": n.get("n_obs"),
                              "fraction": n.get("fraction"),
                              "top_labels": ";".join(f"{t['label']} ({t['fraction']:.3f})" for t in n["top"])})
            if r.get("pairs") is not None:
                prow.append(r["pairs"].assign(label_set=r["label"]))
        save("niches.csv", pd.DataFrame(nrows))
        if prow:
            p = pd.concat(prow)
            save("neighbourhood_enrichment.csv", p[["label_set"] + [c for c in p.columns if c != "label_set"]])
    if cc:
        top, sp = [], []
        for r in cc["runs"]:
            top += [{"label_set": r["label"], **t} for t in r.get("top") or []]
            if r.get("spatial") is not None:
                sp.append(r["spatial"].assign(label_set=r["label"]))
        save("ccc_top.csv", pd.DataFrame(top))
        if sp:
            s = pd.concat(sp)
            save("ccc_spatial.csv", s[["label_set"] + [c for c in s.columns if c != "label_set"]])
    meth = ctx["methods"]
    save("methods.csv", pd.DataFrame([{"app": g["app"], "step": g["step"], "version": run["version"],
                                       "tools": run["tools_text"], "parameters": run["parameters_text"],
                                       "citations": ";".join(map(str, g["cite_numbers"]))}
                                      for g in meth["apps"] for run in g["runs"]]))
    save("citations.csv", pd.DataFrame([{"number": i + 1, "authors": r["short_authors"], "title": r["title"],
                                         "journal": r["journal_display"], "year": r["year"], "volume_pages": r["where"],
                                         "doi": r["doi"], "pmid": r.get("pmid", "")}
                                        for i, r in enumerate(meth["references"])]))
    return written


def _texts(sec: dict[str, Any] | None) -> dict[str, Any] | None:
    if not sec:
        return None
    return {"observations": sec.get("observations", []), "caveats": sec.get("caveats", [])}


def summary_json(ctx: dict[str, Any], own: dict[str, Any]) -> dict[str, Any]:
    st, qc, dom, sv, cmp_, ct, dc, ni, cc = (ctx.get(k) for k in ("st", "qc", "dom", "sv", "cmp", "ct", "dc", "ni", "cc"))
    out: dict[str, Any] = {"app": own["app"], "appVersion": own["version"], "generatedAt": ctx["generated_at"],
                           "title": ctx["title"], "parameters": own["parameters"], "workflows": ctx["workflows"],
                           "sectionErrors": ctx["errors"]}
    if st:
        out["study"] = {k: st.get(k) for k in ("n_samples", "platforms", "kind", "n_obs", "n_genes", "condition_column",
                                               "conditions", "subsampled")} | _texts(st)
    if qc:
        out["qc"] = {k: qc.get(k) for k in ("n_before", "n_after", "n_removed", "by_reason", "filter_mode",
                                            "depth_fold")} | _texts(qc)
    if dom:
        out["domains"] = {"n_domains": dom["n_domains"], "domains": [{k: r[k] for k in ("domain", "n_obs", "fraction",
                                                                                         "n_samples", "conditions", "markers")}
                                                                      for r in dom["rows"]],
                          "agreement": {k: dom["agreement"][k] for k in ("ari", "nmi")} if dom.get("agreement") else None,
                          "best_match": dom.get("best_match")} | _texts(dom)
    if sv:
        out["svg"] = {k: sv.get(k) for k in ("fdr", "n_genes_tested", "n_significant_any", "n_significant_all",
                                             "small_effect_fraction")} | \
            {"per_sample": [{k: r[k] for k in ("sample_id", "n_tested", "n_significant", "median_morans_i")}
                            for r in sv["per_sample"]],
             "consensus_top": [c["gene"] for c in sv["consensus"]]} | _texts(sv)
    if cmp_:
        out["comparison"] = [{k: r.get(k) for k in ("label", "label_column", "condition_column", "reference", "levels",
                                                    "per_level", "fdr", "composition_n_tested",
                                                    "composition_n_significant", "de_n_tested", "de_n_significant")}
                             | {"significant_labels": [c["label"] for c in r.get("composition") or [] if c["significant"]]}
                             | _texts(r) for r in cmp_["runs"]]
    if ct:
        out["cell_types"] = {"types": ct.get("types"), "clusters": [{k: c[k] for k in ("cluster", "cell_type", "n_obs",
                                                                                         "score", "second_best", "margin")}
                                                                    for c in ct.get("clusters") or []],
                             "marker_sets": ct.get("marker_sets")} | _texts(ct)
    if dc:
        out["deconvolution"] = {k: dc.get(k) for k in ("n_obs", "n_cell_types", "types", "dominant_counts", "n_mixed",
                                                       "colocation", "mapping")} | _texts(dc)
    if ni:
        out["niches"] = [{"label": r["label"], "label_column": r["label_column"], "n_obs": r.get("n_obs"),
                          "niches": r.get("niches"), "enriched": r.get("enriched"), "depleted": r.get("depleted")}
                         | _texts(r) for r in ni["runs"]]
    if cc:
        out["cell_cell_communication"] = [{"label": r["label"], "resource": r.get("resource"),
                                           "coverage": r.get("coverage"), "n_consensus": r.get("n_consensus"),
                                           "top": r.get("top"), "anti_correlated": r.get("anti_correlated")}
                                          | _texts(r) for r in cc["runs"]]
    out["methods"] = [{"app": g["app"], "versions": [r["version"] for r in g["runs"]], "citations": g["cite_numbers"]}
                      for g in ctx["methods"]["apps"]]
    out["citations"] = [{"number": i + 1, "doi": r["doi"], "pmid": r.get("pmid") or None, "text": r["text"]}
                        for i, r in enumerate(ctx["methods"]["references"])]
    out["inputData"] = [{"dataset": d["text"], "doi": d.get("doi"), "licence": d["licence"], "n_samples": d["n_samples"]}
                        for d in ctx["credit"]]
    return out
