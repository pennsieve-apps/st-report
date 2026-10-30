"""
A small synthetic set of upstream outputs, following the workflow contract
(conventions.md and each app's README), for the tests and `make run`.

Workflow A (Visium-like spots, 4 samples, ctrl x2 / treat x2): metadata, qc, domains,
svg, compare/domain and compare/annotation. Workflow B (the same spots): metadata, qc,
deconvolution (cell2location), niches/colocation-group. Workflow C (Xenium-like cells, 2
samples): metadata, qc, celltype, celltype-figures, niches/cell-type, ccc/cell-type.
Only the summary tables, summary.json, provenance.json and tiny PNGs are written; none of
the large per-observation files the report never reads.

  python -m tests.make_fixture data/fixture --layout both [--workflow all|a|b|c]

Layouts: linked (INPUT_DIR/<port>/...) and flat (everything under arbitrary nested folder
names, found by file name).
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

TINY_PNG = bytes.fromhex(
    "89504e470d0a1a0a0000000d49484452000000010000000108060000001f15c489"
    "0000000d49444154789c6360000002000154a24f5d0000000049454e44ae426082")

SPOTS = ["s1", "s2", "s3", "s4"]
COND = {"s1": "ctrl", "s2": "ctrl", "s3": "treat", "s4": "treat"}
DONOR = {"s1": "d1", "s2": "d2", "s3": "d3", "s4": "d4"}
CELLS = ["x1", "x2"]
DOMAINS = ["D01", "D02", "D03", "D04"]
TYPES = ["T cell", "B cell", "Epithelial", "Fibroblast"]
DECONV_TYPES = ["T_cell", "B_cell", "Stromal", "Myeloid"]
WORKFLOW_PORTS = {
    "a": ["metadata", "qc", "domains", "svg", "compare"],
    "b": ["metadata", "qc", "deconvolution", "niches"],
    "c": ["metadata", "qc", "celltype", "celltype-figures", "niches", "ccc"],
}


def _prov(app: str, version_key: str = "version", **extra) -> dict:
    return {"app": app, version_key: "0.1.0", "tools": {"python": "3.12", "scanpy": "1.12.5"},
            "parameters": {"SEED": 42, **extra.pop("parameters", {})}, **extra}


def metadata(kind: str) -> dict[str, object]:
    ids = SPOTS if kind == "spot" else CELLS
    rows = []
    for i, s in enumerate(ids):
        rows.append({"sample_id": s, "path": f"input/{s.upper()}", "platform": "visium" if kind == "spot" else "xenium",
                     "kind": kind, "n_obs": 500 + 50 * i if kind == "spot" else 3000 + 500 * i,
                     "n_genes": 2000 if kind == "spot" else 313, "n_obs_in_tissue": 500 + 50 * i if kind == "spot" else 3000 + 500 * i,
                     **({"condition": COND[s], "donor": DONOR[s]} if kind == "spot" else {"condition": "tumour", "replicate": i + 1})})
    prep = {"app": "st-prep", "appVersion": "0.1.0", "parameters": {"MAX_OBS_PER_SAMPLE": 0},
            "samples": {r["sample_id"]: {"platform": r["platform"], "kind": kind,
                                         "coordinate_unit": "pixel" if kind == "spot" else "micron",
                                         "features": {"dropped": {} if kind == "spot" else
                                                      {"Negative Control Probe": {"n_features": 20, "total_counts": 99}}},
                                         "observations": {"n_subsampled_to": None},
                                         "annotation": {"n_annotated": 300} if kind == "spot" else None,
                                         "image": {"lowres": kind == "spot"}, "warnings": [], "n_obs": r["n_obs"],
                                         "n_genes": r["n_genes"], "source": r["path"]} for r in rows},
            "skipped": []}
    return {"samples.csv": pd.DataFrame(rows), "prep.json": prep,
            "provenance.json": _prov("st-prep", "appVersion", parameters={"MAX_OBS_PER_SAMPLE": 0})}


def qc(kind: str) -> dict[str, object]:
    ids = SPOTS if kind == "spot" else CELLS
    rows = []
    for i, s in enumerate(ids):
        before = 500 + 50 * i if kind == "spot" else 3000 + 500 * i
        removed = [5, 2, 40, 3][i] if kind == "spot" else [120, 30][i]
        rows.append({"sample_id": s, "kind": kind, "platform": "visium" if kind == "spot" else "xenium",
                     "n_obs_before": before, "n_obs_after": before - removed, "pct_removed": 100 * removed / before,
                     "n_genes_before": 2000, "n_genes_after": 1900 - 10 * i, "n_genes_dropped": 100 + 10 * i,
                     "median_counts": [500, 4000, 800, 1200][i] if kind == "spot" else [170, 165][i],
                     "median_genes": 300, "median_pct_mt": 0.0 if kind == "cell" else 2.5,
                     "median_control_fraction": 0.01 if kind == "cell" else None, "filter_mode": "mad",
                     "min_counts": 100.0, "min_genes": 80.0, "max_pct_mt": 20.0, "max_control_fraction": 0.1 if kind == "cell" else None,
                     "n_removed_low_counts": removed - 1, "n_removed_low_genes": removed, "n_removed_mt": 0,
                     "n_removed_control": 1 if kind == "cell" else 0, "n_hvg": 2000, "n_pcs": 50, "target_sum": 500, "status": "ok"})
    th = pd.DataFrame([{"sample_id": s, "kind": kind, "mode": "mad", "metric": "total_counts", "direction": "below",
                        "threshold": 100.0, "median": 6.2, "mad": 0.3, "mads": 5, "source": "expm1(...)"} for s in ids])
    files: dict[str, object] = {"qc_summary.csv": pd.DataFrame(rows), "thresholds.csv": th,
                                "summary.json": {"n_samples": len(ids)},
                                "provenance.json": _prov("st-qc", parameters={"MADS": 5.0, "FILTER_MODE": "mad"},
                                                         method={"filter": "MAD outliers"})}
    for s in ids:
        files[f"{s}_qc.png"] = TINY_PNG
        files[f"{s}_spatial_counts.png"] = TINY_PNG
    return files


def domains() -> dict[str, object]:
    rng = np.random.default_rng(1)
    rows = []
    for s in SPOTS:
        n = {"D01": 200, "D02": 150, "D03": 100, "D04": 40 if s == "s3" else 0}
        tot = sum(n.values())
        for d in DOMAINS:
            rows.append({"sample_id": s, "domain": d, "n_obs": n[d], "fraction": n[d] / tot,
                         "domain_n_samples": 1 if d == "D04" else 4})
    markers = []
    genes = {"D01": ["SFTPB", "SFTPC", "NAPSA"], "D02": ["IGHG1", "IGKC", "JCHAIN"], "D03": ["COL1A1", "DCN", "LUM"],
             "D04": ["KRT5", "KRT14", "S100A2"]}
    for d, gs in genes.items():
        for r, g in enumerate(gs, 1):
            markers.append({"domain": d, "gene": g, "logfoldchange": 3 - r * 0.5, "pval_adj": 1e-10 * r,
                            "pct_in": 0.8, "pct_out": 0.2, "rank": r})
    agreement = pd.DataFrame([{"sample_id": s, "ari": float(rng.uniform(0.2, 0.6)), "nmi": 0.4, "n_obs": 300,
                               "n_classes": 3, "n_domains": 4} for s in SPOTS] +
                             [{"sample_id": "all", "ari": 0.35, "nmi": 0.42, "n_obs": 1200, "n_classes": 3, "n_domains": 4}])
    match = pd.DataFrame([
        {"annotation": "TLS", "domain": "D02", "n_annotation": 200, "n_domain": 600, "n_overlap": 150,
         "precision": 0.25, "recall": 0.75, "f1": 0.375},
        {"annotation": "Tumor", "domain": "D04", "n_annotation": 400, "n_domain": 40, "n_overlap": 38,
         "precision": 0.95, "recall": 0.095, "f1": 0.17}])
    summary = {"app": "st-domains", "n_samples": 4, "n_domains": 4, "domains": DOMAINS, "resolution": 0.2,
               "lambda": 0.8, "integration": "harmony", "min_obs_present": 20, "n_genes_banksy": 1800,
               "gene_selection": "1800 highly variable genes from st-qc",
               "annotation": {"ignored": ["Unassigned"], "n_obs_used": 1200, "ari": 0.35, "nmi": 0.42}}
    files: dict[str, object] = {"domain_summary.csv": pd.DataFrame(rows), "domain_markers.csv": pd.DataFrame(markers),
                                "agreement.csv": agreement, "annotation_match.csv": match, "summary.json": summary,
                                "provenance.json": {"app": "st-domains", "version": "0.1.0",
                                                    "tool_versions": {"pybanksy": "1.3.5"},
                                                    "parameters": {"LAMBDA": 0.8, "RESOLUTION": 0.2},
                                                    "method": ["BANKSY features", "Harmony", "Leiden"]}}
    for n in ("umap", "markers_dotplot", "confusion", "domain_composition") + tuple(f"{s}_domains" for s in SPOTS):
        files[f"figures/{n}.png"] = TINY_PNG
    return files


def svg() -> dict[str, object]:
    rng = np.random.default_rng(2)
    genes = ["IGHG1", "SFTPB", "COL1A1", "KRT5"] + [f"G{i:03d}" for i in range(60)]
    rows = []
    for s in SPOTS:
        for g in genes:
            planted = g in ("IGHG1", "SFTPB", "COL1A1") or (g == "KRT5" and s == "s3")
            i = float(rng.uniform(0.3, 0.6)) if planted else float(rng.uniform(-0.01, 0.04))
            p = 1e-20 if planted else (0.01 if i > 0.03 else 0.5)
            rows.append({"sample_id": s, "gene": g, "morans_i": i, "pval": p, "pval_adj": p * 2,
                         "significant": p * 2 < 0.05})
    df = pd.DataFrame(rows)
    df["rank"] = df.groupby("sample_id")["morans_i"].rank(ascending=False).astype(int)
    g = df.groupby("gene")
    cons = pd.DataFrame({"n_samples": g.size(), "n_significant": g["significant"].sum(),
                         "mean_morans_i": g["morans_i"].mean(), "median_rank": g["rank"].median()}).reset_index()
    cons = cons.sort_values(["n_significant", "mean_morans_i"], ascending=False)
    cons["rank"] = range(1, len(cons) + 1)
    summary = {"n_samples": 4, "statistic": "moran", "fdr": 0.05,
               "per_sample": {s: {"n_obs": 500, "graph_mean_degree": 5.8, "condition": COND[s]} for s in SPOTS}}
    files: dict[str, object] = {"svg.csv": df, "svg_consensus.csv": cons, "summary.json": summary,
                                "provenance.json": _prov("st-svg", parameters={"FDR": 0.05}),
                                "consensus_top.png": TINY_PNG, "morans_distribution.png": TINY_PNG}
    for s in SPOTS:
        files[f"{s}_svg_top.png"] = TINY_PNG
    return files


def compare(label_column: str, labels: list[str], significant: str | None) -> dict[str, object]:
    comp_tests = []
    for lab in labels:
        sig = lab == significant
        comp_tests.append({"label": lab, "comparison": "treat vs ctrl", "n_samples": 4, "n_ctrl": 2, "n_treat": 2,
                           "mean_fraction_ctrl": 0.30, "mean_fraction_treat": 0.45 if sig else 0.29,
                           "mean_clr_ctrl": 1.0, "mean_clr_treat": 1.6 if sig else 0.98,
                           "clr_difference": 0.6 if sig else -0.02, "parametric_test": "welch_t", "parametric_stat": 9.0,
                           "parametric_p": 0.001 if sig else 0.6, "parametric_p_adj": 0.004 if sig else 0.8,
                           "rank_test": "mann_whitney_u", "rank_stat": 4, "rank_p": 0.33, "rank_p_adj": 1.0,
                           "tested": True, "significant": sig})
    de, dsum = [], []
    for i, lab in enumerate(labels):
        if i == len(labels) - 1:
            dsum.append({"label": lab, "comparison": "treat vs ctrl", "status": "skipped",
                         "reason": "reference level 'ctrl' has 1 pseudobulks (< MIN_SAMPLES_PER_GROUP=2)", "design": "",
                         "note": "", "n_samples_ctrl": 1, "n_samples_treat": 2, "n_genes_tested": 0,
                         "n_significant": None, "n_up": None, "n_down": None, "top_up": "", "top_down": ""})
            continue
        genes = [("IFI6", 2.5, 1e-8), ("ISG15", 2.0, 1e-6), ("CRYAB", -1.5, 0.02), (f"L{i}G", 1.0, 0.03), ("ACTB", 0.1, 0.9)]
        for g, lfc, p in genes:
            de.append({"label": lab, "comparison": "treat vs ctrl", "gene": g, "base_mean": 100.0, "log2_fold_change": lfc,
                       "lfc_se": 0.3, "stat": lfc / 0.3, "p_val": p / 10, "p_adj": p, "significant": p < 0.05,
                       "n_samples_level": 2, "n_samples_reference": 2})
        dsum.append({"label": lab, "comparison": "treat vs ctrl", "status": "tested", "reason": "", "design": "~ condition",
                     "note": "", "n_samples_ctrl": 2, "n_samples_treat": 2, "n_genes_tested": 5, "n_significant": 4,
                     "n_up": 3, "n_down": 1, "top_up": f"IFI6;ISG15;L{i}G", "top_down": "CRYAB"})
    summary = {"label": label_column.replace("_", "-"), "label_column": label_column, "condition_column": "condition",
               "block_column": None, "reference": "ctrl", "levels": ["ctrl", "treat"],
               "n_samples_per_level": {"ctrl": 2, "treat": 2}, "labels": labels}
    files: dict[str, object] = {"composition_tests.csv": pd.DataFrame(comp_tests), "de.csv": pd.DataFrame(de),
                                "de_summary.csv": pd.DataFrame(dsum), "summary.json": summary,
                                "provenance.json": _prov("st-compare", parameters={"LABEL_COLUMN": label_column,
                                                                                   "CONDITION_COLUMN": "condition", "FDR": 0.05}),
                                "composition.png": TINY_PNG, "de_heatmap.png": TINY_PNG}
    for lab in labels[:-1]:
        files[f"volcano_{lab}.png"] = TINY_PNG
    return files


def deconvolution() -> dict[str, object]:
    rows = []
    for s in SPOTS:
        fr = {"T_cell": 0.4, "B_cell": 0.3, "Stromal": 0.22, "Myeloid": 0.08}
        dom = {"T_cell": 300, "B_cell": 150, "Stromal": 40, "Myeloid": 0}
        for t in DECONV_TYPES:
            rows.append({"sample_id": s, "cell_type": t, "total_abundance": fr[t] * 3000, "mean_abundance": fr[t] * 6,
                         "fraction_of_cells": fr[t], "n_dominant": dom[t] // 4, "fraction_dominant": dom[t] / 2000})
    coloc = pd.DataFrame([{"group": g, "cell_type": t, "weight": w} for g, ws in
                          (("G01", {"T_cell": 1.0, "Stromal": 0.5, "B_cell": 0.1, "Myeloid": 0.0}),
                           ("G02", {"B_cell": 1.0, "Myeloid": 0.3, "T_cell": 0.0, "Stromal": 0.1}),
                           ("G03", {"Myeloid": 1.0, "Stromal": 0.1, "T_cell": 0.0, "B_cell": 0.0}))
                          for t, w in ws.items()])
    by_obs = pd.DataFrame([{"sample_id": s, "obs_id": f"{s}:{i}", "G01": 3.0 if i % 2 else 0.5,
                            "G02": 0.5 if i % 2 else 2.0, "G03": 4.0 if i % 10 == 0 else 0.1}
                           for s in SPOTS for i in range(100)])
    summary = {"n_samples": 4, "n_obs": 2000, "cell_types": DECONV_TYPES, "n_cell_types": 4,
               "reference": {"n_cells_used": 480, "n_types": 4, "dropped_types": ["Rare"]},
               "genes": {"n_matched": 154}, "mapping": {"epochs": 1000, "log1p_correlation": 0.86, "mean_cells_per_obs": 6.2},
               "dominant_type_counts": {"T_cell": 300, "B_cell": 150, "Stromal": 40, "Mixed": 60},
               "n_mixed": 60, "min_dominant_fraction": 0.4}
    files: dict[str, object] = {"type_summary.csv": pd.DataFrame(rows), "colocation.csv": coloc,
                                "colocation_by_obs.csv": by_obs,
                                "summary.json": summary,
                                "provenance.json": {"app": "cell2location", "version": "0.1.0", "tools": {"cell2location": "0.1.5"},
                                                    "parameters": {"SPATIAL_EPOCHS": 1000},
                                                    "citations": [{"name": "cell2location", "doi": "10.1038/s41587-021-01139-4"}]}}
    for n in ("type_composition", "colocation_weights", "reconstruction_qc", "training_history"):
        files[f"figures/{n}.png"] = TINY_PNG
    for s in SPOTS:
        files[f"figures/{s}_abundance.png"] = TINY_PNG
        files[f"figures/{s}_dominant_type.png"] = TINY_PNG
    return files


def niches(label_column: str, labels: list[str], samples: list[str], n_obs: int, degree: float) -> dict[str, object]:
    ne = []
    for s in samples + ["pooled"]:
        for a in labels:
            for b in labels:
                z = 50.0 if a == b else (12.0 if {a, b} == set(labels[:2]) else -8.0)
                ne.append({"sample_id": s, "label_a": a, "label_b": b, "zscore": z, "count": 100})
    comp = []
    for n, main in (("N01", labels[0]), ("N02", labels[2])):
        for lab in labels:
            f = 0.7 if lab == main else 0.1
            comp.append({"niche": n, "label": lab, "fraction": f, "enrichment": f / 0.25,
                         "neighbourhood_fraction": f, "n_obs": n_obs // 2})
    ns = pd.DataFrame([{"sample_id": s, "niche": n, "n_obs": 100, "fraction": f}
                       for i, s in enumerate(samples) for n, f in (("N01", 0.4 + 0.1 * i), ("N02", 0.6 - 0.1 * i))])
    het = [{"label_a": labels[0], "label_b": labels[1], "zscore": 12.0, "count": 100, "homotypic": False}]
    summary = {"label": label_column.replace("_", "-"), "label_column": label_column, "n_samples": len(samples),
               "n_obs": n_obs, "labels": labels,
               "per_sample": {s: {"n_obs": n_obs // len(samples), "graph": "delaunay", "mean_degree": degree} for s in samples},
               "strongest_enriched": [{"label_a": l, "label_b": l, "zscore": 50.0, "count": 100, "homotypic": True} for l in labels],
               "strongest_depleted": [{"label_a": labels[0], "label_b": labels[3], "zscore": -8.0, "count": 5, "homotypic": False}],
               "strongest_enriched_heterotypic": het, "n_pairs_enriched": 5, "n_pairs_depleted": 5, "n_niches": 2,
               "niche_selection": {"requested": 0, "silhouette": {"2": 0.5, "3": 0.45}, "chosen_by": "silhouette"}}
    files: dict[str, object] = {"nhood_enrichment.csv": pd.DataFrame(ne), "niche_composition.csv": pd.DataFrame(comp),
                                "niche_summary.csv": ns, "summary.json": summary,
                                "provenance.json": _prov("st-niches", parameters={"LABEL_COLUMN": label_column}),
                                "nhood_enrichment.png": TINY_PNG, "niche_composition.png": TINY_PNG}
    for s in samples:
        files[f"{s}_niches.png"] = TINY_PNG
    return files


def celltype() -> dict[str, object]:
    ann = pd.DataFrame([
        {"cluster": "C01", "cell_type": "Epithelial", "score": 1.2, "n_obs": 2500, "second_best": "Luminal", "second_score": 1.1},
        {"cluster": "C02", "cell_type": "T cell", "score": 1.8, "n_obs": 1500, "second_best": "CD4 T cell", "second_score": 1.3},
        {"cluster": "C03", "cell_type": "Fibroblast", "score": 1.6, "n_obs": 1200, "second_best": "Smooth muscle", "second_score": 0.8},
        {"cluster": "C04", "cell_type": "B cell", "score": 3.0, "n_obs": 600, "second_best": "Plasma cell", "second_score": 1.0},
        {"cluster": "C05", "cell_type": "Epithelial", "score": 1.0, "n_obs": 500, "second_best": "Luminal", "second_score": 0.9}])
    rows = []
    for i, s in enumerate(CELLS):
        fr = {"Epithelial": 0.45 - 0.1 * i, "T cell": 0.25, "Fibroblast": 0.2 + 0.1 * i, "B cell": 0.1}
        n = 3000 + 500 * i
        rows += [{"sample_id": s, "cell_type": t, "n_obs": int(f * n), "fraction": f} for t, f in fr.items()]
    markers = pd.DataFrame([{"cluster": c, "gene": g, "logfoldchange": 2.0, "pval_adj": 0.0, "pct_in": 0.9, "pct_out": 0.2,
                             "rank": r, "score": 50.0}
                            for c, gs in (("C01", ["EPCAM", "KRT8"]), ("C02", ["CD3E", "TRAC"]), ("C03", ["LUM", "DCN"]),
                                          ("C04", ["MS4A1", "CD79A"]), ("C05", ["KRT7", "ERBB2"]))
                            for r, g in enumerate(gs, 1)])
    summary = {"n_samples": 2, "n_obs": 6300, "n_clusters": 5, "n_cell_types": 4, "fraction_unassigned": 0.0,
               "marker_sets": {"source": "builtin", "n_given": 10, "n_used": 7,
                               "skipped": {"Hepatocyte": {}, "Neuron": {}, "Podocyte": {}}}}
    return {"cluster_annotation.csv": ann, "cell_type_summary.csv": pd.DataFrame(rows), "cluster_markers.csv": markers,
            "summary.json": summary,
            "provenance.json": _prov("st-celltype", parameters={"INTEGRATION": "harmony", "RESOLUTION": 0.5})}


def celltype_figures() -> dict[str, object]:
    files: dict[str, object] = {f"{n}.png": TINY_PNG for n in ("umap_cell_type", "umap_cluster", "cell_type_by_sample",
                                                               "dotplot_marker_sets")}
    for s in CELLS:
        files[f"{s}_cell_type.png"] = TINY_PNG
    files["summary.json"] = {"figures": list(files)}
    files["provenance.json"] = _prov("st-celltype", parameters={"INTEGRATION": "harmony", "RESOLUTION": 0.5})
    return files


def ccc() -> dict[str, object]:
    cons = pd.DataFrame([
        {"rank": 1, "source": "Fibroblast", "target": "T cell", "ligand_complex": "CXCL12", "receptor_complex": "CXCR4",
         "n_samples": 2, "n_samples_top": 2, "mean_magnitude_rank": 1e-4, "mean_specificity_rank": 0.01,
         "min_magnitude_rank": 1e-4, "mean_lr_means": 1.9, "samples": "x1;x2"},
        {"rank": 2, "source": "T cell", "target": "Epithelial", "ligand_complex": "S100A4", "receptor_complex": "ERBB2",
         "n_samples": 2, "n_samples_top": 2, "mean_magnitude_rank": 3e-4, "mean_specificity_rank": 0.02,
         "min_magnitude_rank": 2e-4, "mean_lr_means": 1.8, "samples": "x1;x2"},
        {"rank": 3, "source": "B cell", "target": "T cell", "ligand_complex": "CD70", "receptor_complex": "CD27",
         "n_samples": 2, "n_samples_top": 1, "mean_magnitude_rank": 5e-3, "mean_specificity_rank": 0.1,
         "min_magnitude_rank": 1e-3, "mean_lr_means": 0.9, "samples": "x1;x2"}])
    sp = pd.DataFrame([{"sample_id": s, "ligand": lig, "receptor": rec, "global_morans": r, "global_pval": p,
                        "n_obs_positive": 500, "n_obs_significant": 200, "mean_local": 0.1, "ligand_prop": 0.4,
                        "receptor_prop": 0.4, "rank": k}
                       for s in CELLS for k, (lig, rec, r, p) in enumerate(
                           (("CXCL12", "CXCR4", 0.22, 0.0), ("S100A4", "ERBB2", -0.3, 1.0), ("CD70", "CD27", 0.01, 0.4)), 1)])
    summary = {"label": "cell-type", "label_column": "cell_type", "resource": "consensus", "n_samples": 2,
               "n_interactions": 120, "n_consensus_interactions": 60, "n_consensus_in_every_sample": 50,
               "samples": {s: {"n_lr_pairs_in_data": 23} for s in CELLS}, "parameters": {"TOP_N": 20}}
    files: dict[str, object] = {"lr_consensus.csv": cons, "lr_spatial.csv": sp, "summary.json": summary,
                                "provenance.json": _prov("st-ccc", parameters={"LABEL_COLUMN": "cell_type"})}
    for s in CELLS:
        files[f"{s}_dotplot.png"] = TINY_PNG
        files[f"{s}_source_target.png"] = TINY_PNG
        files[f"{s}_local_CXCL12_CXCR4.png"] = TINY_PNG
    return files


def port_files(workflow: str) -> dict[str, dict[str, object]]:
    """{port-relative path: content}, keyed by port."""
    kind = "cell" if workflow == "c" else "spot"
    out: dict[str, dict[str, object]] = {"metadata": metadata(kind), "qc": qc(kind)}
    if workflow == "a":
        out["domains"] = domains()
        out["svg"] = svg()
        out["compare"] = {f"domain/{k}": v for k, v in compare("domain", DOMAINS, "D02").items()}
        out["compare"].update({f"annotation/{k}": v for k, v in compare("annotation", ["TLS", "Tumor", "Normal"], None).items()})
    if workflow == "b":
        out["deconvolution"] = deconvolution()
        out["niches"] = {f"colocation-group/{k}": v for k, v in
                         niches("colocation_group", ["G01", "G02", "G03", "G04"], SPOTS, 1900, 5.9).items()}
    if workflow == "c":
        out["celltype"] = celltype()
        out["celltype-figures"] = celltype_figures()
        out["niches"] = {f"cell-type/{k}": v for k, v in niches("cell_type", TYPES, CELLS, 6300, 3.2).items()}
        out["ccc"] = {f"cell-type/{k}": v for k, v in ccc().items()}
    return out


def _write(path: Path, content) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(content, pd.DataFrame):
        content.to_csv(path, index=False)
    elif isinstance(content, (dict, list)):
        path.write_text(json.dumps(content, indent=1))
    elif isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(str(content))


FLAT_NAMES = {"metadata": "prep_out/metadata", "qc": "qc_out", "domains": "steps/banksy", "svg": "steps/svg",
              "compare": "steps/compare", "deconvolution": "c2l", "celltype": "typing/labels",
              "celltype-figures": "typing/figures", "niches": "steps/niches", "ccc": "steps/ccc"}


def make_fixture(root: Path, layout: str = "linked", workflow: str = "a") -> dict[str, Path]:
    """Writes root/linked and/or root/flat. workflow: a, b, c, or a combination such as 'ac'."""
    root = Path(root)
    ports: dict[str, dict[str, object]] = {}
    for w in workflow:
        for port, files in port_files(w).items():
            if port in ("metadata", "qc") and port in ports:
                continue  # A and B share the head; C's head would clash with theirs
            ports.setdefault(port, {}).update(files)
    out = {}
    for lay in (("linked", "flat") if layout == "both" else (layout,)):
        base = root / lay
        if base.exists():
            shutil.rmtree(base)
        base.mkdir(parents=True)
        for port, files in ports.items():
            folder = base / (port if lay == "linked" else FLAT_NAMES[port])
            for rel, content in files.items():
                _write(folder / rel, content)
        out[lay] = base
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--layout", choices=["linked", "flat", "both"], default="both")
    ap.add_argument("--workflow", default="a", help="a, b, c or a combination, e.g. abc")
    args = ap.parse_args()
    for lay, path in make_fixture(Path(args.root), args.layout, args.workflow).items():
        print(f"{lay}: {path}")
