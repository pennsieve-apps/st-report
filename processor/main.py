"""
st-report: the last step of the three spatial transcriptomics workflows
(A. Spatial Domain Mapping & Condition Comparison, B. Cell-type Deconvolution & Niches,
C. Single-cell Spatial Profiling). It reads whatever upstream outputs it is given and writes

  OUTPUT_DIR/report/index.html        one self-contained page (no scripts, no external resources)
  OUTPUT_DIR/report/figures/*.png     its figures (its own and copies of upstream ones)
  OUTPUT_DIR/report/provenance.json
  OUTPUT_DIR/tables/*.csv             the tables the page shows
  OUTPUT_DIR/tables/summary.json      every headline number, observation and caveat of the page
  OUTPUT_DIR/tables/provenance.json

Every input is optional; a section is rendered when its inputs are present. The run fails
(exit 1) only when none of the inputs is found, and with exit 2 on a bad parameter.
"""
from __future__ import annotations

import json
import logging
import os
import resource
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
log = logging.getLogger("st-report")
# The Lambda runtime configures the root logger at WARNING before basicConfig can.
log.setLevel(logging.INFO)
for noisy in ("matplotlib", "PIL", "fontTools"):
    logging.getLogger(noisy).setLevel(logging.WARNING)

APP_NAME = "st-report"
PARAMETERS = ("REPORT_TITLE", "FDR", "MAX_ITEMS", "MAX_FIGURES")
DEFAULT_TITLE = "Spatial Transcriptomics Report"
REPO_ROOT = Path(__file__).resolve().parent.parent

WORKFLOWS = {
    "A": "Spatial Domain Mapping & Condition Comparison",
    "B": "Cell-type Deconvolution & Niches",
    "C": "Single-cell Spatial Profiling",
}


class BadParameter(ValueError):
    pass


def _now() -> datetime:
    return datetime.now(timezone.utc)


def version_of(root: Path) -> str:
    for p in (root / "VERSION", Path("/processor/VERSION")):
        if p.is_file():
            return p.read_text().strip()
    return "unknown"


def tool_versions() -> dict[str, str]:
    import jinja2
    import matplotlib
    import numpy
    import pandas

    return {"python": sys.version.split()[0], "pandas": pandas.__version__, "numpy": numpy.__version__,
            "matplotlib": matplotlib.__version__, "jinja2": jinja2.__version__}


def parameters() -> dict[str, Any]:
    def env(name: str, default: str) -> str:
        v = os.environ.get(name)
        return default if v is None or not str(v).strip() else str(v).strip()

    title = env("REPORT_TITLE", DEFAULT_TITLE)
    try:
        fdr = float(env("FDR", "0.05"))
    except ValueError:
        raise BadParameter(f"FDR must be a number, got {os.environ.get('FDR')!r}")
    if not 0 < fdr < 1:
        raise BadParameter(f"FDR must be between 0 and 1, got {fdr}")
    out: dict[str, Any] = {"REPORT_TITLE": title, "FDR": fdr}
    for name, default in (("MAX_ITEMS", "20"), ("MAX_FIGURES", "12")):
        try:
            v = int(env(name, default))
        except ValueError:
            raise BadParameter(f"{name} must be an integer, got {os.environ.get(name)!r}")
        if v < 1:
            raise BadParameter(f"{name} must be at least 1, got {v}")
        out[name] = v
    return out


def peak_rss_mb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return round(r / (1024 * 1024) if sys.platform == "darwin" else r / 1024, 1)


class Figures:
    """Copies upstream PNGs into report/figures/ under a section prefix, at most MAX_FIGURES per section."""

    def __init__(self, fig_dir: Path, max_figures: int):
        self.dir = fig_dir
        self.max = max_figures
        self.used: dict[str, int] = {}
        self.omitted: dict[str, int] = {}

    def copy(self, section: str, src: Path | None, prefix: str = "") -> str | None:
        if src is None or not Path(src).is_file():
            return None
        if self.used.get(section, 0) >= self.max:
            self.omitted[section] = self.omitted.get(section, 0) + 1
            return None
        from processor.figures import safe

        name = safe(f"{section}_{prefix + '_' if prefix else ''}{Path(src).stem}") + ".png"
        shutil.copyfile(src, self.dir / name)
        self.used[section] = self.used.get(section, 0) + 1
        return name

    def many(self, section: str, srcs: list[Path], prefix: str = "") -> list[dict[str, str]]:
        out = []
        for p in srcs:
            name = self.copy(section, p, prefix)
            if name:
                out.append({"file": name, "caption": caption(p.stem)})
        return out

    def own(self, section: str, name: str | None) -> str | None:
        if name:
            self.used[section] = self.used.get(section, 0) + 1
        return name


def caption(stem: str) -> str:
    """'KC1_spatial_counts' -> 'KC1: spatial counts' for figure captions."""
    for suffix, text in (("_spatial_counts", "counts on the tissue, kept and removed"), ("_qc", "QC metrics"),
                         ("_domains", "domains"), ("_svg_top", "top spatially variable genes"),
                         ("_cell_type", "cell types"), ("_niches", "niches"),
                         ("_nhood_enrichment", "neighbourhood enrichment"), ("_dotplot", "top interactions"),
                         ("_source_target", "top interactions per source → target"),
                         ("_abundance", "cell-type abundance"), ("_dominant_type", "dominant type"),
                         ("_colocation", "co-location groups")):
        if stem.endswith(suffix):
            return f"{stem[: -len(suffix)]}: {text}"
    if "_local_" in stem:
        sid, pair = stem.split("_local_", 1)
        return f"{sid}: where {pair.replace('_', ' → ', 1)} is co-expressed (local score)"
    if stem.startswith("volcano_"):
        return f"Volcano: {stem[len('volcano_'):]}"
    if stem.startswith("co_occurrence_"):
        return f"Co-occurrence with {stem[len('co_occurrence_'):].replace('_', ' ')} vs distance"
    return stem.replace("_", " ")


def _section(name: str, fn, errors: list[dict[str, str]], *args):
    """Run one section; a failing section is logged, reported on the page, and the rest goes on."""
    t0 = time.monotonic()
    try:
        result = fn(*args)
        log.info("section %s: %s (%.1fs)", name, "rendered" if result else "no inputs", time.monotonic() - t0)
        return result
    except Exception as exc:  # noqa: BLE001 - one broken input must not sink the report
        log.error("section %s FAILED: %s\n%s", name, exc, traceback.format_exc())
        errors.append({"section": name, "error": f"{type(exc).__name__}: {exc}"})
        return None


def main() -> int:
    started = _now()
    t_start = time.monotonic()
    input_dir, output_dir = os.environ.get("INPUT_DIR"), os.environ.get("OUTPUT_DIR")
    if not input_dir or not output_dir:
        log.error("INPUT_DIR and OUTPUT_DIR must be set")
        return 2
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    try:
        params = parameters()
    except BadParameter as exc:
        log.error("bad parameter: %s", exc)
        return 2
    if not input_dir.is_dir():
        log.error("INPUT_DIR %s is not a folder", input_dir)
        return 1
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    fdr, max_items, max_figs = params["FDR"], params["MAX_ITEMS"], params["MAX_FIGURES"]

    from processor import analysis as A, figures as F, inputs as I, report as R

    inp = I.gather(input_dir)
    for port in I.PORTS:
        found = inp.found.get(port, [])
        log.info("input %s: %s", port, ", ".join(f or "." for f in found) if found else "none")
    if not inp.has_any:
        listing = [inp.rel(p) for p in inp.files]
        log.error("none of the expected inputs was found (samples.csv/prep.json, qc_summary.csv, domain_summary.csv, "
                  "svg_consensus.csv, composition_tests.csv, type_summary.csv, cell_type_summary.csv, celltype figures, "
                  "niche_composition.csv, lr_consensus.csv); INPUT_DIR holds %d file(s): %s", len(listing),
                  ", ".join(listing[:40]) + (" …" if len(listing) > 40 else ""))
        return 1

    report_dir, tables_dir = output_dir / "report", output_dir / "tables"
    fig_dir = report_dir / "figures"
    for d in (report_dir, tables_dir, fig_dir):
        d.mkdir(parents=True, exist_ok=True)
    figs = Figures(fig_dir, max_figs)
    errors: list[dict[str, str]] = []
    cond = A.condition_of(inp)

    st = _section("study design", A.study, errors, inp)
    qc = _section("qc", A.qc, errors, inp)
    if qc:
        qc["figure"] = figs.own("qc", _section("qc figure", F.qc_overview, errors, qc["rows"], qc["unit"], fig_dir))
        qc["figures"] = figs.many("qc", qc.pop("plots"))
    dom = _section("domains", A.domains, errors, inp, max_items)
    if dom:
        dom["figure"] = figs.own("domains", _section("domain coverage figure", F.domain_coverage, errors,
                                                     dom.get("coverage"), cond, dom["unit"], fig_dir))
        dom["figures"] = [{"file": figs.copy("domains", p), "caption": cap} for n, cap in (
            ("umap", "UMAP of the BANKSY embedding, by domain, sample and annotation"),
            ("markers_dotplot", "Top marker genes per domain"),
            ("confusion", "Annotation × domain"),
            ("domain_composition", "Domain fractions per sample")) if (p := dom["upstream_figures"].get(n))]
        dom["figures"] = [f for f in dom["figures"] if f["file"]] + figs.many("domains", dom.pop("maps"))
    sv = _section("svg", A.svg, errors, inp, fdr, max_items)
    if sv:
        sv["figures"] = [{"file": figs.copy("svg", p), "caption": cap} for n, cap in (
            ("consensus_top", "Consensus genes × samples: Moran's I, dot = significant"),
            ("morans_distribution", "Moran's I of every tested gene per sample")) if (p := sv["upstream_figures"].get(n))]
        sv["figures"] = [f for f in sv["figures"] if f["file"]] + figs.many("svg", sv.pop("maps"))
    cmp_ = _section("comparison", A.compare, errors, inp, fdr, max_items)
    if cmp_:
        for r in cmp_["runs"]:
            sec = r["section_key"] = f"compare_{F.safe(r['label'])}"
            others = [l for l in r["levels"] if l != r["reference"]]
            r["figure"] = figs.own(sec, _section("DE count figure", F.de_counts, errors, r.get("de_summary") or [],
                                                 r["label_column"], (others[0], r["reference"]) if len(others) == 1 else None,
                                                 fig_dir))
            r["figures"] = [{"file": figs.copy(sec, p, ""), "caption": cap} for n, cap in (
                ("composition", "Share of each label per sample and condition (adjusted p in the titles)"),
                ("de_heatmap", "Top genes of every DE test: log2 fold change, dot = significant"))
                if (p := r["figures_upstream"].get(n))]
            r["figures"] = [f for f in r["figures"] if f["file"]] + figs.many(sec, r.pop("volcanoes"))
            r.pop("figures_upstream", None)
    ct = _section("cell types", A.celltypes, errors, inp, max_items)
    if ct:
        up = ct.get("upstream_figures") or {}
        ct["figures"] = [{"file": figs.copy("celltype", p), "caption": cap} for n, cap in (
            ("cell_type_by_sample", "Cell-type composition per sample"),
            ("umap_cell_type", "UMAP by cell type"), ("umap_cluster", "UMAP by cluster"),
            ("dotplot_marker_sets", "Marker-set scores per cluster")) if (p := up.get(n))]
        ct["figures"] = [f for f in ct["figures"] if f["file"]]
        if "cell_type_by_sample" not in up and ct.get("composition") is not None:
            name = _section("cell-type composition figure", F.composition_bars, errors, ct["composition"],
                            "Cell-type composition per sample", fig_dir, "celltype_composition.png")
            if name:
                ct["figures"].insert(0, {"file": figs.own("celltype", name), "caption": "Cell-type composition per sample"})
        ct["figures"] += figs.many("celltype", ct.pop("maps", []))
    dc = _section("deconvolution", A.deconvolution, errors, inp, max_items)
    if dc:
        up = dc.get("upstream_figures") or {}
        dc["figures"] = [{"file": figs.copy("deconvolution", p), "caption": cap} for n, cap in (
            ("type_composition", "Estimated cell-type composition per sample"),
            ("colocation_weights", "Co-location groups: cell-type loadings"),
            ("reconstruction_qc", "Observed vs expected counts"),
            ("training_history", "Training history (ELBO)")) if (p := up.get(n))]
        dc["figures"] = [f for f in dc["figures"] if f["file"]]
        if "type_composition" not in up and dc.get("composition") is not None:
            name = _section("deconvolution composition figure", F.composition_bars, errors, dc["composition"],
                            "Estimated cell-type composition per sample", fig_dir, "deconvolution_composition.png")
            if name:
                dc["figures"].insert(0, {"file": figs.own("deconvolution", name), "caption": "Estimated composition"})
        dc["figures"] += figs.many("deconvolution", dc.pop("maps"))
    ni = _section("niches", A.niches, errors, inp, max_items)
    if ni:
        for r in ni["runs"]:
            sec = r["section_key"] = f"niches_{F.safe(r['label'])}"
            r["figures"] = [{"file": figs.copy(sec, p), "caption": cap} for n, cap in (
                ("nhood_enrichment", "Neighbourhood enrichment z-scores, pooled over samples"),
                ("niche_composition", "Composition of each niche")) if (p := r["upstream_figures"].get(n))]
            r["figures"] = [f for f in r["figures"] if f["file"]] + figs.many(sec, r.pop("maps")) + \
                figs.many(sec, r.pop("cooccurrence"))
    cc = _section("cell-cell communication", A.ccc, errors, inp, max_items)
    if cc:
        for r in cc["runs"]:
            sec = r["section_key"] = f"ccc_{F.safe(r['label'])}"
            r["figure"] = figs.own(sec, _section("ccc spatial figure", F.ccc_spatial, errors, r.get("spatial"),
                                                 r["label"], fig_dir))
            r["figures"] = figs.many(sec, r.pop("dotplots")) + figs.many(sec, r.pop("source_target")) + \
                figs.many(sec, r.pop("local"))

    app_version = version_of(REPO_ROOT)
    tools = tool_versions()
    own = {"app": APP_NAME, "version": app_version, "tools": tools, "parameters": params}
    meth = _section("methods", A.methods, errors, inp, own) or {"apps": [], "references": []}
    credit = _section("data credit", A.demo_datasets, errors, inp) or []

    workflows = []
    if dom or cmp_ or (sv and not (ct or dc or ni or cc)):
        workflows.append(WORKFLOWS["A"])
    if dc:
        workflows.append(WORKFLOWS["B"])
    if ct:
        workflows.append(WORKFLOWS["C"])
    if ni and not (dc or ct):
        lc = {str(r["label_column"]) for r in ni["runs"]}
        workflows.append(WORKFLOWS["B"] if "dominant_type" in lc else WORKFLOWS["C"])
    if cc and not ct and WORKFLOWS["C"] not in workflows:
        workflows.append(WORKFLOWS["C"])

    generated_at = started.isoformat(timespec="seconds")
    context = {
        "title": params["REPORT_TITLE"], "generated_at": generated_at, "app_version": app_version, "tools": tools,
        "params": params, "st": st, "qc": qc, "dom": dom, "sv": sv, "cmp": cmp_, "ct": ct, "dc": dc, "ni": ni,
        "cc": cc, "methods": meth, "credit": credit, "errors": errors, "workflows": workflows,
        "omitted": figs.omitted, "max_items": max_items,
    }
    R.write_index(report_dir, context)
    written = R.write_tables(tables_dir, context)
    summary = R.summary_json(context, own)
    with open(tables_dir / "summary.json", "w") as fh:
        json.dump(A.clean(summary), fh, indent=2, ensure_ascii=False)
    written.append("summary.json")

    finished = _now()
    sections = {"study_design": st, "qc": qc, "domains": dom, "svg": sv, "comparison": cmp_, "cell_types": ct,
                "deconvolution": dc, "niches": ni, "cell_cell_communication": cc}
    provenance = {
        "app": APP_NAME, "version": app_version, "tools": tools, "parameters": params,
        "startedAt": generated_at, "finishedAt": finished.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - t_start, 2), "peakRssMb": peak_rss_mb(),
        "runtime": "lambda" if os.environ.get("AWS_LAMBDA_FUNCTION_NAME") else "container",
        "inputsFound": inp.found, "provenanceFiles": [inp.rel(p) for p, _ in inp.provenances],
        "workflows": workflows, "sections": {k: bool(v) for k, v in sections.items()}, "sectionErrors": errors,
        "figuresOmitted": figs.omitted,
        "outputs": {"report": sorted(str(p.relative_to(report_dir)) for p in report_dir.rglob("*") if p.is_file())
                    + ["provenance.json"], "tables": written + ["provenance.json"]},
    }
    for d in (report_dir, tables_dir):
        with open(d / "provenance.json", "w") as fh:
            json.dump(A.clean(provenance), fh, indent=2, ensure_ascii=False)
    log.info("wrote %s (%d section(s), %d error(s), %d figure(s)) and %d table(s) in %.1fs, peak RSS %.0f MB",
             report_dir / "index.html", sum(bool(v) for v in sections.values()), len(errors),
             len(list(fig_dir.glob("*.png"))), len(written), time.monotonic() - t_start, peak_rss_mb())
    return 0


if __name__ == "__main__":
    sys.exit(main())
