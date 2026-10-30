"""
Finding the upstream outputs.

Linked layout: each input port arrives in INPUT_DIR/<port>/ (placement: subfolder).
Flat layout (no links, everything merged side by side): files are found by name anywhere
under INPUT_DIR. Either way an upstream output folder is recognised by the files that only
that app writes (svg_consensus.csv is st-svg's, niche_composition.csv is st-niches', ...),
and its summary.json / provenance.json / figures are read from the same folder. INPUT_DIR
entries are often symlinks, so the walk follows them.

st-compare, st-niches and st-ccc write under a <LABEL>/ folder and can be linked more than
once into the same port: every folder holding their files is one run, keyed by its label.
The large per-observation files (labels.csv, abundance.csv, scores.csv, removed.csv,
lr_cluster.csv, co_occurrence.csv) are never read: the report works from the summaries.
"""
from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

log = logging.getLogger("st-report")

PORTS = ["metadata", "qc", "domains", "svg", "compare", "deconvolution", "celltype", "celltype-figures",
         "niches", "ccc"]

# Files that identify a port's upstream folder (any one of them is enough).
MARKERS: dict[str, tuple[str, ...]] = {
    "qc": ("qc_summary.csv", "thresholds.csv"),
    "domains": ("domain_summary.csv", "domain_markers.csv"),
    "svg": ("svg_consensus.csv", "svg.csv"),
    "compare": ("composition_tests.csv", "de_summary.csv"),
    "deconvolution": ("type_summary.csv", "colocation.csv", "abundance.csv"),
    "celltype": ("cell_type_summary.csv", "cluster_annotation.csv"),
    "celltype-figures": ("umap_cell_type.png", "cell_type_by_sample.png", "dotplot_marker_sets.png",
                         "umap_cluster.png"),
    "niches": ("niche_composition.csv", "nhood_enrichment.csv"),
    "ccc": ("lr_consensus.csv", "lr_spatial.csv"),
}
# Ports whose app writes under <LABEL>/ and may be linked several times.
MULTI = ("compare", "niches", "ccc")


def walk(root: Path) -> list[Path]:
    """Every file under root, following symlinked folders, in a stable order."""
    out: list[Path] = []
    seen: set[str] = set()
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        real = os.path.realpath(dirpath)
        if real in seen:  # a symlink loop
            dirnames[:] = []
            continue
        seen.add(real)
        dirnames.sort()
        for name in sorted(filenames):
            p = Path(dirpath) / name
            if p.is_file():
                out.append(p)
    return out


def read_json(path: Path | None) -> dict[str, Any] | None:
    if path is None or not Path(path).is_file():
        return None
    try:
        with open(path) as fh:
            value = json.load(fh)
        return value if isinstance(value, dict) else None
    except (OSError, ValueError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None


def read_csv(path: Path | None, *cols: str, **kwargs) -> pd.DataFrame | None:
    """A CSV, or None when it is missing, unreadable or lacks one of `cols`."""
    if path is None or not Path(path).is_file():
        return None
    try:
        df = pd.read_csv(path, **kwargs)
    except (OSError, ValueError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
        log.warning("cannot read %s: %s", path, exc)
        return None
    missing = [c for c in cols if c not in df.columns]
    if missing:
        log.warning("%s: not the expected table (no %s); ignored", path, ", ".join(missing))
        return None
    return df


def _under(path: Path, folder: Path) -> bool:
    try:
        path.relative_to(folder)
        return True
    except ValueError:
        return False


@dataclass
class Run:
    """One upstream output folder."""
    port: str
    dir: Path
    label: str = ""
    summary: dict[str, Any] | None = None
    provenance: dict[str, Any] | None = None

    def file(self, name: str) -> Path | None:
        p = self.dir / name
        return p if p.is_file() else None

    def csv(self, name: str, *cols: str, **kwargs) -> pd.DataFrame | None:
        return read_csv(self.file(name), *cols, **kwargs)

    def pngs(self, sub: str = "") -> list[Path]:
        d = self.dir / sub if sub else self.dir
        return sorted(p for p in d.glob("*.png") if p.is_file()) if d.is_dir() else []


@dataclass
class Inputs:
    input_dir: Path
    files: list[Path] = field(default_factory=list)
    samples: pd.DataFrame | None = None
    prep: dict[str, Any] | None = None
    metadata_dir: Path | None = None
    runs: dict[str, list[Run]] = field(default_factory=dict)
    provenances: list[tuple[Path, dict[str, Any]]] = field(default_factory=list)
    found: dict[str, list[str]] = field(default_factory=dict)

    def rel(self, path: Path | None) -> str | None:
        if path is None:
            return None
        try:
            return str(Path(path).relative_to(self.input_dir))
        except ValueError:
            return str(path)

    def one(self, port: str) -> Run | None:
        runs = self.runs.get(port) or []
        return runs[0] if runs else None

    def many(self, port: str) -> list[Run]:
        return list(self.runs.get(port) or [])

    @property
    def has_any(self) -> bool:
        return self.samples is not None or bool(self.prep) or any(self.runs.get(p) for p in PORTS)


def _port_files(inp: Inputs, port: str) -> list[Path]:
    """Files in the port's own folder when it exists and holds any, else every file (flat layout)."""
    d = inp.input_dir / port
    if d.is_dir():
        inside = [f for f in inp.files if _under(f, d)]
        if inside:
            return inside
    return inp.files


def _dirs(inp: Inputs, port: str) -> list[Path]:
    names = set(MARKERS[port])
    files = _port_files(inp, port)
    dirs: list[Path] = []
    for f in files:
        if f.name in names and f.parent not in dirs:
            dirs.append(f.parent)
    if port == "celltype-figures":
        # st-celltype's labels folder can hold figures in the flat layout; never take it for figures.
        dirs = [d for d in dirs if not (d / "cell_type_summary.csv").is_file()] or dirs
    return dirs


def gather(input_dir: Path) -> Inputs:
    inp = Inputs(input_dir=input_dir, files=walk(input_dir))

    # metadata (st-prep): samples.csv with the contract's columns, prep.json
    for f in _port_files(inp, "metadata"):
        if f.name == "samples.csv" and inp.samples is None:
            df = read_csv(f, "sample_id")
            if df is not None and ({"platform", "n_obs"} & set(df.columns)):
                df["sample_id"] = df["sample_id"].astype(str)
                inp.samples, inp.metadata_dir = df, f.parent
        elif f.name == "prep.json" and inp.prep is None:
            prep = read_json(f)
            if prep is not None and isinstance(prep.get("samples"), dict):
                inp.prep = prep
                inp.metadata_dir = inp.metadata_dir or f.parent
    inp.found["metadata"] = [inp.rel(p) for p in ((inp.metadata_dir / "samples.csv"), (inp.metadata_dir / "prep.json"))
                             if p.is_file()] if inp.metadata_dir else []

    for port in MARKERS:
        runs: list[Run] = []
        for d in _dirs(inp, port):
            run = Run(port=port, dir=d, summary=read_json(d / "summary.json"),
                      provenance=read_json(d / "provenance.json"))
            label = ""
            if run.summary:
                label = str(run.summary.get("label") or "")
            if not label and run.provenance:
                label = str((run.provenance.get("parameters") or {}).get("LABEL") or "")
            if not label and d != inp.input_dir / port:
                label = d.name
            run.label = label or port
            if port in MULTI and any(r.label == run.label for r in runs):
                run.label = f"{run.label} ({inp.rel(d)})"
            runs.append(run)
        if len(runs) > 1 and port not in MULTI:
            log.warning("%s: %d folders found (%s); using the first", port, len(runs),
                        ", ".join(inp.rel(r.dir) for r in runs))
            runs = runs[:1]
        inp.runs[port] = runs
        inp.found[port] = [inp.rel(r.dir) or "." for r in runs]

    # every provenance.json (methods & citations), once per distinct file content
    for f in inp.files:
        if f.name == "provenance.json":
            prov = read_json(f)
            if prov and prov.get("app"):
                inp.provenances.append((f, prov))
    return inp
