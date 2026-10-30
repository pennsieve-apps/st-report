"""End-to-end: the synthetic upstream outputs of each workflow (and combinations, and none)
through the processor, in both input layouts, plus parameters, the Lambda handler and an
HTML snapshot."""
from __future__ import annotations

import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pytest

from processor import citations, handler, inputs
from processor import main as app
from tests.make_fixture import make_fixture

SNAPSHOT = Path(__file__).parent / "snapshots" / "index.html"
FIXED_TIME = datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc)
ANCHORS = {"design", "qc", "domains", "svg", "compare", "celltypes", "organisation", "ccc"}


@pytest.fixture(scope="session")
def fixtures(tmp_path_factory) -> dict[str, dict[str, Path]]:
    root = tmp_path_factory.mktemp("fixture")
    return {w: make_fixture(root / w, "both", w) for w in ("a", "b", "c", "ab", "abc")}


def run(input_dir: Path, output_dir: Path, monkeypatch, **params) -> int:
    monkeypatch.setenv("INPUT_DIR", str(input_dir))
    monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    for name in app.PARAMETERS:
        monkeypatch.delenv(name, raising=False)
    for name, value in params.items():
        monkeypatch.setenv(name, value)
    return app.main()


def page(out: Path) -> str:
    return (out / "report" / "index.html").read_text()


def anchors(html: str) -> set[str]:
    return set(re.findall(r'<h2 id="([^"]+)"', html)) & ANCHORS


def summary(out: Path) -> dict:
    return json.loads((out / "tables" / "summary.json").read_text())


def assert_self_contained(out: Path) -> None:
    html = page(out)
    assert "<script" not in html.lower()
    assert not re.search(r'(src|href)="(?!https://doi\.org/|#|figures/)', html), "only local figures and DOI links"
    assert "<link" not in html and "@import" not in html
    for src in re.findall(r'(?:src|href)="(figures/[^"]+)"', html):
        assert (out / "report" / src).is_file(), src


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_workflow_a(fixtures, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixtures["a"][layout], out, monkeypatch) == 0
    html = page(out)
    assert anchors(html) == {"design", "qc", "domains", "svg", "compare"}
    assert_self_contained(out)
    s = summary(out)
    assert s["workflows"] == ["Spatial Domain Mapping & Condition Comparison"]
    assert s["study"]["conditions"] == {"ctrl": 2, "treat": 2}
    assert any("donor" in c and "confounded" in c for c in s["study"]["caveats"])
    assert s["qc"]["n_removed"] == 50 and any("s3 lost more" in o for o in s["qc"]["observations"])
    assert any("8-fold" in o for o in s["qc"]["observations"])  # depth 500 .. 4000
    dom = s["domains"]
    assert dom["n_domains"] == 4 and any("Confined to a single sample: D04" in o for o in dom["observations"])
    assert any("D04 (treat: s3)" in o for o in dom["observations"])
    assert dom["agreement"]["ari"] == pytest.approx(0.35)
    assert s["svg"]["consensus_top"][:3] == ["IGHG1", "SFTPB", "COL1A1"] or set(s["svg"]["consensus_top"][:3]) == {"IGHG1", "SFTPB", "COL1A1"}
    assert any("top-10 markers" in o for o in s["svg"]["observations"])
    runs = {r["label"]: r for r in s["comparison"]}
    assert set(runs) == {"domain", "annotation"}
    assert runs["domain"]["significant_labels"] == ["D02"]
    assert runs["annotation"]["composition_n_significant"] == 0
    assert any("no label's share differs significantly" in o and "2 vs 2" in o for o in runs["annotation"]["observations"])
    assert any("Not tested" in o and "D04" in o for o in runs["domain"]["observations"])
    assert any("change the same way" in o and "IFI6" in o for o in runs["domain"]["observations"])
    tables = {p.name for p in (out / "tables").iterdir()}
    assert {"samples.csv", "qc.csv", "domains.csv", "domain_coverage.csv", "svg_consensus.csv", "composition_tests.csv",
            "de_summary.csv", "de_top_genes.csv", "methods.csv", "citations.csv", "summary.json"} <= tables
    comp = pd.read_csv(out / "tables" / "composition_tests.csv")
    assert set(comp["label_set"]) == {"domain", "annotation"}
    cov = pd.read_csv(out / "tables" / "domain_coverage.csv")
    assert len(cov) == 16 and cov.groupby("sample_id")["fraction"].sum().round(6).eq(1).all()
    for f in ("qc_overview.png", "domain_coverage.png", "de_counts_domain.png"):
        assert (out / "report" / "figures" / f).is_file(), f
    # Methods & citations: every step's provenance, deduplicated references, the report's own last.
    apps = [m["app"] for m in s["methods"]]
    assert apps == ["st-prep", "st-qc", "st-domains", "st-svg", "st-compare", "st-report"]
    dois = [c["doi"] for c in s["citations"]]
    assert len(dois) == len(set(dois)) and "10.1038/s41588-024-01664-3" in dois and "10.1093/bioinformatics/btad547" in dois
    assert "10.1038/s41587-021-01139-4" not in dois  # no cell2location in this workflow
    assert "PMID 38413725" in html  # BANKSY


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_workflow_b(fixtures, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixtures["b"][layout], out, monkeypatch) == 0
    html = page(out)
    assert anchors(html) == {"design", "qc", "celltypes", "organisation"}
    assert_self_contained(out)
    s = summary(out)
    assert s["workflows"] == ["Cell-type Deconvolution & Niches"]
    dc = s["deconvolution"]
    assert dc["n_cell_types"] == 4 and dc["types"][0]["cell_type"] == "T_cell"
    assert any("Myeloid" in o and "dominate no" in o for o in dc["observations"])
    assert any("1,000 epochs" in c for c in dc["caveats"])
    assert any("Mixed" in o for o in dc["observations"])
    groups = {g["group"]: g for g in dc["colocation"]}
    assert groups["G01"]["fraction_strongest"] == pytest.approx(0.5) and groups["G03"]["fraction_strongest"] == pytest.approx(0.1)
    ni = s["niches"][0]
    assert ni["label"] == "colocation-group" and any("co-location groups" in c for c in ni["caveats"])
    assert any("G01 (T_cell / Stromal)" in o for o in ni["observations"])  # groups named by their main types
    assert "10.1038/s41587-021-01139-4" in [c["doi"] for c in s["citations"]]


def test_mostly_mixed_dominant_type(fixtures, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixtures["b"]["linked"], src)
    sm = json.loads((src / "deconvolution" / "summary.json").read_text())
    sm["dominant_type_counts"] = {"B_cell": 20, "Mixed": 1980}
    (src / "deconvolution" / "summary.json").write_text(json.dumps(sm))
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    dc = summary(out)["deconvolution"]
    assert "Mixed" in dc["caveats"][0] and "99%" in dc["caveats"][0]
    assert not any("Dominant type" in o for o in dc["observations"])


@pytest.mark.parametrize("layout", ["linked", "flat"])
def test_workflow_c(fixtures, tmp_path, monkeypatch, layout):
    out = tmp_path / "out"
    assert run(fixtures["c"][layout], out, monkeypatch) == 0
    html = page(out)
    assert anchors(html) == {"design", "qc", "celltypes", "organisation", "ccc"}
    assert_self_contained(out)
    s = summary(out)
    assert s["workflows"] == ["Single-cell Spatial Profiling"]
    ct = s["cell_types"]
    assert any("Close calls" in o and "C01" in o and "C05" in o for o in ct["observations"])
    assert any("Epithelial (C01, C05)" in o for o in ct["observations"])
    assert any("7 of 10" in o for o in ct["observations"])
    assert any("targeted panel" in c for c in s["study"]["caveats"])
    ni = s["niches"][0]
    assert any("sparse" in c for c in ni["caveats"])  # mean degree 3.2
    cc = s["cell_cell_communication"][0]
    assert cc["anti_correlated"] == ["S100A4 → ERBB2"]
    assert any("Only 23 ligand–receptor pairs" in c for c in cc["caveats"])
    assert (out / "report" / "figures" / "ccc_spatial_cell-type.png").is_file()
    tables = {p.name for p in (out / "tables").iterdir()}
    assert {"cell_types.csv", "clusters.csv", "niches.csv", "neighbourhood_enrichment.csv", "ccc_top.csv",
            "ccc_spatial.csv"} <= tables


def test_all_inputs(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["abc"]["linked"], out, monkeypatch) == 0
    html = page(out)
    assert anchors(html) == ANCHORS
    s = summary(out)
    assert len(s["workflows"]) == 3 and not s["sectionErrors"]
    assert {n["label"] for n in s["niches"]} == {"colocation-group", "cell-type"}
    prov = json.loads((out / "report" / "provenance.json").read_text())
    assert all(prov["sections"].values())
    assert_self_contained(out)


def test_two_workflows(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["ab"]["flat"], out, monkeypatch) == 0
    assert summary(out)["workflows"] == ["Spatial Domain Mapping & Condition Comparison",
                                         "Cell-type Deconvolution & Niches"]


@pytest.mark.parametrize("keep", ["metadata", "qc", "domains", "svg", "compare", "deconvolution", "celltype",
                                  "celltype-figures", "niches", "ccc"])
def test_each_port_alone(fixtures, tmp_path, monkeypatch, keep):
    src = tmp_path / "in"
    shutil.copytree(fixtures["abc"]["linked"] / keep, src / keep)
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    s = summary(out)
    assert not s["sectionErrors"], s["sectionErrors"]
    assert_self_contained(out)


def test_all_empty(tmp_path, monkeypatch):
    (tmp_path / "in" / "other").mkdir(parents=True)
    (tmp_path / "in" / "other" / "readme.txt").write_text("nothing here")
    assert run(tmp_path / "in", tmp_path / "out", monkeypatch) == 1
    (tmp_path / "empty").mkdir()
    assert run(tmp_path / "empty", tmp_path / "out2", monkeypatch) == 1


def test_unreadable_input_is_reported(fixtures, tmp_path, monkeypatch):
    src = tmp_path / "in"
    shutil.copytree(fixtures["a"]["linked"], src)
    (src / "svg" / "summary.json").write_text("{not json")
    (src / "domains" / "domain_summary.csv").write_text("a,b\n1,2\n")  # wrong columns: falls back to summary.json
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    s = summary(out)
    assert "svg" in s and "domains" in s


def test_parameters(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    assert run(fixtures["a"]["linked"], out, monkeypatch, FDR="0.001", REPORT_TITLE="My Study", MAX_ITEMS="2",
               MAX_FIGURES="1") == 0
    html = page(out)
    assert "<title>My Study</title>" in html
    s = summary(out)
    runs = {r["label"]: r for r in s["comparison"]}
    assert runs["domain"]["significant_labels"] == []  # adj. p 0.004 > 0.001
    de = pd.read_csv(out / "tables" / "de_summary.csv")
    assert de[de["status"] == "tested"]["n_significant"].max() == 2  # only the 1e-8 and 1e-6 genes
    assert len(pd.read_csv(out / "tables" / "svg_consensus.csv")) == 2
    prov = json.loads((out / "report" / "provenance.json").read_text())
    assert prov["figuresOmitted"].get("qc", 0) > 0 and "more figure(s) not shown" in html


@pytest.mark.parametrize("name,value", [("FDR", "2"), ("FDR", "x"), ("MAX_ITEMS", "0"), ("MAX_FIGURES", "a")])
def test_bad_parameters(fixtures, tmp_path, monkeypatch, name, value):
    assert run(fixtures["a"]["linked"], tmp_path / "out", monkeypatch, **{name: value}) == 2


def test_symlinked_inputs(fixtures, tmp_path, monkeypatch):
    src = tmp_path / "in"
    src.mkdir()
    for port in fixtures["c"]["linked"].iterdir():
        os.symlink(port, src / port.name)
    out = tmp_path / "out"
    assert run(src, out, monkeypatch) == 0
    assert 'id="ccc"' in page(out)


def test_two_compare_links_with_one_label(fixtures, tmp_path):
    src = tmp_path / "in"
    shutil.copytree(fixtures["a"]["linked"], src)
    shutil.copytree(src / "compare" / "domain", src / "compare" / "domain-again")
    found = inputs.gather(src)
    assert [r.label for r in found.many("compare")] == ["annotation", "domain", "domain (compare/domain-again)"]


def test_handler(fixtures, tmp_path, monkeypatch):
    out = tmp_path / "out"
    monkeypatch.setenv("FDR", "0.5")  # left over from a warm start: must be cleared
    result = handler.handler({"inputDir": str(fixtures["c"]["linked"]), "outputDir": str(out),
                              "REPORT_TITLE": "Lambda run"}, None)
    assert result["status"] == "ok" and "FDR" not in os.environ
    assert "<title>Lambda run</title>" in page(out)
    with pytest.raises(RuntimeError):
        handler.handler({"inputDir": str(tmp_path / "missing"), "outputDir": str(out)}, None)


def test_references_complete():
    refs = citations.references()
    for app_name, entry in citations.APPS.items():
        for key in entry["cite"]:
            r = refs[key]
            assert r["doi"] and r["title"] and r["authors"] and r["year"], key
    for ds in citations.DATASETS:
        for k in (ds.get("ref"), ds.get("extra_ref")):
            assert k is None or k in refs


def test_html_snapshot(fixtures, tmp_path, monkeypatch):
    """The page for workflow A's fixture, byte for byte. After an intended change: UPDATE_SNAPSHOTS=1 make test."""
    monkeypatch.setattr(app, "_now", lambda: FIXED_TIME)
    monkeypatch.setattr(app, "tool_versions", lambda: {"python": "3.12", "pandas": "x", "numpy": "x",
                                                        "matplotlib": "x", "jinja2": "x"})
    out = tmp_path / "out"
    assert run(fixtures["a"]["linked"], out, monkeypatch) == 0
    html = page(out)
    if os.environ.get("UPDATE_SNAPSHOTS") or not SNAPSHOT.exists():
        SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
        SNAPSHOT.write_text(html)
    assert html == SNAPSHOT.read_text(), "index.html changed; inspect the diff, then UPDATE_SNAPSHOTS=1 make test"
