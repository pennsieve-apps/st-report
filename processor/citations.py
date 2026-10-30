"""
The papers behind each step, keyed by app name, and the report's own.

The bibliographic data (full author lists, journal, volume, pages, PMID) is in
references.json, written by scripts/verify_citations.py from Crossref / DataCite and
checked against PubMed (title, first author and year). Each app's list follows the
methods papers its own CITATION.cff cites (the demo datasets go to the Input data section
instead). The Citations section lists, deduplicated, the papers of every app whose
provenance.json is among the inputs, plus the report's own.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

REFERENCES = Path(__file__).parent / "references.json"

# What each step is, and the papers it rests on (most important first).
APPS: dict[str, dict[str, Any]] = {
    "st-prep": {"step": "Import (Space Ranger, Xenium, MERSCOPE, CosMx, AnnData, SpatialData)",
                "cite": ["scanpy", "anndata", "spatialdata", "visium"]},
    "st-qc": {"step": "Quality control and normalization",
              "cite": ["heumos", "pipecomp", "seurat_v3", "scanpy", "anndata"]},
    "st-domains": {"step": "Spatial domains (BANKSY)",
                   "cite": ["banksy", "harmony", "leiden", "umap", "sdmbench", "scanpy", "anndata"]},
    "st-svg": {"step": "Spatially variable genes (Moran's I)",
               "cite": ["moran", "geary", "bh", "spatialde", "squidpy", "scanpy", "anndata"]},
    "st-compare": {"step": "Condition comparison (composition, pseudobulk DE)",
                   "cite": ["pydeseq2", "deseq2", "squair", "sccoda", "aitchison", "bh", "scipy", "anndata"]},
    "cell2location": {"step": "Cell-type deconvolution (cell2location)",
                      "cite": ["cell2location", "scvi_tools", "scanpy", "anndata"]},
    "st-celltype": {"step": "Clustering and marker-based cell types",
                    "cite": ["harmony", "leiden", "seurat_2015", "panglaodb", "omnipath", "umap", "scanpy", "anndata"]},
    "st-niches": {"step": "Neighbourhood enrichment and niches",
                  "cite": ["squidpy", "schurch", "cellcharter", "lloyd", "rousseeuw", "scanpy", "anndata"]},
    "st-ccc": {"step": "Cell-cell communication (LIANA+)",
               "cite": ["liana_plus", "liana", "rra", "cellphonedb", "connectome", "natmi", "singlecellsignalr",
                        "moran", "scanpy", "anndata"]},
    "st-report": {"step": "This report", "cite": ["moses_pachter", "spotless", "matplotlib"]},
}

# Pipeline order for the methods section.
ORDER = list(APPS)

# Demo datasets, recognised by the sample paths st-prep recorded.
DATASETS = [
    {"key": "tls", "match": ("TLS_VISIUM_USZ",), "platform": "visium",
     "text": "Dawo S, Nonchev K, Silina K. 10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung (5) "
             "Cancer with Tertiary Lymphoid Structures. Zenodo (2025).",
     "doi": "10.5281/zenodo.14620362", "licence": "CC BY 4.0", "ref": "tls_data"},
    {"key": "lymph", "match": ("V1_Human_Lymph_Node", "Human_Lymph_Node"), "platform": "visium",
     "text": "10x Genomics. Human Lymph Node (Visium Spatial Gene Expression, Space Ranger 1.1.0), "
             "10x Genomics dataset page.",
     "doi": "", "licence": "CC BY 4.0", "ref": None,
     "extra": "With cell2location, the scRNA-seq reference is typically the Tabula Sapiens lymph node "
              "dataset (Tabula Sapiens Consortium, Science 2022), CC BY 4.0.", "extra_ref": "tabula_sapiens"},
    {"key": "xenium", "match": ("Rep1", "Rep2"), "platform": "xenium",
     "text": "10x Genomics. Xenium In Situ: FFPE Human Breast Cancer, replicates 1 and 2 (Xenium Onboard "
             "Analysis 1.0.1), 10x Genomics dataset page; described in Janesick A, et al. Nat Commun 2023.",
     "doi": "", "licence": "CC BY 4.0", "ref": "janesick"},
]


@lru_cache(maxsize=1)
def references() -> dict[str, dict[str, Any]]:
    with open(REFERENCES, encoding="utf-8") as fh:
        return json.load(fh)


def reference(key: str) -> dict[str, Any]:
    r = dict(references()[key])
    r["key"] = key
    r["url"] = f"https://doi.org/{r['doi']}"
    vol = str(r.get("volume") or "")
    issue = str(r.get("issue") or "")
    pages = str(r.get("pages") or "")
    where = vol + (f"({issue})" if issue and vol else "") + (f":{pages}" if pages else "")
    r["where"] = where
    r["journal_display"] = r.get("journal_short") or r.get("journal") or ""
    r["text"] = (f"{r['short_authors']}. {r['title']}. {r['journal_display']} {r['year']}"
                 + (f";{where}" if where else "") + ".")
    return r


def keys_for(app: str) -> list[str]:
    return list(APPS.get(app, {}).get("cite", []))
