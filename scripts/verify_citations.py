#!/usr/bin/env python3
"""
Development tool (not run by the processor): fetches every DOI the report can cite from
Crossref (journal articles) or DataCite (Zenodo, arXiv), looks up its PMID with PubMed
E-utilities (esearch by DOI, then esummary: title, first author and year must match), and
writes processor/references.json with full author lists. CITATION.cff is generated from the
same file (--cff).

  python scripts/verify_citations.py            # refresh processor/references.json
  python scripts/verify_citations.py --cff      # also rewrite CITATION.cff's references

Needs network. The DOIs are in DOIS below: the methods papers each upstream app cites in its
own CITATION.cff, the report's own references, and the demo datasets.
"""
from __future__ import annotations

import html
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "processor" / "references.json"

# key -> DOI. Keys are what processor/citations.py refers to.
DOIS = {
    # shared
    "scanpy": "10.1186/s13059-017-1382-0",
    "anndata": "10.21105/joss.04371",
    "squidpy": "10.1038/s41592-021-01358-2",
    "spatialdata": "10.1038/s41592-024-02212-x",
    "visium": "10.1126/science.aaf2403",
    "scipy": "10.1038/s41592-019-0686-2",
    "bh": "10.1111/j.2517-6161.1995.tb02031.x",
    # st-qc
    "heumos": "10.1038/s41576-023-00586-w",
    "pipecomp": "10.1186/s13059-020-02136-7",
    "seurat_v3": "10.1016/j.cell.2019.05.031",
    # st-domains / st-celltype
    "banksy": "10.1038/s41588-024-01664-3",
    "harmony": "10.1038/s41592-019-0619-0",
    "leiden": "10.1038/s41598-019-41695-z",
    "umap": "10.21105/joss.00861",
    "sdmbench": "10.1038/s41592-024-02215-8",
    "seurat_2015": "10.1038/nbt.3192",
    "panglaodb": "10.1093/database/baz046",
    "omnipath": "10.15252/msb.20209923",
    # st-svg
    "moran": "10.1093/biomet/37.1-2.17",
    "geary": "10.2307/2986645",
    "spatialde": "10.1038/nmeth.4636",
    # st-compare
    "deseq2": "10.1186/s13059-014-0550-8",
    "pydeseq2": "10.1093/bioinformatics/btad547",
    "squair": "10.1038/s41467-021-25960-2",
    "sccoda": "10.1038/s41467-021-27150-6",
    "aitchison": "10.1111/j.2517-6161.1982.tb01195.x",
    # st-ccc
    "liana_plus": "10.1038/s41556-024-01469-w",
    "liana": "10.1038/s41467-022-30755-0",
    "rra": "10.1093/bioinformatics/btr709",
    "cellphonedb": "10.1038/s41596-020-0292-x",
    "connectome": "10.1038/s41598-022-07959-x",
    "natmi": "10.1038/s41467-020-18873-z",
    "singlecellsignalr": "10.1093/nar/gkaa183",
    # st-niches
    "schurch": "10.1016/j.cell.2020.07.005",
    "cellcharter": "10.1038/s41588-023-01588-4",
    "lloyd": "10.1109/TIT.1982.1056489",
    "rousseeuw": "10.1016/0377-0427(87)90125-7",
    # cell2location
    "cell2location": "10.1038/s41587-021-01139-4",
    "scvi_tools": "10.1038/s41587-021-01206-w",
    # the report's own
    "moses_pachter": "10.1038/s41592-022-01409-2",
    "spotless": "10.7554/eLife.88431",
    "matplotlib": "10.1109/MCSE.2007.55",
    # demo data
    "tls_data": "10.5281/zenodo.14620362",
    "janesick": "10.1038/s41467-023-43458-x",
    "tabula_sapiens": "10.1126/science.abl4896",
}


def get(url: str) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "st-report citation check (mailto:support@pennsieve.io)"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as fh:
                return json.load(fh)
        except Exception:  # noqa: BLE001
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"cannot fetch {url}")


def plain(text: str) -> str:
    text = html.unescape(re.sub(r"<[^>]+>", "", text or ""))
    return re.sub(r"\s+", " ", text).strip()


def norm(text: str) -> str:
    text = unicodedata.normalize("NFKD", plain(text)).encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9]", "", text)


def tidy_name(name: str) -> str:
    """Crossref keeps some names as printed: 'MORAN', 'The Tabula Sapiens Consortium*'."""
    name = plain(name).rstrip("*").strip()
    return name.title() if name.isupper() and len(name) > 3 else name


def tidy_title(title: str) -> str:
    """An all-caps title (as printed in 1950) in sentence case."""
    title = plain(title)
    return title[0] + title[1:].lower() if title.isupper() else title


def crossref(doi: str) -> dict:
    m = get("https://api.crossref.org/works/" + urllib.parse.quote(doi))["message"]
    authors = []
    for a in m.get("author", []):
        if a.get("family"):
            authors.append({"family": tidy_name(a["family"]), "given": tidy_name(a.get("given", ""))})
        elif a.get("name"):
            authors.append({"name": tidy_name(a["name"])})
    year = None
    for k in ("published-print", "issued", "published-online", "published"):
        dp = (m.get(k) or {}).get("date-parts")
        if dp and dp[0] and dp[0][0]:
            year = dp[0][0]
            break
    pages = m.get("page") or m.get("article-number") or ""
    journal = (m.get("container-title") or [""])[0]
    short = (m.get("short-container-title") or [journal])[0] if m.get("short-container-title") else journal
    return {"doi": doi, "source": "crossref", "title": tidy_title((m.get("title") or [""])[0]), "authors": authors,
            "journal": plain(journal), "journal_short": plain(short), "year": year,
            "volume": m.get("volume", ""), "issue": m.get("issue", ""), "pages": pages, "type": "article"}


def datacite(doi: str) -> dict:
    a = get("https://api.datacite.org/dois/" + urllib.parse.quote(doi))["data"]["attributes"]
    authors = []
    for c in a.get("creators", []):
        if c.get("familyName"):
            authors.append({"family": c["familyName"], "given": c.get("givenName", "")})
        else:
            authors.append({"name": c.get("name")})
    return {"doi": doi, "source": "datacite", "title": plain(a["titles"][0]["title"]), "authors": authors,
            "journal": a.get("publisher") if isinstance(a.get("publisher"), str) else (a.get("publisher") or {}).get("name", ""),
            "journal_short": "", "year": a.get("publicationYear"), "volume": "", "issue": "", "pages": "",
            "type": "data" if (a.get("types") or {}).get("resourceTypeGeneral") == "Dataset" else "article"}


def pubmed(ref: dict) -> str:
    q = urllib.parse.quote(f"{ref['doi']}[doi]")
    ids = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json")["esearchresult"]["idlist"]
    time.sleep(0.4)
    if len(ids) != 1:
        return ""
    s = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={ids[0]}&retmode=json")["result"][ids[0]]
    time.sleep(0.4)
    title_ok = norm(s.get("title", ""))[:60] == norm(ref["title"])[:60]
    year_ok = str(ref["year"]) in s.get("pubdate", "") or str(ref["year"]) in s.get("epubdate", "")
    first = ref["authors"][0]
    first_name = first.get("family") or first.get("name") or ""
    pm_first = (s.get("authors") or [{}])[0].get("name", "")
    first_name = re.sub(r"^The ", "", first_name)
    author_ok = norm(first_name)[:6] == norm(pm_first.split(" ")[0])[:6] or norm(first_name)[:8] in norm(pm_first)
    if not (title_ok and year_ok and author_ok):
        print(f"  PubMed {ids[0]} does not match {ref['doi']}: title {title_ok}, year {year_ok}, first author {author_ok} "
              f"({pm_first!r} vs {first_name!r}; {s.get('pubdate')})", file=sys.stderr)
        return ""
    return ids[0]


def short_authors(authors: list[dict]) -> str:
    def one(a):
        if a.get("name"):
            return a["name"]
        initials = "".join(p[0] for p in re.split(r"[\s\-\.]+", a.get("given", "")) if p)
        return f"{a['family']} {initials}".strip()
    names = [one(a) for a in authors]
    return ", ".join(names) if len(names) <= 6 else ", ".join(names[:3]) + ", et al."


def verify() -> dict:
    refs = {}
    for key, doi in DOIS.items():
        ref = datacite(doi) if doi.startswith("10.5281/") or doi.startswith("10.48550/") else crossref(doi)
        ref["pmid"] = pubmed(ref) if ref["source"] == "crossref" else ""
        ref["short_authors"] = short_authors(ref["authors"])
        ref["checked"] = time.strftime("%Y-%m-%d")
        refs[key] = ref
        print(f"{key:18} {ref['source']:8} {ref['year']} {ref['short_authors'][:40]:40} PMID {ref['pmid'] or '-':9} {ref['title'][:60]}")
    return refs


def cff_refs(refs: dict) -> str:
    def q(s):
        return json.dumps(str(s), ensure_ascii=False)
    lines = ["references:"]
    for key, r in refs.items():
        lines.append(f"  - type: {'data' if r['type'] == 'data' else 'article'}")
        lines.append(f"    title: {q(r['title'])}")
        lines.append("    authors:")
        for a in r["authors"]:
            if a.get("name"):
                lines.append(f"      - name: {q(a['name'])}")
            else:
                lines.append(f"      - family-names: {q(a['family'])}")
                if a.get("given"):
                    lines.append(f"        given-names: {q(a['given'])}")
        if r.get("journal"):
            lines.append(f"    {'publisher' if r['type'] == 'data' else 'journal'}: " +
                         (f"\n      name: {q(r['journal'])}" if r['type'] == 'data' else q(r['journal'])))
        lines.append(f"    year: {r['year']}")
        if r.get("volume"):
            lines.append(f"    volume: {q(r['volume'])}")
        if r.get("issue"):
            lines.append(f"    issue: {q(r['issue'])}")
        if r.get("pages"):
            p = str(r["pages"]).split("-")
            lines.append(f"    start: {q(p[0])}")
            if len(p) > 1:
                lines.append(f"    end: {q(p[1])}")
        lines.append(f"    doi: {r['doi']}")
        if r.get("pmid"):
            lines.append(f"    notes: {q('PMID ' + r['pmid'])}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    refs = verify()
    OUT.write_text(json.dumps(refs, indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {OUT}")
    if "--cff" in sys.argv:
        cff = ROOT / "CITATION.cff"
        head = cff.read_text().split("references:")[0]
        cff.write_text(head + cff_refs(refs))
        print(f"wrote {cff}")
