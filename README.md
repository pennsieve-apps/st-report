# st-report

The last step of the three spatial transcriptomics workflows, and the only one whose output is
written for a person:

- **A. Spatial Domain Mapping & Condition Comparison** (Visium, Visium HD, Slide-seq):
  `st-prep → st-qc → st-domains → {st-svg, st-compare} → st-report`
- **B. Cell-type Deconvolution & Niches** (Visium + a scRNA-seq reference):
  `st-prep → st-qc → cell2location → st-niches → st-report`
- **C. Single-cell Spatial Profiling** (Xenium, MERSCOPE, CosMx):
  `st-prep → st-qc → st-celltype → {st-niches, st-ccc} → st-report`

It reads whatever upstream outputs it is given and writes one HTML page and the tables behind
it. **Every input is optional**: the report detects what it has and renders the matching
sections, so the same app closes all three workflows (and any mix of them). The page works
offline: no scripts, no external resources, figures in a folder next to it.

**For a biologist.** One place to read what the workflow found: how many spots or cells each
sample had and how many QC removed (and why), which tissue domains or cell types were found and
how much of each section they cover, which genes have a spatial pattern, what differs between
conditions, which cell types sit together and in which niches, and which ligand–receptor pairs
are both expressed and co-located. Each finding is a plain sentence with its numbers ("D04 covers
24% of spots in kidney vs 7.5% in lung; adj. p 0.66"), followed by the caveats that come with
the method: too few sections per group to call a composition change, marker-based labels being
heuristic, Moran's I significant with tiny effects, a 313-gene panel covering 23 ligand–receptor
pairs. It describes the data; it does not give advice.

**For a computational analyst.** The report fits nothing. It reads the summary tables,
`summary.json` and `provenance.json` each upstream app writes (never per-observation files such
as `labels.csv`, `abundance.csv` or `lr_cluster.csv`), recomputes a few cross-step quantities
(domains present in one condition only; spatially variable genes that are also domain or cluster
markers; genes that change in the same direction in most labels; interactions ranked high by
expression but spatially anti-correlated; cell types with a sizeable share that dominate no
spot), re-thresholds significance at its own `FDR` from the upstream adjusted p-values, copies
the upstream figures and adds a few of its own. The Methods section lists every step's version,
tool versions, parameters and method text from its `provenance.json`, and the Citations section
the deduplicated, verified papers of every step present.

## Sections

| Section | Shown when | What it says |
|---|---|---|
| Study design | `metadata` | samples, platform, observations and genes, conditions (st-compare's `CONDITION_COLUMN`, else `condition`) and samples per level; sample-sheet columns that group the samples like the condition (identical, or donor/slide/batch nested in it); annotation coverage, control features dropped, subsampling, st-prep warnings |
| Quality control | `qc` | observations before/after per sample, removals per reason, thresholds, median counts and genes; samples that lost more than the others, depth differences between samples, panels without mitochondrial genes; a per-sample overview figure and st-qc's histograms and tissue maps |
| Spatial domains | `domains` | domains with size, sample coverage (≥ 20 observations), conditions and top markers; a domains × samples heat map; domains confined to one sample or one condition; agreement with the annotation (ARI/NMI per sample and pooled, best-matching domain per class); st-domains' maps, UMAP, dot plot, confusion |
| Spatially variable genes | `svg` | per sample genes tested and significant, median Moran's I; consensus genes with the domains/clusters they mark; share of significant results with I < 0.05; st-svg's figures |
| Condition comparison | `compare` (one block per `<LABEL>/`) | composition tests (mean share per level, CLR difference, adjusted p) and, when nothing is significant, why (samples per group, the smallest p the rank test can reach); pseudobulk DE per label (tested/skipped with reasons, counts up/down, top genes), genes shared by most labels; a DE-count figure, st-compare's composition plot, heat map and volcanoes |
| Cell types | `deconvolution` and/or `celltype` (+ `celltype-figures`) | cell2location: estimated shares, dominant types and Mixed spots, types hidden by `dominant_type`, co-location groups, reference and reconstruction check; st-celltype: clusters with call, score, runner-up and margin (close calls flagged), clusters sharing a call, composition per sample, marker sets not scorable on the panel |
| Spatial organisation | `niches` (one block per `<LABEL>/`) | homotypic and mixed neighbourhood enrichment, segregated pairs, niches and their composition, niche shares per sample, sparse-graph warning; st-niches' heat maps, niche maps, co-occurrence curves |
| Cell–cell communication | `ccc` (one block per `<LABEL>/`) | consensus top interactions with their spatial R, spatially supported pairs, anti-correlated top pairs, panel coverage of the resource; a spatial-R figure and st-ccc's dot plots, source–target heat maps and local maps |
| Methods | any `provenance.json` | per app: version, tools, parameters, method as run, citations |
| Citations | always | deduplicated over the apps present plus the report's own; DOI and PMID |
| Input data | `metadata` | source paths; a data credit when a demo dataset is recognised (TLS Visium, 10x lymph node, 10x Xenium breast) |

A section whose inputs are present but unreadable is reported at the top of the page and in
`provenance.json` (`sectionErrors`); the rest is still rendered.

## Ports

Inputs (all `placement: subfolder`, none required). Linked from a workflow, each arrives in
`INPUT_DIR/<port>/`; without links (everything side by side, nested, symlinks followed) the files
are found by name anywhere under `INPUT_DIR`, and an upstream folder is recognised by a file only
its app writes.

| Input | Files read | From | Recognised by |
|---|---|---|---|
| `metadata` | `samples.csv`, `prep.json`, `provenance.json` | st-prep `metadata` | `samples.csv` with `sample_id` and `platform`/`n_obs`; `prep.json` with `samples` |
| `qc` | `qc_summary.csv`, `thresholds.csv`, `<sample_id>_qc.png`, `<sample_id>_spatial_counts.png`, `provenance.json` | st-qc `qc` | `qc_summary.csv` |
| `domains` | `domain_summary.csv`, `domain_markers.csv`, `agreement.csv`, `annotation_match.csv`, `summary.json`, `figures/*.png`, `provenance.json` | st-domains `labels` | `domain_summary.csv`, `domain_markers.csv` |
| `svg` | `svg.csv`, `svg_consensus.csv`, `summary.json`, `*.png`, `provenance.json` | st-svg | `svg_consensus.csv`, `svg.csv` |
| `compare` | per `<LABEL>/`: `composition_tests.csv`, `de_summary.csv`, `de.csv`, `summary.json`, `composition.png`, `de_heatmap.png`, `volcano_*.png`, `provenance.json` | st-compare (several links side by side) | `composition_tests.csv`, `de_summary.csv` |
| `deconvolution` | `type_summary.csv`, `colocation.csv`, `summary.json`, `figures/*.png`, `provenance.json` | cell2location `labels` | `type_summary.csv`, `colocation.csv`, `abundance.csv` |
| `celltype` | `cluster_annotation.csv`, `cell_type_summary.csv`, `cluster_markers.csv`, `summary.json`, `provenance.json` | st-celltype `labels` | `cell_type_summary.csv`, `cluster_annotation.csv` |
| `celltype-figures` | `umap_cell_type.png`, `umap_cluster.png`, `cell_type_by_sample.png`, `dotplot_marker_sets.png`, `<sample_id>_cell_type.png` | st-celltype `figures` | those PNGs |
| `niches` | per `<LABEL>/`: `nhood_enrichment.csv`, `niche_composition.csv`, `niche_summary.csv`, `summary.json`, `*.png`, `provenance.json` | st-niches (several links side by side) | `niche_composition.csv`, `nhood_enrichment.csv` |
| `ccc` | per `<LABEL>/`: `lr_consensus.csv`, `lr_spatial.csv`, `summary.json`, `*.png`, `provenance.json` | st-ccc (several links side by side) | `lr_consensus.csv`, `lr_spatial.csv` |

A `<LABEL>` is read from the run's `summary.json` (`label`), else its `LABEL` parameter, else the
folder name; two runs with the same label are both kept (the second named after its folder).

Outputs (two, so each in `OUTPUT_DIR/<name>/`):

| Output | Files |
|---|---|
| `report` | `index.html`, `figures/*.png` (the report's own: `qc_overview.png`, `domain_coverage.png`, `de_counts_<LABEL>.png`, `ccc_spatial_<LABEL>.png`, a composition bar chart when the upstream one is missing; and copies of upstream figures, prefixed by section), `provenance.json` |
| `tables` | the CSVs below (each only when its section has inputs), `summary.json`, `provenance.json` |

| Table | Columns |
|---|---|
| `samples.csv` | `samples.csv` of st-prep as read |
| `qc.csv` | st-qc's `qc_summary.csv` + `condition` |
| `domains.csv` | `domain, n_obs, fraction, n_samples` (≥ 20 observations), `samples, conditions, markers` (`;`-joined) |
| `domain_coverage.csv` | `sample_id, domain, n_obs, fraction` |
| `domain_agreement.csv` | `sample_id, ari, nmi, n_obs, ...` per sample and `all` |
| `svg_samples.csv` | `sample_id, n_obs, n_tested, n_significant, fraction_significant, median_morans_i, max_morans_i, top_genes, graph_mean_degree, condition` |
| `svg_consensus.csv` | top `MAX_ITEMS`: `gene, n_samples, n_significant, mean_morans_i, median_rank, rank, top10_marker_of` |
| `composition_tests.csv` | `label_set, label_column, label, comparison, clr_difference, p, p_adj, rank_p_adj, tested, significant` (at `FDR`), `mean_fraction_<level>` |
| `de_summary.csv` | `label_set` + st-compare's `de_summary.csv` (counts recomputed at `FDR` when it differs from st-compare's) |
| `de_top_genes.csv` | top `MAX_ITEMS` per `label_set` by adjusted p: `label, comparison, gene, log2_fold_change, lfc_se, base_mean, p_adj` |
| `cell_types.csv` | `sample_id, cell_type, fraction` (st-celltype) |
| `clusters.csv` | `cluster, cell_type, n_obs, score, second_best, second_score, margin, close, markers` |
| `deconvolution_composition.csv`, `deconvolution_types.csv` | `sample_id, cell_type, fraction_of_cells`; `cell_type, fraction, min, max, n_dominant` |
| `niches.csv` | `label_set, niche, n_obs, fraction, top_labels` |
| `neighbourhood_enrichment.csv` | `label_set` + st-niches' pooled `nhood_enrichment.csv` rows (each pair once) |
| `ccc_top.csv` | top `MAX_ITEMS` consensus interactions: `label_set, rank, source, target, ligand, receptor, n_samples, n_samples_top, mean_magnitude_rank, mean_specificity_rank, spatial_r` |
| `ccc_spatial.csv` | per ligand–receptor pair: `label_set, ligand, receptor, mean_global_morans, n_samples, n_samples_significant` (R > 0, p ≤ 0.05), `n_obs_significant` |
| `methods.csv`, `citations.csv` | per app run: `app, step, version, tools, parameters, citations`; per reference: `number, authors, title, journal, year, volume_pages, doi, pmid` |
| `summary.json` | every headline number, observation and caveat of the page, per section; workflows detected; methods; citations; data credit |

The step fails (exit 1) only when none of the inputs is found; exit 2 on a bad parameter.

## Parameters

| Name | Default | What it does |
|---|---|---|
| `REPORT_TITLE` | `Spatial Transcriptomics Report` | Title of the page. |
| `FDR` | `0.05` | Threshold the report applies to the upstream adjusted p-values when it counts and lists significant composition changes, DE genes and spatially variable genes. Display only; the upstream tables are unchanged. The page says when it differs from the threshold the step used. |
| `MAX_ITEMS` | `20` | Rows per top-N table (consensus genes, DE genes, interactions, enrichment pairs). |
| `MAX_FIGURES` | `12` | Figures per section (per `<LABEL>` for comparison, niches and communication), the report's own included; the page says how many more the upstream output holds. |

## Results on the demo data

Upstream outputs of the workflow repo's `data/work/<demo>/`, assembled as `run-local.sh` does
(each upstream output port copied into `INPUT_DIR/<port>/`); reports in
`data/work/<demo>/report/{report,tables}/`.

**TLS Visium, workflow A** (8 FFPE sections, kidney ×3 and lung ×5; `compare/domain` and
`compare/annotation`): QC removed 218 of 27,715 spots (0.79%; LC5 3.3%); median counts per spot
differ 13-fold between sections (480 in LC1 to 6,032 in LC2), and the probe set has no
mitochondrial genes. BANKSY's 14 domains: 9 present in ≥ 2 sections, 6 in all 8; D10 and D11
only in LC5, D12 only in KC2, D08 only in lung (LC2, LC4). Pooled ARI 0.18 / NMI 0.29 with the
pathologist annotation; TLS → D03 (precision 0.41, recall 0.76). 4,512 of 9,260 genes spatially
variable in at least one section, 6 in all 8, led by immunoglobulins (IGHG1, IGHG2, IGHA1, IGHM),
which are also D03/D07 markers; 36% of significant results have I < 0.05. Composition: no domain
or annotation class differs between lung and kidney (closest D04, 24% vs 7.5% of spots, adj.
p 0.66; with 3 vs 5 sections the rank test cannot go below p = 0.036 before correction over 14
labels). Pseudobulk DE: 3,626 significant gene × domain results in 8 tested domains (D04 most,
1,174), 6 domains not testable (absent from kidney or too small), and 191 genes change the same
way in ≥ 5 of the 8 domains (SFTPB, FASN up in lung; CRYAB, FXYD2 down): an organ-wide shift. The
page says that `tissue` and `condition` group the sections identically and that kidney vs lung
also differs in tumour type and patient.

**Xenium breast cancer, workflow C** (2 replicates, 313-gene panel, full slides; plus st-svg):
QC removed 9,057 of 286,532 cells (3.2%). st-celltype: 10 clusters named as 8 cell types, three
epithelial clusters (C01, C05, C07), one close call (C05 Epithelial 1.04 vs Luminal
hormone-sensing 0.88), 16 of 37 marker sets not scorable on the panel. Niches: same-type
contacts strongest (Epithelial z 630), B cell–T cell the strongest mixed pair (141), epithelium
segregated from T cells (−361); 3 niches by silhouette (0.42; k = 5 0.41). Communication: top
Fibroblast → B/T cell CXCL12 → CXCR4 (spatial R 0.22 in both replicates); S100A4 → ERBB2 ranks
third by expression but is spatially anti-correlated (R −0.32); only 23 ligand–receptor pairs of
the resource have all their genes on the panel.

**Lymph node, workflow B**: the cell2location and st-niches outputs of the lymph-node demo were not
yet available; workflow B is covered by the synthetic fixture.

## Runtime

`cpu: 1024`, `memory: 2048`, `computeTypes: [standard, lambda]` (the image serves both; see
`processor/entry.py` and `processor/handler.py`), timeout 10 min. Measured with
`docker run --memory 3g --cpus 1 --network none` under **amd64 emulation on an M-series Mac**,
peak from the container cgroup's `memory.peak`:

| Inputs | Sections | Peak memory | Wall time (incl. container start) |
|---|---|---|---|
| TLS Visium, workflow A (8 samples, 2 comparisons, 53 figures) | 5 | 198 MB | 8.5 s |
| Xenium breast, workflow C + st-svg (2 full slides, 277 k cells) | 6 | 152 MB | 4.2 s |

The report reads summary tables only, so memory does not grow with the number of spots or cells
(the largest files it reads are st-compare's `de.csv` and st-svg's `svg.csv`, genes × labels or
genes × samples). 2048 MB is the smallest Fargate size for 1 vCPU; as a Lambda function 512 MB
would do. Image: 541 MB on disk (`python:3.12.11-slim-bookworm` + pandas, matplotlib, jinja2,
awslambdaric); the report folder is 19–37 MB on the demos, almost all copied figures.

## Running it

```sh
python3.12 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
make test PYTHON=.venv/bin/python           # synthetic fixture: workflows A, B, C, combinations, each port alone, none
make run PYTHON=.venv/bin/python            # Docker over data/fixture/linked/ -> data/output/
INPUT_HOST=/path/to/inputs make run         # Docker over your own folder (one subfolder per port)
```

`tests/make_fixture.py` writes contract-following upstream outputs for every port (planted: a
domain present in one sample only, a significant composition change, a sample losing more to QC,
an 8-fold depth difference, close cell-type calls, an anti-correlated top interaction, a sparse
neighbour graph, a 23-pair panel). `tests/snapshots/index.html` is the page for workflow A's
fixture; after an intended change to the page, `UPDATE_SNAPSHOTS=1 make test` rewrites it.
`make validate` checks `app.yml` against the manifest schema and `CITATION.cff` with cffconvert.

## Limitations

- **The report is only as good as the upstream outputs.** It re-thresholds and cross-references,
  but takes every statistic (domains, Moran's I, DESeq2 results, abundances, enrichment z-scores,
  LIANA+ ranks) as the steps wrote them; a caveat on the page is the app's own limitation
  restated, not a new test.
- **The condition** is st-compare's `CONDITION_COLUMN` when a comparison is linked, else a
  `condition` column of `samples.csv`; confounding is detected only from other `samples.csv`
  columns (identical groupings, or donor/patient/slide/batch-like columns nested in the
  condition).
- **Demo-dataset credit** is recognised from st-prep's recorded source paths; other datasets get
  no credit line (cite your data yourself).
- **Figures are static PNGs** (no hover, no zoom beyond opening the file); every figure has a
  table on the page or a CSV beside it. At most `MAX_FIGURES` per section are copied.
- The page is light-themed only and in English.
- Several `compare`, `niches` or `ccc` links are reported side by side; nothing is compared
  between them.

## Licences

Our code: Apache-2.0 (`LICENSE`). Bundled components and their licences: `THIRD_PARTY.md` (all
permissive; no GPL code in the image).

## Citations

The report cites, per step, the papers in `processor/citations.py` (the methods papers each
upstream app's own `CITATION.cff` lists), deduplicated, with full bibliographic data in
`processor/references.json`. `scripts/verify_citations.py` fetched every DOI from Crossref
(DataCite for the Zenodo dataset) and its PMID from PubMed E-utilities (`esearch` by DOI, then
`esummary`: title, first author and year must match) on 2026-10-30; all 45 are in
`CITATION.cff` with the author lists exactly as Crossref/DataCite give them. The report's own:

- Moses L, Pachter L. Museum of spatial transcriptomics. *Nat Methods* 19(5):534-546 (2022).
  doi:10.1038/s41592-022-01409-2. PMID 35273392.
- Sang-aram C, Browaeys R, Seurinck R, Saeys Y. Spotless, a reproducible pipeline for
  benchmarking cell type deconvolution in spatial transcriptomics. *eLife* 12:RP88431 (2024).
  doi:10.7554/eLife.88431. PMID 38787371. (Cited with the deconvolution caveats.)
- Hunter JD. Matplotlib: A 2D Graphics Environment. *Comput Sci Eng* 9(3):90-95 (2007).
  doi:10.1109/MCSE.2007.55. (Not in PubMed; checked on Crossref.)

Per step (shown when that step's provenance is among the inputs): scanpy (PMID 29409532),
anndata (doi:10.21105/joss.04371), SpatialData (PMID 38509327), Visium/Ståhl 2016 (PMID 27365449)
(st-prep); Heumos 2023 (PMID 37002403), pipeComp (PMID 32873325), Stuart 2019 (PMID 31178118)
(st-qc); BANKSY (PMID 38413725), Harmony (PMID 31740819), Leiden (PMID 30914743), UMAP
(doi:10.21105/joss.00861), SDMBench (PMID 38491270) (st-domains); Moran 1950
(doi:10.1093/biomet/37.1-2.17), Geary 1954 (doi:10.2307/2986645), Benjamini–Hochberg 1995
(doi:10.1111/j.2517-6161.1995.tb02031.x), SpatialDE (PMID 29553579), squidpy (PMID 35102346)
(st-svg); PyDESeq2 (PMID 37669147), DESeq2 (PMID 25516281), Squair 2021 (PMID 34584091), scCODA
(PMID 34824236), Aitchison 1982 (doi:10.1111/j.2517-6161.1982.tb01195.x), SciPy (PMID 32015543)
(st-compare); cell2location (PMID 35027729), scvi-tools (PMID 35132262) (cell2location);
Satija 2015 (PMID 25867923), PanglaoDB (PMID 30951143), OmniPath (PMID 33749993) (st-celltype);
Schürch 2020 (PMID 32763154), CellCharter (PMID 38066188), Lloyd 1982
(doi:10.1109/TIT.1982.1056489), Rousseeuw 1987 (doi:10.1016/0377-0427(87)90125-7) (st-niches);
LIANA+ (PMID 39223377), LIANA (PMID 35680885), RRA (PMID 22247279), CellPhoneDB (PMID 32103204),
Connectome (PMID 35264704), NATMI (PMID 33024107), SingleCellSignalR (PMID 32196115) (st-ccc).

Demo data: Dawo S, Nonchev K, Silina K. 10x Visium Spatial Transcriptomics Dataset: Kidney (3)
and Lung (5) Cancer with Tertiary Lymphoid Structures. Zenodo (2025).
doi:10.5281/zenodo.14620362 (CC BY 4.0; checked on DataCite). 10x Genomics Human Lymph Node
(Visium) and Xenium FFPE Human Breast Cancer replicates 1 and 2 (CC BY 4.0, dataset pages; the
Xenium data are described in Janesick A, et al. *Nat Commun* 14:8353 (2023),
doi:10.1038/s41467-023-43458-x, PMID 38114474). Tabula Sapiens Consortium. *Science*
376:eabl4896 (2022), doi:10.1126/science.abl4896, PMID 35549404 (the lymph-node scRNA-seq
reference).
