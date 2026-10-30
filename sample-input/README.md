# sample-input

st-report consumes the *outputs* of the other spatial transcriptomics steps (st-prep, st-qc,
st-domains, st-svg, st-compare, cell2location, st-celltype, st-niches, st-ccc), so there is no
public file to download for it directly.

For a local run, generate the synthetic fixture (every port of the three workflows, a few
summary tables and tiny PNGs; no per-observation data):

```sh
make fixture                    # data/fixture/{linked,flat}/ with workflows A, B and C together
make fixture WORKFLOW=a         # workflow A only (metadata, qc, domains, svg, compare)
make run                        # the Docker image over data/fixture/linked/
```

For real data, run the upstream steps on a demo dataset (see
`workflow-spatial-transcriptomics/scripts/fetch-sample-data.sh` and `run-local.sh`) and give the
report one folder per input port, each a copy of the linked upstream output port:

| st-report input | Upstream output | Workflow |
|---|---|---|
| `metadata` | st-prep `metadata` | all |
| `qc` | st-qc `qc` | all |
| `domains` | st-domains `labels` (its whole output) | A |
| `svg` | st-svg (its whole output) | A (also usable in C) |
| `compare` | st-compare (whole output: `<LABEL>/` folders; several runs side by side) | A |
| `deconvolution` | cell2location `labels` | B |
| `celltype` | st-celltype `labels` | C |
| `celltype-figures` | st-celltype `figures` | C |
| `niches` | st-niches (whole output: `<LABEL>/`) | B, C |
| `ccc` | st-ccc (whole output: `<LABEL>/`) | C |

With the workflow repo's `data/work/` (e.g. the TLS Visium demo):

```sh
W=../workflow-spatial-transcriptomics/data/work/tls
mkdir -p in && for p in metadata:prep/metadata qc:qc/qc domains:domains svg:svg compare:compare; do
  mkdir -p in/${p%%:*} && cp -R $W/${p#*:}/. in/${p%%:*}/; done
INPUT_HOST=$PWD/in make run
```

Mount real copies, not symlinks to host paths (they do not resolve inside the container).
