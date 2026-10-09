# Spatial and single-cell analysis of mouse lung after AD5-nCoV inhalation

This repository contains R and Python scripts for analysing mouse lung after AD5-nCoV inhalation using 10x Genomics Visium HD spatial transcriptomics and droplet-based single-cell RNA sequencing. The analyses connect cell-type annotation with cellular neighbourhoods (CNs), spatial relationships, ligand–receptor signals, and plasma-cell and T-cell microenvironments.

## Data availability

Processed spatial single-cell data in AnnData (`.h5ad`) format and an interactive visualisation are publicly available through CNGB STOmics:

**[STDS0000420 — dataset, visualisation and downloads](https://db.cngb.org/stomics/datasets/STDS0000420/summary)**

Download the processed objects from this page before using the spatial analysis scripts. Large data files are hosted separately from this code repository.

The spatial sample identifiers used in the scripts are `d001`, `d011`, `d101`, `d102`, `d211`, `d212`, `d451`, and `d452`. The public export uses the following annotation fields:

| Public AnnData field | Description | Earlier script field |
| --- | --- | --- |
| `min_type` | Detailed cell-type annotation | `type` |
| `major_type` | Broad cell-type group | `type2` |
| `CN` | Final cellular neighbourhood annotation | `nich5` in the export procedure |

Several earlier analysis scripts use `nich4` for CN labels. Check the label definitions before adapting those scripts to the public `CN` field; column renaming alone does not establish that an earlier intermediate annotation is identical to the final annotation.

To inspect a downloaded object locally:

```python
from pathlib import Path
import anndata as ad

input_file = Path("downloaded_data/d211.h5ad")
adata = ad.read_h5ad(input_file)
print(adata.shape)
print(adata.obs.columns.tolist())
print(list(adata.obsm.keys()))
```

## Repository structure

```text
scripts/          Core single-cell, spatial and experimental analyses
space_python/     Neighbourhood discovery using HNSW and Leiden
figures/          Plotting scripts and figure-specific analyses
supplementary/    Additional B-cell and niche analyses
config/           Cell-type order and plotting palettes
FILE_MANIFEST.csv Original-to-archived filename mapping and checksums
```

Script filenames describe their analysis purpose. Numeric prefixes within folders organise the files; they do not define an automated execution order across folders.

## Analysis workflow

| Analysis | Relevant scripts |
| --- | --- |
| Single-cell preprocessing, ambient RNA and doublet assessment, sample tags and targeted S-transcript counts | `scripts/01_*` |
| Single-cell merging, quality filtering, clustering and annotation | `scripts/02_*`, `scripts/03_*` |
| RCTD cell-type mapping and spatial label assignment | `scripts/04_*`, `scripts/05_*` |
| Multi-sample neighbourhood discovery | `space_python/01_multisample_cn_discovery_hnsw_leiden.py` |
| Neighbourhood composition, label standardisation, BALCs refinement and spatial merging | `scripts/06_*`–`scripts/09_*` |
| Distance to bronchovascular structures | `scripts/10_*` |
| Spatial adjacency and label permutation tests | `scripts/11_*` |
| Spatial ligand–receptor scoring using the NicheNet mouse interaction network | `scripts/12_*` |
| Plasma-cell niches, immunoglobulin phenotypes and spatial expression gradients | `scripts/13_*`–`scripts/15_*` |
| T-cell ring-based gene gradients and Hallmark signals | `scripts/16_*`, `scripts/17_*` |
| ELISA and ELISpot statistics and plots | `scripts/18_*` |

Figure scripts cover CN composition and temporal changes, observed-to-expected enrichment (Ro/E), adjacency summaries, ligand–receptor plots, spatial signalling flows, receptor expression and selected tissue regions.

### Cellular neighbourhood discovery

The main discovery script is [`space_python/01_multisample_cn_discovery_hnsw_leiden.py`](space_python/01_multisample_cn_discovery_hnsw_leiden.py). It:

1. Counts neighbouring cell types within a 60-pixel radius, excluding the focal cell.
2. Combines the neighbourhood count matrices across samples.
3. Builds a cosine-distance HNSW index and a symmetric KNN graph with `k=30`.
4. Runs Leiden clustering with resolution `2`, `15` iterations and seed `123`.
5. Merges candidate clusters with cell-type composition cosine similarity of at least `0.9` and assigns size-ranked `niche1`, `niche2`, etc. labels.

The counts used for the HNSW index and the cluster-level proportions used for merging are distinct quantities. The merged numerical labels require biological and anatomical interpretation before assignment to the final 11 CNs. Readers using the public processed objects can start with the supplied `CN` annotation.

`space_python/02_*` preserves a separate single-sample implementation with different parameters (`k=20`, L2 distance and Leiden resolution `3`). `space_python/03_*` operates on existing merged labels, and `space_python/04_*` evaluates their spatial enrichment.

### Spatial interaction analyses

`scripts/11_spatial_adjacency_permutation.py` defines neighbours within 120 pixels and performs 400 label permutations with seed `123`. Its current entry point uses CN labels (`nich4`); analyses of detailed cell types require selecting the appropriate label column.

`scripts/12_spatial_lr_scoring_nichenet.R` uses a NicheNet mouse ligand–receptor network, smooths ligand expression over spatial neighbours, and scores interactions using the square root of smoothed ligand expression multiplied by receptor expression. The archived scoring script uses a 100-pixel radius and at most 100 neighbours. Some downstream signalling-flow figure scripts use a 50-pixel radius for assigning signal sources.

Distances are expressed in the coordinate units used by each input object. Verify pixel-to-physical-distance conversion before interpreting a radius as micrometres.

## Software and inputs

Dependencies vary by analysis. The main packages include:

- **R:** Seurat, spacexr (RCTD), decontX, scDblFinder, tidyverse, Matrix, RANN, future.apply and msigdbr, plus the plotting and spatial packages imported by individual scripts.
- **Python:** AnnData, Scanpy, NumPy, pandas, SciPy, hnswlib, igraph, leidenalg, scikit-learn, matplotlib and seaborn. Selected spatial plots use Trackcell; permutation analyses also use joblib and statsmodels.

Install the packages required by the selected script and consult its imports for additional dependencies. A version-pinned environment is not supplied in this archive.

RCTD and some signalling analyses additionally require an appropriate single-cell reference. The ligand–receptor scoring script expects `nichenet/lr_network_mouse_21122021.rds`. Hallmark analyses require the relevant gene sets, and ELISA/ELISpot scripts require their input tables. These additional inputs are separate from the public spatial h5ad objects.

## Using the archived scripts

Start by choosing the analysis of interest and reviewing its input objects, paths, labels and parameters. The scripts preserve the research workflow and include interactive analysis sections. Some R scripts expect objects such as `spa` (a list of spatial Seurat objects), `spam` (a merged spatial object) or `pbmc` (a single-cell reference) to be available in the session. A downloaded h5ad file is not a substitute for those R objects without conversion and preparation.

Update local paths before running a script. References to the original `配色.R` and `original_order2.txt` correspond to `config/plot_palettes.R` and `config/celltype_order.txt` in this archive. The CN discovery script selects the first eight files from a glob; replace this with an explicit sample list when reproducing the analysis so that sample identity and order are preserved.

Published h5ad objects allow readers to inspect the processed expression data and final annotations without repeating cell reconstruction. Cell reconstruction was performed externally, and its code is not included here. Plots of cell boundaries additionally require the geometry information expected by their plotting scripts.

Some sections require intermediate scores or summary tables. For example, ligand–receptor figures need the relevant LR assay, and gradient plots need previously generated gradient summaries. Check these inputs before running figure scripts. This repository is a collection of analysis scripts, rather than a single command that regenerates every result from raw sequencing data.

## Questions

For questions about code, required inputs or label conventions, please open an [issue](https://github.com/danny273/10XHD-analysis-of-mouse-lung-after-AD5_nCoV/issues). Include the script name, package versions and a concise description of the error or result.

