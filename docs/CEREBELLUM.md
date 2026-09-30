# Developing Cerebellum snRNA-seq integration

The native `snRNA-seq Developing Cerebellum` tab appears between ProteinAtlas.org Profile and TF Targets. It follows the existing Run Analysis gene event. Only Human and Mouse are included, as requested. There is no embedded external app or runtime dependency on its web service.

## Implementation files

Created for this integration: `R/cerebellum.R`, `shiny/cerebellum_interface.R`, `scripts/download-cerebellum.py`, `scripts/download-cerebellum-orthology.py`, `scripts/install-cerebellum-dependencies.R`, `scripts/prepare-cerebellum.R`, `scripts/inspect-cerebellum.R`, `scripts/check-cerebellum.R`, `tests/testthat/test-cerebellum.R`, and this document.

Updated: `R/load_core.R`, `shiny/app.R`, `shiny/interface.R`, `shiny/dashboard.css`, `DESCRIPTION`, `renv.lock`, `scripts/install-dependencies.R`, `scripts/build-deployment-bundle.R`, `scripts/check-depmap-browser.cjs`, `.Rbuildignore`, and `README.md`. These changes extend the existing working tree without replacing earlier app development.

## Sources and attribution

- Sepp M., Leiss K. et al., *Cellular development and evolution of the mammalian cerebellum*, Nature, published online 29 November 2023: https://doi.org/10.1038/s41586-023-06884-x
- Official explorer and released SCE objects: https://apps.kaessmannlab.org/sc-cerebellum-transcriptome/
- Source inspected: https://gitlab.com/kaessmannlab/shiny-mammalian-cerebellum, commit `ed949f66eadb5fdd4ed016f4352df541f10bc3df`.
- Source README declares CC BY 4.0. Retain author/source links, the licence link, and identification of integration changes when redistributing. The UI and exported plots carry attribution. This integration is independently implemented using the documented structure and published data.

The official `final_processed_data/` files used are `hum_sce_final.rds`, `hum_exonic_sce_final.rds`, `mou_sce_final.rds`, and `mou_exonic_sce_final.rds`. Their server Last-Modified dates are 28 October 2023. Local download manifests record full URLs, byte lengths, dates, ETags, and SHA-256 checksums. These are downloaded-byte checksums, not checksums independently signed by the publisher.

## Data integrity and display

Coordinates come directly from the pre-mRNA objects' `reducedDim(sce, "umap2d")`; these are never recomputed. The annotation fields verified in the source app are:

| Control | Source field |
|---|---|
| Broad lineage | `broad_lineage` |
| Cell type | `cell_type` |
| Cell subtype | `subtype` |
| Cell state | `dev_state` |

Annotations are retained from the published embedding cells. Missing annotation is displayed as Unclassified. Identical labels receive identical colours across species. Legends use 13-pixel text and constant-size colour symbols independent of the small plotted points.

Expression comes exclusively from `assay(exonic_sce, "umi")`. Unique exonic cell IDs are matched explicitly to unique embedding cell IDs, validated before indexing. Missing and duplicate IDs fail preparation. Every exonic cell must have a published embedding. The original embedding includes additional nuclei: their expression is unavailable, not zero.

| Species | Published nuclei | Exonic measured nuclei | Without exonic measurement | Exonic genes |
|---|---:|---:|---:|---:|
| Human | 180,956 | 169,280 | 11,676 | 58,301 |
| Mouse | 115,282 | 107,145 | 8,137 | 54,277 |

The display uses raw counts on log(1 + UMI) colour positions, with count-labelled continuous legends and separate species ranges. It does not normalize species to look alike. Grey marks measured zero, pale grey marks no exonic measurement. No cross-species quantitative equivalence is implied.

Interactive plots use Plotly WebGL with pan, zoom, reset, and compact hover. Overviews/backgrounds show a fixed sample stratified by cell type and subtype (50,095 Human; 50,056 Mouse). Every positive-expression cell is included, preserving rare expression signals. All displayed points retain exact original coordinates. PNG/PDF exports include all nuclei, at 9 x 7 inches and 300 dpi for PNG.

## Human-to-mouse matching

`download-cerebellum-orthology.py` saves the official Ensembl BioMart human/mouse table and registry (Ensembl Genes 116 in this installation), complete XML query, retrieval date and checksum. Only unique `ortholog_one2one` mappings are accepted; duplicated one-to-many or many-to-one mappings are excluded. Matching does not use other species or capitalize a human symbol to guess a mouse gene.

Human symbols are matched exactly, case-insensitively. Ensembl IDs (with optional version suffix) are accepted within this module; the other app analyses retain their existing input support. If a modern symbol maps to multiple IDs, it is resolved only when one exact ID occurs in the published human matrix; otherwise it remains unresolved. This handles alternate-locus SLC16A1 records without guessing. A unique gene absent from the human matrix can still display its mouse orthologue. Old stable IDs absent from the saved mapping may require a future curated mapping update.

## Fresh installation

Run from the repository root, using Rscript and Python 3.8+ on PATH (Python 3.8+ supports the download script's assignment expression):

```powershell
Rscript --vanilla scripts/install-dependencies.R --shiny --tests
Rscript --vanilla scripts/install-cerebellum-dependencies.R
python scripts/download-cerebellum.py
python scripts/download-cerebellum-orthology.py
Rscript --vanilla scripts/prepare-cerebellum.R
Rscript --vanilla scripts/check-cerebellum.R
Rscript --vanilla tests/run-tests.R
Rscript --vanilla scripts/serve-shiny.R 3886
```

On this Windows machine, substitute `& 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'` for `Rscript` if R is not on PATH. Preparation requires SingleCellExperiment and its Bioconductor dependencies, Matrix, digest, jsonlite, data.table and filelock. Plotly is the only new web-runtime dependency (with its crosstalk dependency); SCE is unnecessary at runtime.

Default source directory: `data/cerebellum/source`. Default processed directory: `data/cerebellum/processed`. Both are ignored by Git. The two Python download scripts accept an optional destination argument; the R preparation script accepts source and destination arguments. Preparation refuses to overwrite a completed atlas. For a new version, build a new destination, then set `CEREBELLUM_DATA_DIR` to its absolute path and restart the app. The app must receive the **entire processed directory**, including the manifest, index files, gene map and binary count files.

No complete count matrix is loaded at runtime. Each species has a sparse binary file indexed by gene and byte offset; each record stores mapped cell indices and integer UMI counts. Shared metadata/index data are loaded once per process. Gene workers read only the requested gene records, and annotation changes do not start new workers. Rapid queries cancel prior workers and reject stale results. Missing data failures are contained within the module.

## Storage and memory

- Official downloads: 1.660 GB (1.55 GiB).
- Processed assets: approximately 1.589 GB (1.48 GiB).
- Total with retained downloads: approximately 3.25 GB, plus installed packages and exports.
- Shared runtime atlas/index object: measured at 97.1 MiB in R. Per-query count vectors add a few MiB; Plotly serialization and each active browser add overhead. Budget several hundred MiB beyond the existing app and measure concurrent sessions on the deployment host.
- One-time preparation loads one species at a time and transposes its sparse exonic matrix. Allow approximately 8-12 GB available RAM as a conservative provisioning estimate; peak resident memory was not instrumented. No 220-GB upstream SQLite database is constructed.
- Four tested gene queries took approximately 0.04-0.08 seconds each for local lookup after metadata was loaded, excluding worker startup and browser rendering.

## Checks and limitations

Unit checks cover ID joins with shuffled cells, missing/duplicate IDs, zero versus unavailable, missing human/mouse genes, no 1:1 orthologue, rejected ambiguous orthology, four annotation modes, navigation placement, shared gene events and annotation-only updates. Live checks cover HLX/Hlx, MYC/Myc, SLC16A1/Slc16a1, PAX6/Pax6 and an absent query; results and full-data exports are saved under `outputs/cerebellum-acceptance`.

Preparation compares selected reconstructed gene vectors with original sparse assay rows and records embedding hashes. Data files are checked against the download manifest during preparation. Runtime checks verify schema, IDs, sizes and gene-record bounds; they do not rehash 1.6 GB on each query. The prepared manifest records hashes for administrators to audit after transfer.

The display sample is not a cell-abundance estimate. Hover is intentionally compact. Source annotations remain abbreviated. Large full-cell PDF exports may take longer than PNG. The integrated report generator retains its existing contents; cerebellum figures have dedicated downloads in the new tab.
