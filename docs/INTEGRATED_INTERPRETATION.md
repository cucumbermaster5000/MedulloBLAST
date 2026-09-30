# Integrated Interpretation and selected report figures

The Overview now contains a live Integrated Interpretation, a source-linked evidence ledger, cross-dataset pattern cards, and expandable previews of the selected report figures. Existing tabs and analyses remain available. The HTML/PDF report and ZIP preserve the same interpretation, selected figures, captions, sources and settings.

## Architecture

`R/report_evidence.R` operates on existing structured session results. It performs no remote requests, statistical tests, language-model calls or gene-specific branching. Each observation has a stable ID, gene, source, destination tab, text, evidence type, qualitative strength, source statistics and explicit rule flags. Pattern records reference their supporting observation IDs. Every summary sentence retains those IDs. Exported links target the report evidence ledger and source citations; live links navigate to the corresponding application tab.

The rule version is `1.0.0`. A raw R2 omnibus p-value below 0.05 enables a subgroup-association flag. The largest subgroup median is explicitly descriptive: the omnibus test does not establish pairwise enrichment or subgroup specificity. A DepMap gene effect <= -0.5 is a disclosed descriptive screening threshold, not a p-value. Detection in <=1% of measured human developmental nuclei is a disclosed sparse-detection threshold, not a claim of biological absence. Missing measurements never count as zero or as negative evidence. No numerical confidence score is produced.

Supported cross-dataset patterns:

- R2 subgroup association with sparse human developmental detection.
- R2 subgroup association with dependency-screen signals in measured medulloblastoma models.
- R2 subgroup association without any measured medulloblastoma model meeting the dependency threshold.
- Human/mouse developmental detection falling on different sides of the 1% threshold, with explicit species/stage/composition/depth caveats.

All combined patterns are labelled hypothesis-generating. Direct observations, statistical associations, correlations and regulatory evidence remain separate. Source-reported curated target records are distinguished from binding evidence; neither proves regulation in the analysed tumours. Clinical observations retain their exploratory, unadjusted status. Multiple R2 cohorts retain separate provenance and are not treated as replicates merely because both exist.

## Figure selection and reproducibility

Each major plotted output has an **Include in Gene Report** checkbox. Initial defaults are selected; user choices persist during the session. The figure registry uses stable IDs and existing plot objects. Clinical plots retain the current cohort and cutoff. DepMap retains the selected tumour/all-model view. HPA uses one shared plot function for the tab and report, retaining the chosen assay and units. Cerebellum previews reuse the original interactive plotting function; full-data exports reuse its existing PNG/PDF export function and current annotation.

The registry orders tumour expression, subgroups/subtypes, clinical associations, dependency, HPA, cerebellum and functional/regulatory plots. Missing figures are omitted; excluded figures are absent from the report package. Figure exclusion does not erase supporting scientific results from the interpretation. Captions describe the dataset, analysis, gene and available statistics. Survival legends show **Low (n = x)** and **High (n = x)** using the complete cases actually included in the curves, with ties and cutoff changes handled by the existing analysis.

The live report responds to completed module results and figure/settings changes. Opening it performs no retrieval or expensive analysis. Collapsed previews are rendered on demand by Shiny. Generating an export freezes a snapshot in a background worker. Later changes invalidate the old download controls and require generating a fresh export, avoiding stale selections. Export figures remain 300-dpi PNG plus vector PDF companions. The assembled PDF embeds those high-resolution images. The ZIP contains tables, snapshot, `interpretation.json`, `selected_figures.json`, source metadata and reproduction instructions.

## Citations

`report_citations()` centrally maps the actual available sources to R2/cohort provenance, DepMap release metadata, HPA profile/release, the Sepp/Leiss publication and cerebellum manifest, TF source metadata and PubMed query provenance. It reuses existing metadata and does not invent publications. The interpretation JSON preserves the supporting statistics and source references, including individual TF record references in the evidence tables.

## Files

Created: `R/report_evidence.R`, `tests/testthat/test-report-evidence.R`, `scripts/check-report-interpretation.R`, and this document.

Modified: `R/load_core.R`, `R/integrated_report.R`, `R/profile_statistics.R`, `shiny/report_interface.R`, `shiny/interface.R`, `shiny/hpa_interface.R`, `shiny/depmap_interface.R`, `shiny/biology_interface.R`, `shiny/cerebellum_interface.R`, `report/integrated_gene_report.Rmd`, `tests/testthat/test-cerebellum.R`, and `tests/testthat/test-survival-cutoff.R`.

No new package installation or external download was required. The already-installed `base64enc` package is now explicitly declared in `DESCRIPTION` for self-contained figure embedding. `DESCRIPTION`, `scripts/install-dependencies.R`, `scripts/build-deployment-bundle.R` and `docs/INTEGRATED_GENE_REPORT.md` were also updated.

## Validation

Run `Rscript --vanilla tests/run-tests.R`. Tests cover evidence references, statistics, missing HPA/DepMap, sparse development, contradictory dependency, separate R2 cohorts, unsupported protein/pathway claims, figure ordering/exclusion, captions, citation mapping, frozen evidence export, stale-gene rejection and selection-driven download invalidation. Survival tests verify counts for groups with cutoff ties.

`Rscript --vanilla scripts/check-report-interpretation.R` is a local acceptance check using the existing saved MYC report snapshot and installed cerebellum assets. It builds a new timestamped HTML/PDF/ZIP output, verifies evidence anchors and excludes selected mouse plots. It requires that saved acceptance snapshot; unit tests use independent synthetic fixtures.

Validation completed on 26 September 2026: the full suite passed with one existing opt-in live R2 test skipped and three existing package-serialization warnings. Targeted report tests passed after final formatting changes. The MYC HTML/PDF/ZIP acceptance check passed, including embedded-image verification and excluded figures. Browser checks confirmed evidence navigation, figure exclusion and inclusion, selected survival previews and the Low/High sample-count labels. PDF pages containing the interpretation and selected figures were visually reviewed.

## Limits and next improvements

The available HPA profile is not a matched medulloblastoma protein-abundance cohort, so RNA/protein concordance is deliberately not inferred. Target-set enrichment is not independent tumour pathway enrichment, so no regulatory/pathway convergence is claimed. Supporting those patterns requires genuinely comparable protein measurements or independently calculated tumour pathway results, with tissue/context and assay metadata.

Rules use descriptive cutoffs; sensitivity analysis across these cutoffs and explicit cohort/model subgroup matching are useful next improvements. No cell-of-origin, oncogene, clinical-benefit or causality claims are generated. Sparse evidence produces a shorter summary and the explicit no-pattern message rather than padded speculative prose. Large full-data cerebellum PDFs and snapshots increase export time and file size; this does not reload the full count matrix.

## Simplified export — 27 September 2026

The exported report now uses one biological findings section rather than repeating summary, evidence ledger and key findings. It retains evidence anchors, effect estimates, caveats, selected figures and source citations. Long result tables, complete abstracts and software/version inventories remain in the companion ZIP. The reading list contains up to five papers in the selected order, with retraction flags. Live Overview interpretation and underlying analyses are unchanged.
