# Integrated report export. This file consumes a frozen Shiny-session snapshot;
# it performs no retrieval and no statistical analysis.

report_value <- function(x, fallback = "Unavailable") {
  if (is.null(x) || !length(x) || is.na(x[[1]]) || !nzchar(as.character(x[[1]]))) fallback else as.character(x[[1]])
}

report_safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

report_write_table <- function(x, path) {
  if (is.null(x) || !is.data.frame(x) || !ncol(x)) return(FALSE)
  data.table::fwrite(x, path, na = "NA")
  TRUE
}

report_save_plot <- function(plot, stem, figures_dir, width = 9, height = 6) {
  if (is.null(plot)) return(character())
  pdf <- file.path(figures_dir, paste0(stem, ".pdf"))
  png <- file.path(figures_dir, paste0(stem, ".png"))
  ggplot2::ggsave(pdf, plot, device = grDevices::pdf, width = width, height = height)
  ggplot2::ggsave(png, plot, device = "png", width = width, height = height, dpi = 300)
  c(pdf = pdf, png = png)
}

report_hpa_plot <- function(profile, section = "brain") {
  report_hpa_selected_plot(profile,section)
}

report_publication_table <- function(publications) {
  if (is.null(publications) || !length(publications$articles)) return(data.frame())
  do.call(rbind, lapply(publications$articles, function(a) data.frame(
    rank = a$rank, pmid = a$pmid, title = report_value(a$title, ""),
    authors = paste(a$authors, collapse = "; "), journal = report_value(a$journal, ""),
    publication_date = report_value(a$publication_date, ""), doi = report_value(a$doi, ""),
    publication_types = paste(a$publication_types, collapse = "; "), retracted = isTRUE(a$retracted),
    pubmed_url = report_value(a$pubmed_url, ""), abstract = report_value(a$abstract, ""),
    stringsAsFactors = FALSE)))
}

report_key_findings <- function(snapshot) {
  x <- snapshot$main; findings <- character()
  if (!is.null(x$analysis$statistics)) {
    s <- x$analysis$statistics
    findings <- c(findings, sprintf("Broad MB subgroup expression comparison returned Kruskal-Wallis p = %.3g (epsilon-squared %.3g; n = %d). This is an unadjusted distribution comparison.", s$p_value, s$epsilon_squared, s$n))
  }
  if (!is.null(x$subtypes$overall)) findings <- c(findings, sprintf("The molecular-subtype comparison used %d samples and returned p = %.3g; pairwise results use BH adjustment.", x$subtypes$overall$n, x$subtypes$overall$p_value))
  if (!is.null(x$clinical$metastasis$all) && identical(x$clinical$metastasis$all$status, "complete")) {
    s <- x$clinical$metastasis$all$statistics
    findings <- c(findings, sprintf("In complete MB, metastatic versus M0 expression had raw p = %.3g and rank-biserial effect %.3g; this exploratory result is unadjusted for clinical covariates.", s$p_value, s$rank_biserial_metastatic_vs_m0))
  }
  if (!is.null(x$clinical$survival$all) && identical(x$clinical$survival$all$status, "complete")) {
    s <- x$clinical$survival$all$statistics
    findings <- c(findings, sprintf("Complete-MB overall survival used the current %s cutoff (%.5g): log-rank p = %.3g and HR High/Low = %.3g (95%% CI %.3g-%.3g). This is exploratory and unadjusted.", s$cutoff_method, s$cutoff, s$logrank_p, s$hazard_ratio_high_vs_low, s$ci_lower, s$ci_upper))
  }
  d <- snapshot$depmap
  if (!is.null(d$statistics) && is.data.frame(d$statistics) && nrow(d$statistics)) findings <- c(findings, "DepMap findings use the same release and matched cell-line objects shown in the dashboard; dependency does not establish a therapeutic window.")
  if (!length(findings)) findings <- "No inferential result was available in the frozen session snapshot."
  findings
}

report_export_snapshot <- function(snapshot, out_dir) {
  save_selected_plot<-report_save_plot
  # Legacy table export below remains intact; the registry alone decides figures.
  report_save_plot<-function(...)invisible(NULL)
  tables <- file.path(out_dir, "tables"); figures <- file.path(out_dir, "figures"); metadata <- file.path(out_dir, "metadata")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(figures, recursive = TRUE, showWarnings = FALSE); dir.create(metadata, recursive = TRUE, showWarnings = FALSE)
  x <- snapshot$main
  report_write_table(x$retrieval$data, file.path(tables, "r2_patient_expression.csv"))
  report_write_table(x$analysis$statistics, file.path(tables, "r2_subgroup_statistics.csv"))
  report_write_table(x$analysis$summaries, file.path(tables, "r2_subgroup_summaries.csv"))
  report_save_plot(x$analysis$plot, "r2_subgroups", figures, 8, 5.8)
  if (!is.null(x$subtypes)) {
    for (nm in c("data", "summaries", "overall", "pairwise")) report_write_table(x$subtypes[[nm]], file.path(tables, paste0("r2_subtype_", nm, ".csv")))
    report_save_plot(x$subtypes$plot, "r2_subtypes", figures, 13, 7.5)
  }
  if (!is.null(snapshot$pfister)) {
    for (nm in c("data", "summaries", "counts")) report_write_table(snapshot$pfister[[nm]], file.path(tables, paste0("pfister_", nm, ".csv")))
    report_save_plot(snapshot$pfister$plot, "pfister_pan_cancer", figures, 10, 6)
  }
  if (!is.null(x$clinical)) for (ep in c("metastasis", "survival")) for (co in names(x$clinical[[ep]])) {
    a <- x$clinical[[ep]][[co]]; stem <- paste(ep, co, sep = "_")
    report_write_table(a$data, file.path(tables, paste0(stem, "_patients.csv")))
    report_write_table(a$statistics, file.path(tables, paste0(stem, "_statistics.csv")))
    report_write_table(if (ep == "metastasis") a$summaries else a$risk_table, file.path(tables, paste0(stem, "_details.csv")))
    report_save_plot(a$plot, stem, figures, 9, 6)
  }
  b <- snapshot$biology
  if (!is.null(b)) {
    report_write_table(b$annotation$rows, file.path(tables, "functional_annotation.csv"))
    report_write_table(b$targets$evidence, file.path(tables, "tf_target_evidence.csv"))
    report_write_table(b$targets$summary, file.path(tables, "tf_target_summary.csv"))
    set <- snapshot$settings$biology_target_set %||% b$default_target_set %||% "curated"; a <- b$enrichment[[set]]
    if (!is.null(a)) { report_write_table(a$rows, file.path(tables, paste0("tf_target_enrichment_", set, ".csv"))); report_save_plot(a$plot, paste0("tf_target_enrichment_", set), figures, 12, 9) }
  }
  h <- snapshot$hpa
  if (!is.null(h)) {
    for (nm in c("summary", "tissue", "brain", "single_cell", "subcellular", "cancer", "blood", "cell_line", "structure", "interaction")) {
      tab <- tryCatch(hpa_section_table(h, nm), error = function(e) data.frame())
      report_write_table(tab, file.path(tables, paste0("hpa_", nm, ".csv")))
    }
    report_save_plot(report_hpa_plot(h), "hpa_brain", figures, 9, 6)
  }
  d <- snapshot$depmap
  if (!is.null(d)) {
    for (nm in c("mb_models", "all_models", "mb_context", "matched", "statistics", "correlation")) report_write_table(d[[nm]], file.path(tables, paste0("depmap_", nm, ".csv")))
    if (length(d$plots)) for (nm in names(d$plots)) report_save_plot(d$plots[[nm]], paste0("depmap_", report_safe_name(nm)), figures, 9, 6)
  }
  report_write_table(report_publication_table(snapshot$publications), file.path(tables, "pubmed_results.csv"))
  catalog<-report_figure_catalog(snapshot,materialize=TRUE)
  selected<-Filter(function(f)f$selected,catalog)
  for(f in selected)save_selected_plot(f$plot,f$id,figures,
    attr(f$plot,'hpa_width') %||% if(f$id=='r2_subtypes')13 else 9,
    if(grepl('cerebellum|^hpa_',f$id))7 else 6)
  snapshot$interpretation<-report_interpretation(snapshot)
  snapshot$figure_manifest<-lapply(selected,function(f)f[setdiff(names(f),'plot')])
  jsonlite::write_json(snapshot$interpretation,file.path(metadata,'interpretation.json'),auto_unbox=TRUE,pretty=TRUE,null='null',na='null')
  jsonlite::write_json(snapshot$figure_manifest,file.path(metadata,'selected_figures.json'),auto_unbox=TRUE,pretty=TRUE,null='null',na='null')
  if(!is.null(snapshot$cerebellum))for(sp in c('human','mouse')) {
    z<-snapshot$cerebellum$species[[sp]]
    if(identical(z$status,'complete'))report_write_table(data.frame(species=sp,gene=z$gene$symbol,measured=z$n_measured,expressing=z$n_expressing),file.path(tables,paste0('cerebellum_',sp,'_summary.csv')))
  }
  saveRDS(snapshot, file.path(metadata, "report_snapshot.rds"))
  meta <- list(gene = snapshot$gene, generated_at = snapshot$generated_at, settings = snapshot$settings,
    app_version = snapshot$app_version, git_commit = snapshot$git_commit, R = R.version.string,
    sources = snapshot$sources, availability = snapshot$availability,interpretation_version=REPORT_EVIDENCE_VERSION,
    citations=snapshot$interpretation$evidence$citations,figure_selection=snapshot$figure_selection)
  jsonlite::write_json(meta, file.path(metadata, "report_metadata.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  writeLines(c("# Reproduce the integrated report", "", "Open R in the application repository, then run:", "", "```r", "source('R/load_core.R')", "core <- load_r2_core('.')", "snapshot <- readRDS('metadata/report_snapshot.rds')", "core$build_integrated_report(snapshot, 'reproduced-report', '.', include_pdf = TRUE)", "```"), file.path(out_dir, "README.md"))
  invisible(list(tables = tables, figures = figures, metadata = metadata))
}

report_find_pandoc <- function() {
  configured <- Sys.getenv("RSTUDIO_PANDOC", unset = "")
  candidates <- c(configured,
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    "C:/Program Files/RStudio/bin/pandoc")
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(file.path(candidates, if (.Platform$OS.type == "windows") "pandoc.exe" else "pandoc"))]
  if (length(hit)) hit[[1]] else ""
}

report_find_browser <- function() {
  configured <- Sys.getenv("MB_REPORT_BROWSER", unset = "")
  found <- unname(Sys.which(c("msedge", "google-chrome", "chromium", "chromium-browser")))
  candidates <- c(configured, found[nzchar(found)], "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe", "C:/Program Files/Microsoft/Edge/Application/msedge.exe")
  hit <- unique(candidates[nzchar(candidates) & file.exists(candidates)])
  if (length(hit)) hit[[1]] else ""
}

report_print_pdf <- function(html, pdf) {
  browser <- report_find_browser()
  if (!nzchar(browser)) return(FALSE)
  profile <- tempfile("report-browser-"); dir.create(profile)
  on.exit(unlink(profile, recursive = TRUE, force = TRUE), add = TRUE)
  uri <- paste0("file:///", gsub("\\\\", "/", normalizePath(html, winslash = "/", mustWork = TRUE)))
  args <- c("--headless", "--disable-gpu", "--no-pdf-header-footer", shQuote(paste0("--user-data-dir=", profile)), shQuote(paste0("--print-to-pdf=", normalizePath(dirname(pdf), winslash = "/", mustWork = TRUE), "/", basename(pdf))), shQuote(uri))
  tryCatch({
    system2(browser, args, stdout = TRUE, stderr = TRUE, wait = TRUE)
    # Edge is a GUI executable on Windows and can return before its headless
    # child finishes writing the file.
    for(i in seq_len(160L)) {
      if(file.exists(pdf) && is.finite(file.info(pdf)$size) && file.info(pdf)$size > 1000) return(TRUE)
      Sys.sleep(.25)
    }
    FALSE
  }, error = function(e) FALSE)
}

report_zip_dir <- function(report_dir, zip_file) {
  old <- setwd(dirname(report_dir)); on.exit(setwd(old), add = TRUE)
  base <- basename(report_dir)
  if (nzchar(Sys.which("zip"))) {
    status <- system2(Sys.which("zip"), c("-r", shQuote(normalizePath(zip_file, winslash = "/", mustWork = FALSE)), shQuote(base)))
  } else if (nzchar(Sys.which("tar"))) {
    # Windows ships bsdtar; -a selects ZIP from the destination extension and
    # is substantially faster than Compress-Archive for large RDS snapshots.
    status <- system2(Sys.which("tar"), c("-a", "-c", "-f", shQuote(normalizePath(zip_file, winslash = "/", mustWork = FALSE)), shQuote(base)))
  } else if (.Platform$OS.type == "windows") {
    ps <- Sys.which("powershell"); if (!nzchar(ps)) stop("No ZIP utility is available")
    cmd <- sprintf("Compress-Archive -LiteralPath '%s' -DestinationPath '%s' -Force", gsub("'", "''", normalizePath(report_dir, winslash = "\\", mustWork = TRUE)), gsub("'", "''", normalizePath(zip_file, winslash = "\\", mustWork = FALSE)))
    status <- system2(ps, c("-NoProfile", "-Command", shQuote(cmd)))
  } else stop("No ZIP utility is available")
  if (!identical(status, 0L) || !file.exists(zip_file)) stop("ZIP package creation failed")
  zip_file
}

build_integrated_report <- function(snapshot, output_dir, root = ".", include_pdf = TRUE) {
  stopifnot(is.list(snapshot), !is.null(snapshot$main), nzchar(snapshot$gene))
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  report_export_snapshot(snapshot, output_dir)
  snapshot_path <- file.path(output_dir, "metadata", "report_snapshot.rds")
  pandoc <- report_find_pandoc(); if (nzchar(pandoc)) Sys.setenv(RSTUDIO_PANDOC = pandoc)
  if (!requireNamespace("rmarkdown", quietly = TRUE) || !rmarkdown::pandoc_available()) stop("R Markdown and Pandoc are required for the integrated HTML report")
  html_name <- paste0(snapshot$gene, "_integrated_gene_report.html")
  html <- rmarkdown::render(file.path(root, "report", "integrated_gene_report.Rmd"),
    output_file = html_name, output_dir = output_dir, knit_root_dir = output_dir,
    params = list(snapshot_path = snapshot_path, root = root), envir = new.env(parent = globalenv()), quiet = TRUE)
  pdf <- file.path(output_dir, paste0(snapshot$gene, "_integrated_gene_report.pdf"))
  pdf_ready <- isTRUE(include_pdf) && report_print_pdf(html, pdf)
  zip_file <- file.path(dirname(output_dir), paste0(basename(output_dir), ".zip"))
  report_zip_dir(output_dir, zip_file)
  list(html = html, pdf = if (pdf_ready) pdf else NULL, zip = zip_file, directory = output_dir,
    pdf_status = if (pdf_ready) "ready" else "unavailable")
}
