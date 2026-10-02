# Open this file in RStudio and click Run App, or shiny::runApp("shiny").
root <- if (file.exists(file.path("..", "R", "load_core.R"))) ".." else "."
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "load_core.R"), local = TRUE)
core <- load_r2_core(root)
for (package in c("shiny", "bslib")) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install web dependencies with: Rscript scripts/install-dependencies.R --shiny")
}
# The Cloud entry point runs from the repository root, not shiny/.
shiny::addResourcePath("medulloblast-assets", file.path(root, "shiny", "www"))
cfg <- core$read_server_config(file.path(root, "config", "defaults.json"))
cfg$output_dir <- file.path(root, "outputs")
cfg$app_root <- root
cfg$public_app <- TRUE
cfg <- core$provisioning_config(cfg)
options(shiny.sanitize.errors=TRUE)
source(file.path(root,"shiny","display.R"),local=TRUE)
source(file.path(root,"shiny","capacity.R"),local=TRUE)
source(file.path(root,"shiny","runtime.R"),local=TRUE)
source(file.path(root,"shiny","provisioning.R"),local=TRUE)
reference_service <- explorer_provisioning_service(cfg,core)
explorer_log(paste("startup R",getRversion()))
source(file.path(root,"shiny","biology_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","publications_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","depmap_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","hpa_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","cerebellum_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","report_interface.R"),local=TRUE,encoding="UTF-8")
source(file.path(root,"shiny","dashboard.R"),local=TRUE,encoding="UTF-8")
source(file.path(root, "shiny", "interface.R"), local = TRUE, encoding = "UTF-8")
shiny::shinyApp(explorer_ui(), explorer_server(core, cfg, reference_service))
