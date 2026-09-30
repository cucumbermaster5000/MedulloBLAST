# Root entry point for hosting; deploy the project bundle, not shiny/ alone.
source(file.path('config','cloud-defaults.R'),local=TRUE)
apply_cloud_defaults()
start_cloud_resource_diagnostics()
source(file.path("shiny","app.R"),local=TRUE)$value

