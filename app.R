# Root entry point for hosting; deploy the project bundle, not shiny/ alone.
source(file.path('config','cloud-defaults.R'),local=TRUE)
apply_cloud_defaults()
source(file.path("shiny","app.R"),local=TRUE)$value
