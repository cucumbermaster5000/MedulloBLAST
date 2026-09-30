# MedulloBLAST

An independent Shiny application for medulloblastoma gene exploration, developing
cerebellum expression, functional evidence, DepMap and integrated biological reports.
Maintainer: Dave Ng <davengis@gmail.com>.

This source-only release is prepared locally. It is not yet approved for public
deployment. See docs/CLOUD_RELEASE.md for the release gates and exact Cloud settings.

Use repository-root app.R. The production dependency lock is renv.lock; Connect
Cloud consumes the accompanying manifest.json. Reference data are downloaded and
verified in the background when MB_AUTO_PROVISION_DATA=true. No Codex or developer
computer is required by the deployed app. Upstream data services remain required.

Current Shiny DepMap analyses retain the official current-release Breadbox API.
The separately pinned 24Q4 archive supports legacy analyses; it is never silently
substituted for current results. Source release labels and provenance are retained.

Reference downloads, usage limitations, attribution and provisioning behavior are
documented in docs/CLOUD_RELEASE.md and docs/CEREBELLUM.md. The GitHub Release URLs
in config/datasets.json are reserved for the proposed release and are not live yet.

The Free-tier profile sets MB_INTEGRATED_REPORT=false to retain memory headroom.
Biological modules and individual analysis downloads remain available. Full reports
remain available in local deployments. See docs/FREE_TIER_RESOURCE_CHECK.md.
