# MedulloBLAST

An independent Shiny application for medulloblastoma gene exploration, developing
cerebellum expression, functional evidence, DepMap and integrated biological reports.
Maintainer: Dave Ng <davengis@gmail.com>.

Live app: https://cucumbermaster5000-medulloblast.share.connect.posit.cloud/

Published with maintainer approval. Cloud startup and hosted biological queries
have passed under Linux/R 4.6.0. See docs/CLOUD_RELEASE.md for settings and
the scope of resource validation.

Use repository-root app.R. The production dependency lock is renv.lock; Connect
Cloud consumes the accompanying manifest.json. Reference data are downloaded and
verified in the background when MB_AUTO_PROVISION_DATA=true. No Codex or developer
computer is required by the deployed app. Upstream data services remain required.

Current Shiny DepMap analyses retain the official current-release Breadbox API.
The separately pinned 24Q4 archive supports legacy analyses; it is never silently
substituted for current results. Source release labels and provenance are retained.

Reference downloads, usage limitations, attribution and provisioning behavior are
documented in docs/CLOUD_RELEASE.md and docs/CEREBELLUM.md. The GitHub Release URLs
in config/datasets.json point to the published, versioned data release.

The Free-tier profile sets MB_INTEGRATED_REPORT=false to retain memory headroom.
Biological modules and individual analysis downloads remain available. Full reports
remain available in local deployments. See docs/FREE_TIER_RESOURCE_CHECK.md.

Shared free-hosting capacity: up to five connected sessions, one heavy analysis
task at a time. Inactive sessions disconnect after 15 minutes; the app warns
before disconnection and provides a reconnect button. Download results before leaving.

The live counter shows connected browser sessions out of five; multiple tabs may use multiple slots. Release preparation stamps the UTC update date and unique build ID into DESCRIPTION automatically. The displayed semantic version also comes from DESCRIPTION.

Cavalli plots offer on-demand patient CSVs and reproduction ZIPs beside each plot. ZIPs include inclusion/exclusion audits, available statistics, source and method metadata, and a standalone R plotting script. Downloads use completed results without refitting analyses. Cerebellum exports are unchanged.

The Pediatric Pan-Cancer panel includes a sourced legend for all 15 Pfister cancer-type abbreviations and the FPKM expression unit. Original dataset labels are preserved.
