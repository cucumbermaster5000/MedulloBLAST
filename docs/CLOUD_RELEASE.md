# MedulloBLAST Cloud release

Maintainer: Dave Ng <davengis@gmail.com>.
Published dedicated repository: https://github.com/cucumbermaster5000/MedulloBLAST.
Published data release: cerebellum-2026-09-25-v1.

## Validation status and runtime

The Free account has been verified through the Cloud API: 4096 MB memory, one CPU,
five apps, 20 monthly compute credits, and a 1 GiB source bundle limit. Automatic
republishing is disabled. The live app uses one Shiny process and two connections.

Cloud installed all 80 locked packages and started the app under Ubuntu 22.04 / R
4.6.0. Hosted HLX, MYC and SLC16A1 cold/warm queries, cancellation and two sessions
passed. Individual CSV, PNG and PDF exports passed; PDFs were rendered for review.
Separate local Linux installation was skipped at the maintainer's request. Local
preparation still uses Windows R 4.6.1. RELEASE_VALIDATION.json is a preparation
snapshot, not a claim that every future source revision has been tested on Cloud.

Production staging maps config/runtime-renv.lock to renv.lock, README.cloud.md to
README.md, config/cloud.env.example to .env.example, and config/cloud.gitignore to
.gitignore. Test dependencies, atlas preparation scripts, caches, credentials, local
libraries, old reports and development Git history are excluded. The prepared atlas
must never be built on Cloud. The complete app and local report implementation remain.

## Exact initial Cloud settings

- Framework: Shiny for R. Primary file: app.R. Dedicated app repository.
- RAM: 4 GB. CPU: 1. Max Worker Processes: 1. Max Connections: 2 initially.
- Automatically publish on push: off.
- MB_MAX_WORKERS=1; MB_INTEGRATED_REPORT=false; MB_REPORT_PDF=false.
- OMP_NUM_THREADS=1; OPENBLAS_NUM_THREADS=1; MKL_NUM_THREADS=1.
- MB_AUTO_PROVISION_DATA=true; MB_ALLOW_DATA_DOWNLOADS=false.
- R2_DEPMAP_CACHE_DIR=/tmp/medulloblast-v1/cache.
- MB_DATA_ROOT=/tmp/medulloblast-v1/references.
- MB_WORKER_TIMEOUT_SECONDS=600; MB_PROVISION_TIMEOUT_SECONDS=1800.
- NCBI_EMAIL=davengis@gmail.com. NCBI_API_KEY is optional and must stay in Cloud Variables.

Cloud sets R_CONFIG_ACTIVE=connect_cloud. The root app applies these environment defaults
when they are unset; explicit Cloud Variables take precedence. .env.example is documentation,
not a file automatically loaded by the app. Cloud CPU/RAM/process/connection settings must
are configured and verified through the Cloud API. Integrated reports and interpretation overview are
unavailable in this profile; biological modules and individual downloads remain.

Local report-disabled testing peaked at 1.827 GB including active and waiting children,
leaving 54.3% of a conservative 4 GB budget. The 20% margin target is 3.2 GB. See
FREE_TIER_RESOURCE_CHECK.md for methods, timings and limitations. Windows RSS is not Cloud acceptance. The first hosted run reached almost 4 GiB
because Linux retained about 2.5 GB of inactive reference-file cache. The Cloud
provisioner now flushes each verified large file and requests file-specific cache
release using GNU sync/dd, without deleting or modifying the data. The versioned
runtime directory starts empty for the final recovery and capacity check. Final
peak measurements must be taken from the corrected revision; the initial failure
must not be described as passing the 20% headroom requirement.

## Reference recovery

A background worker streams the prepared datasets using config/datasets.json, with fixed
URLs, sizes and SHA-256 hashes. It writes temporary files, checks them, and exposes a
complete dataset only after validation. Interrupted downloads can retry; corrupt complete
directories lose their ready marker and are quarantined. No atlas preparation runs here.
The UI has loading/failure/retry states, and a pending cerebellum gene query starts when
the data become ready. Provisioning has its own 1800-second deadline, separate from the
600-second analysis deadline. Other analyses remain usable while downloads run.

Current Shiny DepMap analyses still use the current-release Breadbox API. The separately
pinned 24Q4-v1 archive supports legacy consumers; it never silently substitutes for current
results. Original archive metadata, source labels and attribution are retained.

Cloud files are ephemeral. The proposed cerebellum URLs are not live until the reviewed
release is published. Independent restarts therefore remain unverified. Dataset recovery
fixture tests pass, but do not establish production download speed, disk use or memory.
https://docs.posit.co/connect-cloud/user/platform/system.html

## Publication review and first-host acceptance

scripts/prepare-cloud-release.R builds a new outputs/cloud-release-*/app folder without
changing the development history. Its parent contains the exact inventory and SHA256SUMS.
Review that candidate and the seven prepared cerebellum assets plus ATTRIBUTION.txt and
SHA256SUMS before public upload. Source authors, DOI, CC BY 4.0 notice, Ensembl provenance,
processing changes and original hashes must accompany the prepared data.

After publication is authorized, publish the fixed data release and app repository, verify
asset URLs/hashes, connect the repository to Cloud, and apply the settings above. Check build
logs for package/R compatibility; test empty-cache start, retry/restart, HLX, MYC, SLC16A1,
cancellation, two sessions and individual downloads. Inspect Cloud memory and CPU charts.
Repeat from another computer with the local app stopped. Do not announce the site until
these checks pass. No GitHub or Cloud publication occurs during local preparation.