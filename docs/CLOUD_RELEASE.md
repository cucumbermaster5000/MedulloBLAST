# MedulloBLAST Cloud release candidate

Maintainer: Dave Ng <davengis@gmail.com>.
Proposed dedicated repository: cucumbermaster5000/MedulloBLAST.
Proposed data release: cerebellum-2026-09-25-v1.

## Validation status and runtime

The maintainer has created a Free Connect Cloud account. Its settings and entitlements
have not yet been inspected. At the maintainer's request, separate Linux testing is
deferred: the first Cloud build and hosted acceptance run will check Linux compatibility.
No local Linux installation is needed, and no Linux pass is claimed.

The source-only candidate targets R 4.6.0, the highest version currently documented by
Connect Cloud. Its manifest and production lock declare that target. Preparation and
local tests use Windows R 4.6.1; RELEASE_VALIDATION.json records both versions honestly.
Declared package R constraints and exact installed versions are checked, but this does
not prove package installation or execution on Linux/R 4.6.0. The development DESCRIPTION
and lockfiles retain their original runtime. Any failure on Cloud must be fixed before
announcing the site. https://docs.posit.co/connect-cloud/user/platform/r.html

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
- R2_DEPMAP_CACHE_DIR=/tmp/medulloblast/cache.
- MB_DATA_ROOT=/tmp/medulloblast/references.
- MB_WORKER_TIMEOUT_SECONDS=600; MB_PROVISION_TIMEOUT_SECONDS=1800.
- NCBI_EMAIL=davengis@gmail.com. NCBI_API_KEY is optional and must stay in Cloud Variables.

Cloud sets R_CONFIG_ACTIVE=connect_cloud. The root app applies these environment defaults
when they are unset; explicit Cloud Variables take precedence. .env.example is documentation,
not a file automatically loaded by the app. Cloud CPU/RAM/process/connection settings must
still be selected in the service UI. Integrated reports and interpretation overview are
unavailable in this profile; biological modules and individual downloads remain.

Local report-disabled testing peaked at 1.827 GB including active and waiting children,
leaving 54.3% of a conservative 4 GB budget. The 20% margin target is 3.2 GB. See
FREE_TIER_RESOURCE_CHECK.md for methods, timings and limitations. Dataset provisioning
and Linux resource use remain to be measured; Windows RSS is not Cloud acceptance.

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