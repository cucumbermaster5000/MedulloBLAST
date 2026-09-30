# Free-tier resource check — 30 September 2026

The analysis workload fits the measured Windows memory budget, but integrated
reports do not meet the 20% headroom target. Disabling only server PDF rendering is
insufficient: simultaneous HTML/ZIP reports peaked at 3.535 GB. The proposed Cloud
profile disables the integrated report with `MB_INTEGRATED_REPORT=false`, preserving
all biological analysis modules and their individual downloads. Full reports remain
available locally by default. This is a local resource screen, not Cloud certification.

## Published target and method

Posit documents 4 GB RAM and 1 CPU for Connect Cloud Free:
https://posit.co/blog/migrating-connect-cloud-posits-unified-publishing-solution
Use 3.2 billion bytes as the conservative acceptance threshold (20% headroom).
The Free account is now created; its actual settings and usage credits remain unconfirmed.

`scripts/check-resources.R` starts an isolated server on port 3888 and sets its CPU
affinity to one allowed logical processor before application startup (helper-process exception documented below). OMP, OpenBLAS and MKL thread counts are one, and MB_MAX_WORKERS=1.
The browser driver exercises actual Shiny controls and downloads. It samples the sum
of RSS for the server and all recursive descendants, including waiting R workers,
supervisors and PDF browser processes. The original R sampler targeted 200 ms but
achieved a median 3.75 seconds; subsequent runs use the native 100 ms sampler. Client browsers
and the measurement process are excluded. Shared pages can be counted more than once;
short spikes between samples can be missed. This is not a hard memory-capped container.

Tests use Windows R 4.6.1, a new empty query cache and the already prepared cerebellum
atlas on disk. They do not establish Ubuntu/R 4.6.0 compatibility, dependency-build
memory, or memory/disk usage during production dataset downloads. Cloud CPU speed,
Linux page-cache accounting and network latency can differ. Confirm final settings
using Cloud's own utilization charts before describing this as Cloud acceptance.

## Full-PDF baseline

Artifacts: `outputs/resource-check-20260929-171337`.

| Workload | Peak server-tree RAM (GB) | Time (seconds) |
|---|---:|---:|
| First HLX query | 1.14 | 82 |
| Repeat HLX | 1.25 | 45 |
| First MYC query | 1.87 | 199 |
| Repeat MYC | 1.36 | 47 |
| First SLC16A1 query | 1.62 | 74 |
| Repeat SLC16A1 | 1.34 | 46 |
| Cancel MYC and switch to HLX | 1.84 | 63 |
| Two sessions: HLX and MYC | 1.73 | 111 |
| First complete report with PDF | 3.77 | 127 |
| Repeated query and complete report | 3.42 | 127 |

All ten workload checks passed and there were no JavaScript errors. Both HTML, PDF
and ZIP downloads completed. The report peak leaves only 5.8% headroom, failing the
20% target despite staying below 4 GB in this run. Reports run on demand; the analysis
peak alone gives no reason to remove a biological module.

An earlier run (`resource-check-20260929-165220`) was stopped after detecting a
repeated UI refresh. Disabled provisioning polled every second and invalidated the
cerebellum-dependent report snapshot. Stable polling now invalidates only when the
dataset state changes. Regression tests cover disabled and unchanged polling.

## Proposed Cloud configuration

- RAM 4 GB; CPU 1; Max Worker Processes 1; Max Connections 2 initially.
- MB_MAX_WORKERS=1; MB_INTEGRATED_REPORT=false; MB_REPORT_PDF=false; OMP_NUM_THREADS=1;
  OPENBLAS_NUM_THREADS=1; MKL_NUM_THREADS=1.
- Keep automatic republishing disabled initially.
- Integrated report generation and overview are disabled; individual analysis downloads remain available.
- Do not raise connection limits above the tested workload without another test.

Reduced-profile measurements and artifact checks must pass before accepting this
configuration. Source release publication, restart recovery from public assets,
supported-R testing and deployment remain separate unfinished gates.

The regression suite passes after fixing Windows quarantine path separators and
revoking the ready marker for corrupted datasets. One opt-in live-R2 test is skipped;
three existing isolated-report `package:stats` warnings remain. The browser workload
does exercise live R2 retrieval separately.

Nothing has been published. A Linux installation on the maintainer's computer is not
required for the eventual host check: it can be performed on Connect Cloud or a
reviewed Linux CI runner once the deployment/runtime prerequisites are resolved.

## Fast-sampled report comparison

Artifacts: `outputs/resource-check-20260929-183952/verified-results.json`.
Sampling: 5,285 samples; median interval 94 ms, maximum 156 ms. Unlike the older
baseline sampler (median 3.75 seconds), this captures brief memory peaks.

| Workload, integrated PDF disabled | Peak GB | Headroom | Seconds |
|---|---:|---:|---:|
| Two warm sessions | 1.761 | 56.0% | 82 |
| First HTML/ZIP report | 2.791 | 30.2% | 107 |
| Repeated query and HTML/ZIP report | 3.204 | 19.9% | 112 |
| Two simultaneous report requests | 3.535 | 11.6% | 206 |

All workload checks passed, with no JavaScript errors. All four ZIPs passed CRC
checks, and their embedded HTML exactly matched the separately downloaded HTML.
The ZIPs contain the expected snapshot and tables. The reports work, but the last
two workloads fail the conservative 3.2 GB threshold. No integrated PDF is present
in this profile, as expected. Individual figure files remain in the packages.

Affinity sampling found all measured R processes on mask 1. One Windows conhost
helper used mask 255 and only 0.0625 CPU-seconds total. This is not an OS-enforced
one-CPU quota for the whole process tree and is not a substitute for Linux limits.

A separate fast-sampled run (`resource-check-20260929-181819`) passed all six gene
queries, cancellation and two-session analysis. Peaks were 1.184-1.889 GB for single
queries and 2.014 GB for two sessions. First queries took 99-205 seconds; repeats
46-77 seconds. It used an empty query cache with prepared reference data on disk,
not an empty reference-data cache. Its report checks incorrectly expected PDF
because the driver did not receive the PDF-off flag; that harness problem was fixed
for the successful focused report run above. Do not treat the earlier report checks
as acceptance results.

## Remaining release gates

- Report-disabled local resource screen passed; results below. Cloud validation remains pending.
- Restore and test supported Ubuntu 22.04 / R 4.6.0 before changing the current
  R 4.6.1 requirement or regenerating final deployment metadata.
- Publish reviewed, checksummed reference assets only after authorization; their
  proposed GitHub Release URLs are not yet live. Fixture recovery tests passed,
  but independent Cloud restart recovery and download-time resource use are untested.
- Stage the dedicated source-only repository and regenerate/check the final manifest
  and bundle. Existing deployment artifacts are drafts after these changes.
- Confirm the new account's actual Free allowances, then verify Cloud resource charts,
  cold starts and individual exports from another computer with the local server off.

Nothing has been published. A Linux environment on this computer is not required:
Linux verification can use a reviewed CI runner or the eventual Cloud environment.
## Final report-disabled profile: PASS (local Windows screen)

Artifacts: `outputs/resource-check-20260930-112105/verified-results.json` and
`workload.json`. Settings: MB_INTEGRATED_REPORT=false, MB_REPORT_PDF=false,
MB_MAX_WORKERS=1; one logical CPU for R processes. New empty query cache,
pre-existing prepared references; provisioning disabled for this resource test.

| Workload | Peak GB | Seconds |
|---|---:|---:|
| First HLX | 1.041 | 142 |
| Repeat HLX | 1.292 | 45 |
| First MYC | 1.694 | 200 |
| Repeat MYC | 1.284 | 31 |
| First SLC16A1 | 1.528 | 68 |
| Repeat SLC16A1 | 1.278 | 20 |
| Cancellation | 1.381 | 39 |
| Two sessions: HLX and MYC | 1.827 | 72 |

All eight analysis workloads and the disabled-report UI check passed, with no
JavaScript errors. Peak 1.827 GB leaves 54.3% of the conservative 4 GB limit unused,
exceeding the 20% target. Sampling: 6,336 samples, median 109 ms, maximum 297 ms.
The regression suite also passed (outputs/regression-20260930.log), including a
server-side test that disabled reports cannot start workers even on injected input.
One opt-in live-R2 unit test was skipped; three existing isolated-report stats
warnings remain. The browser workload independently used live source retrieval.

The regression suite ran alongside the early part of this screen; its processes
are excluded from app RSS but could affect timings. These timings are indicative,
not dedicated-machine latency benchmarks. Current user-facing performance is tens
of seconds for cached queries and up to several minutes for first queries.

Recommended first Free deployment: all biological modules, individual downloads,
one Shiny process, one heavy worker, at most two concurrent sessions; integrated
report/interpretation overview and its HTML/PDF/ZIP generation disabled. No report
source code or existing reports were deleted; local mode retains them by default.
This tested profile is a configuration example and is not automatically applied to
an existing local preview or to a Cloud account. Nothing was published.

Update, 30 September: the maintainer requested deferring separate Linux testing to
the first Cloud deployment. A source-only R 4.6.0 target candidate is being prepared;
this does not change the Windows-only status of the measurements above.
