# Flumina — Nextflow conversion handoff

**Status:** Nextflow conversion complete and verified. FluMut integration extended to
low-frequency variants. Command-line interface and container distribution reworked to
match FluPore. Not yet committed.

---

## Session 3 (2026-07-25): FluPore-style CLI and container distribution

Goal: run Flumina and FluPore the same way — pull one image from Docker Hub, run it
directly under Docker or Apptainer, drive it with short flags and defaults, and keep the
configuration file for the many-parameter case.

### 1. `Scripts/flumina` launcher (new)

A bash launcher in FluPore's shape (`Scripts/flumina_2.0.0`, with `Scripts/flumina` as a
version symlink, mirroring `flupore`/`flupore_1.0.0`). It wraps `nextflow run main.nf`.

Argument letters deliberately match FluPore wherever the meaning matches — `-i -o -t -r
-f -d -s -h --version` mean the same thing in both programs.

| Flag | Meaning | Default |
|---|---|---|
| `-i` | input read directory | `raw_reads` |
| `-o` | output directory | `Flumina_results` |
| `-t` | threads | 4 |
| `-M` | max memory per step | `6.GB` |
| `-c` | configuration file | `config.cfg` if present |
| `-n` | rename CSV | `file_rename.csv` |
| `-r` | reference FASTA | bundled `reference.fa` |
| `-a` | curated AA database | bundled `curated_database.csv` |
| `-m` | metadata CSV | `metadata.csv` |
| `-g` | grouping column | `discrete_host` |
| `-d` `-q` `-f` | min depth / quality / allele frequency | 100 / 30 / 0.01 |
| `-l` | low-frequency FluMut screening at this AF | off |
| `-x` | IRMA config file | `irma_config.sh` if present |
| `-p` | execution profile | `standard` in container, else `docker` |
| `-w` | work directory | `work` |
| `-s` | skip steps: `i`=IRMA, `m`=FluMut | `0` |
| `-R` `-N` | resume / dry run | off |

**Precedence: defaults < config file < command line.** CLI values are held aside during
`getopts` and applied after the config file is read, so a saved configuration can be
reused with one-off overrides (`flumina -c config.cfg -t 24 -o otherFolder`).

The config reader understands the legacy Snakemake-era keys: `DISABLE_IRMA` is mapped
onto the inverted `run_irma`, and `FORCE_UNLOCK` / `OVERWRITE` / `TEMPDIR` / `CLUSTER_JOBS`
are reported as no-longer-needed rather than silently ignored. Unrecognised keys are
reported too, so a typo like `MIN_DPETH` does not quietly do nothing.

### 2. Container is now self-sufficient and is the entry point

`docker run chutter/flumina -i raw_reads -o results` and
`apptainer run flumina_latest.sif -i raw_reads -o results` both work, exactly as FluPore does.

- `ENTRYPOINT ["flumina"]`; pipeline installed to `/opt/flumina`, launcher on PATH
- Nextflow itself is installed **into** the image (pinned `NXF_VER=24.10.4`), so the image
  needs nothing on the host but a container runtime
- The Nextflow framework jar is **pre-warmed at build time** into `/opt/nextflow`
  (`chmod 777`, because Apptainer runs as the invoking user and Docker is commonly run
  with `-u`). Without this every container start re-downloads it, and an air-gapped
  cluster fails outright. Verified with `--network none -u 12345:12345`
- `NXF_OFFLINE=true` stops run-time plugin/version checks
- `FLUMINA_CONTAINER=1` tells the launcher it is already inside the image, so it defaults
  to the `standard` profile and runs tasks in-place instead of nesting a container
- Nextflow install and pipeline COPY sit **after** the conda layer so editing the
  pipeline does not rebuild the slow 4 GB conda layer

### 3. Resource ceilings (`--max_cpus`, `--max_memory`)

Found by running the container: process labels hardcoded 12 GB and 32 GB, so **every run
died immediately** on a stock Docker Desktop VM (7.7 GB) with "process requirement exceeds
available memory". Labels are now capped by `params.max_cpus` / `params.max_memory` via
`capCpus()` / `capMem()`, defaulting to 4 CPUs / 6 GB, which fits a laptop. The `slurm`
profile raises them to 8 / 32 GB, so cluster runs are unaffected.

Note the idiom: the cap functions must be `def` **methods** with the requests written as
lazily-evaluated closures (`cpus = { capCpus(4) }`). Assigning a closure to a variable and
calling it inside `withLabel` fails with `Unknown method invocation capCpus`.

### 4. `irma_config.sh` wired into the pipeline

IRMA sources `./irma_config.sh` from its working directory — that is IRMA's own behaviour
(confirmed in the IRMA script, and in the container's copy). So `--irma_config` simply
stages the user's file into each IRMA task under that name. This works identically under
Docker and Apptainer because nothing in the read-only image is modified.

### 5. Container image reference

The `docker` and `apptainer` profiles previously named a locally-built `flumina:latest`,
which nobody following the README would have. Both now resolve `params.container`
(`chutter/flumina:latest`), with the apptainer profile prefixing `docker://` so Apptainer
pulls and converts it itself — no separate `.sif` to keep in step with the Dockerfile.

### 6. README rewritten in FluPore's structure

Same section order and voice as FluPore's manual: test data → Apptainer (recommended) →
Docker → individual dependencies → required input files → full arguments menu →
configuration files → cluster → output table.

---

## Session 2 (2026-07-25): FluMut low-frequency variant screening

### R rewrite of the FluMut helpers
- `Scripts/rename_for_flumut.R` replaces `rename_for_flumut.py`, for consistency with the
  rest of the pipeline, which is R throughout. Same IRMA→FluMut header translation
  (`>A_HA_H5` → `>sample_HA`)
- `Scripts/apply_lofreq_to_consensus.R` (new) applies LoFreq variants above an allele
  frequency threshold to the IRMA consensus, producing mutated pseudo-consensus sequences

### `FLUMUT_LOWFREQ` process
Separate from `FLUMUT` rather than merged into it, so the two result sets can never be
confused: consensus markers land in `flumut/`, low-frequency markers in `flumut_lowfreq/`.
Opt-in via `--flumut_lowfreq` / `flumina -l <freq>`, threshold via
`--flumut_freq_threshold` (default 0.01).

Use cases: minority H5N1 markers within a host, co-infections or revertants below
consensus level, and tracking emerging resistance across replicates.

**Caveat worth checking before publishing on it:** LoFreq calls positions against the
bwa-aligned reference, while IRMA builds its consensus de novo. If a sample's IRMA
consensus carries an indel relative to the reference, variant coordinates will not line up
and mutations will be written to the wrong positions. Fine for the reference-like samples
tested so far; a coordinate-translation step is needed before trusting this on divergent
samples.

---

## Session 1: Nextflow 2.0 DSL2 conversion

- All 16 Snakemake rules converted to Nextflow processes
- Full relocatability: no absolute paths; all inputs staged into task work dirs; generated
  `config.cfg` uses `OUTPUT_DIRECTORY="."`
- Docker containerization pinned to `linux/amd64` (bioconda has no arm64 gatk4/irma/lofreq)
- Head-to-head equivalence verified against the legacy Snakemake pipeline on `test_dataset`
  (LoFreq byte-identical; GATK identical after removing `bwa mem -t`)

### Bug fixes made during the conversion

| File | Bug | Fix |
|---|---|---|
| `outputSummary.R` | `GROUP_NAMES` treated as a file path, not a column name | map column onto results; validate and list available columns |
| `organizeReads.R` / `organizeIRMA.R` | `is.logical()` on a character vector is always FALSE, so `OVERWRITE` never took effect | parse the string properly |
| `main.nf` (R_ANALYSIS) | reference not at outdir root, where findAAChanges.R expects it | second publishDir for `reference.fa` |
| `main.nf` (FILTER_VARIANTS) | `filter_INDEL` never ran in Snakemake | implemented; revealed a dead rule in the original |
| `main.nf` (GATHER_SAMPLE_VCFS) | identical VCF filenames across samples risked silent label shuffling | per-sample directories before staging |
| `bwa mem` | threading changed insert-size estimates, diverging borderline GATK calls | dropped `-t` to match Snakemake |

---

## Verification status

| Check | Result |
|---|---|
| `flumina -h`, `--help`, `--version` | pass, native and in container |
| Config auto-pickup, legacy keys, typo and obsolete-key reporting | pass |
| CLI overrides config file | pass |
| Argument validation and error messages | pass |
| Container entrypoint, offline, arbitrary uid (`--network none -u 12345:12345`) | pass |
| `nextflow config` resolves for standard/docker/apptainer/slurm | pass |
| Full pipeline run inside the container on `test_dataset` | in progress at handoff |
| Apptainer run | **not tested** — no Apptainer on this Mac; needs a cluster check |
| Docker Hub image | **not pushed** — `chutter/flumina` must be built and pushed |

---

## Snakemake removal and cleanup (2026-07-25)

The Nextflow port is validated, so the legacy pipeline and its scaffolding were removed
rather than carried along. Deleted:

- `Flumina` — the Snakemake driver, with hardcoded `/Users/chutter` conda paths
- `Scripts/snakefile_*.smk` — the three rule files
- `snakemake=7.32.4` from `environment.yaml`, which cut the image from 4.17 GB to 3.77 GB
- `Scripts/organizeReads.R` — only the old driver called it; Nextflow stages reads by symlink
- `test_dataset/config.cfg` and `test_dataset/irma_config.sh` — absolute `/Users/chutter`
  paths, superseded by the `test` profile

The trade-off accepted here: the image can no longer run both pipelines, so the
head-to-head output diff that validated the conversion cannot be re-run. That comparison
is recorded in the Session 1 notes above and its result is committed history.

Explanatory comments in `main.nf` and `nextflow.config` still reference Snakemake. Those
are design rationale — why the output layout, the missing `bwa -t`, and the resurrected
`filter_INDEL` rule are the way they are — and are deliberately kept.

`Scripts/makeGTF.R`, `runSNPGenie.R` and `runWFABC.R` are unreferenced by `main.nf` but
were **kept**: they are unfinished features with live params in `nextflow.config`, not
dead code.

Also fixed: the shipped `irma_config.sh` set `TMP=/Flumina_test/Flumina/irma_tmp`, a path
that exists nowhere. Since the README tells users to copy that file, verbatim use sent
IRMA's temp writes somewhere unwritable. Now `TMP=./irma_tmp`, resolving inside the
per-task work directory.

---

## Cluster parallelism (`-p slurm,apptainer`)

Per-step job submission is a core design goal — it is where Flumina's speed on large runs
comes from — so it is the primary documented cluster path, not a footnote.

**It is fully compatible with the container.** Nextflow runs on the host and submits each
step as its own SLURM job; each of those jobs runs inside the Apptainer image. Full
parallelism and one pinned software environment. The only constraint is that Nextflow
itself cannot run *inside* the container in this mode, because it shells out to `sbatch`,
which is not in the image. `module load nextflow` covers this on most clusters.

New launcher flags: `-Q` queue, `-A` scheduler account options, `-j` max jobs in flight.
New params: `slurm_queue`, `slurm_options`, `queue_size`, `apptainer_cache`.

### Three things fixed to make this actually work

1. **The scheduler queue and account were hardcoded** to `priority` / `--qos=vpru -A
   nadc_iav` in the slurm profile. Any other site's run would have been rejected outright.
   Now params with those values as defaults, overridable per run.

2. **Nextflow rejects any `--param` whose value starts with a dash** — and scheduler
   options nearly always do (`-A account`, `--qos=...`). `--slurm_options "-A myproject"`
   fails with `Unknown option`. The launcher therefore writes the scheduler settings into
   a generated config file passed with `-c`, setting `process.queue`,
   `process.clusterOptions` and `executor.queueSize` directly. Setting the directives
   rather than params also sidesteps config ordering: a `-c` file merges last, so it wins
   over the profile defaults. Verified: `-c` precedence confirmed against a known param.

3. **BSD `mktemp` only expands trailing X's.** The first version used
   `mktemp .../flumina_cluster.XXXXXX.config`, which returned that name *literally* — so
   two concurrent runs on the same node would share one file. Now `mktemp -d` with a
   fixed filename inside, cleaned up by an EXIT/INT/TERM trap.

Also added `apptainer.cacheDir` (default `$HOME/.apptainer/flumina`, override with
`--apptainer_cache`). Without it, a hundred task jobs starting at once each race to
convert the same image. It must live on a filesystem every compute node can see.

`executor.submitRateLimit = '20/1min'` drip-feeds submissions instead of firing hundreds
of `sbatch` calls in one burst, which some schedulers rate-limit or drop.

The launcher fails fast on the two ways this can be set up wrongly: `-p slurm` from inside
the container (reported before any file validation, since it is an architectural blocker
rather than a fixable input), and `-p slurm` with no `sbatch` on PATH.

`job_script_example.sh` leads with this mode — a deliberately tiny 2-core, long-wall-time
job that runs Nextflow and waits on what it submits — and documents the
single-allocation alternative at the bottom for sites with no Nextflow module.

---

## Remaining work

1. **Push the image to Docker Hub** as `chutter/flumina`, which every README instruction
   depends on:
   ```bash
   docker push chutter/flumina:2.0.0 && docker push chutter/flumina:latest
   ```
   Built and verified locally; not pushed.
2. **Test under Apptainer on the cluster** — the one distribution path that could not be
   exercised here, and the one `job_script_example.sh` assumes.
3. **Coordinate translation for low-frequency FluMut screening** (see the caveat above).
4. WFABC and SNPGenie modules remain skeletons.
