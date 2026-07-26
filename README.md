# Flumina
[![DOI](https://zenodo.org/badge/794667242.svg)](https://zenodo.org/doi/10.5281/zenodo.12708893)

A pipeline for processing and calling high-frequency and low-frequency variants from Illumina sequence data for Influenza viruses

The pipeline accomplishes the following:

1) Organize raw read data
2) Remove adaptor contamination and trims low quality reads and bases
3) Assemble cleaned reads into consensus contigs using IRMA
4) Calls high-frequency variants using GATK4
5) Calls low-frequency variants using LoFreq
6) Summarizes variant calling data
7) Link variants to curated amino acids of interest
8) Screen consensus genomes, and optionally low-frequency variants, for H5N1 markers using FluMut
9) Optional dN/dS with SNPGenie and per-site selection coefficients with WFABC

Flumina uses Nextflow for program execution and cluster job submission management. With multi-job cluster management, analyzing many samples in bulk can be accomplished quite rapidly. The speed will increase as more threads are available, as individual pieces of the pipeline can be run in tandem.

A companion pipeline for Oxford Nanopore data, [FluPore](https://github.com/flu-crew/FluPore), is run the same way.

# Running Flumina

There are three ways to run Flumina. The container includes everything, so unless you have a reason to manage the dependencies yourself, use one of the container options.

| | Use when |
|---|---|
| [1. Conda](#1-conda) | You want the tools installed directly, or cannot use containers |
| [2. Container, single job](#2-container-single-job) | Simplest. Everything runs in one allocation |
| [3. Container, multi-job](#3-container-multi-job) | Fastest for many samples. Each step becomes its own cluster job |

All three take the same arguments.

## 1) Conda

Clone the repository and create the environment from the included `environment.yaml`:

```bash
git clone https://github.com/flu-crew/Flumina.git
cd Flumina
conda env create -f environment.yaml -n Flumina
conda activate Flumina
```

Flumina also needs Nextflow, which is not in the conda environment. Install it from https://www.nextflow.io, or load your cluster's module. Then add Flumina to your PATH and run it:

```bash
export PATH=$PATH:"$(pwd)/Scripts"
flumina -i raw_reads -o results
```

## 2) Container, single job

Download the image and run it. Nothing else needs installing, and nothing needs to be cloned:

```bash
apptainer pull -F docker://chutter/flumina
apptainer run flumina_latest.sif -i raw_reads -o results
```

Or with a configuration file, which is picked up automatically if it is named `config.cfg` and sits in your working directory:

```bash
apptainer run flumina_latest.sif -c config.cfg
```

Docker works the same way, but only sees folders you share with it:

```bash
docker run -u $(id -u):$(id -g) -v "$(pwd)":/data -w /data chutter/flumina -i raw_reads -o results
```

Every step runs inside your one allocation, so give the job real resources and match `-t` and `-M` to them.

Apptainer mounts your home and working directory automatically. If your reads are somewhere else, bind that location:

```bash
export APPTAINER_BIND=/project
```

## 3) Container, multi-job

Each step is submitted as its own cluster job, so independent samples and steps run at the same time. Each of those jobs still runs inside the container.

Nextflow runs on the host here, so it needs the pipeline files. Copy them out of the image once — nothing is downloaded:

```bash
module load nextflow apptainer
apptainer pull -F docker://chutter/flumina
apptainer run flumina_latest.sif --export ./flumina
```

Then run it, from a small job with a long wall time:

```bash
./flumina/Scripts/flumina -c config.cfg -p slurm,apptainer
```

Set `-Q` and `-A` to your queue and account. Besides `slurm`, the profiles `pbs`, `pbspro`, `sge` and `lsf` are available. Example job scripts for both container methods are included: `job_script_example_config.sh` and `job_script_example_arguments.sh`.

Two things to note. `-t` and `-M` now apply **per step**, not to the whole run, so 8 threads is usually plenty. And Nextflow works out the container binds itself, so `APPTAINER_BIND` is not needed.

# Arguments

All of Flumina's arguments are laid out in the help menu, accessed with `flumina -h`. Only `-i` and `-o` are needed for a basic run; everything else has a default.

```
    Arguments:
        -i : Input directory of raw paired fastq.gz reads. Default is
             'raw_reads'
        -o : Output directory. Default is 'Flumina_results'
        -t : Number of threads. Default is 4
        -M : Maximum memory any single step may use, e.g. '32.GB'. Default is
             6.GB, which fits a laptop or a stock Docker Desktop VM. Raise it
             on a workstation or cluster. The slurm profile raises it for you
        -c : Configuration file holding any of these parameters in KEY=VALUE
             form (see config.cfg). Command-line arguments override it.
             Default is 'config.cfg' when one exists in the working directory
        -n : CSV with the file name that matches both read pairs in the "File"
             column and the sample name in the "Sample" column. Default is
             'file_rename.csv'
        -r : Reference FASTA with each gene segment as a separate entry.
             Default is the bundled 'reference.fa'
        -a : Curated amino acid database CSV with the columns "Gene",
             "Amino_Acid", and "Type". Default is the bundled
             'curated_database.csv'
        -m : Metadata CSV with at least a "Sample" column. Default is
             'metadata.csv'
        -g : Metadata column name used to group the summaries, e.g. cow versus
             bird versus poultry. Use 'NULL' for no grouping. Default is
             'discrete_host'
        -d : Minimum read depth to keep a variant. Default is 100.
             WARNING: below 100 LoFreq and GATK report false fixations from low
             template input (the founder/jackpot effect), so raise this rather
             than lower it unless you know why you are doing so.
        -q : Minimum quality to keep a variant, 0 disables. Default is 30
        -f : Minimum allele frequency to keep a variant (0-1).
             Default is 0.01 (1%)
        -l : Also screen LOW-FREQUENCY variants for H5N1 markers with FluMut,
             at this minimum allele frequency (0-1). Off by default
        -G : Run SNPGenie dN/dS analysis. Off by default
        -W : Run WFABC selection analysis. Off by default. Needs metadata with
             an individual column and a time-point column
        -x : IRMA configuration file, e.g. to set TMP or SINGLE_LOCAL_PROC.
             Default is 'irma_config.sh' when one exists in the working
             directory
        -p : Execution profile. Container: 'standard', 'docker', 'apptainer'.
             Scheduler: 'slurm', 'pbs', 'pbspro', 'sge', 'lsf'. Combine with a
             comma, e.g. '-p slurm,apptainer'
        -Q : Scheduler queue/partition
        -A : Scheduler account options, passed through as given
        -j : Maximum scheduler jobs in flight at once. Default is 100
        -w : Nextflow work directory. Default is './work'. On a cluster point
             this at scratch
        -s : Skip a portion of Flumina. Options: '0' = skip nothing,
             'i' = skip IRMA assembly, 'm' = skip FluMut marker screening.
             These options may be combined. Default is 0
        -R : Resume a previous run, reusing all successfully completed work
        -N : Dry run. Print the Nextflow command that would be run, then exit

        --export [dir] : copy the pipeline out of the container so Nextflow can
             run on the host, needed for '-p slurm'. Defaults to './flumina'

        --version : prints the version number
```

If a run fails partway through, fix the cause and re-run with `-R` to resume rather than starting over. This only works if the work directory is still there, so do not delete it until you are happy with the results.

# Disk use

The work directory holds every intermediate file and ends up roughly the size of the results again, so a run costs about twice what you keep. Three things help:

`-C` deletes the work directory once a run completes successfully. Use it for production runs; leave it off while you are still getting a run to work, because it discards the resume cache too.

`PUBLISH_MODE="link"` in the config publishes results as hard links instead of copies, so results and work share the same data and a run takes about half the space. It only works when the output and work directories are on the same filesystem. Do not use `symlink` — deleting the work directory would leave broken links.

Old runs can be cleared at any time with Nextflow's own command, which understands what is still needed:

```
nextflow clean -f -before <run-name>
```

`nextflow log` lists the run names. Add `-n` to any of these to see what would be removed without removing it.

# Input files

## Create renaming file

Often the case with multiplexed samples in sequence capture projects, you will find that the names of the reads often are not the desired final names for the sample. To create the renaming file, a .csv file is needed with only two columns: "File" and "Sample". An example is included in the main repo ("example_file_rename.csv").

The "File" column: the unique string that is part of the file name for the two read pairs, while excluding read and lane information. Example:

> ``AX1212_L001_R1.fastq.gz``

> ``AX1212_L001_R2.fastq.gz``

Are the two sets of reads for a given sample. Your "File" column value would then be:

> ``AX1212``

The "Sample" column: What you would like your sample name to be. This will be used up in all downstream analyses unless changed. Ensure that your samples all have unique names and are not contained within each other (e.g. Name_0 is contained within Name_01). Also exclude special characters and replace spaces with underscores. Hyphens are also ok.
In this example, the "Sample" Column would be:
>
>Influenza_virus_AX1212
>

Reads are found recursively, so nested per-sample directories are fine, and both compressed and uncompressed fastq are accepted. A sample listed in the CSV with no reads does not stop the run: it is reported and recorded in `logs/missing_samples.log`, and the rest carry on.

## Setting up a reference

A reference sequence is needed to map the reads and compare amino acid changes to. This reference should have each gene as a separate entry in the fasta file, with the header including only the gene name. For now, multiple CDS reading frames should be included as separate fasta entries. There is no standard reference as the reference would depend on the research question. One is bundled with Flumina and used unless you supply your own with `-r`.

## Setting up configuration file

Flumina uses a configuration file to keep track of the parameters and easily add new ones. An example is included in the main repo, "config.cfg". Every setting has a command-line equivalent, and the command line wins, so a saved configuration can be reused with one-off changes:

```bash
flumina -c config.cfg -t 24 -o differentOutputFolder
```

Nothing in the file is required. A file named `config.cfg` in your working directory is picked up automatically.

## Optional metadata and other files

### sample metadata

A CSV of metadata to join with the amino acid data and summary data can be provided. This CSV file must have at least a column titled "Sample" [capital S] to make the join possible. Without it the summaries are simply not grouped. It is required only for WFABC, which needs an individual column and a time-point column to build allele-frequency time series.

### curated amino acid databases

Normally Flumina will output databases of all the amino acid changes and then a reduced set to those that are nonsynonymous. To create a database of known amino acids of interest, these can be matched to the full amino acid database and separated into a more manageable table. The three columns this CSV must include are

"Gene" - The Gene that matches the names of the gene used in the reference

"Amino_Acid" - The amino acid position in the reference

"Type" - A summary of the function of the amino acid change

Without one, the curated-site summary is skipped and everything else is unaffected.

### irma_config.sh

irma_config.sh, is a configuration file with an example provided here for the program IRMA. Any of the standard IRMA parameters can be changed here and included in your working directory. The 3 essential parameters are provided as an example:

```bash
TMP=./irma_tmp
SINGLE_LOCAL_PROC=2
DOUBLE_LOCAL_PROC=1
```

TMP is the temporary directory that IRMA writes files, and it will use a sometimes unwritable or slow location, so setting this somewhere you can write to and has enough space helps.

Note that SINGLE_LOCAL_PROC should match the number of threads in your job or if using multi-job cluster submission it should match "--cpus-per-task"

# Output

Results are written to the output directory given by `-o`:

| Folder | Contents |
|---|---|
| `variant_analysis/` | Variant tables, amino acid changes, and summaries |
| `variant_analysis/flumut/` | H5N1 markers found in the consensus genomes |
| `variant_analysis/flumut_lowfreq/` | H5N1 markers found in low-frequency variants (`-l`) |
| `IRMA-consensus-contigs/` | Per-sample consensus genomes |
| `IRMA_results/` | Full IRMA output per sample |
| `vcf_files/` | Per-sample GATK and LoFreq VCFs |
| `BAM_files/` | Aligned, sorted, duplicate-marked BAMs |
| `processed-reads/` | Trimmed reads |
| `snpGenie_results/` | Per-site dN/dS estimates (`-G`) |
| `wfabc_analysis/` | Selection coefficients and Ne estimates (`-W`) |
| `logs/` | Per-sample tool logs, and `missing_samples.log` if any sample had no reads |
| `pipeline_info/` | Run timeline, resource report, trace, and the exact config used |
