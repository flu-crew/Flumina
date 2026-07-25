# Flumina: a pipeline for calling high- and low-frequency variants from Illumina influenza A data
[![DOI](https://zenodo.org/badge/794667242.svg)](https://zenodo.org/doi/10.5281/zenodo.12708893)

Welcome to the GitHub repository and manual for the software Flumina! Flumina is an easy-to-use pipeline for Illumina influenza A sequence data from any host and all gene segments that uses existing tools to trim and organize reads, assemble consensus genomes, call both high-frequency and low-frequency variants, link those variants to curated amino acids of interest, and screen for H5N1 molecular markers.

The pipeline accomplishes the following:

1) Organizes raw read data
2) Removes adaptor contamination and trims low quality reads and bases
3) Assembles cleaned reads into consensus contigs with IRMA
4) Calls high-frequency variants using GATK4
5) Calls low-frequency variants using LoFreq
6) Summarizes variant calling data
7) Links variants to curated amino acids of interest
8) Screens consensus genomes, and optionally low-frequency variants, for H5N1 markers using FluMut

Flumina uses Nextflow for program execution and cluster job submission management. With multi-job cluster management, analyzing many samples in bulk can be accomplished quite rapidly. The speed will increase as more threads are available, as individual pieces of the pipeline can be run in tandem.

A companion pipeline for Oxford Nanopore data, [FluPore](https://github.com/flu-crew/FluPore), is run the same way and shares this manual's structure.

# How to Run Flumina

Flumina has dependencies. For your convenience we have built a docker image which provides an isolated filesystem containing all required dependencies, including Nextflow itself. We recommend that you take advantage of this by either using Docker directly or using Docker indirectly via Apptainer (most recommended; formerly Singularity).

## 1) Download the Test Data Set

A four-sample H5N1 test data set is included in this repository under `test_dataset/raw_reads`, along with the matching `example_file_rename.csv` and `example_metadata.csv`. We recommend running this first to confirm Flumina is working before moving to your own data. If you wish, you can skip this step and jump to step 2, "Install Dependencies and Run Flumina", and run Flumina on your data first.

## 2) Install Dependencies and Run Flumina

This section will guide you through how to acquire Flumina, install its dependencies, and do a basic run with either the test data set or your data set. See the next section to explore available arguments.

### Using Apptainer (recommended)

Without going into details, Apptainer is a more secure method of running containers on an HPC system and you can learn more about it here: https://apptainer.org/docs/user/main/index.html
In this case we will be using Apptainer to run a Docker container. Docker is a container system that allows you to run Flumina in an environment we created that already contains all the dependencies you need. For more information about Docker check out their website: https://www.docker.com/
You may need to download Apptainer to use it. You will find installation instructions here: https://github.com/apptainer/apptainer/blob/main/INSTALL.md
If you are working on an HPC cluster you may be able to simply load an Apptainer module.

Once you have either installed or loaded Apptainer, download the docker image via Apptainer with a single command. This will download the file "flumina_latest.sif" in the directory where the command was run. Make sure you have room on your machine for this file before downloading it. Here is the command:

```
apptainer pull -F docker://chutter/flumina
```

After the image has been downloaded, create the folder and move your data into it. Here the input folder is named "raw_reads" to match Flumina's default expectation and should contain all the paired fastq.gz files you want to be part of this run.

```
mkdir raw_reads
mv your_data_*.fastq.gz raw_reads
```

You also need a renaming CSV and a metadata CSV, both described in the "Required input files" section below. With `file_rename.csv` and `metadata.csv` in your working directory, run Flumina:

```
apptainer run flumina_latest.sif -i raw_reads -o results
```

This runs Flumina with default arguments. Everything following "flumina_latest.sif" is fed directly to Flumina. A more realistic example is the following:

```
apptainer run flumina_latest.sif -i raw_reads -o fullFluminaRunOct19 -t 12 -M 32.GB
```

In this case "-t 12" tells Flumina to use 12 threads, "-M 32.GB" raises the per-step memory ceiling from its laptop-safe default, and "-o fullFluminaRunOct19" tells Flumina to put the results in a folder called "fullFluminaRunOct19" which will be created if one does not already exist.

### Using Docker

Docker is a container system which allows you to run Flumina in an environment we created that already contains all the dependencies you need. For more information about Docker check out their website: https://www.docker.com/
To use this method you will first need to download docker to the machine you are using to run Flumina. You will find their download options here: https://docs.docker.com/get-docker

If you are working on an HPC cluster you may be able to simply load a Docker module. On the other hand, some HPC systems believe Docker to be a security risk and therefore allow access to docker images only indirectly via Apptainer (formerly Singularity). If you are working on such an HPC system try the "Using Apptainer" method instead of this one.

We have uploaded the docker image to Docker Hub. In order to download it, go to the directory where your data is stored, then pull the image using the command:

```
docker pull chutter/flumina
```

Then create the input folder and move your reads into it. Unlike Apptainer, Docker only sees the folders you explicitly share with it using `-v`, so mount your working directory and run from there:

```
mkdir raw_reads
mv your_data_*.fastq.gz raw_reads

docker run -u $(id -u):$(id -g) -v "$(pwd)":/data -w /data chutter/flumina -i raw_reads -o results
```

The `-v` argument allows docker to read and write in your current directory, and `-u` makes the results belong to you rather than to root. A more realistic example is the following:

```
docker run -u $(id -u):$(id -g) -v "$(pwd)":/data -w /data chutter/flumina \
    -i raw_reads -o fullFluminaRunOct19 -t 12 -M 32.GB
```

### Installing all dependencies individually (not recommended)

While it is possible to install all of Flumina's dependencies without the docker image, it is not recommended, as the number of dependencies is extensive and there are sometimes conflicts between versions of software which can be hard to predict. If you choose to take this route, the outside programs can be installed through the anaconda environment file provided (version numbers are provided in the environment file for reporting and exact replication).

First, clone this repository:

```
git clone https://github.com/flu-crew/Flumina.git
cd Flumina
```

Install the Anaconda package manager from https://anaconda.org (Miniconda is recommended as it has a smaller footprint), then create the environment:

```
conda env create -f environment.yaml -n Flumina
conda activate Flumina
```

If the environment is to be installed in a specific location, like a project directory on a cluster:

```
conda env create -f environment.yaml -p /PATH/TO/Flumina
conda activate /PATH/TO/Flumina
```

You should only need to activate the environment if running locally. On a cluster, the conda activate line should be placed in your job script.

Flumina also needs Nextflow, which is not in the conda environment because the container ships it separately. Install it from https://www.nextflow.io, then add Flumina to your PATH so it can be accessed from any folder:

```
cd Flumina/Scripts
export PATH=$PATH:"$(pwd)"
```

A better but harder solution is to add the same line to the `.bash_profile` file in your home directory, which makes the change permanent rather than lasting only until you restart your command-line interface.

Now test that Flumina was successfully added to your path by going to a folder outside of the Scripts folder and typing:

```
flumina -h
```

You should see the help menu. After that, Flumina is run exactly as in the container examples above, except that you invoke it directly:

```
flumina -i raw_reads -o results -t 12
```

Note that when running outside the container Flumina defaults to the `docker` profile, meaning Nextflow will run each step inside the Flumina image. If you installed everything with conda instead and want the tools run directly, add `-p standard`.

# Required input files

### The renaming CSV

Often the case with multiplexed samples, you will find that the names of the reads are not the desired final names for the sample. To create the renaming file, a .csv file is needed with only two columns: "File" and "Sample". An example is included in this repository ("example_file_rename.csv").

The "File" column: the unique string that is part of the file name for the two read pairs, while excluding read and lane information. Example:

> ``AX1212_L001_R1.fastq.gz``

> ``AX1212_L001_R2.fastq.gz``

Are the two sets of reads for a given sample. Your "File" column value would then be:

> ``AX1212``

The "Sample" column: What you would like your sample name to be. This will be used in all downstream analyses unless changed. Ensure that your samples all have unique names and are not contained within each other (e.g. Name_0 is contained within Name_01). Also exclude special characters and replace spaces with underscores. Hyphens are also ok. In this example, the "Sample" column would be:

> ``Influenza_virus_AX1212``

Flumina looks for `file_rename.csv` in your working directory by default; point it elsewhere with `-n`.

### Sample metadata

A CSV of metadata to join with the amino acid data and summary data. This file must have at least a column titled "Sample" [capital S] to make the join possible. Any other column may be used to group the summaries with `-g`, for example a host column to compare cow versus bird versus poultry. Flumina looks for `metadata.csv` by default; point it elsewhere with `-m`.

### The reference

A reference sequence is needed to map the reads and compare amino acid changes to. This reference should have each gene as a separate entry in the fasta file, with the header including only the gene name. For now, multiple CDS reading frames should be included as separate fasta entries. There is no standard reference as the reference would depend on the research question. A default reference ships with Flumina and is used unless you supply your own with `-r`.

### Curated amino acid database

Normally Flumina will output databases of all the amino acid changes and then a reduced set of those that are nonsynonymous. To create a database of known amino acids of interest, these can be matched to the full amino acid database and separated into a more manageable table. The three columns this CSV must include are:

"Gene" - The gene that matches the names of the gene used in the reference

"Amino_Acid" - The amino acid position in the reference

"Type" - A summary of the function of the amino acid change

A default database ships with Flumina and is used unless you supply your own with `-a`.

# Flumina's Arguments (user options)

All of Flumina's arguments are laid out in the help menu. This menu can be accessed by typing `flumina -h`, `flumina --help`, or making any mistakes while calling Flumina. The help menu looks like so:

```
    ##########################################################################
                       Welcome to Flumina version 2.0.0!
    ##########################################################################
    Dependencies for this program are nextflow, java 17+, and either docker or
    apptainer. Running the Flumina container image satisfies all of them.
    For this help menu use the argument -h

    WARNING: ENSURE THAT YOU ARE NOT RUNNING THIS ON THE HEADNODE!

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
             at this minimum allele frequency (0-1). LoFreq variants at or
             above this frequency are applied to the IRMA consensus and the
             resulting sequences are screened. Off by default. e.g. '-l 0.01'
        -x : IRMA configuration file, e.g. to set TMP or SINGLE_LOCAL_PROC.
             Default is 'irma_config.sh' when one exists in the working
             directory
        -p : Execution profile: 'standard' (run tools locally), 'docker',
             'apptainer', or 'slurm'. Combine with a comma. Default is
             'standard' inside the Flumina container and 'docker' outside it.
             Use '-p slurm,apptainer' on a cluster to submit every step as its
             own job, running each inside the container — this is what makes
             large runs fast, as independent samples and steps run at once
        -Q : Scheduler queue/partition for -p slurm. Default is 'priority'
        -A : Scheduler account options for -p slurm.
             Default is '--qos=vpru -A nadc_iav'
        -j : Maximum scheduler jobs in flight at once with -p slurm.
             Default is 100. This is the main throughput dial for a large run
        -w : Nextflow work directory holding intermediate files. Default is
             './work'. On a cluster point this at scratch, for example
             -w /scratch/$USER/flumina
        -s : Skip a portion of Flumina. Options: '0' = skip nothing,
             'i' = skip IRMA assembly, 'm' = skip FluMut marker screening.
             These options may be combined. Default is 0
        -R : Resume a previous run, reusing all successfully completed work
        -N : Dry run. Print the Nextflow command that would be run, then exit

        --version : prints the version number
```

For many users the arguments `-i`, `-o`, and `-t` are the only important ones. The following example command specifies all three:

```
flumina -i inputSeqs -o testRunOct19 -t 16
```

Flumina can be operated in a modular fashion using the `-s` argument, which allows the user to skip one or more of Flumina's steps. This is useful if you are not interested in a given analysis, or if a run failed partway through. The following example command skips IRMA assembly and FluMut screening:

```
flumina -s im
```

The `-l` argument turns on low-frequency marker screening. Normally FluMut screens only the IRMA consensus, so a marker present in a minority of the viral population is invisible. With `-l`, LoFreq variants at or above the given allele frequency are applied to the consensus and the resulting sequences are screened as well, which surfaces markers segregating within a host. Results go to a separate `flumut_lowfreq` folder so they are never confused with consensus-level calls:

```
flumina -l 0.05
```

If a run fails partway through, fix the cause and re-run with `-R` to resume rather than starting over:

```
flumina -i raw_reads -o results -R
```

# Optional configuration files

### config.cfg

Every argument above can also be set in a configuration file, which is useful when you have more than a handful of non-default settings or want a record of exactly how a run was configured. An example is included in this repository ("config.cfg"). Values are given one per line in KEY=VALUE form:

```bash
OUTPUT_DIRECTORY="Flumina_results"
READ_DIRECTORY="raw_reads"
RENAME_FILE="file_rename.csv"
METADATA="metadata.csv"
GROUP_NAMES="discrete_host"
MIN_DEPTH=100
THREADS=4
```

Pass it with `-c`:

```
flumina -c config.cfg
```

A file named `config.cfg` in your working directory is picked up automatically, so `flumina` on its own is enough once one exists. **Command-line arguments always override the configuration file**, so a saved configuration can be reused with one-off changes:

```
flumina -c config.cfg -t 24 -o differentOutputFolder
```

### irma_config.sh

`irma_config.sh` is a configuration file for the program IRMA, and any of the standard IRMA parameters can be changed here. The three essential parameters are provided as an example:

```bash
TMP=/path/to/irma_tmp
SINGLE_LOCAL_PROC=2
DOUBLE_LOCAL_PROC=1
```

TMP is the temporary directory that IRMA writes files to; it will otherwise use a sometimes unwritable or slow location, so setting this somewhere you can write to and that has enough space helps. Note that SINGLE_LOCAL_PROC should match the number of threads available to each IRMA job.

A file named `irma_config.sh` in your working directory is picked up automatically; point Flumina at a different one with `-x`.

# Running Flumina on a cluster

This is where Flumina is fastest. With `-p slurm,apptainer` every pipeline step is submitted as its own scheduler job, so independent samples and independent steps run at the same time rather than one after another — with 100 samples the per-sample chain runs about 100-wide. Each of those jobs runs inside the Flumina container, so you still get exactly one pinned software environment.

```
flumina -i raw_reads -o results -p slurm,apptainer -w /scratch/$USER/flumina \
        -Q priority -A "--qos=vpru -A nadc_iav" -j 100 -t 8 -M 32.GB
```

An example SLURM job script is included in this repository ("job_script_example.sh"). Edit the account, partition, and email lines at the top, then submit it:

```
sbatch job_script_example.sh
```

Four things matter more than the rest:

- **Nextflow must run on the host, not inside the container.** It submits work with `sbatch`, and that binary is not in the image. Most clusters provide `module load nextflow`; Flumina will tell you if it is missing. The container is still used for every step — it just is not what wraps Nextflow.
- **The job that runs Nextflow should be tiny.** Two cores is plenty, since it spends its life waiting on the jobs it submits. Give it a long wall time though, as it must outlive them all.
- **Point the work directory at scratch with `-w`,** never at home or an NFS share. Every submitted job reads and writes there, so it must be visible from all compute nodes, and a slow filesystem here is the most common cause of a slow run.
- **`-Q` and `-A` are almost certainly wrong for you.** They default to the flu-crew queue and account and will be rejected on any other cluster. `-j` sets how many jobs may be in flight at once, which is the main throughput dial.

If your cluster has no Nextflow module, or you are only running a handful of samples, everything can instead run inside a single allocation with nothing installed but Apptainer:

```
apptainer run flumina_latest.sif -i raw_reads -o results -t 24 -M 100.GB
```

Steps then run in parallel only up to the cores of that one node, so this is considerably slower for large sample counts. Both approaches are laid out in full in the example job script.

# Output

Results are written to the output directory given by `-o`:

| Folder | Contents |
|---|---|
| `variant_analysis/` | Variant tables, amino acid changes, and summaries |
| `IRMA-consensus-contigs/` | Per-sample consensus genomes |
| `IRMA_results/` | Full IRMA output per sample |
| `vcf_files/` | Per-sample GATK and LoFreq VCFs |
| `BAM_files/` | Aligned, sorted, duplicate-marked BAMs |
| `processed-reads/` | Trimmed reads |
| `flumut/` | H5N1 markers found in the consensus genomes |
| `flumut_lowfreq/` | H5N1 markers found in low-frequency variants (with `-l`) |
| `logs/` | Per-sample tool logs |
| `pipeline_info/` | Run timeline, resource report, trace, and the exact config used |

# Upcoming features

1) WFABC population genetics module
2) SNPGenie dN/dS analysis
3) Code refinement and speedups
