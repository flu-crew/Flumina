#!/usr/bin/env nextflow
/*
 * Flumina 2.0 — Nextflow port of the Snakemake pipeline.
 *
 * Output layout is deliberately identical to the Snakemake version
 * (BAM_files/, vcf_files/, logs/, IRMA_results/, processed-reads/) so the
 * downstream R scripts run unmodified and the two pipelines can be diffed
 * output-for-output on the same test dataset.
 *
 *   nextflow run . -profile test,docker
 *   nextflow run . -profile docker --read_directory reads --rename_file r.csv \
 *                  --reference reference.fa --metadata meta.csv --aa_db db.csv
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    ##########################################################################
                       Welcome to Flumina version ${workflow.manifest.version}!
    ##########################################################################

    Most users should launch Flumina through the `flumina` command rather than
    calling this workflow directly, which gives short arguments and defaults:

        flumina -i raw_reads -o results
        flumina -h

    Calling the workflow directly takes the same settings as --parameters:

    Required:
      --read_directory   Directory of raw paired fastq.gz
      --rename_file      CSV with File,Sample columns
      --reference        Reference FASTA

    Optional inputs:
      --metadata         Metadata CSV containing a Sample column. Without it the
                         summaries are simply not grouped
      --aa_db            Curated amino-acid database CSV. Without it the
                         curated-site join is skipped

    Common:
      --outdir           Output directory              [${params.outdir}]
      --max_cpus         Max CPUs used at once         [${params.max_cpus}]
      --min_depth        Minimum depth to keep a call  [${params.min_depth}]
      --min_quality      Minimum quality to keep       [${params.min_quality}]
      --min_allele_frequency  Minimum allele frequency [${params.min_allele_frequency}]
      --group_names      Metadata column to group by   [${params.group_names}]
      --run_irma         Run IRMA assembly             [${params.run_irma}]
      --irma_config      IRMA parameter file           [${params.irma_config ?: 'none'}]
      --flumut           Screen consensus for markers  [${params.flumut}]
      --flumut_lowfreq   Screen low-freq variants      [${params.flumut_lowfreq}]
      --flumut_freq_threshold  Min AF for low-freq     [${params.flumut_freq_threshold}]
      --snpgenie         Run SNPGenie dN/dS            [${params.snpgenie}]
      --wfabc            Run WFABC selection analysis  [${params.wfabc}]
      -profile           standard | docker | apptainer | slurm | test
    """.stripIndent()
}

/* ==========================================================================
 * Reference preparation
 * --------------------------------------------------------------------------
 * The Snakemake version used five separate rules and a whole second snakefile
 * invocation for this. All four indexing steps are seconds-long and always run
 * together, so splitting them only adds scheduling overhead.
 * ========================================================================== */
process PREPARE_REFERENCE {
    tag   "reference"
    label 'process_low'
    publishDir "${params.outdir}/Reference", mode: params.publish_mode
    // findAAChanges.R reads ${OUTPUT_DIRECTORY}/reference.fa — the old bash
    // driver copied the reference to the output root before calling snakemake,
    // so the R scripts depend on it being there. Reproduce that placement.
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'reference.fa'

    input:
    // stageAs gives the input a fixed, distinct name. Without it, a reference
    // already called reference.fa stages onto the output name and the copy
    // fails with "are the same file".
    path reference, stageAs: 'input_reference'

    output:
    tuple path('reference.fa'), path('reference.fa.*'), path('reference.dict'), emit: index

    script:
    """
    cp -L ${reference} reference.fa
    bwa index -a bwtsw reference.fa
    samtools faidx reference.fa
    gatk CreateSequenceDictionary --REFERENCE reference.fa --OUTPUT reference.dict \\
        --USE_JDK_DEFLATER true --USE_JDK_INFLATER true
    # GATK looks for <base>.dict; some tools look for <file>.dict. Provide both.
    cp reference.dict reference.fa.dict
    """
}

/* ==========================================================================
 * Read processing
 * ========================================================================== */
process FASTP {
    tag   "$sample"
    label 'process_medium'
    publishDir { "${params.outdir}/processed-reads/${sample}" }, mode: params.publish_mode, pattern: '*.fastq.gz'
    publishDir "${params.outdir}/logs",                      mode: params.publish_mode, pattern: '*.{html,json}'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz"), emit: reads
    path "fastp-${sample}.{html,json}",                                              emit: reports

    script:
    """
    fastp --in1 ${r1} --in2 ${r2} \\
        --out1 ${sample}_R1.fastq.gz --out2 ${sample}_R2.fastq.gz \\
        --length_required 60 --low_complexity_filter --complexity_threshold 30 \\
        --trim_poly_x --correction --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --html fastp-${sample}.html --json fastp-${sample}.json --compression 8 \\
        --report_title ${sample}
    """
}

/*
 * IRMA reads an OPTIONAL `irma_config.sh` from its current working directory —
 * that is IRMA's own behaviour, not something added here (see the
 * `[ -r "./irma_config.sh" ]` block in the IRMA script itself). Staging the
 * user's file under that exact name into the task work dir is therefore all
 * that is needed to pass through TMP, SINGLE_LOCAL_PROC and any other IRMA
 * parameter, and it works identically under Docker and Apptainer because
 * nothing inside the read-only image has to be modified.
 *
 * When --irma_config is not set, an empty list is staged, which stages nothing,
 * and IRMA falls back to its module defaults.
 */
process IRMA {
    tag   "$sample"
    label 'process_high'
    publishDir "${params.outdir}/IRMA_results", mode: params.publish_mode

    input:
    tuple val(sample), path(r1), path(r2)
    path irma_cfg, stageAs: 'irma_config.sh'

    output:
    tuple val(sample), path("${sample}"), emit: results

    script:
    """
    IRMA FLU ${r1} ${r2} ${sample}
    """
}

/* ==========================================================================
 * Alignment — GATK best-practice style unmapped-BAM route
 * ========================================================================== */
process FASTQ_TO_SAM {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:  tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path('fastqsam.bam'), emit: bam

    script:
    """
    gatk FastqToSam -FASTQ ${r1} -FASTQ2 ${r2} \\
        -OUTPUT fastqsam.bam -SAMPLE_NAME ${sample} \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

process REVERT_SAM {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:  tuple val(sample), path(bam)
    output: tuple val(sample), path('revertsam.bam'), emit: bam

    script:
    """
    gatk RevertSam -I ${bam} -O revertsam.bam \\
        -SANITIZE true -MAX_DISCARD_FRACTION 0.005 \\
        -ATTRIBUTE_TO_CLEAR XT -ATTRIBUTE_TO_CLEAR XN -ATTRIBUTE_TO_CLEAR AS \\
        -ATTRIBUTE_TO_CLEAR OP -SORT_ORDER queryname \\
        -RESTORE_ORIGINAL_QUALITIES true -REMOVE_DUPLICATE_INFORMATION true \\
        -REMOVE_ALIGNMENT_INFORMATION true \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

process ADD_READ_GROUPS {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:  tuple val(sample), path(bam)
    output: tuple val(sample), path('all_reads.bam'), emit: bam

    script:
    """
    gatk AddOrReplaceReadGroups -I ${bam} -O all_reads.bam \\
        -RGSM ${sample} -RGPU FLOWCELL1.LANE1 -RGID FLOWCELL1.LANE1 \\
        -RGLB LIB-${sample} -RGPL ILLUMINA \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

/*
 * bwa mem is invoked WITHOUT -t, deliberately.
 *
 * bwa estimates the insert-size distribution per read batch, and batch
 * composition depends on the thread count — so -t shifts pairing decisions for
 * a handful of marginal reads. That is enough to change a few borderline GATK
 * calls at segment edges and stop the BAM being bit-identical to the Snakemake
 * baseline, which the original ran single-threaded.
 *
 * If throughput ever matters more than matching that baseline, use
 * '-t N -K 100000000': fixing the chunk size makes threaded output
 * deterministic, though still not identical to -t 1.
 */
process BWA_MAP {
    tag "$sample"; label 'process_high'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('mapped_reads_all.bam'), emit: bam

    script:
    """
    gatk SamToFastq -I ${bam} -FASTQ /dev/stdout \\
        -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true \\
      | bwa mem -M -p ${ref} /dev/stdin \\
      | gatk MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM ${bam} \\
        -OUTPUT mapped_reads_all.bam -R ${ref} -CREATE_INDEX true -ADD_MATE_CIGAR true \\
        -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true \\
        -MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant \\
        -ATTRIBUTES_TO_RETAIN XS \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

process SORT_BAM {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:  tuple val(sample), path(bam)
    output: tuple val(sample), path('mapped_reads_sort.bam'), emit: bam

    script:
    """
    gatk SortSam -INPUT ${bam} -OUTPUT mapped_reads_sort.bam \\
        -CREATE_INDEX true -SORT_ORDER coordinate \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

process MARK_DUPLICATES {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, pattern: '*.bam'
    publishDir { "${params.outdir}/logs/${sample}" },      mode: params.publish_mode, pattern: '*.txt'

    input:  tuple val(sample), path(bam)
    output:
    tuple val(sample), path('mapped_reads_md.bam'), emit: bam
    path 'duplicate_metrics.txt',                   emit: metrics

    script:
    """
    gatk MarkDuplicates -INPUT ${bam} -OUTPUT mapped_reads_md.bam \\
        -CREATE_INDEX true -METRICS_FILE duplicate_metrics.txt \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

process SET_TAGS {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('final_mapped_reads.bam'), path('final_mapped_reads.bai'), emit: bam

    script:
    """
    gatk SortSam -INPUT ${bam} -OUTPUT /dev/stdout -SORT_ORDER coordinate \\
      | gatk SetNmAndUqTags -INPUT /dev/stdin -OUTPUT final_mapped_reads.bam \\
        -CREATE_INDEX true -R ${ref} \\
        -USE_JDK_DEFLATER true -USE_JDK_INFLATER true
    """
}

/* ==========================================================================
 * Variant calling
 * ========================================================================== */
process HAPLOTYPE_CALLER {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, pattern: '*.bam'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode, pattern: '*.vcf'

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('gatk4-haplotype-caller.vcf'), emit: vcf
    path 'haplotype_caller.bam',                           emit: bam

    script:
    """
    gatk HaplotypeCaller -I ${bam} -R ${ref} -O gatk4-haplotype-caller.vcf \\
        -ERC GVCF -ploidy 1 -bamout haplotype_caller.bam
    """
}

process GENOTYPE_GVCF {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(vcf)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('gatk4-unfiltered-genotypes.vcf'), emit: vcf

    script:
    // NOTE: the Snakemake version passed --use-new-qual-calculator, which was
    // removed in GATK 4.1+ (it is the default behaviour now). Dropped here.
    """
    gatk GenotypeGVCFs -V ${vcf} -R ${ref} -O gatk4-unfiltered-genotypes.vcf
    """
}

process SELECT_VARIANTS {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode

    input:  tuple val(sample), path(vcf)
    output:
    tuple val(sample), path('gatk4-unfiltered-snps.vcf'),   emit: snps
    tuple val(sample), path('gatk4-unfiltered-indels.vcf'), emit: indels

    script:
    """
    gatk SelectVariants -V ${vcf} -O gatk4-unfiltered-snps.vcf   --select-type SNP
    gatk SelectVariants -V ${vcf} -O gatk4-unfiltered-indels.vcf --select-type INDEL
    """
}

process FILTER_VARIANTS {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(snps), path(indels)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('gatk4-filtered-snps.vcf'), path('gatk4-filtered-indels.vcf'), emit: vcf

    script:
    // The Snakemake filter_INDEL rule existed but never ran: its output was not
    // listed in `rule all`, so Snakemake never requested it. It runs here.
    """
    gatk VariantFiltration -R ${ref} -V ${snps} -O gatk4-filtered-snps.vcf \\
        -filter "QUAL<30.0"            --filter-name "QUAL" \\
        -filter "QD<2.0"               --filter-name "QD" \\
        -filter "SOR<3.0"              --filter-name "SOR" \\
        -filter "FS<60.0"              --filter-name "FS" \\
        -filter "MQ<40.0"              --filter-name "MQ" \\
        -filter "MQRankSum<-12.5"      --filter-name "MQRankSum" \\
        -filter "ReadPosRankSum<-8.0"  --filter-name "ReadPosRankSum"

    gatk VariantFiltration -R ${ref} -V ${indels} -O gatk4-filtered-indels.vcf \\
        -filter "QD<2.0"               --filter-name "QD" \\
        -filter "QUAL<30.0"            --filter-name "QUAL" \\
        -filter "FS<60.0"              --filter-name "FS" \\
        -filter "ReadPosRankSum<-8.0"  --filter-name "ReadPosRankSum"
    """
}

process LOFREQ {
    tag "$sample"; label 'process_medium'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('lofreq-called-variants.vcf'), emit: vcf

    script:
    """
    lofreq call -f ${ref} -o lofreq-called-variants.vcf ${bam}
    """
}

/* ==========================================================================
 * Downstream R analysis
 * --------------------------------------------------------------------------
 * The R scripts read a KEY=VALUE config.cfg. Rather than rewrite them to take
 * CLI arguments, the pipeline writes a config from its own params — the scripts
 * are used byte-identical to the Snakemake version and stay runnable by hand.
 * ========================================================================== */
/*
 * Gather one sample's VCFs into a directory named after the sample.
 *
 * convertVCFtoTable.R derives the sample name from the directory component of
 * the path it finds (gsub("/.*", "", ...)), so the vcf_files/<sample>/<file>.vcf
 * layout is load-bearing. Nextflow flattens staged filenames, and every sample
 * produces an identically-named lofreq-called-variants.vcf, so they must be
 * separated into per-sample directories before staging.
 */
process GATHER_SAMPLE_VCFS {
    tag "$sample"
    label 'process_low'

    input:  tuple val(sample), path(vcfs, stageAs: 'in/*')
    output: path "${sample}", emit: dir

    script:
    """
    mkdir -p '${sample}'
    cp -L in/* '${sample}'/
    """
}

/*
 * Downstream R analysis — fully relocatable.
 *
 * Every input is staged INTO this task's work directory and the generated
 * config.cfg sets OUTPUT_DIRECTORY=".", so the R scripts resolve everything
 * relative to wherever they happen to run. Nothing references a host path,
 * which is what lets this execute unchanged on a laptop, an HPC node, or an
 * AWS Batch container that has never seen the submitting filesystem.
 *
 * The R scripts themselves are unmodified: they already build every path as
 * paste0(OUTPUT_DIRECTORY, "/..."), so a relative root just works.
 */
process R_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'variant_analysis'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'IRMA-consensus-contigs'
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode, pattern: 'config.cfg'

    input:
    path scripts
    path vcf_dirs,  stageAs: 'vcf_files/*'
    path reference
    path aa_db
    path metadata
    path irma_dirs, stageAs: 'IRMA_results/*'

    output:
    path 'variant_analysis',        emit: results
    path 'config.cfg',              emit: config
    path 'IRMA-consensus-contigs',  emit: irma, optional: true

    script:
    def irma_step = params.run_irma
        ? "Rscript ${scripts}/organizeIRMA.R config.cfg"
        : "echo 'IRMA disabled, skipping organizeIRMA.R'"
    // Both of these are optional. When not supplied nothing is staged, so the
    // input variable is an empty list and would render as an empty string —
    // write the literal NULL the R scripts test for instead.
    def aa_db_cfg    = params.aa_db    ? "${aa_db}"    : 'NULL'
    def metadata_cfg = params.metadata ? "${metadata}" : 'NULL'
    // Grouping is a column OF the metadata, so without metadata there is
    // nothing to group by and outputSummary.R would fail looking for it.
    def group_cfg    = params.metadata ? "${params.group_names}" : 'NULL'
    """
    # Paths are relative to this work dir — never the submitting host.
    cat > config.cfg <<'CFG_END'
OUTPUT_DIRECTORY="."
REFERENCE_FILE="${reference}"
AA_DB="${aa_db_cfg}"
METADATA="${metadata_cfg}"
GROUP_NAMES="${group_cfg}"
MIN_DEPTH="${params.min_depth}"
MIN_QUALITY="${params.min_quality}"
MIN_ALLELE_FREQUENCY="${params.min_allele_frequency}"
# runSNPGenie.R reads these two under different names than the rest of the
# pipeline (MIN_ALLELE_FREQ / MIN_COVERAGE, not MIN_ALLELE_FREQUENCY /
# MIN_DEPTH). Writing both spellings is what makes SNPGenie actually honour the
# same depth and frequency thresholds as every other step; without them it
# silently applies no filtering at all.
MIN_ALLELE_FREQ="${params.min_allele_frequency}"
MIN_COVERAGE="${params.min_depth}"
DEDUP_KEYS="${params.dedup_keys}"
DISABLE_IRMA="${params.run_irma ? 'FALSE' : 'TRUE'}"
SNPGENIE="${params.snpgenie ? 'TRUE' : 'FALSE'}"
INDIVIDUAL_COLUMN="${params.individual_column}"
TIME_COLUMN="${params.time_column}"
GENERATIONS_PER_TIME="${params.generations_per_time}"
FIXATION_CUTOFF="${params.fixation_cutoff}"
WFABC_TIMEOUT="${params.wfabc_timeout}"
MAX_REFINE_ITER="${params.max_refine_iter}"
OVERWRITE="FALSE"
THREADS="${task.cpus}"
CFG_END

    ${irma_step}
    Rscript ${scripts}/convertVCFtoTable.R config.cfg
    Rscript ${scripts}/findAAChanges.R     config.cfg
    Rscript ${scripts}/outputSummary.R     config.cfg
    """
}

/*
 * SNPGenie — per-site dN/dS from the pooled LoFreq calls (Nelson & Hughes 2015).
 *
 * Two scripts in sequence: makeGTF.R turns the reference into the per-segment
 * GTF+FASTA pairs SNPGenie needs, then runSNPGenie.R runs SNPGenie itself over
 * every sample's VCF and collects the per-sample output into combined tables.
 *
 * This ran in the old Snakemake driver and was lost in the Nextflow port — the
 * `snpgenie` parameter survived but nothing acted on it, so setting it to true
 * silently did nothing at all.
 *
 * The config.cfg emitted by R_ANALYSIS is reused rather than regenerated, so
 * these steps cannot drift out of step with the main analysis. It is copied
 * before being appended to because the staged original is a symlink into
 * another task's directory.
 */
process SNPGENIE {
    label 'process_medium'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'snpGenie_results'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'reference_gtf'

    input:
    path scripts
    path config
    path vcf_dirs, stageAs: 'vcf_files/*'
    path reference
    path metadata

    output:
    path 'snpGenie_results', emit: results, optional: true
    path 'reference_gtf',    emit: gtf,     optional: true

    script:
    """
    cp ${config} run_config.cfg
    echo 'THREADS="${task.cpus}"' >> run_config.cfg
    # runSNPGenie.R setwd()s into each per-sample output directory and then
    # keeps using paths built from OUTPUT_DIRECTORY, which only works when that
    # is absolute. The pipeline's relocatable OUTPUT_DIRECTORY="." therefore
    # breaks it the moment it changes directory. Appending wins because the R
    # config parsers keep the last value seen for a key, and \$PWD is resolved
    # at run time inside this task's directory, so nothing is hardcoded.
    echo "OUTPUT_DIRECTORY=\\"\$PWD\\"" >> run_config.cfg

    Rscript ${scripts}/makeGTF.R      run_config.cfg
    Rscript ${scripts}/runSNPGenie.R  run_config.cfg
    """
}

/*
 * WFABC — per-site selection coefficients and Ne from allele-frequency time
 * series (Foll et al. 2015).
 *
 * Unlike every other step this one is not per-sample: it needs the SAME
 * individual sampled at two or more time points, which it reconstructs by
 * joining the variant table to METADATA on INDIVIDUAL_COLUMN and TIME_COLUMN.
 * A metadata file without usable values in those columns yields no usable time
 * series, so runWFABC.R stops with an explanatory error rather than emitting
 * an empty result — hence the optional outputs.
 *
 * wfabc_1 and wfabc_2 are found on PATH; the container builds them from source
 * (see the Dockerfile). Outside the container, set WFABC_PATH in config.cfg.
 */
process WFABC {
    label 'process_medium'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'wfabc_analysis'

    input:
    path scripts
    path config
    path vcf_dirs, stageAs: 'vcf_files/*'
    path reference
    path metadata

    output:
    path 'wfabc_analysis', emit: results, optional: true

    script:
    """
    cp ${config} run_config.cfg
    echo 'THREADS="${task.cpus}"' >> run_config.cfg
    # Same reason as SNPGENIE above: runWFABC.R setwd()s into a per-site
    # directory, so OUTPUT_DIRECTORY has to be absolute.
    echo "OUTPUT_DIRECTORY=\\"\$PWD\\"" >> run_config.cfg

    Rscript ${scripts}/runWFABC.R run_config.cfg
    """
}

/*
 * FluMut — screens consensus genomes against FluMutDB for H5N1 molecular
 * markers of host adaptation, virulence, and antiviral resistance
 * (Giussani et al. 2025, Virus Evolution: doi 10.1093/ve/veaf011).
 *
 * FluMut's default --name-regex is (?P<sample>.+)_(?P<segment>.+): it expects
 * the sample name to be PART OF the header (e.g. >mysample_HA). IRMA's
 * per-sample consensus headers carry no sample name at all — just a bare
 * segment code like >A_HA_H5 — so rename_for_flumut.R injects the sample
 * name from each file's basename and normalises the segment code to
 * FluMutDB's vocabulary (confirmed against flumut_db.sql's `annotations`
 * table: PB2/PB1/PA/HA/NP/NA/MP/NS). Verified against real IRMA output
 * (Bird_Flu/IRMA-consensus-contigs) before wiring this in.
 *
 * FluMut runs once as a single batch across every sample's consensus, which
 * is its intended usage — not once per sample.
 *
 * DELIBERATELY NOT using `flumut --update`: FluMutDB is a living database, and
 * calling --update here would mean the SAME pipeline version can report
 * different markers on different days with no record of why — the exact
 * silent-nondeterminism failure this pipeline's reproducibility work (pinned
 * envs, trace.txt, timeline.html) exists to prevent. The DB version actually
 * used is pinned by the flumut=0.6.5 conda package and is bundled in the
 * container image, so it is fixed for as long as that image tag is fixed.
 * `flumut --version` is captured alongside the results as the provenance
 * record. To deliberately pick up new markers, rebuild the image against a
 * newer flumut pin — a conscious, recorded decision, not an implicit one.
 */
process FLUMUT {
    label 'process_low'
    publishDir "${params.outdir}/flumut", mode: params.publish_mode

    input:
    path scripts
    path consensus_dir, stageAs: 'IRMA-consensus-contigs'

    output:
    path 'markers.tsv',        emit: markers,    optional: true
    path 'mutations.tsv',      emit: mutations,  optional: true
    path 'literature.tsv',     emit: literature, optional: true
    path 'flumut_report.xlsm', emit: report,     optional: true
    path 'flumut_version.txt', emit: version

    script:
    """
    flumut --version > flumut_version.txt

    Rscript ${scripts}/rename_for_flumut.R batch.fasta IRMA-consensus-contigs/*.fasta

    if [ -s batch.fasta ]; then
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m markers.tsv -M mutations.tsv -l literature.tsv \\
               -x flumut_report.xlsm \\
               batch.fasta
    else
        echo "no consensus sequences to screen — skipping flumut" >&2
    fi
    """
}

process FLUMUT_LOWFREQ {
    label 'process_low'
    publishDir "${params.outdir}/flumut_lowfreq", mode: params.publish_mode

    input:
    path scripts
    path consensus_dir, stageAs: 'IRMA-consensus-contigs'
    path vcf_dirs,      stageAs: 'vcf_files/*'

    output:
    path 'markers.tsv',        emit: markers,    optional: true
    path 'mutations.tsv',      emit: mutations,  optional: true
    path 'literature.tsv',     emit: literature, optional: true
    path 'flumut_report.xlsm', emit: report,     optional: true
    path 'flumut_version.txt', emit: version

    script:
    freq_pct = (params.flumut_freq_threshold * 100).toInteger()
    """
    flumut --version > flumut_version.txt

    # Build list of paired FASTA/VCF files for apply_lofreq_to_consensus.R
    r_args=""
    for consensus_fa in IRMA-consensus-contigs/*.fasta; do
        [ -f "\$consensus_fa" ] || continue
        sample=\$(basename "\$consensus_fa" .fasta)
        vcf_file="vcf_files/\${sample}/lofreq-called-variants.vcf"
        if [ -f "\$vcf_file" ]; then
            r_args="\$r_args \$consensus_fa \$vcf_file"
        fi
    done

    if [ -z "\$r_args" ]; then
        echo "no matched IRMA consensus + LoFreq VCF pairs found" >&2
        touch markers.tsv mutations.tsv literature.tsv
        exit 0
    fi

    Rscript ${scripts}/apply_lofreq_to_consensus.R mutated.fasta ${params.flumut_freq_threshold} \$r_args

    if [ ! -s mutated.fasta ]; then
        echo "no low-frequency variants (AF >= ${freq_pct}%) found — skipping flumut" >&2
        touch markers.tsv mutations.tsv literature.tsv
        exit 0
    fi

    Rscript ${scripts}/rename_for_flumut.R batch.fasta mutated.fasta

    if [ -s batch.fasta ]; then
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m markers.tsv -M mutations.tsv -l literature.tsv \\
               -x flumut_report.xlsm \\
               batch.fasta
    else
        echo "no mutated sequences to screen — skipping flumut" >&2
        touch markers.tsv mutations.tsv literature.tsv
    fi
    """
}

/* ==========================================================================
 * Workflow
 * ========================================================================== */
workflow {

    if (params.help) {
        helpMessage()
        return
    }

    // metadata and aa_db are deliberately absent from this list: both are
    // optional. Without metadata the summaries simply are not grouped; without
    // aa_db the curated-site join is skipped. The full variant table and
    // amino-acid table — the substantive outputs — need neither.
    ['read_directory','rename_file','reference'].each { req ->
        if (!params[req]) error("Missing required parameter: --${req}  (see --help)")
    }

    if (params.wfabc && !params.metadata) {
        error("--wfabc needs --metadata: selection is estimated from allele-frequency\n" +
              "  time series, which are reconstructed by joining variants to the\n" +
              "  individual and time-point columns of the metadata.")
    }

    /*
     * Build the sample channel straight from the rename CSV.
     *
     * This replaces organizeReads.R, which physically COPIED every raw fastq
     * into organized-reads/ — doubling disk usage before analysis even started.
     * Nextflow stages by symlink instead, so nothing is duplicated.
     */
    read_dir = file(params.read_directory)

    samples = channel
        .fromPath(params.rename_file)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // the shipped CSVs carry a UTF-8 BOM, which corrupts the first header key
            def key  = row.keySet().find { k -> k.toString().replace('﻿','') == 'File' }
            def pref = row[key]?.toString()?.trim()
            def name = row.Sample?.toString()?.trim()
            if (!pref || !name) error "Bad row in ${params.rename_file}: ${row}"

            def r1 = file("${read_dir}/${pref}*_R1_*.fastq.gz")
            def r2 = file("${read_dir}/${pref}*_R2_*.fastq.gz")
            if (!r1 || !r2) error "No read pair found for '${pref}' in ${read_dir}"
            tuple(name, r1 instanceof List ? r1[0] : r1, r2 instanceof List ? r2[0] : r2)
        }

    ref = PREPARE_REFERENCE(file(params.reference)).index

    trimmed = FASTP(samples).reads

    // Optional IRMA parameter file, staged into every IRMA task as irma_config.sh.
    irma_cfg = params.irma_config
        ? channel.value(file(params.irma_config))
        : channel.value([])

    // Collected IRMA output directories, or an empty list when IRMA is off.
    // An empty list stages nothing, which is how an optional input is expressed.
    irma_dirs = params.run_irma
        ? IRMA(trimmed, irma_cfg).results.map { _s, d -> d }.collect()
        : channel.value([])

    ubam     = FASTQ_TO_SAM(trimmed).bam
    reverted = REVERT_SAM(ubam).bam
    grouped  = ADD_READ_GROUPS(reverted).bam
    mapped   = BWA_MAP(grouped, ref).bam
    sorted   = SORT_BAM(mapped).bam
    marked   = MARK_DUPLICATES(sorted).bam
    final_bam = SET_TAGS(marked, ref).bam

    gvcf     = HAPLOTYPE_CALLER(final_bam, ref).vcf
    geno     = GENOTYPE_GVCF(gvcf, ref).vcf
    selected = SELECT_VARIANTS(geno)
    filtered = FILTER_VARIANTS(selected.snps.join(selected.indels), ref).vcf
    lofreq   = LOFREQ(final_bam, ref).vcf

    // The R stage summarises across ALL samples, so it must wait for every one.
    // Group each sample's VCFs into a directory named after it, then collect —
    // this is both the completion gate and the vcf_files/<sample>/ layout that
    // convertVCFtoTable.R parses sample names out of.
    per_sample = filtered
        .join(lofreq)
        .map { sample, snps, indels, lf -> tuple(sample, [snps, indels, lf]) }

    vcf_dirs = GATHER_SAMPLE_VCFS(per_sample).dir.collect()

    // An empty list stages nothing, which is how an absent optional file is
    // expressed — see the AA_DB/METADATA handling in R_ANALYSIS.
    aa_db_ch    = params.aa_db    ? channel.value(file(params.aa_db))    : channel.value([])
    metadata_ch = params.metadata ? channel.value(file(params.metadata)) : channel.value([])

    r = R_ANALYSIS(
        file("${projectDir}/Scripts"),
        vcf_dirs,
        ref.map { r, _idx, _dict -> r },
        aa_db_ch,
        metadata_ch,
        irma_dirs
    )

    // Needs IRMA consensus, so it can only run when IRMA ran.
    if (params.run_irma && params.flumut) {
        FLUMUT(file("${projectDir}/Scripts"), r.irma)
    }

    // Screen low-frequency variants above threshold for H5N1 markers.
    // Applies LoFreq variants to IRMA consensus sequences, creating
    // mutated pseudo-consensus for FluMut marker screening.
    if (params.run_irma && params.flumut_lowfreq) {
        FLUMUT_LOWFREQ(file("${projectDir}/Scripts"), r.irma, vcf_dirs)
    }

    // Optional population-genetics analyses. Both read the config.cfg written by
    // R_ANALYSIS, which is also what sequences them after it.
    if (params.snpgenie) {
        SNPGENIE(
            file("${projectDir}/Scripts"),
            r.config,
            vcf_dirs,
            ref.map { r_fa, _idx, _dict -> r_fa },
            metadata_ch
        )
    }

    if (params.wfabc) {
        WFABC(
            file("${projectDir}/Scripts"),
            r.config,
            vcf_dirs,
            ref.map { r_fa, _idx, _dict -> r_fa },
            metadata_ch
        )
    }
}
