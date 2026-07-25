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
    Flumina ${workflow.manifest.version}

    Required:
      --read_directory   Directory of raw paired fastq.gz
      --rename_file      CSV with File,Sample columns
      --reference        Reference FASTA
      --metadata         Metadata CSV containing a Sample column
      --aa_db            Curated amino-acid database CSV

    Common:
      --outdir           Output directory              [${params.outdir}]
      --min_depth        Minimum depth to keep a call  [${params.min_depth}]
      --run_irma         Run IRMA assembly             [${params.run_irma}]
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

process IRMA {
    tag   "$sample"
    label 'process_high'
    publishDir "${params.outdir}/IRMA_results", mode: params.publish_mode

    input:
    tuple val(sample), path(r1), path(r2)

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
    """
    # Paths are relative to this work dir — never the submitting host.
    cat > config.cfg <<'CFG_END'
OUTPUT_DIRECTORY="."
REFERENCE_FILE="${reference}"
AA_DB="${aa_db}"
METADATA="${metadata}"
GROUP_NAMES="${params.group_names}"
MIN_DEPTH="${params.min_depth}"
MIN_QUALITY="${params.min_quality}"
MIN_ALLELE_FREQUENCY="${params.min_allele_frequency}"
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

/* ==========================================================================
 * Workflow
 * ========================================================================== */
workflow {

    if (params.help) {
        helpMessage()
        return
    }

    ['read_directory','rename_file','reference','metadata','aa_db'].each { req ->
        if (!params[req]) error("Missing required parameter: --${req}  (see --help)")
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

    // Collected IRMA output directories, or an empty list when IRMA is off.
    // An empty list stages nothing, which is how an optional input is expressed.
    irma_dirs = params.run_irma
        ? IRMA(trimmed).results.map { _s, d -> d }.collect()
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

    R_ANALYSIS(
        file("${projectDir}/Scripts"),
        vcf_dirs,
        ref.map { r, _idx, _dict -> r },
        file(params.aa_db),
        file(params.metadata),
        irma_dirs
    )
}
