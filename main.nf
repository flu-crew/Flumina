#!/usr/bin/env nextflow
/*
 * Flumina 2.0 — Nextflow pipeline.
 *
 * The output layout (BAM_files/, vcf_files/, logs/, IRMA_results/,
 * processed-reads/) matches the earlier Snakemake version, so the downstream R
 * scripts run without changes.
 *
 *   nextflow run . -profile test,docker
 *   nextflow run . -profile docker --read_directory reads --rename_file r.csv \
 *                  --reference reference.fa --metadata meta.csv --aa_db db.csv
 */

nextflow.enable.dsl = 2

/*
 * Locate the R1/R2 pair for one prefix from the rename CSV.
 *
 * Sequencing facilities return reads in layouts that one glob cannot cover:
 * nested per-sample or per-project directories from bcl2fastq, and several
 * mate-marker conventions. Gather files recursively. Then try the mate markers
 * in order of decreasing specificity, so `_R1_` wins over a bare `_1`, and a
 * file like SAMPLE_S1_L001_R1_001.fastq.gz cannot be mistaken for its own mate.
 */
def findReadPair(read_dir, prefix, sample_name) {
    // Accept uncompressed FASTQ and gzipped FASTQ.
    def exts = ['fastq.gz', 'fq.gz', 'fastq', 'fq']

    // Try prefixes in this order. The Sample column is a fallback: a re-run over
    // already-renamed reads no longer has the File value in any filename.
    def gather = { String stem ->
        def out = []
        exts.each { ext ->
            // Both globs are needed: `**` does not match zero directories, so
            // the recursive form alone misses reads directly in read_dir.
            ["${read_dir}/${stem}*.${ext}", "${read_dir}/**/${stem}*.${ext}"].each { pattern ->
                def found = file(pattern)
                if (found) out += (found instanceof List ? found : [found])
            }
        }
        out.unique { it.toString() }.sort { it.name }
    }

    // Use findResult, not a for/break loop: Nextflow's strict parser (25.x+)
    // rejects `for` loops in pipeline scripts.
    def hits = [prefix, sample_name].findAll { it }
                                    .findResult { stem -> gather.call(stem) ?: null } ?: []

    // Return a reason for per-sample problems instead of throwing. A library may
    // legitimately fail sequencing; report and skip that sample so other samples
    // can proceed.
    // Problems with the run as a whole still stop the pipeline.
    if (!hits) {
        return [null, "no read files found matching ${prefix}* or ${sample_name}*"]
    }

    // Ordered most specific first. Each entry is [R1 marker, R2 marker].
    def conventions = [
        ['_R1_', '_R2_'], ['_R1.', '_R2.'], ['.R1.', '.R2.'],
        ['-R1-', '-R2-'], ['-R1.', '-R2.'], ['_1.',  '_2.'],
    ]

    def resolved = conventions.findResult { c ->
        def r1 = hits.findAll { it.name.contains(c[0]) }
        def r2 = hits.findAll { it.name.contains(c[1]) }
        if (r1.size() == 1 && r2.size() == 1) return [r1[0], r2[0]]
        if (r1.size() > 1 || r2.size() > 1) {
            return [null, "matched ${r1.size()} R1 and ${r2.size()} R2 files using " +
                          "'${c[0]}'/'${c[1]}' (${hits*.name.join(', ')}); the File value " +
                          "must identify one sample"]
        }
        return null
    }
    if (resolved) return resolved

    // Last resort: with exactly two files and no recognised marker, take them in
    // sorted order. Every convention puts the first mate first alphabetically,
    // so this is usually correct. It is still a guess, so the pipeline warns.
    if (hits.size() == 2) {
        log.warn """Reads for '${prefix}' carry no recognised mate marker; assuming
  sorted order:  R1 = ${hits[0].name}
                 R2 = ${hits[1].name}
  Check that is correct — the pipeline cannot verify it."""
        return [hits[0], hits[1]]
    }

    return [null, "found ${hits.size()} file(s) (${hits*.name.join(', ')}) but no R1/R2 " +
                  "pair; recognised markers are " +
                  "${conventions.collect { "${it[0]}/${it[1]}" }.join(', ')}, and sorted-order " +
                  "fallback needs exactly two files"]
}

/*
 * Read a boolean parameter safely, whatever its type.
 *
 * From Nextflow 25.x, a command-line `--flag false` arrives as the string
 * "false", not a boolean, and every non-empty string is truthy in Groovy. So
 * `if (params.wfabc)` is true even when the flag is set to false. Convert every
 * boolean parameter through this function.
 */
def asBool(value) {
    if (value instanceof Boolean) return value
    return value?.toString()?.trim()?.toLowerCase() in ['true', 't', 'yes', '1']
}

/*
 * Return whether a program is enabled. This must be a top-level function, not a
 * closure. The strict parser resolves `progOn(...)` as a call to a function, so
 * a closure of that name fails to compile as "not defined".
 */
def progOn(String key) {
    return asBool(params[key])
}

/*
 * Numeric counterpart to asBool. Command-line parameters arrive as strings, and
 * Groovy reads String * Integer as repetition ("0.01" * 100 repeats the
 * string), so direct arithmetic on a parameter fails. Convert every numeric
 * parameter through this function before arithmetic.
 */
def asNum(value, fallback = 0) {
    if (value instanceof Number) return value
    try { return new BigDecimal(value?.toString()?.trim()) }
    catch (ignored) { return fallback }
}

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
                         Tested against the depth the callers see, not the raw
                         count; depth_profiles/ publishes both (col 3, col 4)
      --min_quality      Minimum quality to keep       [${params.min_quality}]
      --min_allele_frequency  Minimum allele frequency [${params.min_allele_frequency}]
      --group_names      Metadata column to group by   [${params.group_names}]
      --ivar             Call variants with iVar too   [${params.ivar}]
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
 * The four indexing steps each take seconds and always run together, so one
 * process runs them all and avoids extra scheduling overhead.
 * ========================================================================== */
process PREPARE_REFERENCE {
    tag   "reference"
    label 'process_low'
    publishDir "${params.outdir}/Reference", mode: params.publish_mode
    /* findAAChanges.R reads ${OUTPUT_DIRECTORY}/reference.fa. The old Bash
     * driver copied the reference to the output root before calling snakemake,
     * so the R scripts depend on it being there. Reproduce that placement.
     *
     * saveAs returns null when the user-supplied reference is already that file.
     * Without the guard this process
     * republishes its own input on top of itself: same bytes, new mtime, and
     * Nextflow's default cache hash includes an input's last-modified time. So
     * PREPARE_REFERENCE could never be cached, every downstream task saw a
     * changed input, and `-resume` re-ran the entire alignment and calling
     * chain. Measured on the swine WGS run: a resume that should have restarted
     * at FLUMUT began re-running BWA_MAP across all 143 samples.
     *
     * Keeping the reference inside the output directory is a natural thing to
     * do — it makes the run self-contained — so this is guarded rather than
     * merely documented.
     */
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'reference.fa',
        saveAs: { fn ->
            def src = file(params.reference).toAbsolutePath().normalize()
            def dst = file("${params.outdir}/${fn}").toAbsolutePath().normalize()
            src == dst ? null : fn
        }

    input:
    // stageAs gives the input a fixed, distinct name. Without it, a reference
    // already called reference.fa stages onto the output name and the copy
    // fails with "are the same file".
    path reference, stageAs: 'input_reference'

    output:
    tuple path('reference.fa'), path('reference.fa.*'), path('reference.dict'), emit: index

    script:
    """
    # A reference from Windows or GISAID may have CRLF line endings. samtools
    # faidx, Biostrings, and the python scripts drop the CR, but bwa keeps it as
    # a sequence character. The shift then makes the callers report most of the
    # segment as fixed variants. Normalise here, because this reference.fa is the
    # copy that is indexed, mapped against, and published. Numbers in HANDOFF.md.
    if [ -n "\$(tr -cd '\\r' < ${reference} | head -c1)" ]; then
        echo "NOTE: reference has CRLF line endings - using a normalised copy"
    fi
    tr -d '\\r' < ${reference} > reference.fa
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
 * IRMA parameters (TMP, SINGLE_LOCAL_PROC, ...) come from an optional user file
 * passed with `--external-config`. IRMA 1.3.5 needs this flag; a staged file
 * alone is not read.
 *
 * When --irma_config is not set, an empty list stages nothing, the flag is
 * omitted, and IRMA uses its module defaults.
 */
process IRMA {
    tag   "$sample"
    label 'process_high'
    // CLEANUP=TRUE (default) publishes only tables/ and logs/. The consensus is
    // saved in IRMA-consensus-contigs/ and IRMA-amended-contigs/, so IRMA_results/
    // does not repeat it. CLEANUP=FALSE (publish_full_irma) publishes the full tree.
    //
    // publishDir cannot select part of a directory output; it publishes the whole
    // dir. So the script builds a small copy in irma_publish/, and this line
    // publishes it with the prefix removed to keep the native
    // IRMA_results/<sample>/ path. Downstream always stages the full <sample>/
    // dir, so the trim does not affect it.
    publishDir "${params.outdir}/IRMA_results", mode: params.publish_mode,
               enabled: params.publish_full_irma
    publishDir "${params.outdir}/IRMA_results", mode: params.publish_mode,
               enabled: !params.publish_full_irma,
               saveAs: { fn -> fn.startsWith('irma_publish/') ? fn.replaceFirst('irma_publish/', '') : null }

    input:
    tuple val(sample), path(r1), path(r2)
    path irma_cfg, stageAs: 'irma_config.sh'

    output:
    tuple val(sample), path("${sample}"), emit: results
    // Small copy (tables/ and logs/) to publish on CLEANUP=TRUE. Absent for the
    // full tree.
    path "irma_publish/${sample}", emit: pub, optional: true

    script:
    def cfg_arg = irma_cfg ? "--external-config run_irma_config.sh" : ""
    def clean_irma = !params.publish_full_irma
    """
    # IRMA 1.3.5 sizes itself with `irma-core num-procs --cap-cores-using-env`,
    # which reads the CPU affinity mask. Under Slurm the mask reflects
    # --cpus-per-task. On PBS without cgroups, or under a Docker CPU quota, the
    # mask shows the whole machine and IRMA oversubscribes the node. Set the CPU
    # count that Nextflow requested instead.
    export LOCAL_PROCS_OVERRIDE=${task.cpus}

    # TMP in an IRMA configuration must be an absolute path that already exists.
    # IRMA builds its working path from it —
    #     ppath="\$TMP"/<user>/IRMAv<version>/<run>-<token>
    # — never creates it, and changes directory as it works, so a relative TMP
    # stops resolving and IRMA produces empty output while it still exits 0.
    #
    # The staged configuration is a symlink into another directory, so rewrite a copy.
    if [ -f irma_config.sh ]; then
        cp irma_config.sh run_irma_config.sh
        irma_tmp=\$(sed -n 's/^[[:space:]]*TMP=//p' run_irma_config.sh | tail -1 | tr -d "\\"'")
        if [ -n "\$irma_tmp" ]; then
            case "\$irma_tmp" in
                /*) ;;
                *)  irma_tmp="\$PWD/\$irma_tmp" ;;
            esac
            mkdir -p "\$irma_tmp" || {
                echo "ERROR: TMP from irma_config.sh cannot be created: \$irma_tmp" >&2
                echo "       IRMA needs this directory and will not create it itself." >&2
                exit 1
            }
            sed -i "s|^[[:space:]]*TMP=.*|TMP=\$irma_tmp|" run_irma_config.sh
        fi
    fi

    IRMA FLU ${cfg_arg} ${r1} ${r2} ${sample}

    # IRMA exits 0 even when it produces no assembly: it writes an empty output
    # skeleton (amended_consensus/, tables/, logs/ ...) and returns success. So
    # check that a consensus exists. An empty result for one sample can be
    # legitimate (too few reads); an empty result for every sample means IRMA
    # failed — see .command.err.
    if ! ls ${sample}/*.fasta >/dev/null 2>&1; then
        echo "WARNING: IRMA produced no consensus sequence for ${sample}." >&2
        echo "         A sample with too few influenza reads can do this legitimately." >&2
        echo "         If it happens for every sample, IRMA itself failed — check this" >&2
        echo "         task's .command.err for repeated 'exec failed' lines." >&2
    fi

    # On CLEANUP=TRUE, clean this sample as soon as IRMA finishes, not at the end
    # of the run. Remove the large intermediate files. Copy the tables and logs
    # into irma_publish/ to publish them. Keep the consensus (*.fasta), the
    # amended consensus, the tables, and the logs in the work dir so the
    # downstream steps can still stage them.
    if [ "${clean_irma}" = "true" ]; then
        find ${sample}/ -mindepth 1 -maxdepth 1 \\
            ! -name '*.fasta' ! -name amended_consensus ! -name tables ! -name logs \\
            -exec rm -rf {} +

        mkdir -p irma_publish/${sample}
        [ -d ${sample}/tables ] && cp -r ${sample}/tables irma_publish/${sample}/tables
        [ -d ${sample}/logs ]   && cp -r ${sample}/logs   irma_publish/${sample}/logs
    fi
    """
}

/* ==========================================================================
 * Alignment — GATK best-practice style unmapped-BAM route
 * ========================================================================== */
process FASTQ_TO_SAM {
    tag "$sample"; label 'gatk'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, enabled: params.publish_intermediate_bams

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
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, enabled: params.publish_intermediate_bams

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
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, enabled: params.publish_intermediate_bams

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
 * Run bwa mem without -t.
 *
 * bwa estimates the insert-size distribution per read batch, and batch
 * composition depends on the thread count. So -t shifts pairing for a few
 * marginal reads, which changes a few borderline GATK calls at segment edges
 * and makes the BAM non-deterministic.
 *
 * If throughput matters more than a reproducible BAM, use '-t N -K 100000000':
 * a fixed chunk size makes threaded output deterministic, though not identical
 * to -t 1.
 */
process BWA_MAP {
    tag "$sample"; label 'process_high'
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, enabled: params.publish_intermediate_bams

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
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, enabled: params.publish_intermediate_bams

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
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, pattern: '*.bam', enabled: params.publish_intermediate_bams
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
    publishDir { "${params.outdir}/BAM_files/${sample}" }, mode: params.publish_mode, pattern: '*.bam', enabled: params.publish_intermediate_bams
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
    // GATK 4.1+ uses the new qual calculator by default, so
    // --use-new-qual-calculator is not needed.
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
    // VariantFiltration flags a record when the expression is true, so each
    // expression below names the condition for a FAILING record. High FS and
    // high SOR indicate strand bias, so those use "greater than"; the rest use
    // "less than", as in GATK's recommendations.
    //
    // Thresholds are GATK's published germline hard-filter recommendations per
    // variant type, with two changes for this data:
    //   - SNP QD uses 15.0, not GATK's 2.0. The QD distribution here is bimodal,
    //     so 15.0 separates thin calls from the dense high-QD peak. Indels keep
    //     QD<2.0.
    //   - Indels take a looser FS and ReadPosRankSum than SNPs.
    //
    // The RankSum filters have no value on a haploid hom-alt call, which has no
    // reference reads to rank. No records are dropped here: the FILTER column
    // annotates the calls, and FluLens decides which to hide.
    """
    gatk VariantFiltration -R ${ref} -V ${snps} -O gatk4-filtered-snps.vcf \\
        -filter "QUAL<30.0"            --filter-name "QUAL" \\
        -filter "QD<15.0"              --filter-name "QD" \\
        -filter "SOR>3.0"              --filter-name "SOR" \\
        -filter "FS>60.0"              --filter-name "FS" \\
        -filter "MQ<40.0"              --filter-name "MQ" \\
        -filter "MQRankSum<-12.5"      --filter-name "MQRankSum" \\
        -filter "ReadPosRankSum<-8.0"  --filter-name "ReadPosRankSum"

    gatk VariantFiltration -R ${ref} -V ${indels} -O gatk4-filtered-indels.vcf \\
        -filter "QD<2.0"               --filter-name "QD" \\
        -filter "QUAL<30.0"            --filter-name "QUAL" \\
        -filter "FS>200.0"             --filter-name "FS" \\
        -filter "ReadPosRankSum<-20.0" --filter-name "ReadPosRankSum"
    """
}

/*
 * LoFreq.
 *
 * -B disables BAQ, matching the iVar path. BAQ suppresses real low-frequency
 * SNPs next to indels. Those extra calls are kept and flagged in the variant
 * table as `dist_to_indel`, rather than lost. Numbers in HANDOFF.md.
 */
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
    lofreq call -B -f ${ref} -o lofreq-called-variants.vcf ${bam}
    """
}

/*
 * iVar — the third caller.
 *
 * Do not provide a GFF. `ivar variants -g` translates amino acids from the
 * start of each reference sequence in frame 1. That is correct for the eight
 * primary ORFs but wrong for M2, NEP, PA-X, and PB1-F2, and it returns a
 * plausible residue for a codon that does not exist rather than an error.
 * Flumina annotates through Scripts/fluORFs.R, which walks the real coding
 * intervals, so iVar contributes nucleotide calls only.
 *
 * mpileup flags and their purposes:
 *   -aa       every position, including zero-coverage ones, so iVar's own
 *             depth filter removes them rather than mpileup's silence
 *   -A        count anomalous read pairs. LoFreq's DP4 requires proper pairs;
 *             iVar is left inclusive so the two callers are not tuned to agree
 *             by construction
 *   -B        no BAQ. On by default, it re-scores base qualities around indels
 *             and suppresses real low-frequency SNPs next to them
 *   -d 0      no depth cap. The default 8000 would truncate deep libraries
 *   -Q 0      let iVar apply the quality threshold via -q, rather than filter
 *             twice at two different values
 *
 * Thresholds come from the same parameters as LoFreq and the R stage, so the
 * three callers use one set of numbers.
 */
process IVAR {
    tag "$sample"; label 'process_medium'
    publishDir { "${params.outdir}/vcf_files/${sample}" }, mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref), path(ref_idx), path(ref_dict)

    output:
    tuple val(sample), path('ivar-called-variants.tsv'), emit: tsv

    script:
    """
    samtools mpileup -aa -A -B -d 0 -Q 0 --reference ${ref} ${bam} \\
      | ivar variants -p ivar-called-variants -r ${ref} \\
          -m ${asNum(params.min_depth)} \\
          -q ${asNum(params.min_quality)} \\
          -t ${asNum(params.min_allele_frequency)}

    # ivar exits 0 with only a header when nothing passes. An absent file would
    # fail the output declaration, so write the header. An empty result is valid
    # for a thin library.
    if [ ! -s ivar-called-variants.tsv ]; then
      printf 'REGION\\tPOS\\tREF\\tALT\\tREF_DP\\tREF_RV\\tREF_QUAL\\tALT_DP\\tALT_RV\\tALT_QUAL\\tALT_FREQ\\tTOTAL_DP\\tPVAL\\tPASS\\n' \\
        > ivar-called-variants.tsv
    fi
    """
}

/* ==========================================================================
 * Downstream R analysis
 * --------------------------------------------------------------------------
 * The R scripts read a KEY=VALUE config.cfg. The pipeline writes that config
 * from its own params, so the scripts run unchanged and stay runnable by hand.
 * ========================================================================== */
/*
 * Gather one sample's VCFs into a directory named after the sample.
 *
 * convertVCFtoTable.R derives the sample name from the directory component of
 * the path (gsub("/.*", "", ...)), so the vcf_files/<sample>/<file>.vcf layout
 * is required. Nextflow flattens staged filenames, and every sample produces an
 * identically-named lofreq-called-variants.vcf, so each sample needs its own
 * directory before staging.
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
 * Per-position depth in reference coordinates. Run for every sample of every
 * run — see the workflow body for why it stopped being conditional.
 *
 * It was built for the low-frequency FluMut screen and now serves two consumers
 * that ask different questions of it: that screen asks "was there any read
 * here at all" (column 3), and the MIN_DEPTH accounting asks "did any caller
 * get to evaluate this position" (column 4). The second is why it is no longer
 * gated on the first.
 *
 * apply_lofreq_to_consensus.R paints LoFreq calls onto the reference, so every
 * position the sample has no reads at was silently taking a REFERENCE base and
 * standing in for real data. FluMut then screened it and reported no marker
 * there — absence of evidence arriving as evidence of absence. That was the
 * documented residual of the 2026-08-01 reference-painting fix, left open
 * because "this process stages no depth source". This is that source.
 *
 * It has to come from the BAM rather than from IRMA's own coverage tables:
 * IRMA's are in ITS consensus coordinates, and reconciling those back to the
 * reference is the exact coordinate problem the reference-painting approach
 * exists to avoid. `samtools depth -a` on the reference-aligned BAM is already
 * in the coordinate system the mask has to apply to.
 *
 * -a emits every position, including zero-coverage positions, which is the
 * point — without it the uncovered positions are simply absent and
 * indistinguishable from a truncated file. -Q 0 because the mask is asking
 * "was there any read here at all", not "was there a confident read".
 *
 * The whole reference is 13,133 bp, so this is ~13k lines per sample and costs
 * nothing to keep.
 *
 * Four columns are emitted: contig, position, raw depth, and the depth the
 * variant callers can actually see. MIN_DEPTH is tested against the second one
 * — iVar gets it as -m against a pileup with overlapping mates zeroed and -q
 * applied — so publishing the raw count alone stated the floor against a number
 * no output carried. Column 4 comes from the same mpileup IVAR runs, making it
 * the caller's own quantity rather than a conversion factor. Columns 1-3 are
 * unchanged and byte-identical to earlier runs. Measurement in HANDOFF.md.
 */
process DEPTH_PROFILE {
    tag "$sample"
    label 'process_low'
    // Publish this file because both consumers are outside this pipeline:
    // FluLens is the only component that maps FluMut's numbering back to
    // reference positions, and the MIN_DEPTH gap can only be read from the file
    // itself, since nothing in the run reports it.
    publishDir "${params.outdir}/depth_profiles", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bai)
    tuple path(ref), path(ref_idx), path(ref_dict)
    path scripts

    output: path "${sample}.depth", emit: depth

    script:
    """
    # Named off the sample so a sample called "raw" cannot collide with it.
    samtools depth -a -Q 0 ${bam} > '${sample}.raw-depth.tmp'

    # Keep this command byte-for-byte identical to the iVar mpileup command. If
    # the commands diverge, column 4
    # stops being the quantity -m is tested against. Change them together.
    samtools mpileup -aa -A -B -d 0 -Q 0 --reference ${ref} ${bam} \\
      | Rscript ${scripts}/visible_depth.R \\
          --raw '${sample}.raw-depth.tmp' \\
          --min-quality ${asNum(params.min_quality)} \\
          --out '${sample}.depth'
    """
}

/*
 * Downstream R analysis — fully relocatable.
 *
 * Every input is staged into this task's work directory, and the generated
 * config.cfg sets OUTPUT_DIRECTORY=".", so the R scripts resolve every path
 * relative to the run directory. Nothing references a host path, so this runs
 * unchanged on a laptop, an HPC node, or a cloud container.
 *
 * The R scripts build every path as paste0(OUTPUT_DIRECTORY, "/..."), so a
 * relative root works.
 */
process R_ANALYSIS {
    label 'process_medium'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'variant_analysis'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'IRMA-consensus-contigs'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: 'IRMA-amended-contigs'
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
    path 'IRMA-consensus-contigs',  emit: irma,     optional: true
    path 'IRMA-amended-contigs',    emit: amended,  optional: true

    script:
    def irma_step = asBool(params.run_irma)
        ? "Rscript ${scripts}/organizeIRMA.R config.cfg"
        : "echo 'IRMA disabled, skipping organizeIRMA.R'"
    // Both inputs are optional. When not supplied, nothing is staged, so the
    // input variable is an empty list and would render as an empty string —
    // write the literal NULL the R scripts test for instead.
    def aa_db_cfg    = params.aa_db    ? "${aa_db}"    : 'NULL'
    def metadata_cfg = params.metadata ? "${metadata}" : 'NULL'
    // Grouping uses a metadata column, so without metadata there is
    // nothing to group by and outputSummary.R would fail looking for it.
    def group_cfg    = params.metadata ? "${params.group_names}" : 'NULL'
    """
    # Paths are relative to this work directory, never to the submitting host.
    cat > config.cfg <<'CFG_END'
OUTPUT_DIRECTORY="."
REFERENCE_FILE="${reference}"
AA_DB="${aa_db_cfg}"
METADATA="${metadata_cfg}"
GROUP_NAMES="${group_cfg}"
MIN_DEPTH="${params.min_depth}"
MIN_QUALITY="${params.min_quality}"
MIN_ALLELE_FREQUENCY="${params.min_allele_frequency}"
MIN_ALT="${params.min_alt}"
# runSNPGenie.R reads these thresholds under different names (MIN_ALLELE_FREQ /
# MIN_COVERAGE, not MIN_ALLELE_FREQUENCY / MIN_DEPTH). Write both spellings so
# SNPGenie uses the same depth and frequency thresholds as every other step;
# without them it applies no filtering.
MIN_ALLELE_FREQ="${params.min_allele_frequency}"
MIN_COVERAGE="${params.min_depth}"
DEDUP_KEYS="${params.dedup_keys}"
# Record which programs ran. The R stage reads the thresholds above; these
# switches document how this run was produced.
IRMA="${asBool(params.run_irma) ? 'TRUE' : 'FALSE'}"
LOFREQ="${asBool(params.run_lofreq) ? 'TRUE' : 'FALSE'}"
GATK4="${asBool(params.run_gatk4) ? 'TRUE' : 'FALSE'}"
IVAR="${asBool(params.ivar) ? 'TRUE' : 'FALSE'}"
FLUMUT="${asBool(params.flumut) ? 'TRUE' : 'FALSE'}"
WFABC="${asBool(params.wfabc) ? 'TRUE' : 'FALSE'}"
SNPGENIE="${asBool(params.snpgenie) ? 'TRUE' : 'FALSE'}"
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
 * GTF+FASTA pairs SNPGenie needs, then runSNPGenie.R runs SNPGenie over every
 * sample's VCF and collects the per-sample output into combined tables.
 *
 * The config.cfg from R_ANALYSIS is reused, not regenerated, so these steps
 * cannot drift from the main analysis. It is copied before being appended to,
 * because the staged original is a symlink into another task's directory.
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
    # runSNPGenie.R setwd()s into each per-sample directory, then keeps using
    # paths built from OUTPUT_DIRECTORY, so OUTPUT_DIRECTORY must be absolute.
    # The relocatable OUTPUT_DIRECTORY="." breaks once it changes directory. The
    # R config parser keeps the last value for a key, so this appended \$PWD
    # (resolved in this task's directory) overrides it without a hardcoded path.
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
 * the sample name in the header (e.g. >mysample_HA). IRMA's per-sample
 * consensus headers carry only a bare segment code like >A_HA_H5, so
 * rename_for_flumut.R injects the sample name from each file's basename and
 * normalises the segment code to FluMutDB's vocabulary (PB2, PB1, PA, HA, NP,
 * NA, MP, NS).
 *
 * FluMut runs once as a single batch across every sample's consensus, which is
 * its intended usage — not once per sample.
 *
 * Do not use `flumut --update`: FluMutDB is a living database, so --update would
 * let the same pipeline version report different markers on different days with
 * no record of why.
 *
 * Pinning the tool is not enough. flumutdb is a separate conda package, and
 * flumut=0.6.5 resolves with database 6.5 or 6.7 depending on the build date,
 * so the database is pinned explicitly too. The database version fixes the
 * markers, not the image tag.
 *
 * `flumut --version` is captured with the results as the provenance record; it
 * reports the tool, not the database. To pick up new markers, move the flumutdb
 * pin — a conscious, recorded decision.
 */
process FLUMUT {
    label 'process_low'
    // Published under variant_analysis/ with the rest of the interpretation:
    // these are marker calls read alongside the amino-acid tables, not a
    // separate kind of output.
    publishDir "${params.outdir}/variant_analysis/flumut", mode: params.publish_mode

    input:
    path scripts
    path consensus_dir, stageAs: 'IRMA-consensus-contigs'
    path reference

    output:
    path 'markers.tsv',            emit: markers,    optional: true
    path 'mutations.tsv',          emit: mutations,  optional: true
    path 'literature.tsv',         emit: literature, optional: true
    path 'flumut_report.xlsm',     emit: report,     optional: true
    path '*_all.tsv',              emit: unfiltered, optional: true
    path 'reference_*.tsv',        emit: refcalls,   optional: true
    // The subtype decision, so consumers read it rather than re-derive it from
    // segment names. See filter_flumut_subtype.R.
    path 'subtype.tsv',            emit: subtype,    optional: true
    path 'flumut_version.txt',     emit: version

    script:
    /*
     * Screen the reference as well, then remove its own findings.
     *
     * FluMut reports every marker a sequence carries, and the reference carries
     * many. Those appear in every sample by construction and say nothing about
     * any of them.
     *
     * markers.tsv drops the rows the reference also has. mutations.tsv is a wide
     * matrix and instead drops only the columns where every sample carries the
     * reference residue — dropping by "the reference has this marker" would
     * discard reversions, and a sample that LOSES a reference marker emits no
     * marker row, so that table is the only place the signal exists.
     *
     * Nothing is discarded: raw output is kept as *_all.tsv and the reference's
     * own calls as reference_*.tsv, so a removal can be explained.
     */
    def subtract = asBool(params.flumut_subtract_reference)
    def keep_hana = asBool(params.flumut_keep_mismatched_ha_na) ? 'TRUE' : 'FALSE'
    """
    flumut --version > flumut_version.txt

    touch batch.fasta
    for f in IRMA-consensus-contigs/*.fasta; do
        if [ -f "\$f" ]; then
            Rscript ${scripts}/rename_for_flumut.R batch.fasta IRMA-consensus-contigs/*.fasta
            break
        fi
    done
    if [ -s batch.fasta ]; then
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m markers.tsv -M mutations.tsv -l literature.tsv \\
               -x flumut_report.xlsm \\
               batch.fasta
    else
        echo "no consensus sequences to screen — skipping flumut" >&2
    fi

    if ${subtract} && { [ -s markers.tsv ] || [ -s mutations.tsv ]; }; then
        Rscript ${scripts}/rename_for_flumut.R ref_batch.fasta ${reference}
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m ref_markers.tsv -M ref_mutations.tsv -l ref_literature.tsv \\
               ref_batch.fasta || true
        Rscript ${scripts}/filter_flumut_reference.R \\
            ref_markers.tsv ref_mutations.tsv \\
            markers.tsv mutations.tsv literature.tsv .
    else
        echo "flumut_subtract_reference disabled or no findings — raw flumut output kept" >&2
    fi

    # Run regardless of the subtraction above because it answers a different
    # question: HA/NA markers are numbered for H5/N1 specifically, so off-subtype
    # they are read against the wrong ruler AND the proteins have diverged too
    # far for equivalence to be assumed. Internal-gene markers are unaffected.
    # IRMA-consensus-contigs is passed as a second subtype source. IRMA writes the
    # subtype it assigned into each consensus header (>A_HA_H5, >A_NA_N1), so a
    # reference with bare A_HA / A_NA names — which is the repo's own reference,
    # and the bundled H5N1 test_dataset — can still be resolved instead of being
    # treated as unconfirmable. The reference still wins where it states one.
    Rscript ${scripts}/filter_flumut_subtype.R ${reference} markers.tsv ${keep_hana} . IRMA-consensus-contigs
    """
}

process FLUMUT_LOWFREQ {
    label 'process_low'
    publishDir "${params.outdir}/variant_analysis/flumut_lowfreq", mode: params.publish_mode

    input:
    path scripts
    path consensus_dir, stageAs: 'IRMA-consensus-contigs'
    path vcf_dirs,      stageAs: 'vcf_files/*'
    path reference
    path depth_files,   stageAs: 'depth/*'

    output:
    path 'markers.tsv',        emit: markers,    optional: true
    path 'mutations.tsv',      emit: mutations,  optional: true
    path 'literature.tsv',     emit: literature, optional: true
    path 'flumut_report.xlsm', emit: report,     optional: true
    path '*_all.tsv',          emit: unfiltered, optional: true
    path 'reference_*.tsv',    emit: refcalls,   optional: true
    // Same as FLUMUT: publish the subtype decision rather than leaving it in a
    // log for every consumer to re-derive.
    path 'subtype.tsv',        emit: subtype,    optional: true
    path 'flumut_version.txt', emit: version

    script:
    freq_pct = (asNum(params.flumut_freq_threshold, 0.01) * 100).toInteger()
    def subtract = asBool(params.flumut_subtract_reference)
    def keep_hana = asBool(params.flumut_keep_mismatched_ha_na) ? 'TRUE' : 'FALSE'
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

    # Reference passed explicitly: the low-frequency screen paints LoFreq calls
    # onto the REFERENCE, not the consensus, because LoFreq's coordinates are the
    # reference's and IRMA's consensus is built de novo. See the script header.
    # depth/ carries one <sample>.depth per sample from DEPTH_PROFILE. Zero-
    # coverage positions are counted and reported, never written as N, because
    # masking manufactures markers; see the script header. The directory is
    # passed rather than the files, so an absent depth source degrades instead of
    # failing. --min-depth reports the exposure against the run's own floor.
    Rscript ${scripts}/apply_lofreq_to_consensus.R mutated.fasta ${params.flumut_freq_threshold} \\
        ${reference} depth --min-depth=${asNum(params.min_depth)} \$r_args

    if [ ! -s mutated.fasta ]; then
        echo "no low-frequency variants (AF >= ${freq_pct}%) found — skipping flumut" >&2
        touch markers.tsv mutations.tsv literature.tsv
        exit 0
    fi

    # No rename step here: apply_lofreq_to_consensus.R already writes
    # >sample_SEGMENT headers. rename_for_flumut.R takes the sample name from the
    # filename, so one combined FASTA would rename every sample to "mutated".
    cp mutated.fasta batch.fasta

    if [ -s batch.fasta ]; then
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m markers.tsv -M mutations.tsv -l literature.tsv \\
               -x flumut_report.xlsm \\
               batch.fasta
    else
        echo "no mutated sequences to screen — skipping flumut" >&2
        touch markers.tsv mutations.tsv literature.tsv
    fi

    # Same reference subtraction as FLUMUT — see the comment there. It matters
    # more here: low-frequency variants only add findings on top of the reference
    # background, so without this the novel calls are harder to see.
    if ${subtract} && { [ -s markers.tsv ] || [ -s mutations.tsv ]; }; then
        Rscript ${scripts}/rename_for_flumut.R ref_batch.fasta ${reference}
        flumut --skip-unmatch-names --skip-unknown-segments \\
               -m ref_markers.tsv -M ref_mutations.tsv -l ref_literature.tsv \\
               ref_batch.fasta || true
        Rscript ${scripts}/filter_flumut_reference.R \\
            ref_markers.tsv ref_mutations.tsv \\
            markers.tsv mutations.tsv literature.tsv .
    else
        echo "flumut_subtract_reference disabled or no findings — raw flumut output kept" >&2
    fi

    # Run regardless of the subtraction above, because it answers a different
    # question: HA/NA markers are numbered for H5/N1, so off-subtype they use the
    # wrong numbering and the proteins have diverged too far to assume
    # equivalence. Internal-gene markers are unaffected. IRMA-consensus-contigs
    # is a SECOND subtype source: IRMA writes the subtype into each consensus
    # header (>A_HA_H5, >A_NA_N1), so a reference with bare A_HA / A_NA names can
    # still be resolved. The reference still wins where it states a subtype.
    Rscript ${scripts}/filter_flumut_subtype.R ${reference} markers.tsv ${keep_hana} . IRMA-consensus-contigs
    """
}

/*
 * Map FluMut marker numbering onto this run's reference.
 *
 * A marker reads `HA1-5:G224S` — residue 224 of HA1, in the numbering of
 * FluMut's own reference. That string does not say where residue 224 falls in
 * the reference this run was mapped against. A constant offset cannot recover
 * it: small proteins carry too few markers to pin an offset, and HA/NA differ
 * by INDELS between subtypes, so no single offset is correct across the protein.
 *
 * FluMut ships its reference sequences and CDS annotations in flumut_db.sqlite,
 * so this is a fact to read and align, not a parameter to estimate. A
 * protein-to-protein alignment absorbs the indels, so HA and NA work on any
 * subtype, not only H5N1.
 *
 * makeGTF.R runs here rather than reuse SNPGENIE's copy, because that process is
 * optional and this must not inherit its switch. It is the same script, so there
 * is still one definition of a CDS.
 */
process FLUMUT_POSITION_MAP {
    label 'process_low'
    // Published beside the marker tables it explains, so a reader finds the two
    // together.
    publishDir "${params.outdir}/variant_analysis/flumut", mode: params.publish_mode

    input:
    path scripts
    path config
    path reference

    output:
    path 'flumut_position_map.tsv', emit: map, optional: true

    script:
    """
    cp ${config} run_config.cfg
    # Absolute, and resolved inside this task's directory, for the same reason
    # SNPGENIE appends it: makeGTF.R builds paths from OUTPUT_DIRECTORY and the
    # relocatable "." does not survive a setwd().
    echo "OUTPUT_DIRECTORY=\\"\$PWD\\"" >> run_config.cfg
    # Point makeGTF.R at the STAGED reference, not the absolute path the config
    # carries. Both work in the container, but only this guarantees the GTF and
    # the translation below come from the same bytes; a numbering derived from a
    # different file than the residues it numbers gives a silent off-by-a-few.
    # The R config parser keeps the last value for a key, so this appended line
    # overrides the config.
    echo "REFERENCE_FILE=\\"\$PWD/${reference}\\"" >> run_config.cfg

    Rscript ${scripts}/makeGTF.R run_config.cfg

    Rscript ${scripts}/flumut_position_map.R \\
        --reference ${reference} \\
        --gtf reference_gtf \\
        --out flumut_position_map.tsv
    """
}

/*
 * Map IRMA's consensus onto this reference, using the same rationale as
 * FLUMUT_POSITION_MAP, one level down.
 *
 * IRMA's consensus is the one FluMut screens. Calls above 50% painted onto the
 * reference are a second, independent derivation of the same quantity, and the
 * two are not interchangeable: IRMA maps to its own refined contig and recruits
 * reads that BWA soft-clips at the segment termini.
 *
 * Reading it needs a coordinate frame, and a contig does not automatically
 * share the reference's: a contig can be truncated, or a base longer, and a
 * single inserted base shifts every later codon, so a position-1 assumption
 * returns confidently wrong residues rather than failing. The frame is aligned
 * per sample per segment, once, here, rather than re-inferred by every reader.
 *
 * IRMA_results is staged for the coverage tables, not the assemblies. The depth
 * floor must be applied to the alignment that PRODUCED the call, because IRMA's
 * depth and BWA's depth disagree exactly where IRMA is interesting.
 */
process IRMA_POSITION_MAP {
    label 'process_low'
    // Beside the amino-acid tables it shares a coordinate frame with, not under
    // flumut/ — this describes IRMA's consensus, not the marker screen.
    publishDir "${params.outdir}/variant_analysis/irma", mode: params.publish_mode

    input:
    path scripts
    path config
    path reference
    path consensus_dir, stageAs: 'IRMA-consensus-contigs'
    path irma_dirs, stageAs: 'IRMA_results/*'

    output:
    path 'irma_position_map.tsv', emit: map,      optional: true
    path 'irma_consensus_aa.tsv', emit: residues, optional: true
    // IRMA's minority calls, placed and stated against the reference base.
    // Corroboration only — IRMA is the one alignment that does not come from
    // the shared BWA BAM, so it is the only caller whose agreement says
    // anything about alignment error.
    path 'irma_variants.tsv',     emit: variants, optional: true

    script:
    """
    cp ${config} run_config.cfg
    echo "OUTPUT_DIRECTORY=\\"\$PWD\\"" >> run_config.cfg
    echo "REFERENCE_FILE=\\"\$PWD/${reference}\\"" >> run_config.cfg

    # Gracefully skip if IRMA assembled nothing, rather than throwing cryptic errors
    shopt -s nullglob
    fastas=(IRMA-consensus-contigs/*.fasta)
    shopt -u nullglob
    if [ \${#fastas[@]} -eq 0 ]; then
        echo "No IRMA consensus contigs assembled. Skipping IRMA_POSITION_MAP." >&2
        touch irma_position_map.tsv irma_consensus_aa.tsv irma_variants.tsv
        exit 0
    fi

    Rscript ${scripts}/makeGTF.R run_config.cfg

    Rscript ${scripts}/irma_position_map.R \\
        --reference ${reference} \\
        --gtf reference_gtf \\
        --contigs IRMA-consensus-contigs \\
        --irma IRMA_results \\
        --min-depth ${asNum(params.min_depth)} \\
        --out irma_position_map.tsv \\
        --out-aa irma_consensus_aa.tsv \\
        --out-var irma_variants.tsv
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

    if (asBool(params.wfabc) && !params.metadata) {
        error("--wfabc needs --metadata: selection is estimated from allele-frequency\n" +
              "  time series, which are reconstructed by joining variants to the\n" +
              "  individual and time-point columns of the metadata.")
    }

    /*
     * Program dependencies. An enabled program must not depend on a disabled
     * program,
     * and at least one variant caller must run. These are the same requirements
     * the config.cfg "Programs to run" section documents.
     */
    if (!(progOn('run_lofreq') || progOn('run_gatk4') || progOn('ivar'))) {
        error("No variant caller is enabled. Set at least one of LOFREQ, GATK4, or IVAR\n" +
              "  to TRUE — with all three FALSE there is nothing to call.")
    }
    if (progOn('flumut') && !progOn('run_irma')) {
        error("FLUMUT=TRUE needs IRMA=TRUE: FluMut screens the IRMA consensus, so without\n" +
              "  IRMA there is no consensus to screen.")
    }
    if (progOn('flumut_lowfreq') && !(progOn('run_irma') && progOn('run_lofreq'))) {
        error("FLUMUT_LOWFREQ=TRUE needs IRMA=TRUE and LOFREQ=TRUE: it applies LoFreq\n" +
              "  variants to the IRMA consensus before it screens.")
    }
    if (progOn('snpgenie') && !progOn('run_lofreq')) {
        error("SNPGENIE=TRUE needs LOFREQ=TRUE: SNPGenie runs on the pooled LoFreq calls.")
    }
    if (progOn('wfabc') && !(progOn('run_lofreq') || progOn('ivar'))) {
        error("WFABC=TRUE needs LOFREQ=TRUE or IVAR=TRUE: selection is estimated from\n" +
              "  allele-frequency time series, and GATK4 reports genotypes, not fractions.")
    }

    /*
     * Build the sample channel straight from the rename CSV. Nextflow stages
     * reads by symlink, so no raw FASTQ is copied.
     */
    read_dir = file(params.read_directory)

    parsed = channel
        .fromPath(params.rename_file)
        .splitCsv(header: true, strip: true)
        .map { row ->
            // the shipped CSVs carry a UTF-8 BOM, which corrupts the first header key
            def key  = row.keySet().find { k -> k.toString().replace('﻿','') == 'File' }
            def pref = row[key]?.toString()?.trim()
            def name = row.Sample?.toString()?.trim()
            if (!pref || !name) error "Bad row in ${params.rename_file}: ${row}"

            def (r1, r2) = findReadPair(read_dir, pref, name)
            // findReadPair returns [null, reason] for a sample it cannot pair
            if (r1 == null) {
                log.warn "Skipping sample '${name}' (${pref}): ${r2}"
                return [name, pref, null, r2]
            }
            [name, pref, r1, r2]
        }

    /*
     * Samples with no usable reads are dropped, not fatal. A failed library is
     * common, and it should not cost the run every other sample. Each one is
     * recorded in logs/missing_samples.log.
     */
    by_status = parsed.branch { _name, _pref, r1, _info ->
        found:   r1 != null
        missing: r1 == null
    }

    by_status.missing
        .map { name, pref, _r1, reason -> "${name}\t${pref}\t${reason}" }
        .collectFile(
            name:     'missing_samples.log',
            storeDir: "${params.outdir}/logs",
            // No trailing newline: newLine:true supplies the separator, and a
            // seed ending in one leaves a blank line under the header.
            seed:     "# Samples in ${params.rename_file} with no usable read pair\n" +
                      "# Written by Flumina ${workflow.manifest.version}\n" +
                      "Sample\tFile\tReason",
            newLine:  true,
            sort:     true
        )

    samples = by_status.found
        .map { name, _pref, r1, r2 -> tuple(name, r1, r2) }
        .ifEmpty {
            error """No samples had a usable read pair, so there is nothing to process.
  Every row in ${params.rename_file} was skipped — see the warnings above and
  ${params.outdir}/logs/missing_samples.log
  A whole-run failure like this usually means --read_directory points at the
  wrong place, or the 'File' column does not match the actual file names."""
        }

    /*
     * One row per sample. Flumina does not merge lanes, so duplicate Sample
     * values would produce two channel entries with the same name. They would
     * then collide in every publishDir and in GATHER_SAMPLE_VCFS and mix two
     * samples' results, so refuse the input rather than produce a corrupt run.
     */
    samples
        .map { sample, _r1, _r2 -> sample }
        .toList()
        .map { names ->
            def dupes = names.countBy { it }.findAll { _k, v -> v > 1 }.keySet()
            if (dupes) {
                error """Duplicate Sample name(s) in ${params.rename_file}: ${dupes.join(', ')}
  Each Sample must appear once. Flumina does not merge lanes: if one sample was
  sequenced across several lanes or runs, concatenate those fastq files first
  and give the result a single row."""
            }
            names
        }

    ref = PREPARE_REFERENCE(file(params.reference)).index

    trimmed = FASTP(samples).reads

    // Optional IRMA parameter file, staged into every IRMA task as irma_config.sh.
    irma_cfg = params.irma_config
        ? channel.value(file(params.irma_config))
        : channel.value([])

    // Collected IRMA output directories, or an empty list when IRMA is off.
    // An empty list stages nothing, which is how an optional input is expressed.
    irma_dirs = asBool(params.run_irma)
        ? IRMA(trimmed, irma_cfg).results.map { _s, d -> d }.collect()
        : channel.value([])

    ubam     = FASTQ_TO_SAM(trimmed).bam
    reverted = REVERT_SAM(ubam).bam
    grouped  = ADD_READ_GROUPS(reverted).bam
    mapped   = BWA_MAP(grouped, ref).bam
    sorted   = SORT_BAM(mapped).bam
    marked   = MARK_DUPLICATES(sorted).bam
    final_bam = SET_TAGS(marked, ref).bam

    // GATK4 variant calling. The GATK tools used for ALIGNMENT above
    // (FASTQ_TO_SAM, MARK_DUPLICATES, ...) are a separate toolchain and always
    // run; run_gatk4 gates only the caller — HaplotypeCaller through
    // VariantFiltration. Each enabled caller is normalised to (sample, [files]).
    def gatk4_ch = null
    if (asBool(params.run_gatk4)) {
        gvcf     = HAPLOTYPE_CALLER(final_bam, ref).vcf
        geno     = GENOTYPE_GVCF(gvcf, ref).vcf
        selected = SELECT_VARIANTS(geno)
        gatk4_ch = FILTER_VARIANTS(selected.snps.join(selected.indels), ref).vcf
                       .map { sample, snps, indels -> tuple(sample, [snps, indels]) }
    }

    def lofreq_ch = asBool(params.run_lofreq)
        ? LOFREQ(final_bam, ref).vcf.map { sample, vcf -> tuple(sample, [vcf]) }
        : null

    def ivar_ch = asBool(params.ivar)
        ? IVAR(final_bam, ref).tsv.map { sample, tsv -> tuple(sample, [tsv]) }
        : null

    // Per-position depth for every sample. Runs unconditionally: the file
    // carries the caller-visible depth column that makes the depth floor
    // auditable, so it is useful even when the low-frequency screen is off. Sits
    // with the callers because that is what it describes; ref and Scripts are
    // for column 4.
    depth_files = DEPTH_PROFILE(final_bam, ref,
                                file("${projectDir}/Scripts")).depth.collect()

    // The R stage summarizes across all samples, so it must wait for every one.
    // Group each sample's VCFs into a directory named after it, then collect.
    // This is both the completion gate and the vcf_files/<sample>/ layout that
    // convertVCFtoTable.R parses sample names from.
    //
    // Fold every enabled caller into one (sample, [files]) tuple, in GATK4,
    // LoFreq, iVar order. `join`, not `mix`, so a sample must be present in every
    // enabled caller; the callers all derive from final_bam, so their sample
    // sets match. GATHER_SAMPLE_VCFS copies whatever it is handed, and
    // convertVCFtoTable.R globs each caller's file by name and tolerates any
    // being absent, so a disabled caller drops its rows. Validation above refuses
    // the no-caller case, so this list is never empty.
    def caller_chs = [gatk4_ch, lofreq_ch, ivar_ch].findAll { it != null }
    per_sample = caller_chs.inject(null) { acc, ch ->
        acc == null ? ch : acc.join(ch).map { sample, a, b -> tuple(sample, a + b) }
    }

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
    if (asBool(params.run_irma) && asBool(params.flumut)) {
        FLUMUT(file("${projectDir}/Scripts"), r.irma, file(params.reference))
    }

    // Where FluMut's numbering lands on this reference. Depends only on the
    // reference and the marker database, not on the samples, so it runs once for
    // whichever screens are enabled rather than per screen.
    if (asBool(params.run_irma) && (asBool(params.flumut) || asBool(params.flumut_lowfreq))) {
        FLUMUT_POSITION_MAP(file("${projectDir}/Scripts"), r.config, file(params.reference))
    }

    // Where IRMA's consensus lands on this reference. Not gated on the FluMut
    // switches — it describes IRMA's assembly, which is worth placing whether or
    // not markers are being screened. Needs the contigs for the sequence and
    // IRMA_results for the coverage the depth floor is read from.
    if (asBool(params.run_irma)) {
        IRMA_POSITION_MAP(file("${projectDir}/Scripts"), r.config,
                          file(params.reference), r.irma, irma_dirs)
    }

    // Screen low-frequency variants above threshold for H5N1 markers.
    // Applies LoFreq variants to IRMA consensus sequences, creating
    // mutated pseudo-consensus for FluMut marker screening.
    if (asBool(params.run_irma) && asBool(params.flumut_lowfreq)) {
        // depth_files is computed unconditionally above; this screen consumes it.
        FLUMUT_LOWFREQ(file("${projectDir}/Scripts"), r.irma, vcf_dirs,
                       file(params.reference), depth_files)
    }

    // Optional population-genetics analyses. Both read the config.cfg written by
    // R_ANALYSIS, which is also what sequences them after it.
    if (asBool(params.snpgenie)) {
        SNPGENIE(
            file("${projectDir}/Scripts"),
            r.config,
            vcf_dirs,
            ref.map { r_fa, _idx, _dict -> r_fa },
            metadata_ch
        )
    }

    if (asBool(params.wfabc)) {
        WFABC(
            file("${projectDir}/Scripts"),
            r.config,
            vcf_dirs,
            ref.map { r_fa, _idx, _dict -> r_fa },
            metadata_ch
        )
    }
}
