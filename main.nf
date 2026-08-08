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

/*
 * Locate the R1/R2 pair for one prefix from the rename CSV.
 *
 * Sequencing cores hand back reads in more shapes than one glob can cover:
 * nested per-sample or per-project directories from bcl2fastq, and several
 * mate-marker conventions. The original pattern here was a single
 * non-recursive `<prefix>*_R1_*.fastq.gz`, which silently found nothing for
 * anything nested — the commonest layout of the lot.
 *
 * Files are gathered recursively, then the mate markers are tried in order of
 * decreasing specificity, so `_R1_` wins over a bare `_1` and a file like
 * SAMPLE_S1_L001_R1_001.fastq.gz cannot be mistaken for its own mate.
 */
def findReadPair(read_dir, prefix, sample_name) {
    // Uncompressed fastq is accepted as well as gzipped: the Snakemake-era
    // organizeReads.R took fastq/fq/fastq.gz/fq.gz, and dropping that would
    // quietly break anyone whose reads are not compressed.
    def exts = ['fastq.gz', 'fq.gz', 'fastq', 'fq']

    // Tried in order. The Sample column is a fallback because reads are often
    // already renamed by the time the pipeline is re-run over them, in which
    // case the File value no longer appears in any filename — another
    // behaviour inherited from organizeReads.R.
    def gather = { String stem ->
        def out = []
        exts.each { ext ->
            // Both globs are needed: `**` does not match zero directories, so
            // reads directly in read_dir are missed by the recursive form alone.
            ["${read_dir}/${stem}*.${ext}", "${read_dir}/**/${stem}*.${ext}"].each { pattern ->
                def found = file(pattern)
                if (found) out += (found instanceof List ? found : [found])
            }
        }
        out.unique { it.toString() }.sort { it.name }
    }

    // findResult rather than a for/break loop: Nextflow's strict parser (25.x+)
    // rejects `for` loops in pipeline scripts outright.
    def hits = [prefix, sample_name].findAll { it }
                                    .findResult { stem -> gather.call(stem) ?: null } ?: []

    // Per-sample problems return a reason instead of throwing. A library that
    // failed to sequence is routine, and losing a whole run's worth of good
    // samples to one bad row helps nobody — they are reported and skipped.
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

    // Last resort, and what organizeReads.R did for every sample: with exactly
    // two files and no recognised marker, take them in sorted order. Every
    // convention in use puts the first mate first alphabetically, so this is
    // usually right — but it is a guess, so say so rather than let a silent
    // mate swap through.
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
 * Read a boolean parameter without trusting its type.
 *
 * Nextflow used to coerce a command-line `--flag false` to a real boolean by
 * matching the type of the config default. From 25.x it does not: the value
 * arrives as the STRING "false", and every non-empty string is truthy in
 * Groovy. So `if (params.wfabc)` fires when wfabc was explicitly disabled, and
 * worse, `--run_irma false` would silently RUN IRMA — a wrong result rather
 * than an error. Tested directly: 24.10.4 coerces, 26.04.3 does not.
 */
def asBool(value) {
    if (value instanceof Boolean) return value
    return value?.toString()?.trim()?.toLowerCase() in ['true', 't', 'yes', '1']
}

/*
 * The numeric sibling of asBool, and it exists for the same reason — the same
 * 25.x change, one type over.
 *
 * A param given on the command line arrives as a STRING. For booleans that gave
 * the trap above. For numbers it is worse than truthy, because Groovy defines
 * String * Integer as REPETITION: "0.01" * 100 is not 1, it is "0.01" written
 * out one hundred times, and .toInteger() then throws on the 400-character
 * result. That is exactly how FLUMUT_LOWFREQ died on the swine WGS run, and
 * because the launcher ALWAYS passes --flumut_freq_threshold, it would have
 * died the same way for every user on every run.
 *
 * Never do arithmetic on a params value directly. Put it through here.
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
 * The Snakemake version used five separate rules and a whole second snakefile
 * invocation for this. All four indexing steps are seconds-long and always run
 * together, so splitting them only adds scheduling overhead.
 * ========================================================================== */
process PREPARE_REFERENCE {
    tag   "reference"
    label 'process_low'
    publishDir "${params.outdir}/Reference", mode: params.publish_mode
    /* findAAChanges.R reads ${OUTPUT_DIRECTORY}/reference.fa — the old bash
     * driver copied the reference to the output root before calling snakemake,
     * so the R scripts depend on it being there. Reproduce that placement.
     *
     * saveAs returns null — i.e. publish nothing — when the reference the user
     * pointed at IS that file already. Without the guard this process
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
 * IRMA parameters (TMP, SINGLE_LOCAL_PROC, ...) come from an optional user file
 * passed with `--external-config`.
 *
 * IRMA 1.0.3 used to source `./irma_config.sh` from its working directory
 * implicitly; 1.3.5 removed that in favour of the explicit flag, so staging the
 * file alone is no longer enough — it has to be named on the command line.
 *
 * When --irma_config is not set, an empty list is staged, which stages nothing,
 * the flag is omitted, and IRMA falls back to its module defaults.
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
    def cfg_arg = irma_cfg ? "--external-config run_irma_config.sh" : ""
    """
    # IRMA 1.3.5 sizes itself with `irma-core num-procs --cap-cores-using-env`,
    # which reads the CPU affinity mask. Under Slurm that reflects
    # --cpus-per-task and is usually right, but affinity is not set everywhere —
    # PBS without cgroups, or Docker's CPU quota, leave it looking like the whole
    # machine, and IRMA would then oversubscribe the node. Stating the number
    # Nextflow actually asked the scheduler for removes the guesswork.
    export LOCAL_PROCS_OVERRIDE=${task.cpus}

    # TMP in an IRMA config must be an ABSOLUTE path that already exists.
    #
    # IRMA builds its working path straight from it —
    #     ppath="\$TMP"/<user>/IRMAv<version>/<run>-<token>
    # — never creates it, and changes directory as it works, so a relative TMP
    # stops resolving partway through. Either mistake produces the same silent
    # failure: reads are counted, then the match stage emits nothing, every
    # table and consensus is empty, and IRMA still exits 0. Measured on one
    # sample: relative TMP gave 0 segments, absolute gave all 8.
    #
    # The staged config is a symlink into another directory, so rewrite a copy.
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

    # IRMA returns 0 even when it has failed outright: it creates its whole
    # output skeleton (amended_consensus/, tables/, logs/ ...) and exits
    # successfully with every directory empty. Nextflow then sees a green task,
    # organizeIRMA.R finds no fasta to collect, and FluMut is skipped silently
    # for want of input. The first anyone knows is a missing result days later.
    #
    # So check the one thing that matters — that a consensus was actually
    # produced. Not fatal on its own, because a sample with too few flu reads
    # legitimately assembles nothing; but if this fires for EVERY sample the
    # cause is IRMA itself, and .command.err will say so.
    if ! ls ${sample}/*.fasta >/dev/null 2>&1; then
        echo "WARNING: IRMA produced no consensus sequence for ${sample}." >&2
        echo "         A sample with too few influenza reads can do this legitimately." >&2
        echo "         If it happens for every sample, IRMA itself failed — check this" >&2
        echo "         task's .command.err for repeated 'exec failed' lines." >&2
    fi
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
    //
    // VariantFiltration flags a record when the expression is TRUE, so every
    // expression below names the condition for a BAD record.
    //
    // FS and SOR were INVERTED until 2026-08-08 — written "FS<60.0" and
    // "SOR<3.0" when it is HIGH values that indicate strand bias. The two
    // filters therefore fired on precisely the records with clean strand
    // statistics, which is why not one record passed on either dataset here:
    // 0 of 100 on the swine run and 0 of 1,442 on the avian test data. The 96
    // records flagged FS;SOR had FS median 0.693 (max 11.6) and SOR median
    // 0.776 — far below the thresholds they were supposedly failing — and the
    // single record carrying real strand bias (SOR=3.212) was the one the SOR
    // filter let through. The other five expressions are "less than" in GATK's
    // own recommendations too, and were always right.
    //
    // Thresholds are GATK's published germline hard-filter recommendations per
    // variant type. Indels legitimately take a looser FS and ReadPosRankSum
    // than SNPs, which the indel block had not reflected.
    //
    // GATK's thresholds are tuned for DIPLOID GERMLINE data at ~30x, and almost
    // none of those assumptions hold here — haploid calling, ~1,770x median
    // depth, balanced strand, and eight short segments that map uniquely.
    // Measured over 1,542 real records (100 swine + 1,442 avian), FOUR of the
    // seven can never fire on data of this shape:
    //
    //   FS>60             0 records   max observed FS  24.5
    //   MQ<40             0 records   min observed MQ  52
    //   MQRankSum<-12.5   0 records   min observed     -0.93
    //   ReadPosRankSum<-8 0 records   min observed     -3.36
    //
    // They are kept anyway: inert costs nothing, they preserve GATK's own
    // convention, and they would start doing work on a differently prepared
    // library. The RankSums additionally have no value on two thirds of records
    // by construction — a haploid hom-alt call has no reference reads to rank
    // against — and a filter on a mostly-absent annotation is not a filter.
    //
    // QD is the one number tuned to this data rather than inherited. GATK's
    // QD<2.0 flagged 3 of 1,542. The distribution is sharply bimodal: 1,446 sit
    // at QD>=25 and the tail is sparse (12 records in 10-15, 8 in 15-20), so
    // QD<15.0 lands in genuinely empty space and flags 38 (2.5%). It is a real
    // quality signal here and not a depth artefact — the concern was that GATK's
    // QUAL saturates at extreme depth and would depress QD for the DEEPEST
    // sites, but the relationship runs the other way: mean depth is 252 for
    // QD<15 against 1,849 for QD>=25, so low QD selects THIN calls, which is
    // what it should do.
    //
    // Indels keep GATK's QD<2.0. The analysis above is SNP-only and the indel
    // VCF is not consumed by the variant table, so there is no basis here for
    // moving it.
    //
    // No truth set exists for any of this, so these are distribution-based
    // judgements. That is also why nothing is DROPPED on the strength of them:
    // the column annotates, and FluLens is where a reader chooses what to hide.
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

/*
 * iVar — the third caller.
 *
 * Deliberately NOT given a GFF. `ivar variants -g` annotates amino acids by
 * translating from the start of each reference sequence in frame 1, which is
 * the exact `ceiling(POS/3)` mistake convertVCFtoTable.R was fixed for: it is
 * right for the eight primary ORFs and meaningless for M2, NEP, PA-X and
 * PB1-F2, and it does not fail — it returns a plausible residue for a codon
 * that does not exist. Flumina annotates through Scripts/fluORFs.R, which walks
 * the real coding intervals, so iVar contributes NUCLEOTIDE calls only and the
 * amino-acid columns come from the same place they do for every other caller.
 *
 * mpileup flags, and why each one is not a default:
 *   -aa       every position, including zero-coverage ones, so iVar's own
 *             depth filter is what removes them rather than mpileup's silence
 *   -A        count anomalous read pairs. LoFreq's DP4 was pinned to REQUIRE
 *             proper pairs; iVar is left inclusive here because its PASS test
 *             is its own, and the two callers are more useful when they are not
 *             tuned to agree by construction
 *   -B        no BAQ. On by default and it re-scores base qualities around
 *             indels, which suppresses real low-frequency SNPs next to them
 *   -d 0      no depth cap. The default 8000 would silently truncate the deep
 *             libraries here — MC-495 runs to 1.4 million matched reads
 *   -Q 0      let iVar apply the quality threshold via -q, rather than
 *             filtering twice at two different values
 *
 * Thresholds come from the same params LoFreq and the R stage use, so the three
 * callers are held to one set of numbers.
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

    # ivar exits 0 having written only a header when nothing passes, and an
    # absent file would fail the output declaration instead of publishing an
    # honest empty result. Both are valid outcomes for a thin library.
    if [ ! -s ivar-called-variants.tsv ]; then
      printf 'REGION\\tPOS\\tREF\\tALT\\tREF_DP\\tREF_RV\\tREF_QUAL\\tALT_DP\\tALT_RV\\tALT_QUAL\\tALT_FREQ\\tTOTAL_DP\\tPVAL\\tPASS\\n' \\
        > ivar-called-variants.tsv
    fi
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
 * Per-position depth in REFERENCE coordinates, for the low-frequency FluMut
 * screen only.
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
 * -a emits EVERY position including zero-coverage ones, which is the entire
 * point — without it the uncovered positions are simply absent and
 * indistinguishable from a truncated file. -Q 0 because the mask is asking
 * "was there any read here at all", not "was there a confident read".
 *
 * The whole reference is 13,133 bp, so this is ~13k lines per sample and costs
 * nothing to keep.
 */
process DEPTH_PROFILE {
    tag "$sample"
    label 'process_low'
    // Published, because the consumer that matters is FluLens rather than this
    // pipeline: it is the only component that maps FluMut's own numbering back
    // to reference positions, so a per-marker "was there any coverage here"
    // check can only be made there.
    publishDir "${params.outdir}/depth_profiles", mode: params.publish_mode

    input:  tuple val(sample), path(bam), path(bai)
    output: path "${sample}.depth", emit: depth

    script:
    """
    samtools depth -a -Q 0 ${bam} > '${sample}.depth'
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
    def irma_step = asBool(params.run_irma)
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
DISABLE_IRMA="${asBool(params.run_irma) ? 'FALSE' : 'TRUE'}"
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
    // The subtype decision, so consumers read it rather than each re-deriving
    // it from segment names — see filter_flumut_subtype.R for why that matters
    // now that the rule has two sources.
    path 'subtype.tsv',            emit: subtype,    optional: true
    path 'flumut_version.txt',     emit: version

    script:
    /*
     * Screening the REFERENCE too, then removing its own findings.
     *
     * FluMut reports every marker a sequence carries, and the reference
     * carries plenty. Those appear in every sample by construction and say
     * nothing about any of them. Measured on the swine WGS run: the reference
     * alone accounts for 84 marker rows, and across 30 samples 2,514 of 2,515
     * reported rows were identical to it. One row in 2,515 carried
     * information, and it was buried under the other 2,514.
     *
     * markers.tsv drops the rows the reference also has. mutations.tsv is a
     * wide matrix and instead drops only the columns where EVERY sample
     * carries the reference residue — dropping by "the reference has this
     * marker" would throw away reversions, and a sample that LOSES a reference
     * marker emits no marker row at all, so that table is the only place the
     * signal exists.
     *
     * Nothing is discarded: raw output is kept as *_all.tsv and the
     * reference's own calls as reference_*.tsv, so a removal can be explained.
     */
    def subtract = asBool(params.flumut_subtract_reference)
    def keep_hana = asBool(params.flumut_keep_mismatched_ha_na) ? 'TRUE' : 'FALSE'
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

    # Runs regardless of the subtraction above, because it answers a different
    # question: HA/NA markers are numbered for H5/N1 specifically, so off-subtype
    # they are read against the wrong ruler AND the proteins have diverged too
    # far for equivalence to be assumed. Internal-gene markers are unaffected.
    # IRMA-consensus-contigs is passed as a SECOND subtype source. IRMA writes the
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
    # onto the REFERENCE, not the consensus, because LoFreq's coordinates are
    # the reference's and IRMA's consensus is built de novo. See the script header.
    # depth/ carries one <sample>.depth per sample from DEPTH_PROFILE. Positions
    # with zero coverage are written as N rather than taking a reference base —
    # see the script header. The directory is passed rather than the files so an
    # absent depth source degrades to the old behaviour instead of failing.
    Rscript ${scripts}/apply_lofreq_to_consensus.R mutated.fasta ${params.flumut_freq_threshold} \\
        ${reference} depth \$r_args

    if [ ! -s mutated.fasta ]; then
        echo "no low-frequency variants (AF >= ${freq_pct}%) found — skipping flumut" >&2
        touch markers.tsv mutations.tsv literature.tsv
        exit 0
    fi

    # No rename step here: apply_lofreq_to_consensus.R already writes
    # >sample_SEGMENT headers. rename_for_flumut.R takes the sample name from
    # the FILENAME, so handing it one combined FASTA renamed every sample to
    # "mutated" and lost sample identity completely.
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

    # Same reference subtraction as FLUMUT — see the comment there for why.
    # It matters more here, not less: applying low-frequency variants only adds
    # findings on top of the reference background, so without this the novel
    # calls are buried even deeper.
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

    # Runs regardless of the subtraction above, because it answers a different
    # question: HA/NA markers are numbered for H5/N1 specifically, so off-subtype
    # they are read against the wrong ruler AND the proteins have diverged too
    # far for equivalence to be assumed. Internal-gene markers are unaffected.
    # IRMA-consensus-contigs is passed as a SECOND subtype source. IRMA writes the
    # subtype it assigned into each consensus header (>A_HA_H5, >A_NA_N1), so a
    # reference with bare A_HA / A_NA names — which is the repo's own reference,
    # and the bundled H5N1 test_dataset — can still be resolved instead of being
    # treated as unconfirmable. The reference still wins where it states one.
    Rscript ${scripts}/filter_flumut_subtype.R ${reference} markers.tsv ${keep_hana} . IRMA-consensus-contigs
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
     * Build the sample channel straight from the rename CSV.
     *
     * This replaces organizeReads.R, which physically COPIED every raw fastq
     * into organized-reads/ — doubling disk usage before analysis even started.
     * Nextflow stages by symlink instead, so nothing is duplicated.
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
     * Samples with no usable reads are dropped rather than fatal. A library
     * that failed to sequence is an ordinary occurrence, and there is no reason
     * for it to cost the run every other sample. They are recorded in
     * logs/missing_samples.log so the omission is on paper rather than only in
     * the console scrollback.
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
     * One row per sample. organizeReads.R allowed several rows for the same
     * Sample and wrote them out as _L001_, _L002_ ... lanes; nothing here merges
     * lanes, so duplicate Sample values would instead produce two channel
     * entries with the same name, which then collide in every publishDir and in
     * GATHER_SAMPLE_VCFS. That silently mixes two samples' results, so refuse
     * the input rather than produce a corrupt run.
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

    gvcf     = HAPLOTYPE_CALLER(final_bam, ref).vcf
    geno     = GENOTYPE_GVCF(gvcf, ref).vcf
    selected = SELECT_VARIANTS(geno)
    filtered = FILTER_VARIANTS(selected.snps.join(selected.indels), ref).vcf
    lofreq   = LOFREQ(final_bam, ref).vcf

    // The R stage summarises across ALL samples, so it must wait for every one.
    // Group each sample's VCFs into a directory named after it, then collect —
    // this is both the completion gate and the vcf_files/<sample>/ layout that
    // convertVCFtoTable.R parses sample names out of.
    //
    // Two shapes rather than one join made conditional. `join` against an empty
    // channel emits NOTHING, so folding an optional iVar into a single join
    // would not degrade to "no iVar rows" — it would drop every sample from the
    // R stage and produce an empty run that still exits green.
    if (asBool(params.ivar)) {
        per_sample = filtered
            .join(lofreq)
            .join(IVAR(final_bam, ref).tsv)
            .map { sample, snps, indels, lf, iv -> tuple(sample, [snps, indels, lf, iv]) }
    } else {
        per_sample = filtered
            .join(lofreq)
            .map { sample, snps, indels, lf -> tuple(sample, [snps, indels, lf]) }
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

    // Screen low-frequency variants above threshold for H5N1 markers.
    // Applies LoFreq variants to IRMA consensus sequences, creating
    // mutated pseudo-consensus for FluMut marker screening.
    if (asBool(params.run_irma) && asBool(params.flumut_lowfreq)) {
        // Depth is computed only when this screen actually runs — it is the one
        // consumer, and 143 extra tasks are not worth paying for otherwise.
        depth_files = DEPTH_PROFILE(final_bam).depth.collect()
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
