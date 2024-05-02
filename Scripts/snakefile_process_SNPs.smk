# #Read the list of sample names from the directories in the input directory
sample_names = [d for d in os.listdir("organized-reads") if os.path.isdir(os.path.join("organized-reads", d))]
print(sample_names)

CLUSTER_JOBS = os.getenv('CLUSTER_JOBS', 'FALSE')

rule all:
    input:
        expand("processed-reads/{sample}/{sample}_R1.fastq.gz", sample=sample_names),
        expand("processed-reads/{sample}/{sample}_R2.fastq.gz", sample=sample_names),
        expand("BAM_files/{sample}/fastqsam.bam", sample = sample_names),
        expand("BAM_files/{sample}/revertsam.bam", sample = sample_names),
        expand("BAM_files/{sample}/all_reads.bam", sample = sample_names),
        expand("BAM_files/{sample}/mapped_reads_all.bam", sample = sample_names),
        expand("BAM_files/{sample}/mapped_reads_sort.bam", sample = sample_names),
        expand("BAM_files/{sample}/mapped_reads_md.bam", sample = sample_names),
        expand("logs/{sample}/duplicate_metrics.txt" , sample = sample_names),
        expand("BAM_files/{sample}/final_mapped_reads.bam", sample = sample_names),
        expand("BAM_files/{sample}/haplotype_caller.bam", sample = sample_names),
        expand("vcf_files/{sample}/gatk4-haplotype-caller.vcf", sample = sample_names),
        expand("vcf_files/{sample}/gatk4-unfiltered-genotypes.vcf", sample = sample_names),
        expand("vcf_files/{sample}/gatk4-unfiltered-snps.vcf", sample = sample_names),
        expand("vcf_files/{sample}/gatk4-unfiltered-indels.vcf", sample = sample_names),
        expand("vcf_files/{sample}/gatk4-filtered-snps.vcf", sample = sample_names),
        expand("vcf_files/{sample}/lofreq-called-variants.vcf", sample = sample_names),
        expand("IRMA_results/{sample}/", sample = sample_names)

        
#Process reads with fastp
rule fastp_process:
    input:
        R1="organized-reads/{sample}/{sample}_L001_READ1.fastq.gz",
        R2="organized-reads/{sample}/{sample}_L001_READ2.fastq.gz"
    params:
        cluster=CLUSTER_JOBS
    output:
        "processed-reads/{sample}/{sample}_R1.fastq.gz",
        "processed-reads/{sample}/{sample}_R2.fastq.gz"
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2}"
        " --out1 {output[0]} --out2 {output[1]}"
        " --length_required 60 --low_complexity_filter --complexity_threshold 30"
        " --trim_poly_x --correction --detect_adapter_for_pe"
        " --dedup dup_calc_accuracy 5"
        " --html logs/fastp-{wildcards.sample}.html --json logs/fastp-{wildcards.sample}.json --compression 8"
        " report_title {wildcards.sample}"

#Converts fastq to sam/bam
rule run_irma:
    input:
        R1="processed-reads/{sample}/{sample}_R1.fastq.gz",
        R2="processed-reads/{sample}/{sample}_R2.fastq.gz"
    params:
        cluster=CLUSTER_JOBS
    output:
        directory("IRMA_results/{sample}/")     
    shell:
        "IRMA FLU {input.R1} {input.R2} {output}"

#Converts fastq to sam/bam
rule fastq_to_sam:
    input:
        R1="processed-reads/{sample}/{sample}_R1.fastq.gz",
        R2="processed-reads/{sample}/{sample}_R2.fastq.gz"
    output:
        "BAM_files/{sample}/fastqsam.bam"     
    shell:
        "gatk FastqToSam -FASTQ {input.R1} -FASTQ2 {input.R2}"
        " -OUTPUT {output} -SAMPLE_NAME {wildcards.sample}"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#Revert sam stuff to clean bam file
rule revert_sam:
    input:
        "BAM_files/{sample}/fastqsam.bam"     
    output:
        "BAM_files/{sample}/revertsam.bam"     
    shell:
        "gatk RevertSam -I {input} -O {output}"
        " -SANITIZE true -MAX_DISCARD_FRACTION 0.005"
        " -ATTRIBUTE_TO_CLEAR XT -ATTRIBUTE_TO_CLEAR XN -ATTRIBUTE_TO_CLEAR AS"
        " -ATTRIBUTE_TO_CLEAR OP -SORT_ORDER queryname"
        " -RESTORE_ORIGINAL_QUALITIES true -REMOVE_DUPLICATE_INFORMATION true"
        " -REMOVE_ALIGNMENT_INFORMATION true"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#Converts fastq to sam/bam
rule replace_read_group:
    input:
        "BAM_files/{sample}/revertsam.bam"     
    output:
        "BAM_files/{sample}/all_reads.bam"     
    shell:
        "gatk AddOrReplaceReadGroups -I {input} -O {output}"
        " -RGSM {wildcards.sample} -RGPU FLOWCELL1.LANE1 -RGID FLOWCELL1.LANE1"
        " -RGLB LIB-{wildcards.sample} -RGPL ILLUMINA"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#sam to fastq, bam mapping, merge bam
rule bwa_map:
    input:
        reads="BAM_files/{sample}/all_reads.bam",    
        reference="Reference/reference.fa"
    params:
        cluster=CLUSTER_JOBS
    output:
        "BAM_files/{sample}/mapped_reads_all.bam"     
    shell:
        "gatk SamToFastq -I {input.reads} -FASTQ /dev/stdout"
        " -CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true |"
        " bwa mem -M -p {input.reference} /dev/stdin |"
        " gatk MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM {input.reads}"
        " -OUTPUT {output} -R {input.reference} -CREATE_INDEX true -ADD_MATE_CIGAR true"
        " -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true"
        " -MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#sorts bam file
rule sort_bam:
    input:
        "BAM_files/{sample}/mapped_reads_all.bam"     
    output:
        "BAM_files/{sample}/mapped_reads_sort.bam"     
    shell:
        "gatk SortSam -INPUT {input} -OUTPUT {output}"
        " -CREATE_INDEX true -SORT_ORDER coordinate"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#marks duplicate reads
rule mark_duplicates:
    input:
        "BAM_files/{sample}/mapped_reads_sort.bam"     
    output:
        #directory("logs"),
        reads="BAM_files/{sample}/mapped_reads_md.bam",
        metrics="logs/{sample}/duplicate_metrics.txt"      
    shell:
        "gatk MarkDuplicates -INPUT {input} -OUTPUT {output.reads}"
        " -CREATE_INDEX true -METRICS_FILE {output.metrics}"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#sorts bam file and sets some tags to make final set of mapped reds
rule set_tags:
    input:
        reads="BAM_files/{sample}/mapped_reads_md.bam",
        reference="Reference/reference.fa" 
    output:
        "BAM_files/{sample}/final_mapped_reads.bam"
    shell:
        "gatk SortSam -INPUT {input.reads}"
        " -OUTPUT /dev/stdout -SORT_ORDER coordinate |"
        " gatk SetNmAndUqTags -INPUT /dev/stdin -OUTPUT {output}"
        " -CREATE_INDEX true -R {input.reference}"
        " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"

#calls haplotypes on each same
rule haplotype_caller:
    input:
        reads="BAM_files/{sample}/final_mapped_reads.bam",
        reference="Reference/reference.fa"
    params:
        cluster=CLUSTER_JOBS
    output:
        bam="BAM_files/{sample}/haplotype_caller.bam",
        vcf="vcf_files/{sample}/gatk4-haplotype-caller.vcf"
    shell:
        "gatk HaplotypeCaller -I {input.reads}"
        " -R {input.reference} -O {output.vcf}"
        " -ERC GVCF -ploidy 1"
        " -bamout {output.bam}"

#join genotypes haplotype called files
rule genotyping:
    input:
        vcf="vcf_files/{sample}/gatk4-haplotype-caller.vcf",
        reference="Reference/reference.fa"
    output:
        "vcf_files/{sample}/gatk4-unfiltered-genotypes.vcf"
    shell:
        "gatk GenotypeGVCFs -V {input.vcf}"
        " -R {input.reference} -O {output}"
        " --use-new-qual-calculator true"

# Selects only the SNPs from the VCF
rule select_SNP:
    input:
        "vcf_files/{sample}/gatk4-unfiltered-genotypes.vcf"
    output:
        "vcf_files/{sample}/gatk4-unfiltered-snps.vcf"
    shell:
        "gatk SelectVariants -V {input}"
        " -O {output} --select-type SNP"

# Selects only the INDELS from the VCF
rule select_INDEL:
    input:
        "vcf_files/{sample}/gatk4-unfiltered-genotypes.vcf"
    output:
        "vcf_files/{sample}/gatk4-unfiltered-indels.vcf"
    shell:
        "gatk SelectVariants -V {input}"
        " -O {output} --select-type INDEL"

# Selects only the SNPs from the VCF
rule filter_SNP:
    input:
        vcf="vcf_files/{sample}/gatk4-unfiltered-snps.vcf",
        reference="Reference/reference.fa"
    output:
        "vcf_files/{sample}/gatk4-filtered-snps.vcf"
    shell:
        "gatk VariantFiltration"
        " -R {input.reference}" 
        " -V {input.vcf}"
        " -O {output}"
        " -filter \"QUAL<30.0\" --filter-name \"QUAL\""
        " -filter \"QD<2.0\" --filter-name \"QD\""
        " -filter \"SOR<3.0\" --filter-name \"SOR\""
        " -filter \"FS<60.0\" --filter-name \"FS\""
        " -filter \"MQ<40.0\" --filter-name \"MQ\""
        " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum\""
        " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum\""
        

# Selects only the SNPs from the VCF
rule filter_INDEL:
    input:
        vcf="vcf_files/{sample}/gatk4-unfiltered-indels.vcf",
        reference="Reference/reference.fa"
    output:
        "vcf_files/{sample}/gatk4-filtered-indels.vcf"
    shell:
        "gatk VariantFiltration"
        " -R {input.reference}" 
        " -V {input.vcf}"
        " -O {output}"
        " -filter \"QD<2.0\" --filter-name \"QD\""
        " -filter \"QUAL<30.0\" --filter-name \"QUAL\""
        " -filter \"FS<60.0\" --filter-name \"FS\""
        " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum\""

# Runs LoFreq for low frequency variants 
rule run_loFreq:
    input:
        reads="BAM_files/{sample}/final_mapped_reads.bam",
        reference="Reference/reference.fa"
    params:
        cluster=CLUSTER_JOBS
    output:
        "vcf_files/{sample}/lofreq-called-variants.vcf"
    shell:
        "lofreq call -f {input.reference}"
        " -o {output} {input.reads}"
