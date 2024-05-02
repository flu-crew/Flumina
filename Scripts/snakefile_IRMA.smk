# #Read the list of sample names from the directories in the input directory
sample_names = [d for d in os.listdir("processed-reads") if os.path.isdir(os.path.join("processed-reads", d))]
print(sample_names)

rule all:
    input:
        expand("processed-reads/{sample}/{sample}_R1.fastq.gz", sample=sample_names),
        expand("processed-reads/{sample}/{sample}_R2.fastq.gz", sample=sample_names),
	expand("IRMA_results/{sample}/", sample=sample_names)
        
#Process reads with fastp
rule fastp_process:
    input:
        R1="organized-reads/{sample}/{sample}_L001_READ1.fastq.gz",
        R2="organized-reads/{sample}/{sample}_L001_READ2.fastq.gz"
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
        cluster_params="-p priority --qos=vpru -A nadc_iav -D /project/nadc_iav/hpai-2024-variants --mem 40G --cpus-per-task=2"
    output:
        directory("IRMA_results/{sample}/")     
    shell:
        "IRMA FLU {input.R1} {input.R2} {output}"

