import os

#Expected files at the end
rule all:
    input:
        "Reference/reference.fa",
        "Reference/reference.dict",
        "Reference/reference.fa.dict",
        "Reference/reference.fa.amb",
        "Reference/reference.fa.ann",
        "Reference/reference.fa.bwt",
        "Reference/reference.fa.pac",
        "Reference/reference.fa.sa",
        "Reference/reference.fa.fai"

#Copies the reference to a new directory to make a bunch of reference files
isExist = os.path.exists("Reference")

if not isExist:
    os.mkdir("Reference")

rule reference_copy:
    input:
        "reference.fa" 
    output:
        "Reference/reference.fa"   
    shell:
        """
        cp reference.fa Reference/reference.fa
        """

#Indexes reference
rule bwa_index:
    input:
        "Reference/reference.fa" 
    output:
        "Reference/reference.fa.amb",
        "Reference/reference.fa.ann",
        "Reference/reference.fa.bwt",
        "Reference/reference.fa.pac",
        "Reference/reference.fa.sa"        
    shell:
        "bwa index -a bwtsw {input}"

#Creates samtools index file from reference
rule samtools_faidx:
    input:
        "Reference/reference.fa" 
    output:
        "Reference/reference.fa.fai"     
    shell:
        "samtools faidx {input}"

#Creates gatk4 sequence dictionary from reference
rule gatk_dictionary:
    input:
        "Reference/reference.fa" 
    output:
        "Reference/reference.dict"     
    shell:
        "gatk CreateSequenceDictionary --REFERENCE {input} --OUTPUT {output}"
        " --USE_JDK_DEFLATER true --USE_JDK_INFLATER true"

#Copies the dict file with a new name
rule dict_copy:
    input:
        "Reference/reference.dict" 
    output:
        "Reference/reference.fa.dict"   
    shell:
        """
        cp Reference/reference.dict Reference/reference.fa.dict
        """
