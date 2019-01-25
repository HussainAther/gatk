from snakemake import shell
shell("module load GATK")

samples = ["Undetermined_S0_L001_R1_001", 
           "Undetermined_S0_L001_R2_001", 
           "Undetermined_S0_L002_R1_001",
           "Undetermined_S0_L002_R2_001"]

samples_fastq_gz = ["Undetermined_S0_L001_R1_001.fastq.gz", 
           "Undetermined_S0_L001_R2_001.fastq.gz", 
           "Undetermined_S0_L002_R1_001.fastq.gz",
           "Undetermined_S0_L002_R2_001.fastq.gz"]


rule map_reads:
    input:
        reads="{sample}.fastq.gz"
    output:
        temp("mapped/{sample}.sorted.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.27.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        "mapped/{sample}.sorted.bam"
    output:
        bam="dedup/{sample}.bam",
        metrics="qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    wrapper:
        "0.26.1/bio/picard/markduplicates"


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "mapped/{sample}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "dedup/{sample}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f

rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
    output:
        bam=protected("recal/{sample}.bam")
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    wrapper:
        "0.27.1/bio/samtools/index"
