samples = ["Undetermined_S0_L001_R1_001", 
           "Undetermined_S0_L001_R2_001", 
           "Undetermined_S0_L002_R1_002",
           "Undetermined_S0_L001_R2_002"]

if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/regions.bed"
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {input} > {output}"


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions) +
            config["params"]["gatk"]["HaplotypeCaller"])


rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions = []
#        regions="called/regions.bed" if config["processing"].get("restrict-regions") else []
    output:
        gvcf=protected("called/{sample}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra=get_call_variants_params
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref=config["ref"]["genome"],
#        gvcfs="/data/athersh/dna-seq-gatk-variant-calling/00-All.vcf.gz"
        gvcfs=expand("called/{sample}.g.vcf.gz", sample=samples)
    output:
        gvcf="called/all.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
#        gvcf=expand("/data/athersh/dna-seq/gatk-variant-calling/{sample}.g.vcf.gz",sample=samples)
        gvcf="called/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"

