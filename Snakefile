import pandas as pd

from snakemake.utils import report as reporter

shell("module load GATK")

report: "report/workflow.rst"

rule all:
    input:
        ["annotated/all.vcf.gz",
        "qc/multiqc.html",
        "plots/depths.svg",
        "plots/allele-freqs.svg"]

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")
units = pd.read_table(config["units"], dtype=str).set_index(["sample"], drop=False)
validate(units, schema="schemas/units.schema.yaml")
validate(config, schema="schemas/config.schema.yaml")

if config["filtering"]["vqsr"]:
    filter_type="recalibrated"
else:
    filter_type="hardfiltered"

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),

##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample."""
    fastqs = units.loc[(sample), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def is_single_end(sample):
    """Return True if sample is single end."""
    return pd.isnull(units.loc[(sample), "fq2"])

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(platform=units.loc[(sample), "platform"])

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    return "{Batch}".format(**wildcards)

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("recal/{sample}.bam", sample=wildcards.sample)

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default

def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=input.regions, config["params"]["gatk"]["HaplotypeCaller"])

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

if "restrict-regions" in config["processing"]:
    rule compose_regions:
        input:
            config["processing"]["restrict-regions"]
        output:
            "called/regions.bed"
        conda:
            "envs/bedops.yaml"
        shell:
            "bedextract {input} > {output}"

rule snpeff:
    input:
        "filtered/all.vcf.gz",
    output:
        vcf=reporter("annotated/all.vcf.gz", caption="report/vcf.rst", category="Calls"),
        csvstats="snpeff/all.csv"
    log:
        "logs/snpeff.log"
    params:
        reference=config["ref"]["name"],
        extra="-Xmx6g"
    wrapper:
        "0.27.1/bio/snpeff"

rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"],
        regions = []
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
        gvcf="called/all.g.vcf.gz"
    output:
        vcf=temp("genotyped/all.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.log"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format("SNP" if wildcards.vartype == "snvs" else "INDEL")

rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="genotyped/all.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/selectvariants"

def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}

rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="filtered/all.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.hardfiltered.vcf.gz")
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"

rule recalibrate_calls:
    input:
        vcf="filtered/all.{vartype}.vcf.gz"
    output:
        vcf=temp("filtered/all.{vartype}.recalibrated.vcf.gz")
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]
    log:
        "logs/gatk/variantrecalibrator/{vartype}.log"
    wrapper:
        "0.27.1/bio/gatk/variantrecalibrator"

rule merge_calls:
    input:
        vcf=expand("filtered/all.{vartype}.{filtertype}.vcf.gz", vartype=["snvs", "indels"], filtertype=filter_type)
    output:
        vcf="filtered/all.vcf.gz"
    log:
        "logs/picard/merge-filtered.log"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"

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


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}.zip"
    wrapper:
        "0.27.1/bio/fastqc"

rule samtools_stats:
    input:
        "recal/{sample}.bam"
    output:
        "qc/samtools-stats/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}.log"
    wrapper:
        "0.27.1/bio/samtools/stats"

rule multiqc:
    input:
        expand(["qc/samtools-stats/{u.sample}.txt",
        "qc/fastqc/{u.sample}.zip",
        "qc/dedup/{u.sample}.metrics.txt"], u=units.itertuples()), "snpeff/all.csv"
    output:
        reporter("qc/multiqc.html", caption="report/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"

rule vcf_to_tsv:
    input:
        "annotated/all.vcf.gz"
    output:
        reporter("tables/calls.tsv.gz", caption="report/calls.rst", category="Calls")
    conda:
        "envs/rbt.yaml"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | rbt vcf-to-txt -g --fmt DP AD --info ANN | gzip > {output}"

rule plot_stats:
    input:
        "tables/calls.tsv.gz"
    output:
        depths=reporter("plots/depths.svg", caption="report/depths.rst", category="Plots"),
        freqs=reporter("plots/allele-freqs.svg", caption="report/freqs.rst", category="Plots")
    conda:
        "envs/stats.yaml"
    script:
        "scripts/plot-depths.py"
