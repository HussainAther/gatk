samples = ["Undetermined_S0_L001_R1_001", 
           "Undetermined_S0_L001_R2_001", 
           "Undetermined_S0_L002_R1_001",
           "Undetermined_S0_L002_R2_001"]

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
        report("qc/multiqc.html", caption="../report/multiqc.rst", category="Quality control")
    log:
        "logs/multiqc.log"
    wrapper:
        "0.27.1/bio/multiqc"
