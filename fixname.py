o = open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed.fa", "w")

with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fa", "r") as file:
    for line in file:
        if line.startswith(">"):
            gene_name = line.split(" ")[0].replace(">", "")
            o.write(" ".join(line.split(" ")[:3]) + " gene:" + gene_name + " ".join(line.split(" ")[4:]))
        else:
            o.write(line)
