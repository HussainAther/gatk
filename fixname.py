o = open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed.fa", "w")

with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fa", "r") as file:
    for line in file:
        if line.startswith(">"):
#            gene_name = line.split(" ")[3].replace("gene:", "")
#            o.write(">" + gene_name + " " + " ".join(line.split(" ")[1:3]) + " gene:" + gene_name + " " + " ".join(line.split(" ")[4:]))
            gene_name = line.split(" ")[0].replace(">", "")
            o.write(">" + gene_name + " " + " ".join(line.split(" ")[1:3]) + " gene:" + gene_name + " " + " ".join(line.split(" ")[4:]))
        else:
            o.write(line)
