import mygene

o = open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed.fa", "w")
o1 = open("/data/athersh/gatk/10ensnames", "w")

ens = []
gene_name = []
with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fa", "r") as file:
    for line in file:
        if line.startswith(">"):
            gene = line.split(" ")[3].split(":")[1]
            gene_name.append(line.split(" ")[3].split(":")[1])
            ens.append(line.split(" ")[0].replace(">", ""))
            o.write(">" + str(gene) + " " + " ".join(line.split(" ")[1:]))
        else:
            o.write(str(line))

with o1 as file:
    for i in ens:
        file.write(i)
        file.write("\n")


