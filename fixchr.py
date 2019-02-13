o = open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed2.fa", "w")

with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed.fa", "r") as file:
    for line in file:
        if line.startswith(">"):
            chromosome = line.split(":")[1]
            o.write(">" + chromosome + " "+ " ".join(line.split(" ")[1:]))
        else:
            o.write(line)
