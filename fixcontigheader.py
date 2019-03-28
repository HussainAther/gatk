
headerlist = []

o = open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed3.fa", "w") 

with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed2.fa", "r") as file:
    for line in file:
        if line.startswith(">"):
            length = int(line.split(":")[3]) - int(line.split(":")[2]) - 1
            o.write(line.split(" ")[0] + "." + str(length) + " " + " ".join(line.split(" ")[1:]))
        else:
            o.write(line)
