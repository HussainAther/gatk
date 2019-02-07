import mygene 
mg = mygene.MyGeneInfo()

o = open("data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed.fa", "w")

ens = []
with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fa", "r") as file:
    if line.startswith(">"):
        ens.append(line.split(" ")[0].replace(">", ""))

print(ens)
