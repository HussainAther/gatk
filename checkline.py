with open("/data/athersh/gatk/Danio_rerio.GRCz10.cdna.all.fixed2.fa.headers.uniq", "r") as file:
    for line in file:
        if len(line.split(" ")) > 1:
            if (line.split(" ")[6]) != "1":
                print(line.split(" ")[6])
