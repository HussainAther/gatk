chr_list = []

o = open("/data/athersh/gatk/mapped/Undetermined_S0_L001-1.sorted.newheader.sam", "w")

with open("/data/athersh/gatk/mapped/Undetermined_S0_L001-1.sorted.header.sam", "r") as file:
    for line in file:
        if "chr" in line:
            chro = line.split("\t")[1].replace("SN:", "")  
            if chro not in chr_list:
                o.write("\t".join(line.split("\t")[:-1]))
                o.write("\n")
                chr_list.append(chro)
        else:
            o.write(line)
