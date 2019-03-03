chr_list = []

o = output=("/data/athersh/gatk/mapped/Undetermined_S0_L001-1.sorted.newheader.sam", "w")

with open("/data/athersh/gatk/mapped/Undetermined_S0_L001-1.sorted.header.sam", "r") as file:
    for line in file:
        if "chr" in line:
            chro = line.split("\t")[1].replace("SN:", "")  
            if chro not in chr_list: 
                chr_list.append(chro)

print(chr_list)
