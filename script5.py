
#---------------Making final fasta file
import re
with open ('Final_haps.txt') as f1:
    l=[]
    for i in f1:
        j=i.rstrip()
        k=j.split("\t")
        barcode=k[0]
        hap=k[1]
        hap1=hap.split('_')
        hp2=str(hap1)
        res1 = "".join(re.findall("[a-zA-Z]+", hp2))  #picking only alleles
        pos1="".join(re.findall("[0-9]+", hp2))  ##contains information of position of the alleles
        t=barcode,res1,pos1
        l.append(t)
f= open("haplotype.fasta", "a")
for p in l[1:]:
        f.write('>'+p[0]+'\n')
        f.write(p[1]+ '\n')
f.close()

