# -*- coding: utf-8 -*-

fin = open("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/hg38_ref/gencode.v21.annotation.gtf", "r")
fout = open("/Users/liuzhe/Desktop/cityu/ceRNA_network/02_ceRNA-Projects/TCGA_data/hg38_ref/geneInfo.txt", "w")


for line in fin:
    line = line.rstrip()
    if(line.startswith("#")):
        continue
    else:
        s = line.split("\t")
        if(s[2] != "gene"):
            continue
        else:
            ss = s[8]
            sss = ss.split('''"''')[1]
            geneID = sss.split('''.''')[0]
            #print(geneID)
            geneName = ss.split('''"''')[9]
            geneType = ss.split('''"''')[5]
            fout.write(geneID + "\t" + geneName + "\t" + geneType + "\n")

fout.close()

