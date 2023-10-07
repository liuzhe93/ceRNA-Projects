
lncRNA_symbol = open("lncRNAs.txt", "r")
for line in lncRNA_symbol:
    line = line.rstrip()
    lncRNA.append(line)
lncRNA_symbol.close()

sequences = {}
ac = ""
seq = ""
fr = open("gencode.v21.lncRNA_transcripts.fa", "r")
for line in fr:
    if(line.startswith(">") and seq != ""):
        sequences[ac] = seq
        seq = ""
    if(line.startswith(">")):
        ac = line.split("|")[5]
    else:
        seq = seq + line.strip()
#print(sequences)
fr.close()

fout = open("output_seq.txt", "w")
for k in sequences.keys():
    if(k in lncRNA):
        seq_lnc = sequences.get(k)
        fout.write(">"+k+"\n")
        fout.write(seq_lnc+"\n")

fout.close()
