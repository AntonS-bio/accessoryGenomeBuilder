import subprocess
import sys

samplesDir=sys.args[1]

fileName="KpV14"

gff=subprocess.run("grep \"hypothetical protein\" "+samplesDir+"/prokkaExFasta/"+fileName+".gff | grep -v UniProt", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

genes=gff.stdout.decode().split("\n")

for gene in genes:
    values=gene.split("\t")
    fasta=subprocess.run("bedtools getfasta -fi "+samplesDir+fileName+".fasta -bed <(echo -e \""+values[0]+"\t"+values[3]+"\t"+values[4]+"\")", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)    
    fasta=fasta.stdout.decode().strip()
    print(fasta)
    blastx=subprocess.run("blastx -query <(echo -e \""+fasta+"\") -task blastx-fast -db refseq_protein -evalue 0.001 -max_target_seqs 25 -outfmt '6 qaccver saccver pident length evalue stitle sid' -remote", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
    blastx=blastx.stdout.decode().strip().split("\n")
    print(blastx)
    break

