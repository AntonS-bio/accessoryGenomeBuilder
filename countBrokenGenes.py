#take all fastas in specified directory and check all for specific gene
from os import listdir, getcwd, makedirs, remove, chdir, path
from os.path import isfile, join, splitext, exists
import shutil
import pandas as pd
import subprocess
import sys
import multiprocessing as mp
from Bio import SeqIO


# import configs.Kp as KpConfig
# config=KpConfig.configData()

# wd=config.outputDir
# samplesDir=config.samplesDir

def splitBlastID(idValue):
    sampleData=idValue.replace("::",":").split(":")
    sampleFile=sampleData[0]
    sampleContig=sampleData[1]
    starEnd=sampleData[2].split("-")
    return [sampleFile, sampleContig, int(starEnd[0]), int(starEnd[1])]


cliques={} #{key=B#####, values=vertex Sample::Chr:from-to}
samplesPerClique={} ##{key=B#####, values=number of samples with this clique complete or partial gene}
cliqueSeqLenght={} ##{key=B#####, values=lenght of clique sequence}
with open(wd+"betweenCliques.csv") as file:
    for line in file:
        values=line.strip().split("\t")
        if not values[1] in cliques:
            cliques[values[1]]=values[0]
            samplesPerClique[values[1]]=0
            cliqueSeqLenght[values[1]]=0


SamplesToExclude=[]
maxLengthDifference=0.5
minIdentity=90
proteinSearch=False

samplesfiles = [f for f in listdir(samplesDir) if isfile(join(samplesDir, f)) and (splitext(f)[1]==".fasta" or splitext(f)[1]==".fna")]




counter=1
if exists(wd+"temp.fasta"):
    remove(wd+"temp.fasta")
for clique in cliques.keys():
    print(counter/len(cliques))
    #create blast DB
    cliqueCoordinates=splitBlastID(cliques[clique])
    bedFile="\""+cliqueCoordinates[1]+"\t"+str(cliqueCoordinates[2])+"\t"+str(cliqueCoordinates[3])+"\t"+clique+"\""
    subprocess.call("bedtools getfasta -nameOnly -fi "+samplesDir+cliqueCoordinates[0]+".fasta -bed <(echo -e  "+bedFile+") >> temp.fasta", shell=True,  executable="/bin/bash")
    counter+=1

for record in SeqIO.parse(wd+"temp.fasta", "fasta"):
    cliqueSeqLenght[record.id]=len(record.seq)


with open(wd+"cliqueLenghts.tsv", "w") as output:
    for clique in cliqueSeqLenght:
        output.write(clique+"\t"+str(cliqueSeqLenght[clique])+"\n")

if exists(wd+"/tempBlastDB/"):
    shutil.rmtree(wd+"/tempBlastDB/")
makedirs(wd+"/tempBlastDB/")
subprocess.call("makeblastdb -in temp.fasta -title temp -out "+wd+"/tempBlastDB/temp -dbtype nucl \
    -blastdb_version 4 1>/dev/null", shell=True)

counter=0
for sample in samplesfiles:
    cliquesInSample=set()
    print(counter/len(samplesfiles))
    sequence=subprocess.run("blastn -query "+samplesDir+sample+" -task 'megablast' \
            -max_target_seqs 1000000000 -db "+wd+"/tempBlastDB/temp \
            -num_threads 8 -evalue 1.0E-5 -word_size 21 \
            -outfmt \"6 delim=  qseqid qstart qend sseqid sstart send pident evalue qseq\""
            , shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
    blastHits=sequence.stdout.decode().split("\n")
    for blasthit in blastHits:
        blasthit=blasthit.split("\t")
        if len(blasthit)>1 and abs(int(blasthit[4])-int(blasthit[5]))>=100 and float(blasthit[6])>minIdentity:
            cliquesInSample.add(blasthit[3])
    for item in cliquesInSample:
        samplesPerClique[item]=samplesPerClique[item]+1
    counter+=1


with open(wd+"brokenCliques.tsv", "w") as output:
    for clique in samplesPerClique:
        output.write(clique+"\t"+str(samplesPerClique[clique])+"\n")

