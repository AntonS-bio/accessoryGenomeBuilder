import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC, Gapped, generic_dna
from os import listdir, chdir, chdir, remove, mkdir
from os.path import isfile, join, exists
import multiprocessing as mp
import statistics
import shutil

import configs.Kp as KpConfig
config=KpConfig.configData()

wd=config.wd
outputDir=config.outputDir
tempDir=config.tempDir
refSample=config.refSample
samplesDir=config.samplesDir
sampleGffDir=config.sampleGffDir
AMRDir=config.AMRDir
AMRfile=config.AMRfile
postfix=config.postfix
blastIdentityThreshold=config.blastIdentityThreshold
coreGenomeThreshold=config.coreGenomeThreshold
permittedLengthVariance=config.permittedLengthVariance




#set-up temporary directory structure:
if exists(tempDir):
    shutil.rmtree(tempDir)
#shutil.rmtree(outputDir)
#os.mkdir(outputDir)
if not exists(outputDir):
    mkdir(outputDir)
    
mkdir(tempDir)
tempSubDirs=["AllvsAll/", "blastResults/", "genes/", "nonCoreGenes/", "unallignedFasta/"]
for dir in tempSubDirs:
    mkdir(tempDir+dir)

def runBlast(sample):
    sample=sample.replace(postfix,"")
    chdir(tempDir+"/genes/")
    subprocess.call("blastn -query " +samplesDir+sample+postfix+" -task 'megablast' -db "+refSample+"genes -max_target_seqs 1000000000 -num_threads 8 -evalue 0.0000000001 -word_size 20 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident qseq sstrand \" > "+tempDir+"/blastResults/"+sample+".csv", shell=True)

def runMafft(gene):
    gene=gene+".fasta"
    subprocess.call("mafft --auto "+tempDir+"/unallignedFasta/"+gene+" > "+tempDir+"/allignedFasta/"+gene, shell=True)

def splitBlastID(idValue):
    sampleData=idValue.replace("::",":").split(":")
    sampleFile=sampleData[0]
    sampleContig=sampleData[1]
    starEnd=sampleData[2].split("-")
    return [sampleFile, sampleContig, int(starEnd[0]), int(starEnd[1])]

def runNonCoreBlast(file):
    prefix=file.replace(".fasta", "")
    print(prefix)
    subprocess.call("blastn -query "+ tempDir+ "AllvsAll/"+file+" -task 'megablast' -db allNoneCoreGenes -max_target_seqs 1000000000 -num_threads 3 -evalue 0.0000000001 -word_size 20 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident sstrand \" >> "+tempDir+"/AllvsAll/"+prefix+".csv", shell=True)
    return "done"




##Use reference fasta to create a new fasta file where each record is a gene from the reference fasta file
##each gene is based on Prokka (or another) annotation

ExclSamples=[] #open(wd+"exclude.tsv").read().splitlines()

samples=[]
#Get list of samples to use
with open(outputDir+"samples.txt") as file:
    for line in file:
        if line.strip() not in ExclSamples:
            samples.append(line.strip())
refSample=samples[0].replace(postfix,"")


refFasta=SeqIO.to_dict(SeqIO.parse(samplesDir+refSample+postfix, "fasta"))
genes={}

with open(tempDir+"/genes/"+refSample+"genes.fasta", "w") as geneFile:
    for line in open(sampleGffDir+refSample+".gff"):
        if line[0]!="#":
            values=line.split("\t")
            if len(values)>2 and values[2]=="CDS":
                id=values[8].split(";")[0].replace("ID=","").replace("\"","")
                genes[id]=[]
                if values[6]=="+":
                    sequence=refFasta[values[0]][int(values[3]):int(values[4])]
                else:
                    sequence=refFasta[values[0]][int(values[3]):int(values[4])].reverse_complement()
                if len(sequence.seq)>100:
                        geneFile.write(">"+id+"\n")
                        geneFile.write(str(sequence.seq)+"\n")

chdir(tempDir+"/genes/")

##Make the fasta of genes into blast DB
subprocess.call("makeblastdb -in "+ tempDir+"/genes/"+refSample+"genes.fasta -blastdb_version 4 -title "+refSample+"genes -out "+refSample+"genes -dbtype nucl", shell=True)


##run blast of individual samples vs all genes of reference
if __name__ == '__main__':

    pool=mp.Pool(30)
    results=pool.map(runBlast, samples)
    pool.close()
    pool.join()

    #collect gene lengths from results of blast of sample vs reference genes
    print("colleting genes from blast results")
    geneLengths=dict.fromkeys(genes)
    geneSequences={} #key=samplename, value={refGene, sequence} 
    geneCoordinates={} #key=sample name, value={refGene, [Contig, Start, End]} 
    for sample in samples:
        sampleGenes={}
        sample=sample.replace(postfix,"")
        geneSequences[sample]={}
        geneCoordinates[sample]={}
        for line in open(tempDir+"/blastResults/"+sample+".csv"):
            values=line.strip().split("\t")
            allignmentLength=int(values[2])-int(values[1])
            if float(values[5])>blastIdentityThreshold and ( values[3] not in sampleGenes or sampleGenes[values[3]]<allignmentLength ):
                sampleGenes[values[3]]=allignmentLength
                if values[7]=="minus":
                    values[6]=str(Seq(values[6]).reverse_complement())
                geneSequences[sample][values[3]]=values[6]
                geneCoordinates[sample][values[3]]=[values[0], values[1], values[2]]
        for sampleGene in sampleGenes.keys():
            genes[sampleGene].append(sampleGenes[sampleGene])

    #Check which genes are +/- 5% of reference file gene length, this identifies complete genes. 
    #Core genomes are complete genes present in >90% of the samples. 
    print("writing genes that passed quality checks")
    chdir(tempDir)
    genesToAllign=set()
    coreGenes=set()
    with open(outputDir+"samplesPerGene.csv","w") as tempOutput:
        tempOutput.write("Gene name\tTotal genes\tComplete Genes\tMedian Gene Length\n")
        medianLength=0
        for gene in genes.keys():
            complete=0
            if len(genes[gene])>0:
                medianLength=statistics.median(genes[gene])
                geneMedianLB=statistics.median(genes[gene])*(1-permittedLengthVariance)
                geneMedianUB=statistics.median(genes[gene])*(1+permittedLengthVariance)
                for sample in genes[gene]:
                    complete+=1 if (sample>geneMedianLB and sample<geneMedianUB and sample>50) else 0 #100 is min length of the fragment

            if complete>(len(samples)*(1-coreGenomeThreshold)): #this is for composition of core genome and identification of MGEs
                coreGenes.add(gene)

            if complete==len(samples): #this is for phylogenetic tree reconstruction
                genesToAllign.add(gene)
                with open(tempDir+"/unallignedFasta/"+gene+".fasta", "w") as allignedFile:
                    for sample in geneSequences.keys():
                        allignedFile.write(">"+sample+"\n")
                        allignedFile.write(geneSequences[sample][gene]+"\n")
                #extract alligned sequences
            
            tempOutput.write(gene+"\t"+str(len(genes[gene]))+"\t"+str(complete)+"\t"+str(medianLength)+"\n")

    print("collecting non-core genes")
    for sample in geneCoordinates.keys():
        bedData=""
        for gene in geneCoordinates[sample].keys():
            if gene in coreGenes:
                bedData=bedData+'\t'.join(geneCoordinates[sample][gene])+"\n"
        #write genes in to bed file, and merge the entries as there are some overlapping/odd length genes produced by prokka
        with open(tempDir+"temp.bed", "w") as bedFile:
            bedFile.write(bedData)
        subprocess.call("sort -o "+tempDir+"temp.bed -k 1,1 -k2,2n "+tempDir+"temp.bed", shell=True,  executable="/bin/bash")
        subprocess.call("bedtools merge -d 5 -i "+tempDir+"temp.bed > sampleCore.bed", shell=True,  executable="/bin/bash")
        #extract non-core genes from sample's gff
        subprocess.call("bedtools intersect -v -a "+sampleGffDir+sample+".gff -b "+tempDir+"sampleCore.bed | grep CDS > "+tempDir+"/nonCoreGenes/"+sample+".gff", shell=True,  executable="/bin/bash")


    
    #generate aggregate non-Core fasta for blast search
    if exists(tempDir+"/nonCoreGenes/allNoneCoreGenes.fasta"):
        remove(tempDir+"/nonCoreGenes/allNoneCoreGenes.fasta")

    for sample in samples:
        subprocess.call("awk -v OFS=\"\t\" '{print $1,$4,$5,\""+sample+"\"}' "+tempDir+"/nonCoreGenes/"+sample+".gff | grep -v \"#\" > temp.bed", shell=True, executable="/bin/bash")
        subprocess.call("bedtools getfasta -name -fi "+samplesDir+sample+postfix+" -bed temp.bed >> " + tempDir+"/nonCoreGenes/allNoneCoreGenes.fasta", shell=True,  executable="/bin/bash")

    #Run all-vs-all non-core genes blast search.
    chdir(tempDir+"/nonCoreGenes/")
    print("Making non-core DB")
    subprocess.call("makeblastdb -in "+ tempDir+"/nonCoreGenes/allNoneCoreGenes.fasta -blastdb_version 4 -title allNoneCoreGenes -out allNoneCoreGenes -dbtype nucl", shell=True)
    print("Non-core blast")

    #the very large fasta file seems to create a problem for blats, so split it into multiple files N genes
    #generate aggregate non-Core fasta for blast search
    nonCoreSubsetFiles=[]
    subsetFileGenesCount=1000
    counter=0
    SubFilesCount=1
    #totalCounter=0
    records=list(SeqIO.parse(tempDir+"/nonCoreGenes/allNoneCoreGenes.fasta",'fasta'))
    for record in records:
        if (counter==0 and SubFilesCount==1) or counter==subsetFileGenesCount:
            if SubFilesCount!=1:
                tempFastaFile.close()
            SubFilesCount+=1
            tempFileName="allNonCoreSubFile_"+str(SubFilesCount)+".fasta"
            nonCoreSubsetFiles.append(tempFileName)
            tempFastaFile=open(tempDir+"AllvsAll/"+tempFileName, "w")
            counter=0

       
        tempFastaFile.write(record.format("fasta"))
        counter+=1
    if not tempFastaFile.closed:
        tempFastaFile.close()

    pool=mp.Pool(20)
    results=pool.map(runNonCoreBlast, nonCoreSubsetFiles)
    pool.close()
    pool.join()

    if exists(tempDir+"/AllVsAll.csv"):
        remove(tempDir+"/AllVsAll.csv")

    with open(tempDir+"/AllVsAll.csv", "w") as AllVsAll:
        for file in nonCoreSubsetFiles:
            file=file.replace(".fasta",".csv")
            with open(tempDir+"AllvsAll/"+file,"r") as source:
                for line in source:
                    AllVsAll.write(line)

    ##filter out those blast results that are self-vs-self or reverse complements (A-B and B-A)

    processedGenes=set()
    geneCount={}
    linesToWrite=[]
    vertexIDs=set()
    for line in open(tempDir+"AllVsAll.csv"):
        values=line.strip().split("\t")
        source=splitBlastID(values[0]) #[sampleFile, sampleContig, starEnd[0], starEnd[1]]
        target=splitBlastID(values[3]) 
        lengthRatio=abs( source[3]- source[2] )/abs( target[3]- target[2] )
        if values[0]!=values[3] and (values[0]+values[3]) not in processedGenes and lengthRatio > (1-permittedLengthVariance) and lengthRatio < (1+permittedLengthVariance) and float(values[5])>blastIdentityThreshold: # 5 is the pident
            vertexIDs.add(values[0])
            linesToWrite.append(line)
            processedGenes.add(values[0]+values[3]) #this captures that in all vs all each hit is duplicated
            processedGenes.add(values[3]+values[0]) #this captures that in all vs all each hit is duplicated
            for gene in [values[0], values[3]]:
                if gene in geneCount:
                    geneCount[gene]+=1
                else:
                    geneCount[gene]=1

    #this second loop is needed to avoid all-Vs-all full genenome (as opposed non-core genome) blast which would take ages
    filteredBlast=open(outputDir+"betweenEdges.csv", "w")
    for line in linesToWrite:
        values=line.strip().split("\t")
        if geneCount[values[0]]<(len(samples)*(1-coreGenomeThreshold)) and geneCount[values[3]]<(len(samples)*(1-coreGenomeThreshold)):
            filteredBlast.write(values[0]+"\t"+values[3]+"\tBetween\n")


    filteredBlast.close()

    with open(outputDir+"vertices.csv", "w") as verticesFile:
        for id in vertexIDs:
            verticesFile.write(id+"\n")
    
    chdir(wd)
    import groupGenesToMge
    groupGenesToMge.run(wd, sampleGffDir, outputDir)

    import generateCliques
    generateCliques.run(outputDir)

    import cliqueGeneMatrix_v2
    cliqueGeneMatrix_v2.run(outputDir)
    
    import annotateVertices
    annotateVertices.run(sampleGffDir, outputDir)

    import identifyAMRgenes
    identifyAMRgenes.run(AMRDir, AMRfile, outputDir,samplesDir)
