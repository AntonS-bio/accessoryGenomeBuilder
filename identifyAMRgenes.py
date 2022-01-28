import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC, Gapped, generic_dna
from os import listdir, chdir, chdir, remove
from os.path import isfile, join, exists
import multiprocessing as mp
import csv
import pickle

samplesDir=""
withinNodes={}
betweenNodes={}
withinAMRGenes={}
allNodes={}
geneLengths={}
AssemblySequences={}

def splitBlastID(nodeName):
    values=nodeName.replace("::",":").split(":")
    fileName=values[0]
    #print("Node: "+nodeName)
    #print("Node values:"+'\t'.join(str(f) for f in values))
    chr=values[1]
    [start, end]=values[2].split("-")
    return [fileName, chr, start, end]

def LoadFastasToMemory():
    global allNodes, AssemblySequences

    #generate two dictionaries. One mapping nodes to files and one coordinates to node.
    #Can do in one nested dicitonary, but this is easier to work with
    geneCoodinates={} #{B/W nodeID: [[coordinates1],[coordiantes2]]
    assemblyNodes={} #{file: [nodeIDs]}
    for graphNode in allNodes.keys():
        for gene in allNodes[graphNode]:
            if graphNode not in geneCoodinates:
                geneCoodinates[graphNode]=[]
            values=splitBlastID(gene)
            geneCoodinates[graphNode].append(values)
            if values[0] not in assemblyNodes:
                assemblyNodes[values[0]]=[]
            assemblyNodes[values[0]].append(graphNode)

    assemblyFiles=list(assemblyNodes.keys())

    counter=0
    for assembly in assemblyFiles:
        print("Loading fastas progress: " + str(int(counter/len(assemblyFiles)*100)))
        counter+=1
        record_dict = SeqIO.to_dict(SeqIO.parse(samplesDir+assembly+".fasta", "fasta"))
        for node in assemblyNodes[assembly]:
            AssemblySequences[node]=[]
            for gene in geneCoodinates[node]:
                #gene is splitBlastID output: [fileName, chr, start, end]
                if gene[0]==assembly:
                    record = SeqRecord(
                                Seq(str(record_dict[gene[1]].seq[int(gene[2]):int(gene[3])])),
                                id=assembly+"::"+gene[1]+":"+gene[2]+"-"+gene[3],
                                name=node)
                    AssemblySequences[node].append(record)
           
def runBlast(graphNode):
    #print(graphNode)
    nodeGenesAMR={} #key=gene coordinates, values=blast result
    global allNodes, geneLengths, AssemblySequences
    #for gene in allNodes[graphNode]:
    for record in AssemblySequences[graphNode]: #each record in single gene, record name = nodeID
        #geneCoodinates=splitBlastID(gene)
        #geneSequences=subprocess.run("bedtools getfasta -fi "+samplesDir+geneCoodinates[0]+".fasta -bed <(echo -e \""+geneCoodinates[1]+"\t"+geneCoodinates[2]+"\t"+geneCoodinates[3]+"\")", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
        #for sequence in geneSequences.stdout.decode().split("\t"):
        #for sequence in record:
        blastNucleotide=subprocess.run("blastn -query <(printf \'"+str(record.seq)+"\' ) -task 'megablast' -db AMRnucleotide -max_target_seqs 1000000000 -num_threads 4 -evalue 0.0000000001 -word_size 11 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident \" ", executable="/bin/bash", shell=True, stdout=subprocess.PIPE)
        blastProtein=subprocess.run("blastx -query <(printf \'"+str(record.seq)+"\' )  -db AMRprotein -max_target_seqs 1000000000 -num_threads 4 -evalue 0.0000000001 -word_size 4 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident \" ", executable="/bin/bash", shell=True, stdout=subprocess.PIPE)
        blastResult=blastNucleotide.stdout.decode().split("\n")+blastProtein.stdout.decode().split("\n")
        bestMatch="" ##there should be one best match per gene. Defined as at least 95% identity and match length of at least 90% of gene.
        bestMatchIdentity=0
        for value in blastResult:
            #print(value)
            if value!="":
                blastOutputs=value.split("\t")
                if float(blastOutputs[5])>95 and float(blastOutputs[5])>bestMatchIdentity and abs(int(blastOutputs[2])-int(blastOutputs[1]))>=geneLengths[blastOutputs[3]]*0.9:
                    bestMatch=value
                    bestMatchIdentity=float(blastOutputs[5])

        if len(bestMatch)>1:
            nodeGenesAMR[record.id]=bestMatch
        if graphNode in betweenNodes:
            break #there is no need to do blast search for every gene in betweenClique as they are all the same.
    #print("Results:"+str(len(nodeGenesAMR)))
    # if len(nodeGenesAMR)>0:
    #     print(nodeGenesAMR)
    return nodeGenesAMR



def run(AMRDir, AMRfile, outputDir, PassedSamplesDir):
    global samplesDir
    samplesDir=PassedSamplesDir
    global withinNodes
    withinNodes={}
    global betweenNodes
    betweenNodes={}
    global withinAMRGenes
    withinAMRGenes={}

    counter=0
    subprocess.call("pkill -9 identifyAMRgenes", shell=True) #cleans up existing running threads if any.

    ##collect all cross clique vertex IDs
    counter=0
    for line in open(outputDir+"crossCliqies.csv"):
        values=line.strip().split("\t")
        for value in values:
            if value[0]=="B":
                betweenNodes[value]=[]
            elif value[0]=="W":
                withinNodes[value]=[]
        counter+=1


    ##collect gene coordinates for each cross clique vertex ID
    for line in open(outputDir+"betweenCliques.csv"):
        values=line.strip().split("\t")
        if values[1] in betweenNodes:
            betweenNodes[values[1]].append(values[0])

    for line in open(outputDir+"withinCliques.csv"):
        values=line.strip().split("\t")
        if values[1] in withinNodes:
            withinNodes[values[1]].append(values[0])

    #flip betweenNodes and withinNodes key value pairs
    withinGenes={}
    betweenGenes={}
    for node in withinNodes:
        for gene in withinNodes[node]:
            withinGenes[gene]=node
    for node in betweenGenes:
        for gene in betweenGenes[node]:
            betweenGenes[gene]=node


    refGenesMetaData={} #contains type, subtype, etc for AMR genes
    with open(AMRDir+AMRfile) as file:
        file.readline()
        csv_reader=csv.reader(file, delimiter=',')
        for line in csv_reader:
            for id in line[8:11]:
                if id!="":
                    refGenesMetaData[id]=[line[4], line[5], line[6], line[7], line[1], line[2]] #type, subtype, class, subclass, genefamily, genename

    chdir(AMRDir)
    #collect the lengths of all AMR genes and proteins. This is used in determining completeness of blast result

    global geneLengths
    for record in SeqIO.parse(AMRDir+"nucleotide.fasta", "fasta"):
        geneLengths[record.id]=len(record.seq)
    for record in SeqIO.parse(AMRDir+"protein.fasta", "fasta"):
        geneLengths[record.id]=len(record.seq)*3 #blast gives result in nucleotides, but length is in amino acids

    subprocess.call("makeblastdb -in "+ AMRDir+"nucleotide.fasta -blastdb_version 4 -title AMRnucleotide -out AMRnucleotide -dbtype nucl", shell=True)
    subprocess.call("makeblastdb -in "+ AMRDir+"protein.fasta -blastdb_version 4 -title AMRprotein -out AMRprotein -dbtype prot", shell=True)
    
    if __name__ == 'identifyAMRgenes':# or standAloneRun:
    #if True:
        global allNodes
        allNodes={**withinNodes, **betweenNodes}
        LoadFastasToMemory()
        #sys.exit()
        pool=mp.Pool(40)
        withinResults=pool.map(runBlast, withinNodes)
        pool.close()
        pool.join()
        joinedWithinResults={}
        for result in withinResults:
           if len(result)>0:
               joinedWithinResults={**joinedWithinResults, **result}

        pool=mp.Pool(20)
        betweenResults=pool.map(runBlast, betweenNodes)
        pool.close()
        pool.join()
        joinedBetweenResults={}
        for result in betweenResults:
           if len(result)>0:
               joinedBetweenResults={**joinedBetweenResults, **result}

        print("pickling")
        pickle.dump( joinedWithinResults, open(outputDir+"joinedWithinResults.pkl", "wb") )
        pickle.dump( joinedBetweenResults, open(outputDir+"joinedBetweenResults.pkl", "wb") )

        joinedWithinResults=pickle.load( open(outputDir+ "joinedWithinResults.pkl", "rb" ) )
        joinedBetweenResults=pickle.load( open(outputDir+ "joinedBetweenResults.pkl", "rb" ) )

        #sys.exit()
        for key in joinedWithinResults:
            values=joinedWithinResults[key].split("\t")
            #print(key+'\t'.join(refGenesMetaData[values[3]]))

        newMetaData=open(outputDir+"MetaData.tsv","w")
        allResults={**joinedWithinResults, **joinedBetweenResults}
        allCliques= {**betweenNodes, **withinNodes}
        for line in open(outputDir+"vertexMetaData.tsv"):
            cliqueNode=line.strip().split("\t")[0]
            line=line.strip()
            if cliqueNode[0]=="B":
                #there is a inconsistent number of columns in W/B metadata. For the moment, pad B with two tabs
                line=line.strip()+"\t"+"\t"
            types={}
            subtypes={}
            classes={}
            subclasses={}
            geneFamilies={}
            geneNames={}
            geneCount="between"
            fileCount="within"
            fileName="between"
            if cliqueNode not in allCliques:
                continue
            for gene in allCliques[cliqueNode]:
                if gene in allResults:
                    ncbiID=allResults[gene].split("\t")[3]
                    i=0
                    for dictionary in [types, subtypes, classes, subclasses, geneFamilies, geneNames]:
                        if refGenesMetaData[ncbiID][i] not in dictionary: 
                            dictionary[refGenesMetaData[ncbiID][i]]=0
                        dictionary[refGenesMetaData[ncbiID][i]]+=1
                        i+=1
                if cliqueNode in withinGenes: #only withinNodes have more than one gene and unique file name
                    geneCount=len(withinGenes[cliqueNode])
                    fileName=splitBlastID(withinGenes[cliqueNode][0])[0]
                if cliqueNode in betweenNodes: #only betweenNodes has multiple files
                    fileCount==len(betweenNodes[cliqueNode])
        
            additionalMetaValues=[""]*9
            additionalValue=""
            i=0
            for dictionary in [types, subtypes, classes, subclasses, geneFamilies, geneNames]:
                for key in dictionary:
                    additionalMetaValues[i]=additionalMetaValues[i]+key+"x"+str(dictionary[key])+";"
                i+=1
            additionalMetaValues[6]=str(geneCount)
            additionalMetaValues[7]=fileName
            additionalMetaValues[8]=str(fileCount)
            newMetaData.write(line+"\t"+'\t'.join(additionalMetaValues)+"\n")
        newMetaData.close()
        #print(joinedWithinResults)
        #print(joinedBetweenResults)

