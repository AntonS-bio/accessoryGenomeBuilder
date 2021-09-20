import subprocess
from os import listdir, chdir, chdir, remove
import multiprocessing as mp

import pickle

samplesDir=""
withinNodes={}
betweenNodes={}
withinAMRGenes={}
allNodes={}

def splitBlastID(nodeName):
    values=nodeName.replace("::",":").split(":")
    fileName=values[0]
    #print("Node: "+nodeName)
    #print("Node values:"+'\t'.join(str(f) for f in values))
    chr=values[1]
    [start, end]=values[2].split("-")
    return [fileName, chr, start, end]


def runBlast(graphNode):
    print(graphNode)
    nodeGenesAMR={} #key=gene coordinates, values=blast result
    global allNodes
    for gene in allNodes[graphNode]:
        geneCoodinates=splitBlastID(gene)
        geneSequences=subprocess.run("bedtools getfasta -fi "+samplesDir+geneCoodinates[0]+".fasta -bed <(echo -e \""+geneCoodinates[1]+"\t"+geneCoodinates[2]+"\t"+geneCoodinates[3]+"\")", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)
        geneLength=abs(int(geneCoodinates[1])-int(geneCoodinates[2]))          
        for sequence in geneSequences.stdout.decode().split("\t"):
            blastNucleotide=subprocess.run("blastn -query <(printf \'"+sequence+"\' ) -task 'megablast' -db virRefSeqNuc -max_target_seqs 1000000000 -num_threads 4 -evalue 0.0000000001 -word_size 11 -outfmt \" 6 delim=  qseqid qstart qend sseqid evalue pident \" ", executable="/bin/bash", shell=True, stdout=subprocess.PIPE)
            blastResult=blastNucleotide.stdout.decode().split("\n")
            bestMatch="" ##there should be one best match per gene. Defined as longest given identity of over 95%.
            bestMatchIdentity=0
            for value in blastResult:
                #print(value)
                if value!="":
                    blastOutputs=value.split("\t")
                    if float(blastOutputs[5])>95 and float(blastOutputs[5])>bestMatchIdentity and abs(int(blastOutputs[2])-int(blastOutputs[1]))>=geneLength*0.9:
                        bestMatch=value


            if len(bestMatch)>1:
                nodeGenesAMR[gene]=bestMatch
            if graphNode in betweenNodes:
                break #there is no need to do blast search for every gene in betweenClique as they are all the same.
    print("Results:"+str(len(nodeGenesAMR)))
    if len(nodeGenesAMR)>0:
        print(nodeGenesAMR)
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
    subprocess.call("pkill -9 identifyVirusGenes", shell=True) #cleans up existing running threads if any.

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
        for line in file:
            if line[0]==">":
                values=line.strip().split(" ")
                id=values[0].replace(">","")
                refGenesMetaData[id]=' '.join(str(f) for f in values[1:len(values)])

    chdir(AMRDir)

    subprocess.call("makeblastdb -in "+ AMRDir+"virRefSeqNuc.fasta -blastdb_version 4 -title virRefSeqNuc -out virRefSeqNuc -dbtype nucl", shell=True)


    if __name__ == 'identifyVirusGenes':
        global allNodes
        allNodes={**withinNodes, **betweenNodes}
        pool=mp.Pool(30)

        #temp={}
        #for key in list(withinNodes.keys())[1:50]:
        #    temp[key]=withinNodes[key]
        #withinNodes=temp

        withinResults=pool.map(runBlast, withinNodes)
        pool.close()
        pool.join()
        joinedWithinResults={}
        for result in withinResults:
            if len(result)>0:
                joinedWithinResults={**joinedWithinResults, **result}

        pool=mp.Pool(50)

        #temp={}
        #for key in list(betweenNodes.keys())[1:50]:
        #    temp[key]=betweenNodes[key]
        #betweenNodes=temp

        betweenResults=pool.map(runBlast, betweenNodes)
        pool.close()
        pool.join()
        joinedBetweenResults={}
        for result in betweenResults:
            if len(result)>0:
                joinedBetweenResults={**joinedBetweenResults, **result}

        print("pickling")
        pickle.dump( joinedWithinResults, open(outputDir+"joinedWithinResultsVirus.pkl", "wb") )
        pickle.dump( joinedBetweenResults, open(outputDir+"joinedBetweenResultsVirus.pkl", "wb") )

        joinedWithinResults=pickle.load( open(outputDir+ "joinedWithinResultsVirus.pkl", "rb" ) )
        joinedBetweenResults=pickle.load( open(outputDir+ "joinedBetweenResultsVirus.pkl", "rb" ) )

        #sys.exit()
        for key in joinedWithinResults:
            values=joinedWithinResults[key].split("\t")
            #print(key+'\t'.join(refGenesMetaData[values[3]]))

        newMetaData=open(outputDir+"MetaDataVirus.tsv","w")
        allResults={**joinedWithinResults, **joinedBetweenResults}
        allCliques= {**betweenNodes, **withinNodes}
        for line in open(outputDir+"MetaDataIS.tsv"):
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
            virusName={}
            geneCount="between"
            fileCount="within"
            fileName="between"
            if cliqueNode not in allCliques:
                continue
            virusID=""
            virLen=""
            isVir="No"
            #print(allResults.keys())
            for gene in allCliques[cliqueNode]:
                if gene in allResults:
                    blastValues=allResults[gene]
                    blastValues=blastValues.split("\t")
                    virusID=blastValues[3]+" "+refGenesMetaData[blastValues[3]].strip()
                    virLen=len(allCliques[cliqueNode])
                    isVir="Yes"
                    break

            print(virusID)
            newMetaData.write(line.strip()+"\t"+virusID+"\t"+str(virLen)+"\t"+isVir+"\n")
        newMetaData.close()
        #print(joinedWithinResults)
        #print(joinedBetweenResults)