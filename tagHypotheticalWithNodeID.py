
# import configs.Kp as KpConfig
# config=KpConfig.configData()

# wd=config.wd
# outputDir=config.outputDir

# tempDir=config.tempDir
# refSample=config.refSample
# samplesDir=config.samplesDir
# sampleGffDir=config.sampleGffDir

newGffDir="destination for gff files where 'hypothetical' is replaced with graph vertex id."
tempDir=wd+"tempNCBI_Subset/"
debug=True

def splitBlastID(nodeName):
    values=nodeName.replace("::",":").split(":")
    fileName=values[0]
    #print("Node: "+nodeName)
    #print("Node values:"+'\t'.join(str(f) for f in values))
    chr=values[1]
    [start, end]=values[2].split("-")
    return [fileName, chr, start, end]

gffFiles={} #{fileName : [ [nodeID, blastIDsplit], [nodeID, blastIDsplit] ]}
with open(outputDir+"betweenCliques.csv") as file:
    for line in file:
        node=line.strip().split("\t")
        [filename, chr, start, end]=splitBlastID(node[0])
        if filename not in gffFiles:
            gffFiles[filename]=[]
        gffFiles[filename].append([node[1],filename, chr, start, end]) #node[1] is B# node id

for gffName in gffFiles.keys():
    outputfile=open(newGffDir+gffName+".gff", "w")
    with open(sampleGffDir+gffName+".gff") as gffFile:
        for line in gffFile:
            if line[0:2]=="##":
                #this is header line, print as is
                outputfile.write(line)
            else:
                lineValues=line.strip().split("\t")
                geneChr=lineValues[0]
                geneStart=lineValues[3]
                geneEnd=lineValues[4]
                if lineValues[8].find("product=hypothetical protein")>0 and lineValues[8].find("RefSeq")<0 and lineValues[8].find("UniProtKB:")<0 and lineValues[8].find("RM-I.faa")<0:
                    #find the ID of this gene if it exists
                    nodeID=""
                    for node in gffFiles[gffName]:
                        if node[2]==geneChr and node[3]==geneStart and node[4]==geneEnd:
                            nodeID=node[0]
                            break
                    if nodeID!="": #in some cases, hypothetical protein is not part of a accessory genome. Skip those cases
                        commentValues=lineValues[8].split(";")
                        for i in range(len(commentValues)):
                            if commentValues[i].find("inference=")>-1:
                                commentValues[i]="inference=graphNode:"+nodeID
                                lineValues[8]=';'.join(commentValues)
                outputfile.write('\t'.join(lineValues)+"\n")
    outputfile.close()
    #break
