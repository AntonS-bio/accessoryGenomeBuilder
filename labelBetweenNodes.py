# import configs.Kp as KpConfig
# config=KpConfig.configData()
# outputDir=config.outputDir
# sampleGffDir=config.sampleGffDir


def splitBlastID(idValue):
    sampleData=idValue.replace("::",":").split(":")
    sampleFile=sampleData[0]
    sampleContig=sampleData[1]
    starEnd=sampleData[2].split("-")
    return [sampleFile, sampleContig, int(starEnd[0]), int(starEnd[1])]

def getGeneProduct(data): #data=[sampleID, [ {key=W###/B###, value=[sample::chr::from-to] } ]]
    #parses GFF file and for each vertex (samples::chr:from-to) identifies the product of the gene.
    #for speed, get all coordinates into a set
    geneProductDic={}
    for clique in data[1].keys():
        for vertex in data[1][clique]:
            geneProductDic[vertex]=""
    #parse GFF
    with open(sampleGffDir+data[0]+".gff") as gffFile:
        for line in gffFile:
            if line[0]!="#":
                gffColumns=line.split("\t")
                lineVertexID=data[0]+"::"+gffColumns[0]+":"+gffColumns[3]+"-"+gffColumns[4]
                if lineVertexID in geneProductDic:
                    for value in gffColumns[8].split(";"):
                        if value[0:8]=="product=":
                            geneProductDic[lineVertexID]=value.replace("product=","").strip()
    return geneProductDic


print("Gathering between and within cliques")
betweenVertices={} #vertex ID, clique W/BID
samplesNodesVerticesDic={} # {key=sampleID, value={key=node i.e. W### or B###, value= [vertices i.e. samples::chr:from-to] } }

for line in open(outputDir+"betweenCliques.csv"): #file must have ID as first column, second column ID is optional
    values=line.strip().split("\t")
    vertex=splitBlastID(values[0])
    if vertex[0] not in samplesNodesVerticesDic:
        samplesNodesVerticesDic[vertex[0]]={}
    if values[1] not in samplesNodesVerticesDic[vertex[0]]:
        samplesNodesVerticesDic[vertex[0]][values[1]]=[]
    samplesNodesVerticesDic[vertex[0]][values[1]].append(values[0])
    betweenVertices[values[0]]=values[1]


print("Getting genes names")
counter=0
vertexProductsDic={}
for sample in samplesNodesVerticesDic.keys():
    vertexProductsDic.update(getGeneProduct([sample, samplesNodesVerticesDic[sample]]))
    counter+=1
    print(float(counter/len(samplesNodesVerticesDic.keys())))

cliqueProducts={} #{key=clique B###, value=selected product, hypothetical/putative are last to be picked}
for clique in betweenVertices.values():
    cliqueProducts[clique]=""

for vertex in vertexProductsDic.keys():
    if cliqueProducts[betweenVertices[vertex]]=="":
        cliqueProducts[betweenVertices[vertex]]=vertexProductsDic[vertex]
    #decide if current product is better than already assigned
    if cliqueProducts[betweenVertices[vertex]].find("hypothetical protein")>-1 or cliqueProducts[betweenVertices[vertex]].find("putative protein")>-1:
        cliqueProducts[betweenVertices[vertex]]=vertexProductsDic[vertex]

with open(outputDir+"betweenCliquesProducts.tsv","w") as output:
    output.write("Cliques\tProduct\n")
    for key in cliqueProducts:
        output.write(key+"\t"+cliqueProducts[key]+"\n")

