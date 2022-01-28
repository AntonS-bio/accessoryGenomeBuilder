import time
from SampleMetadata import Metadata as sm
from VertexMetadata import Metadata as vm


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

def run(config):
    outputDir=config.outputDir
    global sampleGffDir
    sampleGffDir=config.sampleGffDir
    metaDataDir=config.wd


    print("Loading sample metadata")
    rawMetaData={} 
    for line in open(metaDataDir+"sampleMetaData.txt"):
        values=line.strip().split("\t")
        rawMetaData[values[0]]=sm(values[0], values[1], values[2], values[3])

    crossCliqueNodes=[]
    for line in open(outputDir+"crossCliqies.csv"):
        crossCliqueNodes.append(line.strip().split("\t")[0])
        crossCliqueNodes.append(line.strip().split("\t")[1])
    crossCliqueNodes=set(crossCliqueNodes)
    print("Gathering between and within cliques")
    betweenVertices={} #vertex ID, clique W/BID
    withinVertices={} #vertex ID, clique W/BID
    allVertices=set() #all vertex IDs
    samplesNodesVerticesDic={} # {key=sampleID, value={key=node i.e. W### or B###, value= [vertices i.e. samples::chr:from-to] } }

    for DicFile in zip([betweenVertices,withinVertices],["betweenCliques.csv","withinCliques.csv"]):
        for line in open(outputDir+DicFile[1]): #file must have ID as first column, second column ID is optional
            values=line.strip().split("\t")
            if values[1] in crossCliqueNodes: #if clique not part of crossCliques, there is no point annotating it
                vertex=splitBlastID(values[0])
                if vertex[0] not in samplesNodesVerticesDic:
                    samplesNodesVerticesDic[vertex[0]]={}
                if values[1] not in samplesNodesVerticesDic[vertex[0]]:
                    samplesNodesVerticesDic[vertex[0]][values[1]]=[]
                samplesNodesVerticesDic[vertex[0]][values[1]].append(values[0])
                DicFile[0][values[0]]=values[1]

    #allVertices=set(betweenVertices.keys())|set(withinVertices.keys()) #all samplesID::chr::from-to for which to get product

    print("Getting genes names")
    counter=0
    vertexProductsDic={}
    for sample in samplesNodesVerticesDic.keys():
        vertexProductsDic.update(getGeneProduct([sample, samplesNodesVerticesDic[sample]]))
        counter+=1
        #print(float(counter/len(samplesNodesVerticesDic.keys())))
    
    geneVertices=[] #this is for legacy reason, the annotation speed has been improved a lot, but rest of logic left unchanged.
    for key in vertexProductsDic.keys():
        # if vertexProductsDic[key]=="":
        #     print(key)
        geneVertices.append([key, vertexProductsDic[key]])



    if __name__ == 'annotateVertices':


        geneVerticesProducts={}
        for val in geneVertices:
            geneVerticesProducts[val[0]]=val[1]


        geneVerticesEdgesCount={}
        for collection in [betweenVertices, withinVertices]:
            for vertex in collection:
                if vertex not in geneVerticesEdgesCount:
                    geneVerticesEdgesCount[vertex]=0
                geneVerticesEdgesCount[vertex]+=1


        print("Got all gene names")
        print("Annotating gene vertices")
        verticesMetadata={}
        for geneVertex in geneVerticesProducts:
            metadata=vm("gene",geneVertex)
            metadata.splitBlastID(geneVertex)
            metadata.selectProductName([geneVerticesProducts[geneVertex]])
            verticesMetadata[geneVertex]=metadata

        print("Annotating nonGene verices")
        nonGeneMetadata={} #key=str, value=vm
        originalVertexCount=len(geneVerticesProducts)
        print(len(geneVerticesProducts))
        while len(geneVerticesProducts)>0:
            geneVertex=list(geneVerticesProducts.keys())[0]
            #print(geneVertex+"\t"+str(1-len(geneVerticesProducts)/originalVertexCount))
            verticesToRemove=set()
            if geneVertex in betweenVertices and betweenVertices[geneVertex] not in nonGeneMetadata:
                betweenMetaData=vm("between", betweenVertices[geneVertex])
                #collect all products for this vertex
                betweenMetaData.set_ST("between")
                products=set()
                samples=set()
                instanceCount=0
                for key in geneVerticesProducts:
                        #because multiple samples:chr::start-end nodes can point to same nonGene node, there is no point is looping over them all
                        #so only collect products and samples with specific genes in between nonGene node and mark all the genes pointing to the nonGene node for deletion
                    if key in betweenVertices and betweenVertices[key] == betweenVertices[geneVertex]:
                        verticesToRemove.add(key)
                        products.add(geneVerticesProducts[geneVertex])
                        samples.add(verticesMetadata[geneVertex].get_sampleID())
                        instanceCount+=1
                betweenMetaData.selectProductName(products)
                betweenMetaData.set_geneCount(instanceCount)
                betweenMetaData.set_sampleCount(len(samples))
                #create annotation for between vertex
                nonGeneMetadata[betweenVertices[geneVertex]]=betweenMetaData
            if geneVertex in withinVertices and withinVertices[geneVertex] not in nonGeneMetadata:
                withinMetaData=vm("within", withinVertices[geneVertex])
                #collect all products for this vertex
                withinMetaData.set_ST(rawMetaData[verticesMetadata[geneVertex].get_sampleID()].get_ST())
                withinMetaData.set_year(rawMetaData[verticesMetadata[geneVertex].get_sampleID()].get_year())
                withinMetaData.set_location(rawMetaData[verticesMetadata[geneVertex].get_sampleID()].get_location())
                samples=set()
                instanceCount=0
                for key in geneVerticesProducts:
                    if key in withinVertices and withinVertices[key] == withinVertices[geneVertex]:
                        #because multiple samples:chr::start-end nodes can point to same nonGene node, there is no point is looping over them all
                        #so only collect number of genes the within nonGene node contains and mark all the genes pointing to the nonGene node for deletion
                        verticesToRemove.add(key)
                        samples.add(verticesMetadata[geneVertex].get_sampleID())
                        instanceCount+=1
                withinMetaData.set_product(["within"])
                withinMetaData.set_geneCount(instanceCount)
                withinMetaData.set_sampleCount(len(samples))
                #create annotation for between vertex
                nonGeneMetadata[withinVertices[geneVertex]]=withinMetaData

            #check that tentative removal list is actually safe to remove by verifying that their nodes are present in nonGeneMetadata
            temp=set()
            for item in verticesToRemove:
                if item not in betweenVertices and (item in withinVertices and withinVertices[item] in nonGeneMetadata):
                    #sample::chr:start-end vertex is only present in within node and is already added to nonGeneMetadata
                    temp.add(item)
                elif item not in withinVertices and (item in betweenVertices and betweenVertices[item] in nonGeneMetadata):
                    #sample::chr:start-end vertex is only present in within node and is already added to nonGeneMetadata
                    temp.add(item)
                elif (item in betweenVertices and betweenVertices[item] in nonGeneMetadata) and (item in withinVertices and withinVertices[item] in nonGeneMetadata):
                    temp.add(item)

            verticesToRemove=temp

            #remove those vertices that have already been examined
            for m in verticesToRemove: #this has to be done after both within and between vertices have been processed because they overlap
                del geneVerticesProducts[m]


        print("metadata\t"+str(len(nonGeneMetadata)))
    
        print("Writing output")
        output=open(outputDir+"vertexMetaData.tsv", "w")
        output.write("ID\tNodeType\tST\tProduct\tYear\tLocation\n")
        for value in nonGeneMetadata.values():
            output.write('\t'.join(str(f) for f in [value.get_vertexID(),value.get_vertexType(),value.get_ST(),value.get_product(),value.get_year(),value.get_location()])+"\n")

