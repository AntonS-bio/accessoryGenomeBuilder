

def run(wd, sampleGffDir, outputDir):

    def sortGenomicFeatures(data):#{sample, {chr, [ [start, end] ] }}
        for sample in data.keys():
            for chr in data[sample].keys():
                data[sample][chr].sort(key=lambda e: e[0])


    sampleData={} #{sample, {chr, [ [start, end] ] }}
    for line in open(outputDir+"vertices.csv"): #vertices.csv
        values=line.strip().replace("::",":").split(":")
        sample=values[0]
        chr=values[1]
        start=int(values[2].split("-")[0])
        end=int(values[2].split("-")[1])
        if sample in sampleData:
            if chr in  sampleData[sample]:
                sampleData[sample][chr].append([start,end])
            else:
                sampleData[sample][chr]=[ [start,end] ]
        else:
            sampleData[sample]={chr: [ [start,end] ] }

    sortGenomicFeatures(sampleData)


    print("loading gff")
    
    gffData={} #{sample, {chr, [ [start,end] ] }}

    for sample in sampleData:
        for line in open(sampleGffDir+sample+".gff"):
            if line[0]!="#":
                values=line.strip().split("\t")
                if values[2]=="CDS":
                    chr=values[0]
                    [start, end]=[int(values[3]), int(values[4])]
                    if sample in gffData:
                        if chr in  gffData[sample]:
                            gffData[sample][chr].append([start,end])
                        else:
                            gffData[sample][chr]=[ [start,end] ]
                    else:
                        gffData[sample]={chr: [ [start,end] ] }

    sortGenomicFeatures(gffData)

    print("Merging genomic features into MGEs")
    #find the density of non-core genes.
    lookBackDepth=6
    globalFeatureCounter=1
    mgeFeatures={} #{sample : {chr: {feature: [start, end]}}}
    for sample in gffData.keys():
        for chr in gffData[sample].keys():
            nonCoreGenes=[0]*len(gffData[sample][chr])
            if chr in sampleData[sample]:
                insideInsertion=False
                firstNonCoreFeatureIndex=0
                lastNonCoreFeatureIndex=0
                featureCounter=0

                for feature in gffData[sample][chr]:#this needs to check if feature [region of sample/chr overlaps with any of features in gffData[sample][chr][[start,end]]]
                    #determine if gff feature is core or non-core
                    nonCoreGenes[featureCounter]=0
                    for startEnd in sampleData[sample][chr]:
                        if not( feature[1]<startEnd[0] or feature[0]>startEnd[1]):
                            #this is non-core gene
                            nonCoreGenes[featureCounter]=1
                            break
                    featureCounter+=1
                
                    if featureCounter>=lookBackDepth:
                        nonCoreGenesCount=sum(nonCoreGenes[(featureCounter-lookBackDepth):featureCounter])

                        if nonCoreGenesCount>(lookBackDepth*0.75) and featureCounter!=len(gffData[sample][chr]):
                            if not insideInsertion:
                                #find the last core gene, the insert starts right after it
                                for i in range((featureCounter-lookBackDepth-1),featureCounter):
                                    if nonCoreGenes[i]==1: #first non-core genes in among preceeding genes
                                        firstNonCoreFeatureIndex=0 if i==0 else i #i=0 is corner case, normally, the last examined feature is the first noncore
                                        break

                            insideInsertion=True
                        elif (insideInsertion and nonCoreGenesCount<(lookBackDepth*0.5) or (insideInsertion and featureCounter==len(gffData[sample][chr]))): #checks on non core genes count being higher than x% of look back depth.
                            if insideInsertion or featureCounter==len(gffData[sample][chr]): #the second accomodates the case where contig ends in an insertion.
                                #insertion sequence has ended a few genes earlier, find the last noncore gene
                                insideInsertion=False
                                if featureCounter==len(gffData[sample][chr]):
                                    lastNonCoreFeatureIndex=featureCounter-1
                                else: #find the index of the last non-core gene.
                                    for i in range(featureCounter-1,-1,-1):
                                        if (i-1>0 and nonCoreGenes[i]==1 and nonCoreGenes[i-1]==1) or i==0:
                                            lastNonCoreFeatureIndex=i
                                            break
                                for k in range(firstNonCoreFeatureIndex,lastNonCoreFeatureIndex+1):
                                    if sample not in mgeFeatures: 
                                        mgeFeatures[sample]={chr: {}}
                                    if (chr not in mgeFeatures[sample]):
                                        mgeFeatures[sample][chr]={globalFeatureCounter: []}
                                    if globalFeatureCounter not in mgeFeatures[sample][chr]:
                                        mgeFeatures[sample][chr][globalFeatureCounter]=[]

                                    mgeFeatures[sample][chr][globalFeatureCounter].append(gffData[sample][chr][k])
                            
                                globalFeatureCounter+=1


                                #if globalFeatureCounter>20:
                                #    sys.exit()


            if sum(nonCoreGenes)>len(gffData[sample][chr])*0.5: #whole contig is probably MGE
                if sample not in mgeFeatures: 
                    mgeFeatures[sample]={chr: {}}
                #here, any features already identified on chromomosome are completely replaced
                #so the who choromosome is newly added to dictionary to remove previous info. This avoids dulication of MGEs
                mgeFeatures[sample][chr]={globalFeatureCounter: []}
                mgeFeatures[sample][chr][globalFeatureCounter]=gffData[sample][chr]
                globalFeatureCounter+=1



    print("Generating within MGE edges")
    output=open(outputDir+"within.csv","w")
    processedEdges=set()
    for sample in mgeFeatures:
        for chr in mgeFeatures[sample]:
            for feature in mgeFeatures[sample][chr].keys():
                for gene in mgeFeatures[sample][chr][feature]:
                    for targetGene in mgeFeatures[sample][chr][feature]:
                        if gene[0]!= targetGene[0] and gene[1]!= targetGene[1] and len(mgeFeatures[sample][chr][feature])>5:
                            sourceNode=sample+"::"+chr+":"+str(gene[0])+"-"+str(gene[1])
                            targetNode=sample+"::"+chr+":"+str(targetGene[0])+"-"+str(targetGene[1])
                            if (sourceNode+targetNode) not in processedEdges:
                                output.write(sourceNode+"\t"+targetNode+"\tWithin\t"+str(feature)+"\n")
                                processedEdges.add(sourceNode+targetNode)
                                processedEdges.add(targetNode+sourceNode)

    output.close()


