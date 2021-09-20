def run(outputDir):

    def saveCliqueDic(fileName, vertices, connectionType): #connection type is B or W for between/within
        output=open(outputDir+fileName, "w")
        for vertex in vertices:
            output.write(vertex+"\t"+connectionType+str(vertices[vertex])+"\n")
        output.close()

    cliques={}
    betweenCliques={} #vertex id - between clique


    print("Loading betweenEdges")
    for line in open(outputDir+"betweenEdges.csv"):
        values=line.strip().split("\t")
        if values[0] not in betweenCliques:
            cliques[len(betweenCliques)]=set([values[0]])
            betweenCliques[values[0]]=len(betweenCliques)
        if values[1] not in betweenCliques:
            cliques[len(betweenCliques)]=set([values[1]])        
            betweenCliques[values[1]]=len(betweenCliques)


    #initially, every vertex is it's own clique
    mergedCliques=1
    print("Generating between cliques")
    while mergedCliques!=0:
        print(mergedCliques)
        mergedCliques=0
        for line in open(outputDir+"betweenEdges.csv"):
            values=line.strip().split("\t")
            if betweenCliques[values[0]]!=betweenCliques[values[1]]:
                #copy all elements of clique 2 into clique 1
                cliques[ betweenCliques[values[0]] ] = cliques[betweenCliques[values[0]]]|cliques[betweenCliques[values[1]]]
                cliqueToDelete=betweenCliques[values[1]]
                for vertex in cliques[betweenCliques[values[1]]]:
                    betweenCliques[vertex]=betweenCliques[values[0]]
                del cliques[cliqueToDelete]
                mergedCliques+=1



    crossCliques=set()
    withinCliques={}
    betweenCliqueLinksCount={} #betwenn clique ID B# -> count of connecting W#
    print("Generating within cliques")
    for line in open(outputDir+"within.csv"):
        values=line.strip().split("\t")
        if values[0] in betweenCliques:
            crossCliques.add("W"+values[3]+"\tB"+str(betweenCliques[values[0]]))
            withinCliques[values[0]]=values[3]
        elif values[1] in betweenCliques:
            crossCliques.add("W"+values[3]+"\tB"+str(betweenCliques[values[1]]))
        
    
        withinCliques[values[0]]=values[3]
        withinCliques[values[1]]=values[3] #each values is a pair of two genes, withinCliques is a dictionary key: gene, value=withinClique. 

    saveCliqueDic("betweenCliques.csv",betweenCliques, "B")
    saveCliqueDic("withinCliques.csv",withinCliques,"W")

    #this is not the best way of doing it; the aim is to remove those B# which connect to only one W#
    betweenLinksCount={}
    for edge in crossCliques:
        values=edge.split("\t")
        if values[1] not in betweenCliqueLinksCount:
            betweenCliqueLinksCount[values[1]]=0
        betweenCliqueLinksCount[values[1]]+=1


    output=open(outputDir+"crossCliqies.csv", "w")
    for edge in crossCliques:
        values=edge.split("\t")
        if betweenCliqueLinksCount[values[1]]>1:
            output.write(edge+"\n")

    output.close()
    output=open(outputDir+"allCliques.csv", "w")
    for vertex in set(betweenCliques.keys()):
        output.write("B"+str(betweenCliques[vertex])+"\t"+str(vertex)+"\tBetween\n")
    for vertex in set(withinCliques.keys()):
        output.write("W"+str(withinCliques[vertex])+"\t"+str(vertex)+"\tWithin\n")
    output.close()

