#takes file and generate a matrix for PCA analysis. Matrix consists of rows as clique or file ids and genes RefSeq from prokka as columns

import configs.Kp as KpConfig
config=KpConfig.configData()

dataDir=config.outputDir
cliquesAssemblies={} #clique : set(files)
assemblyCliques={} #assembly : set(cliques)

def run(dataDir):
    
    cliquesFile=dataDir+"betweenCliques.csv"

    #collect files names and gene coordinates for each clique
#    assemblyCliques={}
    for line in open(cliquesFile):
        [nodeName, clique]=line.strip().split("\t")

        # if clique not in cliquesRefSeq:
        #     cliquesRefSeq[clique]=[]

        values=nodeName.replace("::",":").split(":")
        fileName=values[0]

        if clique not in cliquesAssemblies:
            cliquesAssemblies[clique]=set()
        cliquesAssemblies[clique].add(fileName)

        if fileName not in assemblyCliques:
            assemblyCliques[fileName]=set()
        assemblyCliques[fileName].add(clique)


    #geneate the matrix
    totalCliques=0
    for clique in cliquesAssemblies.keys():
        totalCliques+=1 if len(cliquesAssemblies[clique])>1 else 0
    
    outputMatrixAssemblies=[]
    assemblyIndices={} #assembly id : row index in the matrix

    for sample in assemblyCliques.keys():
        assemblyIndices[sample]=len(assemblyIndices)
        outputMatrixAssemblies.append([0]*(totalCliques+1))
        outputMatrixAssemblies[-1][0]=sample

    cliqueIndices={}
    for clique in cliquesAssemblies.keys():
        if len(cliquesAssemblies[clique])>1:
            if clique not in cliqueIndices: #skip genes that occur in single sample
                cliqueIndices[clique]=len(cliqueIndices)+1 #+1 due to clique id/sample id column

            for assembly in cliquesAssemblies[ clique ] :
                outputMatrixAssemblies[ assemblyIndices[ assembly ] ][ cliqueIndices[clique] ]=1

    cliqueNames=[0]*(totalCliques+1)
    cliqueNames[0]="Clique"
    for clique in cliqueIndices:
        cliqueNames[cliqueIndices[clique]]=clique


    output=open(dataDir+"/pcaMatrixAssemblies.tsv", "w")
    output.write('\t'.join(str(f) for f in cliqueNames)+"\n")
    for line in outputMatrixAssemblies:
        output.write('\t'.join(str(f) for f in line)+"\n")

run(dataDir)