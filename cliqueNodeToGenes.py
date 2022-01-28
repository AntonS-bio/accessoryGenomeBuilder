import sys
import subprocess
#from MGE.scripts.accessoryGenome.configs.Kp import configData
import configs.Liver_60 as KpConfig

config=KpConfig.configData()


# config=KpConfig.configData()
# wd=config.wd
# samplesDir=config.samplesDir
# sampleGffDir=config.sampleGffDir


vertexName=sys.argv[1]  # either W## or B####
wd=config.wd
samplesDir=config.samplesDir

if vertexName.find("W")>-1:
    cliquesList=config.outputDir+"/withinCliques.csv"
else:
    cliquesList=config.outputDir+"/betweenCliques.csv"

regionStart=""
regionEnd=""
totalCliques=0
subprocess.call("echo '' > cliqueGenes.fasta ", executable="/bin/bash", shell=True)				
for line in open(cliquesList):
    if line.strip().split("\t")[1]==vertexName:
        totalCliques+=1
        nodeName=line.strip().split("\t")[0]

        values=nodeName.replace("::",":").split(":")
        fileName=values[0]
        chr=values[1]
        [start, end]=values[2].split("-")
        if regionStart=="":
                regionStart=start
        regionEnd=end

        gff=subprocess.run(f'bedtools intersect -a {config.sampleGffDir+fileName}.gff -b <(echo -e \""{chr}"\t"{start}"\t"{end}"\")', shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

        gffLine=gff.stdout.decode().split("\t")
        product=gffLine[8]
        orientation=gffLine[6]
        #start=gffLine[3]
        #end=gffLine[4]
        #print('\t'.join(f for f in [fileName]+gffLine))

        ##to get all B-type nodes fastas into same file
        #subprocess.call(f'echo \'>{fileName}\' >> cliqueGenes.fasta', executable="/bin/bash", shell=True)
        #bedline=f'{chr}\t{str(min(int(start), int(end)))}\t{str(max(int(start), int(end)))}\tNone\t.\t{orientation}'
        #subprocess.call(f'bedtools getfasta -s -fi {samplesDir}{fileName}.fasta -bed <(echo -e "{bedline}") | grep -v ">" >> cliqueGenes.fasta ', executable="/bin/bash", shell=True)				


#print("Total cliques connected: "+str(totalCliques))
#print(">"+vertexName)
#print(product)
if vertexName.find("B")>-1:
    pass
    #subprocess.call("echo '>"+vertexName+ "' > "+wd+"geneFastas/"+vertexName+".fasta", executable="/bin/bash", shell=True)				
    subprocess.call("bedtools getfasta -s -fi "+config.samplesDir+fileName+
        ".fasta -bed <(echo -e \""+chr+"\t"+str(start)+
        "\t"+str(end)+"\t"+vertexName+"\t1\t"+orientation+"\") | grep -v '>' ", 
        executable="/bin/bash", shell=True)				
    #subprocess.call("bedtools getfasta -fi "+samplesDir+fileName+".fasta -bed <(echo -e \""+chr+"\t"+str(min(int(start), int(end)))+"\t"+str(max(int(start), int(end)))+"\") | grep -v '>' >> "+wd+"geneFastas/"+vertexName+".fasta", executable="/bin/bash", shell=True)				
if vertexName.find("W")>-1:
    subprocess.call("bedtools getfasta -fi "+config.samplesDir+fileName+".fasta -bed <(echo -e \""+chr+"\t"+str(min(int(regionStart), int(regionEnd)))+"\t"+str(max(int(regionStart), int(regionEnd)))+"\") | grep -v '>' ", executable="/bin/bash", shell=True)				
#subprocess.call("echo -e '>"+vertexName+"' >> withinClique.fasta", executable="/bin/bash", shell=True)		
#subprocess.call("bedtools getfasta -fi "+samplesDir+fileName+".fasta -bed <(echo -e \""+chr+"\t"+str(min(int(regionStart), int(regionEnd)))+"\t"+str(max(int(regionStart), int(regionEnd)))+"\") | grep -v '>' >> withinClique.fasta", executable="/bin/bash", shell=True)		
        
