import sys
import subprocess


vertexName=sys.argv[1] #either W### or B####
blastSubject=sys.argv[2] #target for search, fasta file

if vertexName.find("W")>-1:
    cliquesList=wd+"/withinCliques.csv"
else:
    cliquesList=wd+"/betweenCliques.csv"

regionStart=""
regionEnd=""
totalCliques=0
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

        gff=subprocess.run("bedtools intersect -a "+sampleGffDir+fileName+".gff -b <(echo -e \""+chr+"\t"+start+"\t"+end+"\")", shell=True, executable="/bin/bash", stdout=subprocess.PIPE)

        gffLine=gff.stdout.decode().split("\t")
        if len(gffLine)>1:
            annotation=gffLine[8].strip().split(";")
            RefSeq=""
            Product=""
            for value in annotation:
                if value.find("RefSeq:")>-1:
                    RefSeq=value.split(":")[-1]
                if value.find("product=")>-1:
                    Product=value.split("=")[1]
            subprocess.call("echo '>"+vertexName+"' >> temp.fasta ", executable="/bin/bash", shell=True)				
            subprocess.call("bedtools getfasta -fi "+samplesDir+fileName+".fasta -bed <(echo -e \""+chr+"\t"+str(min(int(start), int(end)))+
                    "\t"+str(max(int(start), int(end)))+"\") | grep -v '>' >> temp.fasta", executable="/bin/bash", shell=True)
            subprocess.call("~/generalScripts/quickblast.sh "+blastSubject+" temp.fasta", executable="/bin/bash", shell=True)
            subprocess.call("rm temp.fasta", executable="/bin/bash", shell=True)
            print(Product+";"+RefSeq)
            break
            

