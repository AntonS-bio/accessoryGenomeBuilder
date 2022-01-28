# config={
# "wd":"directory in which any work will be done file will be",
# "outputDir": "directory where final output files go",
# "tempDir":"directory for temp files",
# "refSample": "FileName without .fasta. This will not impact results, algorithm simply takes genes from it to start finding core genes",
# "samplesDir": "Directory with fasta files",
# "sampleGffDir": "Directory with gff files",
# "AMRDir": "Directory with AMR genes (optional)",
# "AMRfile": "NCBI's AMR ref genes file ex. NCBI_refgenes_16-Jul-20.csv",
# "postfix": "extension of fasta files (fna/fasta)",
# "blastIdentityThreshold":95, #min required identity between genes to be considered the same
# "coreGenomeThreshold": 0.01, #share of samples that need to have genes for it to be "core"
# "permittedLengthVariance": 0.05 #how much lenght of gene can deviate from median lenght of homologues. More than this mean genes leads to genes being classified as different
# }



class configData:
    def __init__(self):
        self.wd = "some/dir"  #directory in which any work will be done file will be
        self.outputDir = self.wd+"Kp"  #directory where final output files go"
        self.tempDir=self.wd+"Kp"+"temp" #directory for temp files
        self.refSample="some/dir" #FileName without .fasta. This will not impact results, algorithm simply takes genes from it to start finding core genes
        self.samplesDir="some/dir"  #Directory with fasta files
        self.sampleGffDir="some/dir" #Directory with gff files
        self.AMRDir="AMRDir" #Directory with fasta for AMR genes (optional)
        self.AMRfile="NCBI_refgenes_16-Jul-20.csv" #"NCBI's AMR ref genes file (optional)"
        self.postfix=".fasta" #extension of fasta files (fna/fasta)
        self.blastIdentityThreshold=95 #min required identity between genes to be considered "same" genes and grouped together 
        self.coreGenomeThreshold=0.01 #share of samples that need to have genes for it to be "core"
        self.permittedLengthVariance=0.05 #how much lenght of gene can deviate from median lenght of homologues. More than this mean genes leads to genes being classified as different        


    # def get_sampleID(self):
    #     return self.__sampleID
    # def set_sampleID(self, sampleID):
    # 	self.__sampleID = sampleID

    # def get_ST(self):
    #     return self.__ST
    # def set_ST(self, ST):
    # 	self.__ST = ST