import SampleMetadata
class Metadata(SampleMetadata.Metadata):
    def __init__(self, vertexType, idvalue):
        self.set_vertexType(vertexType)
        self.set_vertexID(idvalue)
        super(SampleMetadata.Metadata,self).__init__()
        self.clear()
        self.set_product("NA")


    def splitBlastID(self, idValue):
        sampleData=idValue.replace("::",":").split(":")
        sampleFile=sampleData[0]
        sampleContig=sampleData[1]
        starEnd=sampleData[2].split("-")
        super(Metadata,self).set_sampleID(sampleFile)

    def selectProductName(self, products):
        if self.get_vertexType()=="between" or self.get_vertexType()=="gene":
            if len(products)==1:
                self.set_product(products.pop())
                return

            #this applies to between vertices
            sorted_products=sorted(products, key=len, reverse=True) 
            for product in sorted_products:
                if product not in ["hypothetical protein", "putative protein"]:
                    self.set_product(product)
                    return product
            self.set_product(sorted_products[0]) #case where all products are hypothetical
            return sorted_products[0]
        else:
            self.set_product("within")
            return "within"

    def get_vertexID(self):
    	return self.__vertexID
    def set_vertexID(self, vertexID):
    	self.__vertexID = vertexID

    def get_vertexType(self): #within/between/gene
    	  return self.__vertexType
    def set_vertexType(self, vertexType):
    	self.__vertexType = vertexType

    def get_vertexID(self):
    	return self.__vertexID
    def set_vertexID(self, vertexID):
    	self.__vertexID = vertexID

    def get_product(self):
        if self.__vertexType=="between" or self.__vertexType=="gene":
            return self.__product
        else:
            return "NA"
    def set_product(self, product):
        self.__product = product

    def get_geneCount(self):
    	return self.__geneCount
    def set_geneCount(self, geneCount):
    	self.__geneCount = geneCount

    def get_lactamase(self):
    	return self.__lactamase
    def set_lactamase(self, lactamase):
    	self.__lactamase = lactamase

    def get_sampleCount(self):
    	return self.__sampleCount
    def set_sampleCount(self, sampleCount):
    	self.__sampleCount = sampleCount

