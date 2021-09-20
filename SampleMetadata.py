class Metadata:
    def __init__(self, sampleID, location, ST, year):
        self.__ST = ST
        self.__year = year
        self.__location=location
        self.__sampleID=sampleID

    def clear(self):
        self.__ST = ""
        self.__year = ""
        self.__location=""
        self.__sampleID=""

    def get_sampleID(self):
        return self.__sampleID
    def set_sampleID(self, sampleID):
    	self.__sampleID = sampleID

    def get_ST(self):
        return self.__ST
    def set_ST(self, ST):
    	self.__ST = ST

    def get_year(self):
        return self.__year
    def set_year(self, year):
    	self.__year = year

    def get_location(self):
        return self.__location
    def set_location(self, location):
        self.__location = location

    def get_virulence(self):
    	return self.__virulence
    def set_virulence(self, virulence):
    	self.__virulence = virulence

    def get_kLocus(self):
    	return self.__kLocus
    def set_kLocus(self, kLocus):
    	self.__kLocus = kLocus

    def get_oLocus(self):
    	return self.__oLocus
    def set_oLocus(self, oLocus):
    	self.__oLocus = oLocus