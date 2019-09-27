###################################################
###################################################
###################################################
# classes for use with pipeline_primerdesign.py

###################################################
###################################################
###################################################

import cgatcore.iotools as IOTools

class PrimerSet(object):
    '''
    class for a primer set
    '''
    def __init__(self):

        self.name = None
        self.forwardseq = None
        self.reverseseq = None
        self.forwardtm = 0
        self.reversetm = 0
        self.size = 0
        self.forwardgc = 0
        self.reverse = 0
        self.forwardlength = 0
        self.reverselength = 0


    def readName(self, inf):
        inf = IOTools.open_file(inf)
        return inf.readline()[len("PRIMER PICKING RESULTS FOR "):-1]
        
    def readForward(self, inf):
        inf = IOTools.open_file(inf)
        for line in inf.readlines():
            if line.startswith("LEFT PRIMER"):
                data = line[:-1].split(" ")
                data = [x for x in data if x != ""]
                forward_seq, forward_gc, forward_tm, forward_len = data[10], data[5], data[4], data[3]
        return forward_seq, forward_gc, forward_tm, forward_len

    def readReverse(self, inf):
        inf = IOTools.open_file(inf)
        for line in inf.readlines():
            if line.startswith("RIGHT"):
                data = line[:-1].split(" ")
                data = [x for x in data if x != ""]
                reverse_seq, reverse_gc, reverse_tm, reverse_len = data[10], data[5], data[4], data[3]
        return reverse_seq, reverse_gc, reverse_tm, reverse_len

    def readSize(self, inf):
        inf = IOTools.open_file(inf)
        for line in inf.readlines():
            if line.startswith("PRODUCT SIZE"):
                data = line[:-1].split(" ")
                data = [x for x in data if x != ""]
                size = data[2].strip(",")
        return size

    def parse(self,
              attributes=[]):
        self.name, self.size, self.forwardseq, self.forwardgc, self.forwardtm, self.forwardlength, self.reverseseq, self.reversegc, self.reversetm, self.reverselength = attributes
        return self

