#####################################################
#####################################################
#####################################################
# script to convert the output from the ABI machine
# into nice formatted table
#####################################################
#####################################################
#####################################################


import sys
import string

assert sys.argv[1] and sys.argv[2], "must specify ABI output and plate layout files"

inf = sys.argv[1]
layout = sys.argv[2]

###############################
# function definitions
###############################

def buildWellToSample(layout):
    '''
    build a container mapping
    well to sample id
    '''
    columns = range(1,25,1)
    rows = list(string.ascii_uppercase[0:16])
    well2sample = {}
    c = 0
    for line in layout.readlines():
        data = line[:-1].rstrip().split("\t")
        assert len(data) == 24, "plate layout in wrong format, please review input"
        for i in range(len(data)):
            well = rows[c] + str(i + 1)
            if data[i] == "":
                sample = "NA"
            else:
                sample = data[i]
            well2sample[well] = sample
        c += 1
    return well2sample

################################
################################
################################

def buildWellToValueAndGene(infile):
    '''
    build a mapping of well to Ct value
    '''
    # skip lines starting with the following
    # strings
    skip_starting = ["Block",
                     "Calibration",
                     "Chemistry",
                     "Experiment",
                     "Instrument",
                     "Passive",
                     "Quantification",
                     "Signal",
                     "Stage",
                     "Well"]
    well2ct = {}
    well2gene = {}
    for line in infile.readlines():
        if not line: 
            continue
        if True in [line.startswith(x) for x in skip_starting]:
            continue
        data = line[:-1].split("\t")
        well, gene, value = data[1], data[3], data[4] 
        well2ct[well] = value
        well2gene[well] = gene

    return well2ct, well2gene


################################
################################
################################

well2sample = buildWellToSample(open(layout))
well2ct, well2gene = buildWellToValueAndGene(open(inf))

sys.stdout.write("well\tsample\tCt\tgene\n")
for well, sample in well2sample.iteritems():
    if sample == "NA":
        continue
    
    sys.stdout.write("%s\t%s\t%s\t%s\n" % (well, sample, well2ct[well], well2gene[well]))

