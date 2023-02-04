# Harlan Stevens
import re #regex
import sys

bowtieFilePath = sys.argv[1]  #One the commandline, input your bowtie file name.


#Read through the bowtie file and check for all the reads that did not align perfectly. Count the edit distances for each.
def parseSamMisreads(fileName):
    countDict= {}
    with open(fileName, "r") as file:
        for line in file:
            if line.startswith("CM"): #make sure line was valid read and not a header
                regex = r"NM:i:(\d+)" #return numbers from NM:i:someNumber. According to Bowtie Manual this says:
                # "The edit distance; that is, the minimal number of one-nucleotide edits
                #(substitutions, insertions and deletions) needed to transform the read string into the reference string.
                #Only present if SAM record is for an aligned read.

                NM = re.findall(regex, line) #Should only exist once
                if len(NM)>=1: 
                    totalMisreads = int(NM[0])
                    if totalMisreads in countDict:
                        countDict[totalMisreads]+=1
                    else:
                        countDict[totalMisreads]=1
        return countDict
                
countDict = parseSamMisreads(bowtieFilePath) #Return as dictionary
print(countDict)