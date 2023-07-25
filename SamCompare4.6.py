#Most recent version of the SamCompare file.
#This file reads through an Art file (with the original reads) and a generated sam file from an Aligner (bowtie2,BWA, Chromap, GSnap, or Subread)
#And then compares them to find out which reads are "Correctly aligned", "Incorrectly aligned", "half-correctly aligned", or
#Unmapped. 
#Count files are generated, and incorrect lines are printed out for further downstream analysis.
import sys
import time
import os
import glob
start = time.time()

artFile = sys.argv[1] #Art_hs_GRCh38_*_*_sort.sam 
samFolder = sys.argv[2] 
outputDir = sys.argv[3]
restartFile = sys.argv[4]
countDataFile = sys.argv[5]


UnmappedFlag = {'69', '73','77', '85', '87','89','93', '101', '103', '117'} #First of pair that are unmapped
FirstOfPairFlag = {'65', '67', '81', '83', '97', '99', '113', '115'} #First of pair that aligned

samFile = glob.glob(f"{samFolder}/*.sam")[0] #Only one sam file in this folder
#samFile= glob.glob(f"{samFolder}/*shorter.sam")[0] #For Chromap

print(samFile)

IncorrectLine = []
temp = []

#Count all of these
CorrectCount = 0
IncorrectCount = 0 #Both aligned to the wrong place
matchCount = 0 #total Count
UnmappedCount = 0
halfCorrectCount=0 #one correctly aligned, while the other incorrectly aligned.

#Read through the Aligned file, add any mapped reads into a temp nested list.

with open(samFile, 'r') as a2file:
    for line in a2file:
        if line[0] != '@':
            split = line.split('\t')
            if split[1] in UnmappedFlag:
                UnmappedCount += 1
            elif split[1] in FirstOfPairFlag:
                #read name, chromosome, position1, position2 (end pos)
                temp.append([split[0],split[2],split[3],split[7]])

middle = time.time()
print("Read samFile")
print(middle - start)
print("temp len", len(temp))

#Read the art file, compare the temp list of mapped reads against this.
counter=0
with open(artFile, 'r') as a1file:
    for index, line in enumerate(a1file):
    
        artfile = line.split('\t')
        if  artfile[1] in FirstOfPairFlag:

            #counter is for the temp, after we find all the mapped reads, we don't need to do anything.
            if counter== len(temp):
                continue


            samLine = temp[counter]
            if artfile[0]== samLine[0]:
                matchCount += 1
                #If the name, chromosome, and positions line up, it is correctly mapped.
                if artfile[2]== samLine[1] and artfile[3]== samLine[2] and artfile[7]==samLine[3]:
                    CorrectCount +=1
                #if both positions are wrong, it is incorrect
                elif artfile[3]!= samLine[2] and artfile[7]!=samLine[3]: #the first doesnt line up
                    IncorrectCount += 1
                    IncorrectLine.append(samLine[0])
                #If only one position is wrong, it is half correct (likely means that it is in a repeat region)
                elif artfile[3]!= samLine[2] or artfile[7]!=samLine[3]: 
                    halfCorrectCount+=1
                else:
                    IncorrectCount += 1
                    IncorrectLine.append(samLine[0])
                    print(samLine[0])
                counter+=1

middle1=time.time()
print("Compare files done done")
print(middle1 - middle)
print(len(IncorrectLine))

outputFileEnd = "newFlag_faster"+samFile.split("/")[-1]
outputFile = os.path.join(outputDir, outputFileEnd)


prevLine= "0 0 0"
isSecond= False

#From this IncorrectLine list, we only have the first in pair.
#Go back to the samFile (aligned file) to find the second in pair and write those to a new file.
with open(outputFile, "w") as outFile:
    with open(samFile) as file:
        incorrectLineCount=0
        for line in file:
            if not line.startswith("@"):
                l = line.split("\t")
                if incorrectLineCount== len(IncorrectLine):
                    continue
                
                if l[0]== prevLine.split("\t")[0]:
                    isSecond=True
                else: 
                    isSecond=False
                    
                if isSecond: #the first in pair
                    if l[0] == IncorrectLine[incorrectLineCount]:
                        outFile.write(prevLine)
                        outFile.write(line)
                        incorrectLineCount+=1

                prevLine = line

#Write all count data.
with open(countDataFile, "w") as countFile:
    countFile.write(f"IncorrectCount:\t{IncorrectCount}\n")
    countFile.write(f"HalfCorrectCount:\t{halfCorrectCount}\n")
    countFile.write(f"CorrectCount:\t{CorrectCount}\n")
    countFile.write(f"UpmappedCount:\t{UnmappedCount}\n")
    countFile.write(f"TotalCount:\t{matchCount}\n")


end=time.time()
print("Write files. End Time: ",end - middle1)
print("IncorrectCount: ", IncorrectCount)
print("HalfCorrectCount: ",halfCorrectCount)
print("CorrectCount: ", CorrectCount)
print("Unmapped Count: ", UnmappedCount)
print("TotalCount: ", matchCount)
