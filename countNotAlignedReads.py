# Harlan Stevens
import sys
import re
import os
from collections import defaultdict
import copy

#Inputs from commandline: python3 countNotAlignedReads.py <artFile path> <bowTieFile path> <readLength>
artFile = sys.argv[1]
bowTieFile = sys.argv[2]
readLength= sys.argv[3]
folder= os.path.split(artFile)[0]

#output art file.
newArtFile = os.path.join(folder, "new_art_file.sam")

#Flag scores in the bowtie file that indicate it the read aligned multiple times
FLAGSCORES = [371,435,355,403,339,419,323,387]

#Read through bowtie file and return a set of all reads that aligned multiple times
def findNotAlignedReads(bowTieFile):
    setOfMultipleAligns = set()
    with open(bowTieFile) as file1:
        for line in file1:
            if line[0]!= "@":
                l = line.split("\t")
                if int(l[1]) in FLAGSCORES:
                    setOfMultipleAligns.add(l[0])
    
    return setOfMultipleAligns

#Create a new ART file that only contains the reads that aligned multiple times
def limitArtFile(artFile, newArtFile, setOfMultipleAligns):
    with open(artFile) as file:
        with open(newArtFile, "w+") as outputFile:
            for line in file:
                if line[0]!="@":
                    l = line.split("\t")
                    if l[0] in setOfMultipleAligns:
                        outputFile.write(line)
    
#Read the new art file that only contains the multiple aligned reads and turn it into a dictionary with the relevent data.
def parseNewFile(newArtFile):
    infoDic = defaultdict(list)
    with open(newArtFile) as file:
        for line in file:
            l = line.split("\t")
            infoDic[l[0]].append([l[1], l[2], l[3], l[5]])
            
    return infoDic #dictionaryOfListOfLists
            
    
    
#Read through bowtie file and count how many reads are aligned  correctly, incorrectly, or do not align.
def countNotAlignedReads(bowTieFile, infoDic, setOfMultipleAligns):
    artDic = copy.deepcopy(infoDic)  #dictionary of relevent data from the newArtFile

    countDic = {"correct": 0, "incorrect":0, "noMatch":0}
        
    #Final dictionary will contain info on if it aligned correctly because the sequence was better,
    #if it aligned correctly by chance, or if it had a worse sequence but still aligned correctly.
    #The same is done with the reads that aligned to the wrong sequence.
    #If the read did not align, we check if there is a mismatch in the sequence or not.
    moreInfoDic = {"correct": {"chance":0, "betterSeq":0, "worseSeq":0},
                   "incorrect": {"chance":0, "betterSeq":0, "worseSeq":0},
                   "noMatch":{"hasMismatch":0, "noMismatch":0}
    }
    
    regex = r"NM:i:(\d+)"
    
    
    prevName = ""
    flag = "no match"
    wentThroughNames = False    
    totalNMCorrect = 0
    totalNMIncorrect = 0
    
    mismatchCheck = str(readLength)+"="
    
    #Read through the bowtie file
    with open(bowTieFile) as file1:
        for line in file1:
            if line[0]!= "@":
                l = line.split("\t")
                #Check if read is aligned multiple times
                if l[0] in setOfMultipleAligns:
                    currentName = l[0] #read name
                    if prevName != currentName:
                        #if the read is new, fill in data from the previous read.
                        if flag == "correct": #if the last read was "correct"
                            if totalNMCorrect<totalNMIncorrect:
                                moreInfoDic["correct"]["betterSeq"]+=1
                            elif totalNMCorrect==totalNMIncorrect:
                                moreInfoDic["correct"]["chance"]+=1
                            else:
                                moreInfoDic["correct"]["worseSeq"]+=1
                                
                        elif flag =="incorrect": #if the last read was "incorrect"
                            if totalNMCorrect<totalNMIncorrect:
                                moreInfoDic["incorrect"]["worseSeq"]+=1
                            elif totalNMIncorrect==totalNMCorrect:
                                moreInfoDic["incorrect"]["chance"]+=1
                            else:
                                moreInfoDic["incorrect"]["betterSeq"]+=1

                                
                        elif flag=="no match" and wentThroughNames: #if the last read was "no match" and it existed (no match is default)

                            countDic["noMatch"]+=1
                            #Check if the read that didn't align has mismatches by going through the art dictionary.
                            hasMismatch = False
                            for i in artDic[prevName]:
                                if i[3]!=mismatchCheck:
                                    hasMismatch=True
                                    
                            if hasMismatch:
                               moreInfoDic["noMatch"]["hasMismatch"]+=1
                            else:
                                moreInfoDic["noMatch"]["noMismatch"]+=1
                        
                        #After checking those and updating the dictionary, clear the data.
                        wentThroughNames = False    
                        flag = "no match"
                        totalNMCorrect = 0
                        totalNMIncorrect = 0
                        prevName = currentName
                    
                    #Find if there are any mismatches, which NM:i:(\d+) tells.
                    mismatches = re.findall(regex, line)
                    
                    #l[1] is the flagscore. If it is not there, it should have aligned correct.                    
                    if int(l[1]) not in FLAGSCORES:
                        if len(mismatches)>=1: #should never be more than once!
                            totalNMCorrect+= int(mismatches[0])
                        
                        #To be sure, go through the art dictionary and remove the corresponding reads.
                        for i in artDic[currentName]:
                            if i[2] == l[3] and i[1] == l[2]:
                                flag = "correct"  #should never exist if the flag was "incorrect" before
                                countDic["correct"]+=0.5 #should be there twice!
                                artDic[currentName].remove(i)
                                break
                    #l[1] is the flagscore. If it is  there, it aligned incorrectly.                
                    if int(l[1]) in FLAGSCORES:
                        if len(mismatches)>=1: #should never be more than once!
                            totalNMIncorrect+= int(mismatches[0])
                        
                        #To be sure, go through the art dictionary and remove the corresponding reads.
                        for i in artDic[currentName]:
                            if i[2]== l[3] and i[1] == l[2]:
                                flag = "incorrect" #should never exist if the flag was "correct" before
                                countDic["incorrect"]+=0.5 #should be there twice!
                                artDic[currentName].remove(i)
                                break
        
                    wentThroughNames = True    
    
    #At the end of the file, check for the last read that could have aligned multiple times.
    if flag == "correct":
        if totalNMCorrect<totalNMIncorrect:
            moreInfoDic["correct"]["betterSeq"]+=1
        elif totalNMCorrect==totalNMIncorrect:
            moreInfoDic["correct"]["chance"]+=1
        else:
            moreInfoDic["correct"]["worseSeq"]+=1
                                
    elif flag =="incorrect":
        if totalNMCorrect<totalNMIncorrect:
            moreInfoDic["incorrect"]["worseSeq"]+=1
        elif totalNMIncorrect==totalNMCorrect:
            moreInfoDic["incorrect"]["chance"]+=1
        else:
            moreInfoDic["incorrect"]["betterSeq"]+=1
                                
    elif flag=="no match":
        countDic["noMatch"]+=1
        hasMismatch = False
        for i in artDic[currentName]:
            if i[3]!=mismatchCheck:
                hasMismatch=True
        if hasMismatch:
            moreInfoDic["noMatch"]["hasMismatch"]+=1
        else:
            moreInfoDic["noMatch"]["noMismatch"]+=1

    #return the dictionaries to be printed and analyzed.
    return countDic, moreInfoDic


#Call the four functions, returning the two dictionaries at the end.

setOfMultipleAligns = findNotAlignedReads(bowTieFile)
limitArtFile(artFile, newArtFile, setOfMultipleAligns)
infoDic = parseNewFile(newArtFile)
countDic, moreInfoDic = countNotAlignedReads(bowTieFile, infoDic, setOfMultipleAligns)
print(countDic)
print(moreInfoDic)