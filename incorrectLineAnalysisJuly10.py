#To better understand where incorrect reads are coming from, we segmented those reads into 
#1) Reads that came from the same chromosome
#2) Reads that look like they came from a different chromosome, but actually came from the same one (Art file can report duplicate chromosomes)
#3) Reads from different Chromosomes


import re
import sys
import os

cipherFile  = sys.argv[1] #outputFileJuly17 #Maps Art chromosome names back to their official chromosomes
samFile = sys.argv[2] #newFlag_faster*.sam
artFile = sys.argv[3] #limitedArtFile.sam
outputFolder= sys.argv[4] 
nmNumber=sys.argv[5] #25, 50,100, or 150

print(outputFolder)
nmValue = str(nmNumber)+"M" #150M
cigarScore =str(nmNumber)+"=" #100= is good
chromDic = {}
with open(cipherFile) as file:
    for line in file:
        l = line.strip().split(",")
        chromDic[l[0]]=l[1]

#print(chromDic)
sameChrom = []
sameChromNM= 0
sameChromCigar= 0

sameButDiff= []
sameButDiffNM= 0
sameButDiffCigar= 0

diffChrom = []
diffNM = 0
diffCigar= 0

badChroms = set()

#Read the incorrect Lines file and seperate them into each category.
#NM is not important right now.
with open(samFile) as incorrectFile:
    for line in incorrectFile:
        if line.startswith("@"):
            continue
        
        l = line.strip().split()
        chromMapped = l[0].split("-")[0]
        chromOther = l[2]
        nm = l[5]

        if chromMapped==chromOther:
            sameChrom.append([l[0],l[3]])
            if nm!=nmValue:
                sameChromNM+=1
            continue
        if chromMapped not in chromDic:
            badChroms.add(chromMapped)
            continue
        if chromOther not in chromDic:
            badChroms.add(chromOther)
            continue
        first= chromDic[chromMapped]
        second = chromDic[chromOther]
        
        if first==second:
            sameButDiff.append(l[0])
            if nm!=nmValue:
                sameButDiffNM+=1
            continue
        else:
            diffChrom.append(l[0])
            if nm!=nmValue:
                diffNM+=1 

print(badChroms)
counterSameChrom =0
counterSameButDiff =0
counterDiffChrom=0
over250Dist= 0
under250Dist=0

#Read through the limitedArtFile.sam (Generated earlier) and find the cigar scores+distances from the original position.
with open(artFile) as artF:
    for line in artF:
        l= line.strip().split()

        if counterSameChrom<len(sameChrom):
            if l[0] == sameChrom[counterSameChrom][0]:
                dist =int(sameChrom[counterSameChrom][1])-int(l[3])
                if dist>=-250 and dist<=250:
                    under250Dist+=1
                else:
                    over250Dist+=1
                if l[5]!=cigarScore:
                    sameChromCigar+=1
                counterSameChrom+=1

        if counterSameButDiff<len(sameButDiff):
            if l[0] == sameButDiff[counterSameButDiff]:
                if l[5]!=cigarScore:
                    sameButDiffCigar+=1
                counterSameButDiff+=1

        if counterDiffChrom<len(diffChrom):
            if l[0] == diffChrom[counterDiffChrom]:
                if l[5]!=cigarScore:
                    diffCigar+=1
                counterDiffChrom+=1


outputFile =os.path.join(outputFolder,"newFlagLineAnalysisOutput")
print(outputFile)
with open(outputFile,"w") as output:
    output.write(f"SameCount: {len(sameChrom)} Bad NM: {sameChromNM} Bad Cigar: {sameChromCigar} Over250Dist: {over250Dist} Under250Dist: {under250Dist}\n")
    output.write(f"SameButDiffCount: {len(sameButDiff)} Bad NM: {sameButDiffNM} Bad Cigar: {sameButDiffCigar}\n")
    output.write(f"DiffCount: {len(diffChrom)} Bad NM: {diffNM} Bad Cigar: {diffCigar}\n")

print(f"SameCount: {len(sameChrom)} Bad NM: {sameChromNM} Bad Cigar: {sameChromCigar} Over250Dist: {over250Dist} Under250Dist: {under250Dist}\n")
print(f"SameButDiffCount: {len(sameButDiff)} Bad NM: {sameButDiffNM} Bad Cigar: {sameButDiffCigar}\n")
print(f"DiffCount: {len(diffChrom)} Bad NM: {diffNM} Bad Cigar: {diffCigar}\n")
