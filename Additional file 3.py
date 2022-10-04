#used to make SHRiMP1 output compairable to SAM files
import sys

forwardUnmapped = open(sys.argv[1], 'r')
reverseUnmapped = open(sys.argv[2], 'r')
forwardMapped =  open(sys.argv[3], 'r')

unmapped = []
unmappedCount = 0

for line in forwardUnmapped:
    unmapped.append(line[:-3])
    unmappedCount += 1

for line in reverseUnmapped:
    if line[:-3] not in unmapped:
        unmapped.append(line[:-3])
        unmappedCount +=1

fMapped = forwardMapped.readlines()

with open('SHRiMP1_fix_output.txt', 'w') as f:
    for item in fMapped:
        split = item.split('\t')
        key = split[0]
        chrome = split[2]
        if key[:-3] not in unmapped:
            chromeSplit = chrome.split(' ')
            f.write(key[:-2]+'\t'+chromeSplit[0]+'\t'+split[3]+'\n')
