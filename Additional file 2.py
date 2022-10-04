#used to list all the FLAG scores in a SAM file
import csv
import sys

align = open(sys.argv[1], 'r')
SamCSV_Flag = str(sys.argv[2])+"_SamCSV_Flag.csv"

flag = []
for line in align:
    if line[0] != '@':
        alignSplit = line.split('\t')
        if alignSplit[1] not in flag:
            flag.append(alignSplit[1])

with open(SamCSV_Flag, 'wb') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows([flag])
