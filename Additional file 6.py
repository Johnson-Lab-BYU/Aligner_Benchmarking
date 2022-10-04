#Used to compare SAM files to SAM files
#both files must be sorted with sort command first
import sys

align1 = open(sys.argv[1], 'r')
align2 = open(sys.argv[2], 'r')
#sam files used entered via user input
CorrectCount = 0
IncorrectCount = 0
matchCount = 0
UnmappedCount = 0
nearMatch = 0
IncorrectLine = []

a1lines = align1.readlines()
a2lines = align2.readlines()
FirstOfPair = []
SecondOfPair = []
Unmapped = []
nearMatchList = []
UnmappedFlag = ['4','69','77','89','101','109','117']
FirstOfPairFlag = ['65','67','73','81','83','97','99','113','115']

for line in a2lines:
    if line[0] != '@':
        split = line.split('\t')
        if split[1] in UnmappedFlag:
            UnmappedCount += 1
            Unmapped.append(line)
        elif split[1] in FirstOfPairFlag:
            FirstOfPair.append(line)
        else:
            SecondOfPair.append(line)

for line in a1lines:
    if line[0] != '@':
        Reference = line.split('\t')
        if Reference[1] in FirstOfPairFlag:
            for item in FirstOfPair:
                Test = item.split('\t')
                if Reference[0] == Test[0]:
                    matchCount += 1
                    if Reference[2] == Test[2]:
                        if Reference[3] == Test[3]:
                            CorrectCount += 1
                        elif Reference[3] != Test[3]:
                            IncorrectCount += 1
                            if abs(int(Reference[3]) - int(Test[3])) in range(1,21):
                                nearMatch += 1
                                nearMatchList.append(line)
                            elif abs(int(Reference[3]) - int(Test[3])) > 20:
                                IncorrectLine.append(line)
                    elif Reference[2] != Test[2]:
                        IncorrectCount += 1
                        IncorrectLine.append(line)
                    FirstOfPair.remove(item)
                    break
                else:
                    break

countlist = str(sys.argv[3])+'_count.txt'
f = open(countlist, 'a+')
f.write('correct:'+'\t'+str(CorrectCount)+'\n'+'Incorrect:'+'\t'+str(IncorrectCount)+'\n'+'Match Count:'+'\t'+str(matchCount)+'\n'+'Unmapped Count:'+'\t'+str(UnmappedCount)+'\n'+'Near Match Count:'+'\t'+str(nearMatch)+'\n')
f.close()

IncorrectList = str(sys.argv[3])+'_incorrectReads.sam'
h = open(IncorrectList, 'a+')
for i in IncorrectLine:
    h.write(i)
h.close()

nearMatchDoc = str(sys.argv[3])+'_NearMatchReads.sam'
t = open(nearMatchDoc, 'a+')
for z in nearMatchList:
    t.write(z)
t.close()

UnmappedReads = str(sys.argv[3])+'_unmappedReads.sam'
g = open(UnmappedReads, 'a+')
for j in Unmapped:
    g.write(j)
g.close()
