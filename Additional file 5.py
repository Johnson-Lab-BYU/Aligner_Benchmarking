#used to compare shrimp1 fixed output files or Maq fixed output files to SAM files
#both files must be sorted with sort command first
import sys

align1 = open(sys.argv[1], 'r')
align2 = open(sys.argv[2], 'r')

CorrectCount = 0
IncorrectCount = 0
matchCount = 0
nearMatch = 0
IncorrectLine = []
nearMatchList = []

a1lines = align1.readlines()
a2lines = align2.readlines()
sortedReads = sorted(a2lines)
FirstOfPairFlag = ['65','67','73','81','83','97','99','113','115']

for line in a1lines:
    if line[0] != '@':
        Reference = line.split('\t')
        if Reference[1] in FirstOfPairFlag:
            for item in sortedReads:
                Test = item.split('\t')
                if Reference[0] == Test[0]:
                    matchCount += 1
                    if Reference[2] == Test[1]:
                        if int(Reference[3]) == int(Test[2]):
                            CorrectCount += 1
                        elif int(Reference[3]) != int(Test[2]):
                            IncorrectCount += 1
                            if abs(int(Reference[3]) - int(Test[2])) in range(1,21):
                                nearMatch += 1
                                nearMatchList.append(line)
                            elif abs(int(Reference[3]) - int(Test[2])) > 20:
                                IncorrectLine.append(line)
                    elif Reference[2] != Test[1]:
                        IncorrectCount += 1
                        IncorrectLine.append(line)
                    sortedReads.remove(item)
                    break
                else:
                    break

countlist = str(sys.argv[3])+'_count.txt'
f = open(countlist, 'a+')
f.write('correct:'+'\t'+str(CorrectCount)+'\n'+'Incorrect:'+'\t'+str(IncorrectCount)+'\n'+'Match Count:'+'\t'+str(matchCount)+'\n'+'Near Match:'+'\t'+str(nearMatch)+'\n')
f.close()

IncorrectList = str(sys.argv[3])+'_incorrectReads.txt'
h = open(IncorrectList, 'a+'
for i in IncorrectLine:
    h.write(i)
h.close()

nearMatchDoc = str(sys.argv[3])+'_NearMatchReads.sam'
t = open(nearMatchDoc, 'a+')
for z in nearMatchList:
    t.write(z)
t.close()
