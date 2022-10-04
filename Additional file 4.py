#used to make Maq output compairable to SAM files
import sys

forwardMapped = open(sys.argv[1], 'r')

unmapped = []
pairedRead = []
unmappedCount = 0
pairedReadCount = 0

fMapped = forwardMapped.readlines()

with open('Maq_fix_output.txt', 'w') as f:
    for item in fMapped:
        split = item.split('\t')
        key = split[0]
        if key[-2:] == '/1':
            if split[5] == '18':
                f.write(key[:-2]+'\t'+split[1]+'\t'+split[2]+'\n')
            elif split[5] != '18':
                unmapped.append(item)
                unmappedCount += 1
        elif key[-2:] == '/2':
            pairedRead.append(item)
            pairedReadCount += 1

with open('Maq_unmapped.txt', 'w') as g:
    g.write('Unmapped Count:'+str(unmappedCount)+'\n')
    for line in unmapped:
        g.write(line)

with open('Maq_paired.txt', 'w') as g:
    g.write('Paired Read Count:'+str(pairedReadCount)+'\n')
    for line in pairedRead:
        g.write(line)
