import sys
def samPositions(samfile):
    infile = open(samfile,'r')
    geneDict = {}
    for line in infile:
        if line[0] != '@':
            lineSplit = line.split('\t')
            query = lineSplit[0]
            querySplit = query.split('_')
            gene = querySplit[0]
            scaffold = lineSplit[2]
            start = int(lineSplit[3])
            cigar = lineSplit[5]
            totalLength = cigarCalc(cigar)
            stop = start + totalLength
            if scaffold != '*':
                scaffDict = {scaffold:(start,stop)}
                if gene not in geneDict:
                    geneDict[gene] = scaffDict
                else:
                    oldScaffDict = geneDict[gene]
                    if scaffold not in oldScaffDict:
                        oldScaffDict[scaffold] = (start,stop)
                    else:
                        scaffInfo = oldScaffDict[scaffold]
                        oldStart = scaffInfo[0]
                        oldStop = scaffInfo[1]
                        if start < oldStart:
                            oldStart = start
                        if stop > oldStop:
                            oldStop = stop
                        oldScaffDict[scaffold] = (oldStart,oldStop)
                    geneDict[gene] = oldScaffDict
    infile.close()
    for gene in geneDict:
        finalScaffDict = geneDict[gene]
        for scaffold in finalScaffDict:
            scaffInfo = finalScaffDict[scaffold]
            start = scaffInfo[0]
            stop = scaffInfo[1]
            sys.stdout.write(gene + '\t' + scaffold + '\t' + str(start) + '\t' + str(stop) + '\n')
    

def cigarCalc(cigar):
    totalLength = 0
    currString = ''
    for char in cigar:
        if char == 'M':
            totalLength += int(currString)
            currString = ''
        elif char == 'N':
            totalLength += int(currString)
            currString = ''
        elif char == 'I':
            totalLength += int(currString)
            currString = ''
        elif char == 'D':
            currString = ''
        elif char == 'S':
            totalLength += int(currString)
            currString = ''
        elif char == 'H':
            currString = ''
        elif char == 'P':
            currString = ''
        elif char == 'X':
            totalLength += int(currString)
            currString = ''
        elif char == '=':
            totalLength += int(currString)
            currString = ''
        else:
            currString += char
    return totalLength

samPositions(sys.argv[1])