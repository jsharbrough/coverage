import sys
import pickle
def indexVCF(vcfFile):
    infile = open(vcfFile,'r')
    depthDict = {}
    for line in infile:
        if line[0] != '#':
            lineSplit = line.split('\t')
            info = lineSplit[7]
            infoSplit = info.split(';')
            dp = infoSplit[0]
            depth = int(dp[3:])
            currScaffold = lineSplit[0]
            pos = int(lineSplit[1])
            if currScaffold in depthDict:
                scaffDict = depthDict[currScaffold]
                scaffDict[pos] = depth
            else:
                scaffDict = {pos:depth}
            depthDict[currScaffold] = scaffDict             
    with open(vcfFile + '.pickle','wb') as f:
        pickle.dump(depthDict,f)

def coverage(vcfFile,regionsFile=False):
    try:
        indexFile = open(vcfFile + '.pickle', 'rb')
    except IOError:
        indexVCF(vcfFile)
        indexFile = open(vcfFile + '.pickle', 'rb')
    depthDict = pickle.load(indexFile)
    indexFile.close()
    try:
        infile = open(regionsFile,'r')
        sys.stdout.write('Gene\tScaffold\tStart\tStop\tMean Depth\tMedian Depth\n')
        for line in infile:
            newLine = line
            while newLine[-1] == '\n' or newLine[-1] == '\t' or newLine[-1] == '\r':
                newLine = newLine[0:-1]
            lineSplit = newLine.split('\t')
            gene = lineSplit[0]
            scaffold = lineSplit[1]
            start = int(lineSplit[2])
            stop = int(lineSplit[3])
            depths = []
            i = start
            scaffDepths = depthDict[scaffold]
            totalDepth = 0
            while i < stop:
                posDepth = scaffDepths[i]
                totalDepth += posDepth
                depths.append(posDepth)
                i += 1
            depths.sort()
            meanDepth = totalDepth/(float(stop - start))
            if len(depths)%2 == 0:
                median = float(depths[len(depths)/2] + depths[len(depths)/2 - 1])/2.0
            else:
                median = depths[len(depths)/2]
            sys.stdout.write(gene + '\t' + scaffold + '\t' + str(start) + '\t' + str(stop) + '\t' + str(meanDepth) + '\t' + str(median) + '\n')
    except IOError:
        sys.stdout.write('Mean Depth\tMedian Depth\n')
        depths = []
        totalDepth = 0
        for scaffold in depthDict:
            scaffDict = depthDict[scaffold]
            for pos in scaffDict:
                depth = scaffDict[pos]
                depths.append(depth)
                totalDepth += depth
        depths.sort()
        meanDepth = totalDepth/(float(len(depths)))
        if len(depths)%2 == 0:
            median = float(depths[len(depths)/2] + depths[len(depths)/2 - 1])/2.0
        else:
            median = depths[len(depths)/2]
        sys.stdout.write(str(meanDepth) + '\t' + median + '\n')

def help():
    helpStatement = 'USAGE\n\n\tpython coverage.py foo.vcf regions.txt > coverage.txt\n'
    sys.stdout.write(helpStatement)

if len(sys.argv) == 3:
    coverage(sys.argv[1],sys.argv[2])
elif len(sys.argv) == 2:
    coverage(sys.argv[1])
else:
    help()