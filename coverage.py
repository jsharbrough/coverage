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
            if dp[0:3] == 'DP=':
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

def coverage(vcfFile,regionsFile=False, maxDepth=False):
    try:
        indexFile = open(vcfFile + '.pickle', 'rb')
    except IOError:
        indexVCF(vcfFile)
        indexFile = open(vcfFile + '.pickle', 'rb')
    depthDict = pickle.load(indexFile)
    indexFile.close()
    if regionsFile != False:
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
                    if maxDepth != False:
                        if posDepth < maxDepth:
                            totalDepth += posDepth
                            depths.append(posDepth)
                    else:
                        totalDepth += posDepth
                        depths.append(posDepth)
                    i += 1
                depths.sort()
                meanDepth = totalDepth/(float(len(depths)))
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
                    posDepth = scaffDict[pos]
                    if maxDepth != False:
                        if posDepth < maxDepth:
                            totalDepth += posDepth
                            depths.append(posDepth)
                    else:
                        totalDepth += posDepth
                        depths.append(posDepth)
            depths.sort()
            meanDepth = totalDepth/(float(len(depths)))
            if len(depths)%2 == 0:
                median = float(depths[len(depths)/2] + depths[len(depths)/2 - 1])/2.0
            else:
                median = depths[len(depths)/2]
            sys.stdout.write(str(meanDepth) + '\t' + str(median) + '\n')
    else:
        sys.stdout.write('Mean Depth\tMedian Depth\n')
        depths = []
        totalDepth = 0
        for scaffold in depthDict:
            scaffDict = depthDict[scaffold]
            for pos in scaffDict:
                posDepth = scaffDict[pos]
                if maxDepth != False:
                    if posDepth < maxDepth:
                        totalDepth += posDepth
                        depths.append(posDepth)
                else:
                    totalDepth += posDepth
                    depths.append(posDepth)
        depths.sort()
        meanDepth = totalDepth/(float(len(depths)))
        if len(depths)%2 == 0:
            median = float(depths[len(depths)/2] + depths[len(depths)/2 - 1])/2.0
        else:
            median = depths[len(depths)/2]
        sys.stdout.write(str(meanDepth) + '\t' + median + '\n')


def help():
    helpStatement = '\nUSAGE\n\n\tpython coverage.py -i <vcfFile> <options> > coverage.txt\n\nOPTIONS\n\t-i\t<vcf file>\tInput File [REQUIRED]. If no input file is given, the program will output the \n\t\t\t\thelp menu. File should be in VCF format for a single individual per file.\n\t-r\t<regions file>\tRegions file. Tab-delimited file containing the region name, the scaffold \n\t\t\t\tname, the start, and stop position for each region. One region per line. \n\t\t\t\tsamPositions.py can be used to generate this file or it can be manually \n\t\t\t\tgenerated.\n\t-d\t<int>\t\tMax Depth. Sites with coverage greater than d will be ignored.\n\n'
    sys.stderr.write(helpStatement)

sys.stdout.write(str(sys.argv))
optionDict = {}
i = 1
while i < (len(sys.argv)-1):
    if sys.argv[i] == '-i':
        optionDict['vcfFile'] = sys.argv[i+1]
    elif sys.argv[i] == '-r':
        optionDict['regionsFile'] = sys.argv[i+1]
    elif sys.argv[i] == '-d':
        optionDict['maxDepth'] = int(sys.argv[i+1])
    i += 1

if len(optionDict) > 0 and 'vcfFile' in optionDict:
    if 'regionsFile' in optionDict:
        if 'maxDepth' in optionDict:
            coverage(optionDict['vcfFile'],optionDict['regionsFile'],optionDict['maxDepth'])
        else:
            coverage(optionDict['vcfFile'],optionDict['regionsFile'])
    elif 'maxDepth' in optionDict:
        coverage(optionDict['vcfFile'],False,optionDict['maxDepth'])
    else:
        coverage(optionDict['vcfFile'])
else:
    help()