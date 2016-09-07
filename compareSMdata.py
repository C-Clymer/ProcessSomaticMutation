####################################################################
import sys
import statistics
import os
####################################################################
'''
This function will take two input files and compare them
for ourFile, it should be created by processSomaticMutation.py script
for theirFile, it should be taken from CancerBrowser and be their genomicmatrix
We need to compare every patient/gene+value from each file
We want to count how many were ==, !=, did not exist on theirs, did not exist on ours
'''


def main():
    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return
    (ourFile, theirFile) = args

    compareFiles(ourFile,theirFile)

    return
def compareFiles(ourFile,theirFile):

    notOnOurs=0
    theirTable={}
    ourTable={}

    print("Reading their file...")
    #create table based on their file
    readInFile(theirFile, theirTable)    
    print("Reading our file...")
    #create table based on our file
    readInFile(ourFile, ourTable)    

    #compare tables
    print("Comparing mutations...")
    (equal, notEqual) = compareTables(theirTable,ourTable, notOnOurs)

    #create output string
    output=""
    output+="Mutation match count: " + str(equal) +'\n'
    output+="Mutation nonmatch count: " + str(notEqual) + '\n'
    percentOff=(notEqual/(equal+notEqual))*100
    output+=str(percentOff) + "% off"
    outputFile=os.path.dirname(ourFile)+'/NEWComparison.txt'
    print(output)
    with open(outputFile, 'w') as outFile:
        outFile.write(output)
    return

def findMAFline(tGene, patientName, inputFile):
    #open the file in correct codec so we dont get any errors from DOS files
    with open(inputFile, encoding='iso-8859-15') as inFile:
        for line in inFile:
            temp = line.rstrip().split('\t')
            if temp[0] == tGene and temp[15][8:15]==patientName:

                return line
    print("This shouldnt be happening")
    return


'''
Read in the inputFile and create a 2d array using dicts
Top array will be patient key, 2nd array will be gene key
each gene key points to a value which is taken from the input file

Arg1(inputFile) = file to be read in
Arg2(mainDict) = top level array

Inputfile should be formatted like so:
patients    ID1 ID2 ID3 ... etc
gene1   x   x   x   ...     x
gene2   x   x   x   ...     x
...
etc   x   x   x   ...     x
'''
def readInFile(inputFileName, mainDict):
    #create table based on their file
    with open(inputFileName) as file:
        header = [] #Used as a parallel key array
        for line in file:
            lSplit = line.strip().split("\t")
            gene = lSplit[0]
            if gene == "genes" or gene == "sample":
                header = lSplit #set patient key array
                for i in range(1, len(header)):
                    #Create a subdictionary for every key
                    #Gets only the last 7 characters of patient
                    header[i] = header[i][len(header[i])-7:]
                    if not header[i] in mainDict:
                        mainDict[header[i]] = {}
            else:
                for i in range(1, len(header)):
                    patient = header[i]
                    num = lSplit[i].split(".")[0]
                    mainDict[patient][gene] = int(num)
    return

'''
This function will take two tables, built from the input files, and compare each mutation to eachother
We are working with a 2d array, array[patient][gene]=value
We need to compare all matching patients for genes, then compare all matching gene values
add to correct variable as needed
We will also need to compare for extra genes per patient, if a gene is not present in our patient but is in theirs
we will scan through looking for all values. if 0 then equal++, if >1 notEqual++
We will also need to compare for all patients they have and we dont, same steps as above
'''
def compareTables(theirTable,ourTable, notOnOurs):
    sharedPatients,theirPatients,ourPatients = compareLists(theirTable, ourTable)
    equal=0
    notEqual=0
    sharedGenes,theirGenes,ourGenes=compareLists(theirTable[sharedPatients[0]],ourTable[sharedPatients[0]])

    for sPatient in sharedPatients:
        for sGene in sharedGenes:
            if (theirTable[sPatient][sGene]>=1) == (ourTable[sPatient][sGene]>=1):
                equal+=1
            else:
                notEqual+=1
        if len(theirGenes)>0:
            for tGene in theirGenes:
                if theirTable[sPatient][tGene]>0:
                    notEqual += 1
                else:
                    equal +=1

        if len(ourGenes)>0:
            for oGene in ourGenes:
                if theirTable[sPatient][oGene]>0:
                    notEqual += 1
                else:
                    equal +=1
    if len(theirPatients)>0:
        for tPatient in theirPatients:
            theirGenes=list(theirTable[tPatient])
            for tGene in theirGenes:
                if theirTable[tPatient][tGene]>0:
                    notEqual += 1
                else:
                    equal +=1
    if len(ourPatients )>0:
        for oPatient in ourPatients:
            ourGenes=list(ourTable[oPatient])
            for oGene in ourGenes:
                if ourTable[oPatient][tGene]>0:
                    notEqual += 1
                else:
                    equal +=1


    return (equal, notEqual)




def compareLists(theirDict, ourDict):
    theirList=list(theirDict)
    ourList=list(ourDict)
    theirSharedIndex=[]
    ourSharedIndex=[]
    sharedList=[]
    popCount=0

    for i in range(len(theirList)):

        for j in range(len(ourList)):
            
            if theirList[i]==ourList[j]:
                sharedList.append(theirList[i])
                theirList[i] = False
                ourList[j] = False
                break

    ls = len(sharedList)
    lo = len(ourList)
    lt = len(theirList)
    ourListShortened = [False]*(lo - ls)
    theirListShortened = [False]*(lt - ls)

    j = 0
    for i in range(lo):
        if ourList[i] != False:
            ourListShortened[j] = ourList[i]
            j = j + 1

    j = 0
    for i in range(lt):
        if theirList[i] != False:
            theirListShortened[j] = theirList[i]
            j = j + 1
    return sharedList, theirListShortened, ourListShortened

    '''
    for i in theirList:
        if i in ourList:
            sharedList.append(i)
            ourList.pop(ourList.index(i))
            theirList.pop(theirList.index(i))
    return sharedList, theirList,ourList
    '''
    '''
    for i in range(len(theirList)):
        x=i-popCount
        for j in range(len(ourList)):
            if theirList[x]==ourList[j]:
                sharedList.append(theirList[x])
                theirList.pop((x))
                ourList.pop((j))
                popCount+=1 
                break 
    return sharedList, theirList, ourList 
    '''

    '''
    for i in range(len(theirList)):
        for j in range(len(ourList)):
            if theirList[i]==ourList[j]:
                sharedList.append(theirList[i])
                ourSharedIndex.append(j)
                theirSharedIndex.append(i)         
    theirSharedIndex.sort()
    theirSharedIndex=theirSharedIndex[::-1]
    ourSharedIndex.sort()
    ourSharedIndex=ourSharedIndex[::-1]
    for i in theirSharedIndex:
        theirList.pop(i)
    for i in ourSharedIndex:
        ourList.pop(i)

    return sharedList,theirList,ourList
    '''

'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) != 3:
        print("Usage: OurFile.txt TheirFile")
        return ()

    #otherwise set each var to the respected argument
    ourFile = args[1]
    theirFile = args[2]

    #return all new vars
    return (ourFile, theirFile)
####################################################################
#call main function to start with
if __name__ == "__main__":
    main()
