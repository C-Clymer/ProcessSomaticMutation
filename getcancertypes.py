####################################################################
import glob
import sys
import statistics
import os
####################################################################
'''
compares all.tsv and our patients not in pancan to get their cancer type, one time run..all done
'''

def main():
    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return
    (theirFile, ourFile)=args
    ourSet=set()
    outputFile="./Cancernotinpan.txt"
    output=""
    with open(ourFile) as inFile:
        for line in inFile:
            temp = line.rstrip().split('\t')
            check=temp[1]
            check=check.split('-')
            if check[3][:2] < "09":
                ourSet.add((check[2],temp[2]))
    with open(theirFile) as inFile:
        temp=inFile.readline()
        temp=temp.split('\t')
        temp.pop(len(temp)-1)
        notPAN=set(temp)
    for i in notPAN:
        for j in ourSet:
            if i == j[0]:
                output+=str(j[0])+'\t'+str(j[1])+'\n'
    with open(outputFile, 'w') as outFile:
        outFile.write(output)
    return


'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) != 3:
        print("Usage: theirFile ourFile.tsv")
        return ()

    #otherwise set each var to the respected argument
    theirFile = args[1]
    ourFile = args[2]

    #return all new vars
    return (theirFile, ourFile)
####################################################################
#call main function to start with
if __name__ == "__main__":
    main()
