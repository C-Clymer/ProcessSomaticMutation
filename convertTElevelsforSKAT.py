####################################################################
import sys
####################################################################
'''
this program will format the data we already have into an input 
that will be accepted into the SKAT library for R
'''

def main():
    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return

    inputFile, outputFile = args
    table = readInTable(inputFile)
    printOut(table, outputFile)
    return

def readInTable(inputFile):
    with open(inputFile) as inFile:
        patients=inFile.readline().strip().split('\t')
        patients=patients[1:]
        for i in range(len(patients)):
            patients[i]=patients[i][8:]
        values=inFile.readline().strip().split('\t')
        values=values[1:]
        table=zip(patients,values)
    return table

def printOut(table, outputFile):
    with open(outputFile, 'w') as outFile:
        for i in table:
            outFile.write(i[0] + '\t' + i[1] + '\n')

'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) != 3:
        print("Usage: input output ")
        return ()

    #otherwise set each var to the respected argument
    arg1 = args[1]
    arg2 = args[2]

    #return all new vars
    return (arg1, arg2)
####################################################################
#call main function to start with
if __name__ == "__main__":
    main()
