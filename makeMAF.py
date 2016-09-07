####################################################################
import glob
import sys
import statistics
import os
####################################################################

#globals
folderFlag=False


'''
this program takes input of MAF or folder of MAFs
it will save all the lines from the MAF file that contains a patient that we have in our
all.tsv file(hardcoded, sorry!)
This newly created MAF file will output into the folder that contains the MAFs we are 
working with
'''
def main():
    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return

    #ourPatient file (hardcoded, sorry)
    ourPatientFile="/home/valsh/Work/clymer/patient_files/all.tsv"
    ourPatients=readOurs(ourPatientFile)

    if folderFlag==True:
        inputFolder=args
        globPattern=inputFolder+'/'+'*.maf'
        inputFiles=glob.glob(globPattern)
    else:
        inputFiles=[args]
    finalOutput=""
    outputFile=os.path.dirname(inputFiles[0])+'/ourPatients.maf'
    for inputFile in inputFiles:
        print(inputFile + " ...")
        finalOutput+=readMAF(inputFile, ourPatients)
    with open(outputFile, 'w') as outFile:
        print("Writing output to: " + outputFile)
        outFile.write(finalOutput)
    return


def readMAF(inputFile, ourPatients):
    #open the file in correct codec so we dont get any errors from DOS files
    with open(inputFile, encoding='iso-8859-15') as inFile:
        output=""
        #loop for each line
        for line in inFile:
            #if line starts with # or "Hugo_" skip it (Hugo_ indicates the header)
            if line[0] != "#" and line[0:5] != "Hugo_":
                temp=line.split('\t')
                patient=temp[15].split('-')
                if patient[2] in ourPatients:
                    output+=('\n'+line.strip())
    return output


'''
Read in all our patients into a set and return it
'''
def readOurs(ourFile):
    ourSet=set()
    with open(ourFile) as inFile:
        for line in inFile:
            temp = line.rstrip().split('\t')
            check=temp[1]
            check=check.split('-')
            if check[3][:2] < "09":
                ourSet.add(check[2])
    return ourSet

'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
    global folderFlag
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) != 3 and len(args)!= 2:
        print("Usage: file.maf")
        print("Alt Usage: -f ./FOLDER/")
        return ()

    #otherwise set each var to the respected argument
    if len(args) == 2:
        arg1=args[1]
    else:
        folderFlag=True
        arg1 = args[2]

    #return all new vars
    return arg1
####################################################################
#call main function to start with
if __name__ == "__main__":
    main()