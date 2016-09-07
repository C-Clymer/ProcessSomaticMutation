####################################################################
import glob
import sys
import os
####################################################################
'''
this script will compare all the patients in our files with the patients
in the PANCAN maf files, it will then print out all the matching patients
along with a total count of matching
'''

def main():
    allpatients="/home/valsh/Work/server/patient_files/all.tsv"
    pancanfolder="/home/valsh/Work/server/SomaticMutationData/PANCAN/"
    allset=set()
    pancanset=set()
    globPattern=pancanfolder+'*.maf'
    inputFiles=glob.glob(globPattern)

    with open(allpatients) as infile:
        for line in infile:
            temp = line.rstrip().split('\t')
            check=temp[1]
            check=check.split('-')
            if check[3][:2] < "09":
                allset.add(check[2])
    for i in inputFiles:
        with open(i) as infile:
            mafset=readMAF(i)
            pancanset=set.union(pancanset,mafset)
    print(len(pancanset))
    print(len(allset))
    print(len(set.intersection(pancanset, allset)))
    outputfile="/home/valsh/Work/server/patient_files/inPANCAN.txt"
    with open(outputfile, 'w') as outfile:
        for i in set.intersection(pancanset, allset):
            print( i, file=outfile)
    return


def readMAF(inputFile):
    fullPatientSet=set()
    #open the file in correct codec so we dont get any errors from DOS files
    with open(inputFile, encoding='iso-8859-15') as inFile:
        #loop for each line
        for line in inFile:
            #if line starts with # or "Hugo_" skip it (Hugo_ indicates the header)
            if line[0] != "#" and line[0:5] != "Hugo_" and len(line)>3:
                #split the line up based on tabs and store in list
                temp = line.rstrip().split('\t')
                #split the 15th position in line list into its own list
                #this is the barcode for the tumor sample
                check=temp[15]
                check=check.split('-')
                #fullPatientSet.add(temp[15])
                fullPatientSet.add(check[2])

    return fullPatientSet

####################################################################
#call main function to start with
if __name__ == "__main__":
    main()
