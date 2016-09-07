####################################################################
import sys
import os
import glob
####################################################################

'''
This script will process an input list of patients, the cancer mapping file, 
and the MAFs folder.
It will output a table of the patients and their cancer type along with
all the MAFs the patient can be found in.
This output file will be in the MAFs folder it was sent.
'''
def main():

	#check for arguments by calling parseArgs
	args = parseArgs(sys.argv)
	#if there are no arguments then close program
	if args == ():
		return

	patientsFile, cancerMappingFile, inputFolder = args

	dataTable, patientSet, cancerMappingTable =getCancer(patientsFile,cancerMappingFile)

	globPattern=inputFolder + '**/**/**/**/*.maf'
	inputFiles=glob.glob(globPattern)

	for inputFile in inputFiles:
		print(inputFile + " ...")
		readMAF(inputFile, patientSet, dataTable, cancerMappingTable)

	writeTable(inputFolder,dataTable)

	return

'''
This function will take the input patient file and the cancer mapping file.
It will read in the cancer mapping file and set it up as a dictionary with the
	keys being the 4 digit patient ID.
It will then read in the patients file and get the cancer type from the table
	that was just created.
We will return: dict[cancertype]->list[list[patient,list[]], ...]
	along with a set of all of our patients and the cancer mapping table
'''
def getCancer(patientsFile, cancerMappingFile):

	with open(cancerMappingFile) as inFile:
		cancerMappingTable={}
		for line in inFile:
			temp=line.strip('\n').split('\t')
			patient=temp[0][8:]
			cancer=temp[1]
			cancerMappingTable[patient]=cancer
	cancerPatients=set(cancerMappingTable.keys())

	with open(patientsFile) as inFile:
		inputPatientsTable={}
		patientSet=set()
		patientList=inFile.readline().strip('\n').split('\t')
		for p in patientList:
			patient=p[8:]
			cancer=cancerMappingTable[patient]
			data=[patient,[]]
			patientSet.add(patient)
			if not cancer in inputPatientsTable.keys():
				inputPatientsTable[cancer]={}
			if not patient in inputPatientsTable[cancer].keys():
				inputPatientsTable[cancer][patient]=[]
	# count=0
	# for i in inputPatientsTable.keys():
	# 	count+=len(inputPatientsTable[i])
	# print(count)

	return	inputPatientsTable, patientSet, cancerMappingTable

'''
This function takes input of a MAF file and the patient set
It scans the MAF looking for our patients, when it finds the
patient it will add the MAF file name into the dataTable
'''
def readMAF(inputFile, patientSet, dataTable, cancerMappingTable):
	#open the file in correct codec so we dont get any errors from DOS files
	with open(inputFile, encoding='iso-8859-15') as inFile:
		file=os.path.basename(inputFile)
		#loop for each line
		for line in inFile:
			#if line starts with # or "Hugo_" skip it (Hugo_ indicates the header)
			if line[0] != "#" and line[0:5] != "Hugo_":
				temp=line.split('\t')
				#determine if this line contains our patient, if so then do work
				patient=temp[15][8:12]
				if patient in patientSet:
					cancer=cancerMappingTable[patient]
					if not file in dataTable[cancer][patient]:
						dataTable[cancer][patient].append(os.path.basename(inputFile))

	return 


'''
This function will write out our data table to a .tsv file in our input folder
We are outputting like so:

cancerType	patientID	MAFfile1,	...,	etc
...
etc

'''
def writeTable(folder, dataTable):
	outputFile=folder+'patientMAFs.tsv'
	with open(outputFile, 'w') as outFile:
		for cancerType in dataTable.keys():
			for patient in dataTable[cancerType]:
				outputLine=cancerType + '\t' + patient + '\t'
				for file in dataTable[cancerType][patient]:
					outputLine+=file + '\t'

				outFile.write(outputLine + '\n')
	return


'''
Check arguments passed by user and returns arguments.
'''
def parseArgs(args):
	#check for commandline arguments, change compare value to be 1 more than
	#amount you are expecting, if it is less than expected then
	#print out the expected input
	if len(args) != 4:
		print("Usage: PatientsFile CancerMappingFile /MAFfolder/")
		return ()

	#otherwise set each var to the respected argument
	arg1 = args[1]
	arg2 = args[2]
	arg3 = args[3]

	#return all new vars
	return (arg1, arg2, arg3)
####################################################################
#call main function to start with
if __name__ == "__main__":
	main()
