####################################################################
import glob
import sys
import os
####################################################################
'''
Goal for this program:
Take commandline arguments for cancer type folder taken from TCGA website
and process the file_manifest.txt, FILE_SAMPLE_MAP.txt, and .maf files

Program will take input of .maf files and create an output for the total count
of somatic mutations per patient per gene, also the total count of mutations
minus silent mutations

MAF files organize data with tab spaceing between columns
These are the columns we are interested in:
*(1)hugo symbol - gene
*(9)Variant_Classification - translational effect of variant allele
        Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
        Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site,
        Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank,
        IGR, Intron, RNA, Targeted_Region
*(16)tumor_sample_barcode - BCR aliquot barcode for the tumor sample
'''

def main():
    #commandline argument will modify this to deal with folders of maf files
    global folderFlag
    folderFlag=False
    #commandline argument to determine if we are just comparing two preivously created tables    
    comparisonFile=False

    #check for arguments by calling parseArgs
    args = parseArgs(sys.argv)
    #if there are no arguments then close program
    if args == ():
        return

    #determine what flags need to be dealt with    
    if len(args)==2:
        if args[0]==True:
            (folderFlag, inputFolder)=args
        else:
            (folderFlag,inputFile)=args
    elif len(args)==3:
        if args[0]==False:
            (folderFlag, inputFile, comparisonFile) = args
        else:
            (folderFlag, inputFolder, comparisonFile)=args

    #if we are working with a folder, grab all the maf file paths and place into list for loop
    if folderFlag==True:
        globPattern=inputFolder+'/'+'*.maf'
        inputFiles=glob.glob(globPattern)
    #set output files names into the right directory    
    if folderFlag == False:
        #we want to output file into subfolder (ProcessedData) where inputfile is located
        #check if folder already exists, if it doesnt then create it
        if not os.path.exists(os.path.dirname(inputFile)+"/ProcessedData/"):
            os.makedirs(os.path.dirname(inputFile)+"/ProcessedData/")
        #set output file (same location with .SM_DATA.txt appended to the end)
        if not os.path.exists(os.path.dirname(inputFile)+"/ProcessedData/"+os.path.basename(inputFile)[:-4]+'/'):
            os.makedirs(os.path.dirname(inputFile)+"/ProcessedData/"+os.path.basename(inputFile)[:-4]+'/')
        outputFile=os.path.dirname(inputFile)+ "/ProcessedData/"+os.path.basename(inputFile)[:-4]+'/SM_DATA.txt'
        outputSingleFile=os.path.dirname(inputFile)+ "/ProcessedData/"+os.path.basename(inputFile)[:-4]+'/SINGLE_SM_DATA.txt'
        fullPatientOut=os.path.dirname(inputFile)+ "/ProcessedData/"+os.path.basename(inputFile)[:-4]+'/PATIENT_LIST.txt'
    #set output file names into a main folder when given mutliple maf files
    if folderFlag==True:
        if not os.path.exists(inputFolder+"ProcessedData/"):
            os.makedirs(os.path.dirname(inputFolder+"ProcessedData/Folder"))
        outputFile=inputFolder+'ProcessedData/Folder/SM_DATA.txt'
        outputSingleFile=inputFolder+'ProcessedData/Folder/SINGLE_SM_DATA.txt'
        fullPatientOut=inputFolder+'/ProcessedData/Folder/PATIENT_LIST.txt'

    #gene list for output and patientDict to hold all the data
    geneList=[]
    patientDict={}

    #run on single file
    if folderFlag==False:
        #process the file, call readMAF function
        print('\n' + "Reading: " + os.path.basename(inputFile))
        fullPatientSet = readMAF(inputFile,patientDict,geneList)

    #run on multiple files
    if folderFlag==True:
        fullPatientSet=set()
        for inputFile in inputFiles:
            print('\n' + "Reading: " + os.path.basename(inputFile))
            fullPatientSetNew = readMAF(inputFile,patientDict,geneList)
            fullPatientSet.update(fullPatientSetNew)

    #output the information to output file using writeTable function
    print("Writing table to : " + os.path.basename(outputFile))
    writeTable(outputFile, patientDict, geneList)
    print("Writing single table to : " + os.path.basename(outputSingleFile))
    writeSingleTable(outputSingleFile, patientDict, geneList)
    #output finished and the total count of patients/genes to console
    print("Finished: Total count of patients = " + str(len(patientDict)) + ' Genes = ' + str(len(geneList)))

    #Uncomment if we want a list of all the patients
    ##print("Writing full patient list to: " + os.path.basename(fullPatientOut))
    ##with open(fullPatientOut, 'w') as outFile:
    ##    for i in fullPatientSet:
    ##        outFile.write(i + '\t')

    #compare files if 2 arguments in commandline
    if comparisonFile!=False:
        compareFiles(outputSingleFile, comparisonFile, inputFile)

    return

'''
This function will read in the file, line by line, for every new patient sample we get
we will crate a new dictionary entry for it (e.g. 01A vs 06A)
for each new gene encountered in each sample we will make a new dict entry into
patientDict[patient][sample]. Each patientDict will hold dictionary of samples with
a dictionary of genes, each gene dictionary will hold a list of 2 values,
count of total somatic mutations and count of silent type somatic mutations
'''
def readMAF(inputFile, patientDict,geneList):
    fullPatientSet=set()
    #open the file in correct codec so we dont get any errors from DOS files
    with open(inputFile, encoding='iso-8859-15') as inFile:
        #loop for each line
        for line in inFile:
            #if line starts with # or "Hugo_" skip it (Hugo_ indicates the header)
            if line[0] != "#" and line[0:5] != "Hugo_":
                #split the line up based on tabs and store in list
                temp = line.rstrip().split('\t')
                #split the 15th position in line list into its own list
                #this is the barcode for the tumor sample
                fullPatientSet.add(temp[15])
                patient=temp[15].rsplit('-')
                #we need to grab the patient ID and the tumor sample used
                patientSample=patient[3]
                patient=patient[2]
                #if the gene(temp[0]) is not in the geneList, add it
                if temp[0] not in geneList:
                    geneList.append(temp[0])
                #if the patient is not in patientDict, add them to it
                #and create a dict for samples
                if patient not in patientDict:
                    patientDict[patient]={}
                #if the sample is not yet in patientDict[patient] then add it
                #and make a new dict for gens
                if patientSample not in patientDict[patient]:
                    patientDict[patient][patientSample]={}
                #if the current gene is in the gene dict then add 1 to the total count
                if temp[0] in patientDict[patient][patientSample]:
                    patientDict[patient][patientSample][temp[0]][0]+=1
                #if it is not, add the gene and set its list to [1,0]
                #(1 total mutation since we have just come across it)
                else:
                    patientDict[patient][patientSample][temp[0]]=[1,0]
                #determine if this mutation is silent, if it is then increment the silent count
                if temp[8]=='Silent':
                    patientDict[patient][patientSample][temp[0]][1]+=1
        #send the patientDict to comparePDicts function so we can
        #determine which sample to keep
        comparePDicts(patientDict)
    return fullPatientSet

'''
This function takes the patientDict dictonary
and scans through it finding patients with more than 1 sample
when it comes across a patient it will take the sum of the total mutations
inside each patient sample and determine which has a higher count
the higher count is what we will use to print out
'''
def comparePDicts(patientDict):
    print("Comparing multi-sample patients" + '\n' +'-'*8)
    #get a list of keys from patientDict and step through them
    pKeys=patientDict.keys()
    for pKey in pKeys:
        #get  alist of keys from current patient for their samples
        #if they have more than 1 sample, start processing
        sKeys=sorted(patientDict[pKey].keys())
        if len(sKeys) > 1:
            #create a countSums list so we can store the values of each samples
            #mutation count
            countSums = []
            #print out patient and the keys they hold BEFORE removal
            print(str(pKey))
            print('Before:' + '\t' + str(sKeys), end='\t')
            #loop through the sample based on gene keys and add up all the total
            #counts, append it to countSums so we can compare
            for sKey in sKeys:
                gKeys=patientDict[pKey][sKey].keys()
                sSum=0
                for gKey in gKeys:
                    sSum=sSum+patientDict[pKey][sKey][gKey][0]
                countSums.append(sSum)
            #print the countSums so we can see the values BEFORE removal
            print(str(countSums))
            #set i to -1, this will hold the index position into the list so we
            #know which sample to remove
            #the countSums and sKeys list are related so this will work
            i=-1
            #loop through each value and compare it to the other values
            #store the higher value INDEX into i
            for value1 in countSums:
                for value2 in countSums:
                    if value1 < value2:
                        i=countSums.index(value2)
            #loop through sKey list, if the INDEX of current sKey does not equal
            #the INDEX we found above, delete it from the patientDict
            for sKey in sKeys:
                if sKeys.index(sKey)!=i:
                    patientDict[pKey].pop(sKey,None)
            #print out what is leftover from AFTER the deletion
            print('After:' +'\t' + str(patientDict[pKey].keys()) +'\t' + str(countSums[i]))
            print('-'*8)
    return

'''
This function will output the table created from patientDict into the
outputFile specified at commandline

Should look like this:
genelist  patient1-tCount  patient1-dcount    patient2-tCount patient2-dCount    etc
gene    x   x   x   x   etc
gene    x   x   x   x   etc
'''
def writeTable(outputFile, patientDict, geneList):
    with open(outputFile, 'w') as outFile:
        #get list of patientDict keys so we can print out header and step through
        #the patientDict for final output
        pkeys=sorted(patientDict.keys())
        #print header
        header = "genes"
        for p in pkeys:
            sKeys=patientDict[p].keys()
            for sKey in sKeys:
                #give each patient 2 columns, 1 for total count and 1 for deleterious count
                header= header + '\t' + p +'-'+sKey[:2]+'_t'+'\t'+ p + '-' + sKey[:2] + '_dt'
        print(header,file=outFile)
        #create each gene row by stepping through the geneList
        for g in geneList:
            row=g
            for p in pkeys:
                sKeys=patientDict[p].keys()
                for s in sKeys:
                    #if the patient has no value for this gene, print out 0's
                    if g not in patientDict[p][s]:
                        dCount = 0
                        tCount = 0
                    else:
                        dCount = patientDict[p][s][g][0]-patientDict[p][s][g][1]
                        tCount = patientDict[p][s][g][0]
                    row=row+'\t'+str(tCount)+'\t'+str(dCount)
            outFile.write(row+'\n')
    return


'''
This function will output the table created from patientDict into the
outputFile specified at commandline

Should look like this:
genelist  patient1-dCount   x   x   x   etc
gene    x   x   x   x   etc
'''
def writeSingleTable(outputFile, patientDict, geneList):
    with open(outputFile, 'w') as outFile:
        #get list of patientDict keys so we can print out header and step through
        #the patientDict for final output
        pkeys=sorted(patientDict.keys())
        #print header
        header = "genes"
        for p in pkeys:
            sKeys=patientDict[p].keys()
            for sKey in sKeys:
                #give each patient 2 columns, 1 for total count and 1 for deleterious count
                header= header + '\t' + p +'-'+ sKey[:2]
        print(header,file=outFile)
        #create each gene row by stepping through the geneList
        for g in geneList:
            row=g
            for p in pkeys:
                sKeys=patientDict[p].keys()
                for s in sKeys:
                    #if the patient has no value for this gene, print out 0's
                    if g not in patientDict[p][s]:
                        dCount = 0
                        tCount = 0
                    else:
                        dCount = patientDict[p][s][g][0]-patientDict[p][s][g][1]
                    row=row+'\t'+str(dCount)
            outFile.write(row+'\n')
    return



'''
Check commandline arguments passed by user and returns arguments.
'''
def parseArgs(args):
    #check for commandline arguments, change compare value to be 1 more than
    #amount you are expecting, if it is less than expected then
    #print out the expected input
    if len(args) == 1 or len(args)>4:
        print("Usage: inputFile.maf")
        print("Alt usage: inputFile.maf comparisonFile")
        print("Alt usage: -f inputFolder")
        print("Alt usage: -f inputFolder comparisonFile")
        return ()
    #validate file extension
    if args[1].endswith('.maf') != 1:
        if args[1]=='-f':
            if len(args)==4:
                inputFolder=args[2]
                comparisonFile=args[3]
                folderFlag=True
                return (folderFlag,inputFolder,comparisonFile)
            if len(args)==3:
                inputFolder=args[2]
                folderFlag=True
                return (folderFlag, inputFolder)
        print("Please specify a valid .maf file")
        return ()

    if len(args)==3:
        inputFile=args[1]
        comparisonFile=args[2]
        folderFlag = False
        return (folderFlag, inputFile, comparisonFile)
    else:
        #otherwise set each var to the respected argument
        inputFile = args[1]
        folderFlag=False
        return (folderFlag, inputFile)
    #if it hits this return, something went wrong
    return ()



'''
This function will take two input files and compare them
for ourFile, it should be created by processSomaticMutation.py script
for theirFile, it should be taken from CancerBrowser and be their genomicmatrix
We need to compare every patient/gene+value from each file
We want to count how many were ==, !=, did not exist on theirs, did not exist on ours
'''
def compareFiles(ourFile,theirFile,inputFile):
    equal=0				#holds count of matches
    notEqual=0			#holds count of mismatches
    theirPatients=[]	#list of their patients, used during read in
    theirGenes=set()	#set of their genes
    theirTable={}		#their Dict with keys based on gene
    ourPatients=[]		#list of our patients, used during read in
    ourGenes=set()		#set of our genes
    ourTable={}			#our Dict with keys based on gene

    print("Reading their file...")
    #create table based on their file
    with open(theirFile) as theirIn:        
        temp=theirIn.readline().rstrip().split('\t')[1:]
        for place in temp:
            theirPatients.append(place[8:])
        for line in theirIn:
            tempLine=line.rstrip().split('\t')
            gene=tempLine[0]
            theirGenes.add(gene)
            theirTable[gene]=list(zip(theirPatients, tempLine[1:]))

    print("Reading our file...")
    #create table based on our file
    with open(ourFile) as ourIn:
        temp=ourIn.readline().rstrip().split('\t')[1:]
        ourPatients=temp
        for line in ourIn:
            tempLine=line.rstrip().split('\t')
            gene = tempLine[0]
            ourGenes.add(gene)
            ourTable[gene]=list(zip(ourPatients, tempLine[1:]))

    #convert patients to sets for faster processing from here on
    theirPatients=set(theirPatients)
    ourPatients=set(ourPatients)

    #get sets of all the shared patients and genes
    sharedPatients=set.intersection(theirPatients,ourPatients)
    sharedGenes = set.intersection(theirGenes, ourGenes)

  	#start output for MAFLine output if we arnt working with an entire folder
    if folderFlag is False:
        with open(inputFile) as inFile:
            notEqualOutput=inFile.readline()

    #step through dict based on gene, comparing patients.
    #if they have 1, we need 1+, if they have 0 we need 0
    print("Comparing mutations...")
    for tGene in sharedGenes:
        for x in theirTable[tGene]:
            for y in ourTable[tGene]:
                if y[0]==x[0]:
                    if x[1]=='1' or x[1]==1 or x[1]=='1.0' or x[1]==1.0:
                        if y[1]>='1' or y[1]==1:
                            equal+=1
                        else:
                            notEqual+=1                                
                    else:
                        if y[1]=='0' or y[1]==0:
                            equal+=1
                        else:
                            notEqual+=1
                            if folderFlag is False:
                                notEqualOutput+= findMAFline(tGene, y[0], inputFile)

    #Check all genes they have and we dont, if 1 then it is a mismatch
    onlyTheirGenes=set.difference(theirGenes,ourGenes)                              
    for tGene in onlyTheirGenes:
        for x in theirTable[tGene]:
            if x[1]=='1' or x[1]==1 or x[1]=='1.0' or x[1]==1.0:
                notEqual+=1
            if x[1]=='0' or x[1]==0 or x[1]=='0.0' or x[1]==0:
                equal+=1

    #Check all genes we have and they dont, if 1 then it is a mismatch
    onlyOurGenes=set.difference(ourGenes,theirGenes)                             
    for oGene in onlyOurGenes:
        for x in ourTable[oGene]:
            if x[1]=='1' or x[1]==1:
                notEqual+=1
            if x[1]=='0' or x[1]==0:
                equal+=1

    #Output final results, to console and to outFile
    output="Our patients: " + str(len(ourPatients))+ " Their Patients: " +str(len(theirPatients)) + " Equal count: "+str(len(sharedPatients)) +'\n'
    output+="Our genes: " + str(len(ourGenes))+ " Their genes: " +str(len(theirGenes)) + " Equal count: "+str(len(sharedGenes)) + '\n'
    output+="Mutation match count: " + str(equal) +'\n'
    output+="Mutation nonmatch count: " + str(notEqual) + '\n'
    percentOff=(notEqual/(equal+notEqual))*100
    output+=str(percentOff) + "% off"
    outputFile=os.path.dirname(ourFile)+'/Comparison.txt'
    print(output)
    with open(outputFile, 'w') as outFile:
        outFile.write(output)

	#This code will print out the maf line results of any mismatched genes where they have 1 and we have 0
	#will only execute if the comparison was not run on an entire folder and if there is data to display
    if folderFlag == False:
    	if len(notEqualOutput)>614:
	        outputFile=os.path.dirname(ourFile)+'/NotEqualMAFlines.txt'
	        with open(outputFile, 'w') as outFile:
	            outFile.write(notEqualOutput)
    return


'''
This function will use the given data to find the correct line in the MAF line to return as a string
It is not used when doing a scan of all .maf files in a folder
'''
def findMAFline(tGene, patientName, inputFile):
    #open the file in correct codec so we dont get any errors from DOS files
    with open(inputFile, encoding='iso-8859-15') as inFile:
        for line in inFile:
            temp = line.rstrip().split('\t')
            if temp[0] == tGene and temp[15][8:15]==patientName:
                return line
    print("==========================")
    print("This shouldnt be happening")
    print("==========================")
    return


####################################################################

#call main function to start with
if __name__ == "__main__":
    main()
