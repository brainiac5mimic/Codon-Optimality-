# creator: Mark Holton
# file name: codonExtractor.py
# 
# 
# 

import sys
#import matplotlib
#matplotlib.use('agg')
#import matplotlib.pyplot as plt
import copy
import pickle
from klepto.archives import dir_archive


inFile = sys.argv[1]
outFile = sys.argv[2]
outDict = sys.argv[3]
sequenceExcel = sys.argv[4]
# istance variables 

rnadict ={  }
#form of rnadict example
# "gene":{'geneName':{},'codons':{ }, 'scores':{ } }
#codons = listOutput in findCodon 
#scores = codonScores


codonScores=[]
otherScores={}
listOfScores = []
scoreList =[]

print('Start')

def testPrint():
    print('\n')
    print(str(rnadict.get('F45E6.3').get('scores')))
    print('\n')
    print(str(rnadict.get('F45E6.3').get('sum_of_5_on_each_side')))
    print(str(rnadict.get('F45E6.3').get('change_by')))
    print('\n')
    print(str(rnadict.get('F45E6.3').get('percentiles')))

# Functions
def compareCodons(lscores,allScores,gene):
    #lscores list of the scores for current gene 
    #allScores dictionry 
    otherScores[gene]={}
    findPercentiles(lscores, allScores, gene)
    sumof5(lscores, allScores, gene)
    changeby(lscores, allScores, gene)

def findPercentiles(lscores,allScores,gene):
    centile=[]
    raw_cent=[]
    i=1
    while i < len(lscores):
        temp=i/len(lscores)*100
         
        temp=round(temp,2)
        
        centile.append(str(temp)+'%' )
        raw_cent.append(temp)
        i=i+1    
    rnadict[gene].update({'percentiles':centile})
    rnadict[gene].update({'raw_percentiles':raw_cent})

def changeby(lscores,allScores,gene):
    system="change_by"
    change_by=[0]
    rnadict.get(gene)[system]={}
    i=1
    while i < len(lscores):
        temp=lscores[i]-lscores[i-1]
         
        temp=round(temp,3)
        
        change_by.append(temp)
        i=i+1
    rnadict[gene].update({system: change_by})



def sumof5(lscores,allScores,gene):
    system="sum_of_5_on_each_side"
    sum_of_5_on_each_side=[]
    rnadict.get(gene)[system]={}

    i=5
    while i+5 < len(lscores):
        temp=(lscores[i-5]+lscores[i-4]+lscores[i-3]+lscores[i-2]+lscores[i-1]+lscores[i]+lscores[i+1]+lscores[i+2]+lscores[i+3]+lscores[i+4]+lscores[i+5])/11
         
        temp=round(temp,3)
         
        sum_of_5_on_each_side.append(str(temp))
        i=i+1


    rnadict[gene].update({system: sum_of_5_on_each_side})

    
def getInfo(excelLineArray):
    #format of sequenceInfo.txt Gene_Name	Alias	Number_of_reads	Normalized_expression	Median	PolyA_class	25%_percentile	75%_percentile	MAX_polyA	MIN_polyA	RANGE_polyA	Mean	Geometric_mean	Standard_deviation	Fop	Codon_class
    #Gene_Name    Alias    Number_of_reads    Normalized_expression    Median    PolyA_class    25%_percentile    75%_percentile    MAX_polyA    MIN_polyA    
    #RANGE_polyA    Mean    Geometric_mean    Standard_deviation    Fop    Codon_class
    for line in excelLineArray:
        lineArray = line.split('\t')
        geneName= lineArray[0].strip()
        print(geneName)
        rnadict[geneName]={'geneName': geneName,'Alias':lineArray[1].strip(),    'Number_of_reads':lineArray[2].strip(),    'Normalized_expression':lineArray[3].strip(),    'Median':lineArray[4].strip(),    'PolyA_class':lineArray[5].strip(),    '25%_percentile':lineArray[6].strip(),    '75%_percentile':lineArray[7].strip(),    'MAX_polyA':lineArray[8].strip(),    'MIN_polyA':lineArray[9].strip(),   'RANGE_polyA':lineArray[10].strip(),    'Mean':lineArray[11].strip(),    'Geometric_mean':lineArray[12].strip(),    'Standard_deviation':lineArray[13].strip(),    'Fop':lineArray[14].strip(),    'Codon_class':lineArray[15].strip()    }
    print("\nend")
            


def emptyScores():
    tempcod= len(codonScores)
    x=0
    while x<tempcod:
        del codonScores[0]
        x=x+1
    

def findCodon(strInput,ind):
    listOutput=[]
    index = strInput.find("ATG");
    geneName = strInput[0:index].strip()
    while index+3 < len(strInput):
        codon = strInput[index:index+3]
        listOutput.append(codon);
        codonScores.append(codonScore(codon))
        index = index+3;
    
    score=codonScores
    print(geneName)
    print(len(score))
    if(geneName in rnadict):
      print(geneName)
      rnadict[geneName].update({ 'codons': listOutput, 'scores': copy.deepcopy(score) })
      compareCodons(score, listOutput, geneName)
    elif(geneName[0:len(geneName)-1] in rnadict and 'scores' not in rnadict[geneName[0:len(geneName)-1]]):
      #if(geneName[len(geneName)-1:len(geneName)]=='a'):
      print(geneName +'  trunk  '+geneName[0:len(geneName)-1])
      geneName=geneName[0:len(geneName)-1]
      rnadict[geneName].update({ 'codons': listOutput, 'scores': copy.deepcopy(score) })
      compareCodons(score, listOutput, geneName)
    
    #else:
      #rnadict[geneName]={ 'geneName': geneName, 'codons': listOutput, 'scores': copy.deepcopy(score) }
    return stringCodons(listOutput)

def stringScores(scoreArray):
    strIndex = 0;
    strReturn = "Scores are: "
    while strIndex < len(scoreArray):
        strReturn = strReturn +  str(scoreArray[strIndex]) +" \t"
        strIndex = strIndex +1;
    emptyScores()
    return strReturn 

def findCodonFile( arrayInput):
    returnList=[]
    ind = 1;
    for element in arrayInput:
        returnList.append(findCodon(element,ind)+'\n' +stringScores(codonScores)+'\n')
        ind=ind+1;
    return returnList


def stringCodons(listOutput):
    strIndex = 0;
    strReturn = "Codons are: "
    while strIndex < len(listOutput):
        strReturn = strReturn +  str(listOutput[strIndex]) +" \t"
        strIndex = strIndex +1;
    return strReturn


#def splot(yPoints,ind,name):
#   y = yPoints   # a list of numbers
#    #plt.subplot(1,1,1)
#    plt.plot(y)                     # draw the graph
#    plt.savefig(name.strip()+'.png',bbox_inches='tight')
#    plt.clf()



def codonScore(codon):
    score = 0.0;
    if codon == "TTT":
        score = 0.98
    if codon == "TTC":
        score = 1.02
    if codon == "TTA":
        score = 0.75
    if codon == "TTG":
        score = 1.37
    if codon == "CTT":
        score = 1.73
    if codon == "CTC":
        score = 0.91
    if codon == "CTA":
        score = 0.54
    if codon == "CTG":
        score = 0.69
    if codon == "ATT":
        score = 1.63
    if codon == "ATC":
        score = 0.88
    if codon == "ATA":
        score = 0.49
    if codon == "ATG":
        score =  1.00 #--
    if codon == "GTT":
        score = 1.69
    if codon == "GTC":
        score = 0.77
    if codon == "GTA":
        score = 0.72
    if codon == "GTG":
        score = 0.82
    if codon == "TAT":
        score = .97
    if codon == "TAC":
        score = 1.03
    if codon == "TAA":
        score = 1.87
    if codon == "TAG":
        score = .54
    if codon == "CAT":
        score = 1.13
    if codon == "CAC":
        score = .87
    if codon == "CAA":
        score = 1.39
    if codon == "CAG":
        score = .61
    if codon == "AAT":
        score = 1.1
    if codon == "AAC":
        score = .9
    if codon == "AAA":
        score = .84
    if codon == "AAG":
        score = 1.16
    if codon == "GAT":
        score = 1.36
    if codon == "GAC":
        score = .64
    if codon == "GAA":
        score = 1.15
    if codon == "GAG":
        score = .85
    if codon == "TCT":
        score = 1.47
    if codon == "TCC":
        score = 0.98
    if codon == "TCA":
        score = 1.44
    if codon == "TCG":
        score = .83
    if codon == "CCT":
        score = 0.52
    if codon == "CCC":
        score = 0.23
    if codon == "CCA":
        score = 2.75
    if codon == "CCG":
        score = 0.51
    if codon == "ACT":
        score = 1.34
    if codon == "ACC":
        score = 1.02
    if codon == "ACA":
        score = 1.15
    if codon == "ACG":
        score = 0.49
    if codon == "GCT":
        score = 1.64
    if codon == "GCC":
        score = 1.06
    if codon == "GCA":
        score = 0.99
    if codon == "GCG":
        score = 0.31
    if codon == "TGT":
        score = 1.14
    if codon == "TGC":
        score = .86
    if codon == "TGA":
        score = .59
    if codon == "TGG":
        score = 1.00   # ---
    if codon == "CGT":
        score = 1.84
    if codon == "CGC":
        score = .73
    if codon == "CGA":
        score = 1.07
    if codon == "CGG":
        score = .31
    if codon == "AGT":
        score = .76
    if codon == "AGC":
        score = .52
    if codon == "AGA":
        score = 1.79
    if codon == "AGG":
        score = .26
    if codon == "GGT":
        score = .7
    if codon == "GGC":
        score = .28
    if codon == "GGA":
        score = 2.85
    if codon == "GGG":
        score = .16
    return score

def save_obj(obj, name ):
    #demo = dir_archive(name, obj, serialized=True) 
    #demo.dump()
    #del demo
    with open(name +'.p', 'wb') as fp:
        pickle.dump(obj,fp,2)

#def load_obj(name ):
    #with open(name + '.pickle', 'rb') as f:
        #return pickle.load(f)
#plt.plot(0)                     # draw the graph
with open(sequenceExcel,'r') as t:
    sequencelines = t.readlines()  
    
getInfo(sequencelines)

with open(inFile,'r') as i:
    lines = i.readlines()    #print(lines)
#plt.figure(num=None, figsize=(100, 10), dpi=80, facecolor='w', edgecolor='k')
#plt.set_size_inches(18.5, 10.5)    
processedLines = findCodonFile(lines)
    


    #print(processedLines)
with open(outFile,'w') as o:
   for line in processedLines:
       o.write(line + '\n')
        
        
#print(str(rnadict))
save_obj(rnadict, 'outDict')

#plt.savefig('test.png')
#plt.show()                      # show it to me!

print("\nDone")
    
