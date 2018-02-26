'''
Created on Jun 16, 2017

@author: Mark Holton
'''
import sys
import xlwt
import xlsxwriter
import matplotlib


import matplotlib.pyplot as plt

import pickle
from klepto.archives import dir_archive

from copy import deepcopy


dictFile = sys.argv[1]
rnadict = {  }
lengthOfrawicent=[]
lengthOfrawacent=[]
lengthOffiveicent=[]
lengthOffiveacent=[]
lengthOfchangeicent=[]
lengthOfchangeacent=[]

groupOfgenes=[]




plt.figure(num=None, figsize=(100, 10), dpi=80, facecolor='w', edgecolor='k')

with open(dictFile, 'rb') as fp:
  rnadict = pickle.load(fp)

print("dictionary loaded")
    


def table(label,ist):
    print(label+'\t', end="")

    for element in ist:
        print(str(element)+'\t', end="")
    print(' ')

def vertable(label,ist):
    print(label+':\t'+str(ist))
    
def tablemultytab(label,ist,numTabs):
    print(label+'\t', end="")
    tab=''
    for element in ist:
        i=0
        while(i<numTabs+1):
            tab=tab+'\t'
            i=i+1
        print(str(element)+tab, end="")
    print(' ')

def find_every_ten_percentile(raw,five,change,centile,base_centile):
    ret=[0,0,0,0,0, 0,0,0,0,0 ]
    temp=[10,10,10,10,10, 10,10,10,10,10,]
    indexoftemp=deepcopy(ret)
    counter=0
    for element in base_centile:
        if(element%10==0):
            ret[int(element/10)-1]=counter
        if(element%10 < temp[int(element/10)-1]):
            temp[int(element/10)-1]=element%10
            indexoftemp[int(element/10)-1]=counter
        counter=counter+1
    counter=0  
    for num in temp:
        if num!=10:
            lement = base_centile[int(indexoftemp[counter])]
            ret[int(lement/10)-1]=lement
        counter=counter+1
    ret=indexoftemp
    return ret

def complex_table(raw, five, change, centile, base_centile,percentiles_list):
    imraw=[]
    imfive=[]
    imchange=[]
    imcentile=[]
    for element in percentiles_list:
        imraw.append(raw[element])
        imfive.append(five[element])
        imchange.append(change[element])
        imcentile.append(centile[element])
    table('raw',imraw)
    table('five',imfive)
    table('change',imchange)
    table('centile',imcentile)

def graphScores(found_genes):
    
    for searchTerm in found_genes:
    
        raw= rnadict.get(searchTerm).get('scores')
        five= rnadict.get(searchTerm).get('sum_of_5_on_each_side')
        change= rnadict.get(searchTerm).get('change_by')
        base_centile= rnadict.get(searchTerm).get('raw_percentiles')
        splot(raw, "Raw Scores")
        splot(five, "Sum of 5 on each side")
        splot(change, "Change by")
        splot(base_centile, "Raw Percentiles")
        
        
def splot(yPoints,name):
    y = yPoints   # a list of numbers
    plt.plot(y)                     # draw the graph
    plt.savefig('images/'+name.strip()+'.png',bbox_inches='tight')
    plt.clf()        

def minindex(rawscores,add):     
    miindex =[add]
    index=0
    for score in rawscores:
        if(float(score) == float(rawscores[miindex[0]-add])):
            miindex.append(index+add)
        if(float(score) < float(rawscores[miindex[0]-add])):
            miindex=[index+add]
        
        index=index+1
    return miindex


def maxindex(rawscores,add):     
    maindex = [add]
    index=0
    for score in rawscores:
        if(float(score) == float(rawscores[ maindex[0] -add]) ):
            maindex.append(index+add)
        if(float(score) > float(rawscores[maindex[0]-add])):
            maindex=[index+add]
        
        index=index+1
    return maindex
    
def maxin(rawscores):     
    maxindex = 0
    maxscore=0
    index=0
    for score in rawscores:
        if(float(score) > float(rawscores[maxindex])):
            maxindex=index
            maxscore=score
        index=index+1
    return maxscore

def miax(raw,scoresystem):
    minin=minindex(raw,0)
    maxin=maxindex(raw,0)
    numofscores=len(raw)
    mincent=minin/numofscores
    maxcent=maxin/numofscores
    print('\nThe max of %s is the %dth codon and its percentile is %f' %(scoresystem,maxin,maxcent))
    print(str(otherMax))
    print('The min of the %s is the %dth codon and its percentile is %f' %(scoresystem,minin,mincent))
    print(str(otherMin))   

def search(searchTerm, catagory):
    found_genes = [] 
    if catagory=='Scores':
            found_genes.append(searchTerm)
            raw= rnadict.get(searchTerm).get('scores')
            five= rnadict.get(searchTerm).get('sum_of_5_on_each_side')
            change= rnadict.get(searchTerm).get('change_by')
            centile= rnadict.get(searchTerm).get('percentiles')
            base_centile= rnadict.get(searchTerm).get('raw_percentiles')


            table('raw',raw)
            table('five',five)
            table('change',change)
            table('centile',centile)
            
            miax('raw',raw)
            miax('five',five)
            miax('change',change)
            
            
            important_percentiles=find_every_ten_percentile(raw, five, change, centile, base_centile)
            print('\nImportant percentiles')
            complex_table(raw, five, change, centile, base_centile,important_percentiles)
            
    for gene in rnadict.keys():
        if catagory in rnadict[gene]:
            subdict = rnadict[gene].get(catagory)
            if searchTerm in subdict: 
                outPrint()
                print(rnadict[gene].get("geneName"))
                found_genes.append(rnadict[gene].get("geneName"))
                print(str(rnadict[gene].get("codons")))
                vertable('geneName',rnadict[gene].get("geneName"))
                vertable('Alias',rnadict[gene].get("Alias"))
                vertable('Number_of_reads',rnadict[gene].get("Number_of_reads"))
                vertable('Normalized_expression',rnadict[gene].get("Normalized_expression"))
                vertable('Median',rnadict[gene].get("Median"))
                vertable('PolyA_class',rnadict[gene].get("PolyA_class"))
                vertable('25%_percentile',rnadict[gene].get("25%_percentile"))
                vertable('75%_percentile',rnadict[gene].get("75%_percentile"))
                vertable('MAX_polyA',rnadict[gene].get("MAX_polyA"))
                vertable('MIN_polyA',rnadict[gene].get("MIN_polyA"))
                vertable('RANGE_polyA',rnadict[gene].get("RANGE_polyA"))
                vertable('Mean',rnadict[gene].get("Mean"))
                vertable('Geometric_mean',rnadict[gene].get("Geometric_mean"))
                vertable('Standard_deviation',rnadict[gene].get("Standard_deviation"))
                vertable('Fop',rnadict[gene].get("Fop"))
                vertable('Codon_class',rnadict[gene].get("Codon_class"))
        
                #print(rnadict[gene])
                print("\n")
    if len(found_genes)>2:
        print("Found Genes: "+str(found_genes).strip('[]'))
    input("Press enter to continue")
    if(input("Do you wish to graph the scores?\n")=="Yes"):
        graphScores(found_genes);
    if(input("Do you wish to export the scores to excel?\n")=="Yes"):
        output('CodonScores.xls',found_genes);
        
        
        
        
def output(filename,found_genes):
    book = xlwt.Workbook()
    for gene in found_genes:
      she = book.add_sheet(gene)
      searchTerm=gene
      raw= rnadict.get(searchTerm).get('scores')
      five= rnadict.get(searchTerm).get('sum_of_5_on_each_side')
      change= rnadict.get(searchTerm).get('change_by')
      centile= rnadict.get(searchTerm).get('percentiles')
      base_centile= rnadict.get(searchTerm).get('raw_percentiles')
      #print(str(raw))
      col=0
      she.write(0,col,'raw')
      she.write(1,col,'five')
      she.write(2,col,'change')
      she.write(3,col,'percentiles')
      she.write(4,col,'raw_percentiles')
      
      row=0
      col=1
      for num in raw:
        she.write(row, col, num)
        col+=1
      row=1
      col=6
      for num in five:
        she.write(row, col, num)
        col+=1
      row=2
      col=1
      for num in change:
        she.write(row, col, num)
        col+=1
      row=3
      col=1
      for num in centile:
        she.write(row, col, num)
        col+=1
      row=4
      col=1
      for num in base_centile:
        she.write(row, col, num)
        col+=1
      
      
           
    book.save(filename)
    
def getCentile(scores,total):
    index = 0
    output=[]
    while(index<len(scores)):
        output.append(scores[index]/total)
        index=index+1
    return output
    
def add5tolist(lisT):
    output =[]
    for num in lisT:
        output.append(num+5)
    return output

def retminmaxlist(gene):
    raw= rnadict[gene].get('scores')
    five= rnadict[gene].get('sum_of_5_on_each_side')
    change= rnadict[gene].get('change_by')
    #print(gene)
    minraw=minindex(raw,0)
    maxraw=maxindex(raw,0)
    numofscores=len(raw)
    rawicent=getCentile(minraw,numofscores)
    rawacent=getCentile(maxraw,numofscores)
    
    minfive=minindex(five,5)
    maxfive=maxindex(five,5)
    numofscores=len(raw)
    fiveicent=getCentile(minfive,numofscores)
    fiveacent=getCentile(maxfive,numofscores)
    
    minchange=minindex(change,0)
    maxchange=maxindex(change,0)
    numofscores=len(raw)
    changeicent=getCentile(minchange,numofscores)
    changeacent=getCentile(maxchange,numofscores)
    
    scores=[rawicent,rawacent,fiveicent,fiveacent,changeicent,changeacent]
    
    if(len(rawacent)>150):
      #input(gene+'------------------------------------')
      groupOfgenes.append(gene)
#K04H4.1
#ZK617.1
#C09D1.1
#K07E12.1
#F01G12.5

    else:
      lengthOfrawicent.append(len(rawicent))
      lengthOfrawacent.append(len(rawacent))
      lengthOffiveicent.append(len(fiveicent))
      lengthOffiveacent.append(len(fiveacent))
      lengthOfchangeicent.append(len(changeicent))
      lengthOfchangeacent.append(len(changeacent))
      
    return scores

def allmiaxtoExcel(filename):
    #raw min max five change median PolyA_class
    #file with all min max
    book = xlwt.Workbook()
    she = book.add_sheet('sheet2')
    row=0
    for gene in rnadict:
        if 'Median' in rnadict[gene]:
          retminmaxlist(gene)
    print('\n\n'+str(gene))
    #scores=[rawicent,rawacent,fiveicent,fiveacent,changeicent,changeacent]
    afone=1+maxin( lengthOfrawicent)
    aftwo=afone+maxin(lengthOfrawacent)
    afthree=aftwo+maxin(lengthOffiveicent)
    affour=afthree+maxin(lengthOffiveacent)
    affive=affour+maxin(lengthOfchangeicent)
    afsix=affive+maxin(lengthOfchangeacent)
    
    print(afone)
    print(aftwo)
    print(afthree)
    print(affour)
    print(affive)
    print(afsix)
    
    she.write(row,0,'geneName')
    she.write(row,1,'raw min')
    she.write(row,2,'raw max')
    she.write(row,3,'five min')
    she.write(row,4,'five max')
    she.write(row,5,'change min')
    she.write(row,6,'change max')
    she.write(row,7,'median')
    she.write(row,8,'PolyA_class')
    
    rows=[1,afone,aftwo,afthree,affour,affive]
    row=1
    for gene in rnadict:
        if gene=='K07E12.1':
            sam=0
        elif 'Median' in rnadict[gene]:
          col=0
          she.write(row,col,str(rnadict[gene].get('geneName')))
          col=1
          index=0
          for lis in retminmaxlist(gene):
              col=rows[index]
              if(col>256):
                break
              for num in lis:
                  if(col>256):
                    break
                  she.write(row, col, num)
                  col+=1
              index=index+1
          she.write(row,afsix,rnadict[gene].get('Median'))
          she.write(row,afsix+1,rnadict[gene].get('PolyA_class'))
          row+=1
        
    book.save(filename+'All')
        
        
def outPrint():
    print("\nFound\n")
    
    
def sGraph():
    book = xlsxwriter.Workbook('AllFIVEScoresForGraph.xlsx')
    she = book.add_worksheet('sheet2')
    tooLong=[]
    row=0
    for gene in rnadict:  
          if(row>=65536):
              print(gene)
              break
          col=0
          she.write(row,col,str(rnadict[gene].get('geneName')))
          col=1
          five = rnadict.get(gene).get('sum_of_5_on_each_side')
          raw= rnadict[gene].get('scores')
          minfive=minindex(five,5)
          maxfive=maxindex(five,5)
          numofscores=len(raw)
          fiveicent=getCentile(minfive,numofscores)
          fiveacent=getCentile(maxfive,numofscores)
          she.write(row, col, fiveicent[0])
          she.write(row, 2, fiveacent[0])
          row+=1
    book.close() 
    print('start\t'+str(tooLong))
    print(len(rnadict))



#if(input("Do you wish to export the scores to excel?\n")=="Yes"):
#        allmiaxtoExcel('CodonMinMaxScores.xls');
sGraph()
userStillUsing = True
while(userStillUsing):
    oneInput = input("What would you like to search by?(geneName    Alias    Scores    Number_of_reads    Normalized_expression    Median    PolyA_class    25%_percentile    75%_percentile    MAX_polyA    MIN_polyA    RANGE_polyA    Mean    Geometric_mean    Standard_deviation    Fop    Codon_class): \n")
    if(oneInput == "geneName"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Alias"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Scores"):
        twoInput = input(("What gene do you wish to search for?: \n" ))
        search(twoInput, oneInput)
    if(oneInput == "Number_of_reads"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Normalized_expression"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Median"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "PolyA_class"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "25%_percentile"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "75%_percentile"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "MAX_polyA"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "MIN_polyA"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "RANGE_polyA"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Mean"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Geometric_mean"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Standard_deviation"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Fop"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "Codon_class"):
        twoInput = input(("What %s do you wish to search for?: \n" % oneInput))
        search(twoInput, oneInput)
    if(oneInput == "close"):
        userStillUsing= False
    
    
    
    
