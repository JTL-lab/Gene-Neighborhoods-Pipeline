#!/usr/bin/env python
# coding: utf-8



#get_ipython().system('pip install dna_features_viewer')
import  dna_features_viewer
from dna_features_viewer import GraphicFeature, GraphicRecord




################Imports#################
import numpy as np
import pandas as pd
import os
import glob
#get_ipython().system('pip install biopython')
from Bio import SeqIO
import itertools
#et_ipython().system('pip install reportlab')
import sys


### Specify the path of the RGI files and extract the .txt files from the folder
#path = os.path.join(os.path.expanduser('~'),'Fin_RGI')

path=sys.argv[1]
#print (path)
readfiles=glob.glob(os.path.join(path,"*.txt"))

####Read each rgi file(.txt) using pandas and store them in a dataframe####
dataframelist=[]
filenames=[]
for i in readfiles:
    a=os.path.basename(i)
    a=a.split(".")
    #print(a[0])
    filenames.append(a[0])
    i=pd.DataFrame(pd.read_csv(i,sep="\t"))
    dataframelist.append(i)




datadict={}### dictionary of dataframes with keys as filenames ###
for i in range(len(filenames)):
    a=dataframelist[i]
    b=filenames[i]
    datadict[b]=pd.DataFrame(a)  

def parti(c):
  head,sep,tail=c.partition("_")
  return (head)


def locus_generator(c,genome_name):
    a,b=c.split("_")
    tag=genome_name+"_"+b
    #print(tag)
    return tag


for k,v in datadict.items():
    newcon=[]
    Locus_Temp=[]
    for i in v["Contig"]:
        newcon.append(parti(i))
        Locus_Temp.append(locus_generator(i,k))
    v['req_cont']=newcon    ###adding a new column "req_cont" into dataframe
    v['Locus_Tag']=Locus_Temp

####function to extract the required data from gbk files(entire genomes:) using biopython####    
def extract(infile):
  gene_start=[]
  gene_end=[]
  gene_strand=[]
  gene_name=[]
  loc_tag=[]
  function=[]
  protein_seq=[]
  contig_name=[]
  unique=[]

  for index, record in enumerate(SeqIO.parse(infile, "genbank")):
    #print("index %i, ID = %s, length %i, with %i features"% (index, record.id, len(record.seq), len(record.features)))
    for i in record.features:
        if i.type == "CDS" and "gene" in i.qualifiers:
          locations=i.location
          gene_start.append(locations.start)
          gene_end.append(locations.end)
          gene_strand.append(locations.strand)
          loc_tag.append(i.qualifiers['locus_tag'])
          function.append(i.qualifiers['product'])
          protein_seq.append(str(i.qualifiers['translation']))
          gene_name.append(i.qualifiers['gene'])
          contig_name.append(record.id)
        elif i.type =="CDS":
          locations=i.location
          gene_start.append(locations.start)
          gene_end.append(locations.end)
          gene_strand.append(locations.strand)
          loc_tag.append(i.qualifiers['locus_tag'])
          function.append(i.qualifiers['product'])
          protein_seq.append(str(i.qualifiers['translation']))
          gene_name.append("Unidentfied")
          contig_name.append(record.id)
          
  salmonella_gene_frame=pd.DataFrame()
  salmonella_gene_frame['GeneStart']=gene_start
  salmonella_gene_frame['GeneEnd']=gene_end
  salmonella_gene_frame['GeneStrand']=gene_strand
  salmonella_gene_frame['Locus_Tag']=loc_tag
  salmonella_gene_frame['GeneName']=gene_name
  salmonella_gene_frame['Product']=function
  salmonella_gene_frame['ProteinSequence']=protein_seq
  salmonella_gene_frame['contig_name']=contig_name
  
  #print(contig_name)
  for i in contig_name:
    if i not in unique:
      unique.append(i)


  return salmonella_gene_frame,unique

####specifying the location/path names og gbk files and loading them####
#filepath=os.path.join(os.path.expanduser('~'),'Fin_gbk')
filepath=sys.argv[2]
gbkfiles=glob.glob(os.path.join(filepath,"*.gbk"))


# In[12]:


gbk_names=[]
for i in gbkfiles:
    a=os.path.basename(i)
    a=a.split(".")
    #print(a[0])
    gbk_names.append(a[0])

## generating the keys various dictionaries required  ###

gbknames=[]
uniquenames=[]
datasetslist=[]
carbapene_gene_list=[]
neighboring_10list=[]
for i in gbk_names:
    uniquenames.append(str(i))
    datasetslist.append(str(i))


### extract the data from gbk files and store in a dict of dataframes####
gbkdict={}
uniquedict={}
for i in range(len(gbk_names)):
    a=gbk_names[i]
    b=uniquenames[i]
    gbkdict[a],uniquedict[b]=extract(gbkfiles[i])



#### manipulate the locus tag and protein sequence ###
for j,i in gbkdict.items():
    i['Locus_Tag']=i['Locus_Tag'].apply(lambda i:str(i).replace("[","").replace("]","").replace("'",""))
    i['ProteinSequence']=i['ProteinSequence'].str.strip('[]')

####function to create groups based on contigs ###

def make_groups(frame):
  #print(frame)
  
  group=frame.groupby(frame['contig_name'])
  datasets = {}
  for groups, data in group:
    datasets[groups] = data
  return datasets

#### creating a dict for groups #####
datasetdict={}
for i in range(len(gbk_names)):
    for j,k in gbkdict.items():
        datasetdict[j]=make_groups(k)
        
        
#### find the unique drugclasses of all the genomes and store it in uniquedrugdict####
uniquedrugdict={}

for j,k in datadict.items():
    uniquedrugclasses=[]
    for l in range(len(k)):
        if k['Drug Class'][l]  not in uniquedrugclasses:
                uniquedrugclasses.append(k["Drug Class"][l])
    uniquedrugdict[j]=uniquedrugclasses
        


###### having the minimum has static value extract the genome that has minimum number of drug classes has key##
min=100
minarray=[]
for i,j in uniquedrugdict.items():
    #print(len(j))
    if len(j)<min:
        min=len(j)
    if len(j)==min:
        minarray.append(i)
min_drug_key=minarray[0]
print(min_drug_key)


listofdrugnames_modified=[]
for i,j in uniquedrugdict.items():
    if i==min_drug_key:
        for k in j:
            if len(k.split(";"))>1:
                listofdrugnames_modified.append(k.split("; ")[0].split(" ")[0]+str("_")+ k.split(";")[1].split(" ")[1])
            else:
                listofdrugnames_modified.append(k.split("; ")[0].split(" ")[0])


#creating separate dictionary for all the unique drugclasses of the genome "key"
def createeachdict_drug(drugindex):
    
    temp_dict={}
    for j,k in datadict.items():
        temp_dict[j] =k[k['Drug Class']==drugindex]
    
    return temp_dict
   

main_dictionary={}
for i in range(len(listofdrugnames_modified)):
    main_dictionary[listofdrugnames_modified[i]]=createeachdict_drug(uniquedrugdict[min_drug_key][i])

Dict_multigene_instances={}
for i,j in main_dictionary.items():
    Dict_multigene_instances[i]={k: v for k, v in j.items() if len(v)>1}


Dict_singlegene_instances={}
for i,j in main_dictionary.items():
    #print(i)
    temp_dict={}
    for a,b in j.items():
        if len(b)==1:
            temp_dict[a]=b
            
    Dict_singlegene_instances[i]=temp_dict


emptykeyslist=[]
for i,j in Dict_singlegene_instances.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del Dict_singlegene_instances[i]

emptykeyslist=[]
for i,j in Dict_multigene_instances.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del Dict_multigene_instances[i]




unique_ARO_dict={}
unique_ARO_names=[]
for j,k in Dict_multigene_instances.items():
    for n,m in k.items():   
        m.reset_index(drop=True, inplace=True)
        unique_ARO_names.append(m['Best_Hit_ARO'])
        unique_ARO_dict[j]=m['Best_Hit_ARO']


Required_drug_names=[]
Required_drug_dictkeys=[]

for i in range(len(unique_ARO_names)):
    for j in unique_ARO_names[i]:
        #print(j)
        if j not in Required_drug_names:
            Required_drug_names.append(j)


for i in Required_drug_names:
    Required_drug_dictkeys.append(i.split(" ")[0])



def createARO_drug(drug):
    temp_dict={}
    temp_check_len_aray=[]
            
    for j,k in datadict.items():
        temp_check_len_aray=k[k['Best_Hit_ARO']==drug]
        if len(temp_check_len_aray)>0:
            temp_dict[j]=temp_check_len_aray
                       
    return temp_dict


Dict_each_drug_len_greater2={}
for i in range(len(Required_drug_dictkeys)):
    Dict_each_drug_len_greater2[Required_drug_dictkeys[i]]=createARO_drug(Required_drug_names[i])



def contigend_visualization(contig_array,end_direction,genome):
    
    Features=[]
    temp_array_totrack_length=[]   
      
    contig_array.reset_index(drop=True, inplace=True)
    
    
    if end_direction=="upward":
        b=GraphicFeature(start=contig_array["GeneStart"][0]-1500, end=contig_array["GeneStart"][0], strand=+1,color="#0A090A",label="Contig_Ends")
        Features.append(b)  
        temp_array_totrack_length.append((contig_array["GeneStart"][0]-1500,contig_array["GeneStart"][0]))
    
    
    for i in range(len(contig_array)):
            #print(contig_array["GeneStart"][i])
            a=GraphicFeature(start=contig_array["GeneStart"][i], end=contig_array["GeneEnd"][i], strand=contig_array["Strand"][i],color=contig_array["Genecolor"][i],label=str(contig_array["GeneName"][i]))
            Features.append(a)
            temp_array_totrack_length.append((contig_array["GeneStart"][i],contig_array["GeneEnd"][i]))
    
            
    if end_direction=="downward":
        b=GraphicFeature(start=contig_array["GeneEnd"].iloc[-1]+1500, end=contig_array["GeneEnd"].iloc[-1]+2500, strand=+1,color="#0A090A",label="Contig_Ends")
        Features.append(b) 
        temp_array_totrack_length.append((contig_array["GeneEnd"].iloc[-1]+1500,contig_array["GeneEnd"].iloc[-1]+2500))
    
   
    length=temp_array_totrack_length[-1][1] - temp_array_totrack_length[0][0]

    record = GraphicRecord(first_index=temp_array_totrack_length[0][0],sequence_length=length,features=Features)  
    ax,_=record.plot(figure_width=20,strand_in_label_threshold=7)
    name=str(genome)+".png"
    ax.figure.savefig(name)
            
    


# In[51]:


### Finding neighbours(5 upstream and 5 dowmstream genes for each card gene)
def find_neighbor(rec,uname,data,number_of_genes,drug,genome):
  neighbor_genes=[]
  rec.reset_index(drop=True, inplace=True)

    

  for j in range(len(rec['Start'])):
        m=[]
        n=[]
        upwardgenes=[]
        downwardgenes=[]
        recarray=[]
        i=rec.loc[j].Start
        k=rec.loc[j].req_cont
        #print(i)
        if number_of_genes==10:
            g=10
        elif number_of_genes==14:
            g=14
            
        if k in uname:
            newlist=data[k]
            newlist.reset_index(drop=True, inplace=True)
            for l in range(len(newlist)):
                if newlist['GeneStart'][l]> i:
                    downwardgenes.append((newlist['GeneStart'][l],newlist['GeneEnd'][l],newlist['GeneStrand'][l],newlist['Locus_Tag'][l],newlist['ProteinSequence'][l],newlist['GeneName'][l],"#ffcccc"))
                else:
                    upwardgenes.append((newlist['GeneStart'][l],newlist['GeneEnd'][l],newlist['GeneStrand'][l],newlist['Locus_Tag'][l],newlist['ProteinSequence'][l],newlist['GeneName'][l],"#ccccff"))
                #print(upwardgenes) 
            
            newu=pd.DataFrame(upwardgenes,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor"])
            newd=pd.DataFrame(downwardgenes,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor"])
            m=(newu.iloc[(newu['GeneStart']-i).abs().argsort()[:g]]).sort_values(by="GeneStart")
            #print(m)
            
            n=(newd.iloc[(newd['GeneStart']-i).abs().argsort()[:g]]).sort_values(by="GeneStart")
            #print(n)
            
            
            recarray.append((rec['Start'][j],rec['Stop'][j],rec['Orientation'][j],rec['Locus_Tag'][j],rec['Predicted_Protein'][j],rec['Best_Hit_ARO'][j],"#C70039"))
            o=pd.DataFrame(recarray,columns=["GeneStart","GeneEnd","Strand","Locus_Tag","ProteinSequence","GeneName","Genecolor"])
            
            
            e = pd.concat([m,o,n], ignore_index=True)
            e.reset_index(drop=True, inplace=True)
            neighbor_genes.append(e)
            #print(e)      
                       
            
        else:
            print(k)
            print(uname)
            print("contig does not exist----"+"drugclass:"+drug)
            print(rec['Best_Hit_ARO'][j])
            
  if(len(m)!=0):
    if(len(n)<g ):
        
        #print("contig ends downward------- : "+str(drug)+str("at position -------")+ str(len(n)))
        if number_of_genes==10:
            contigend_visualization(e,"downward",genome)


        
        #contigend_visualization(plotindex=1)
        neighbor_genes=[]
        
   
    elif (len(m)<g):
        #print("contig ends upward-------: "+str(drug)+str("at position ----") + str(len(m)))
        if number_of_genes==10:
            contigend_visualization(e,"upward",genome)
#     
        
        neighbor_genes=[]          

  return neighbor_genes


# In[188]:


multipleinstance_neighboringdict_combined_range_10={}

for i,j in Dict_each_drug_len_greater2.items():
    Dict_neighboring_genes_range10={}
    for k,l in j.items():
        temp=find_neighbor(j[k],uniquedict[k],datasetdict[k],10,i,k)
        if len(temp)>0:
            Dict_neighboring_genes_range10[k]=temp 
    multipleinstance_neighboringdict_combined_range_10[i]=Dict_neighboring_genes_range10
            

emptykeyslist=[]
for i,j in multipleinstance_neighboringdict_combined_range_10.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del multipleinstance_neighboringdict_combined_range_10[i]


multipleinstance_neighboringdict_combined_range_14={}

for i,j in Dict_each_drug_len_greater2.items():
    Dict_neighboring_genes_range14={}
    for k,l in j.items():
        temp=find_neighbor(j[k],uniquedict[k],datasetdict[k],14,i,k)
        if len(temp)>0:
            Dict_neighboring_genes_range14[k]=temp
    multipleinstance_neighboringdict_combined_range_14[i]=Dict_neighboring_genes_range14


emptykeyslist=[]
for i,j in multipleinstance_neighboringdict_combined_range_14.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del multipleinstance_neighboringdict_combined_range_14[i]



singleinstance_neighboringdict_combined_range10={}
for i,j in Dict_singlegene_instances.items():
    Dict_neighboring_genes_range10={}
    for k,l in j.items():
        temp=find_neighbor(j[k],uniquedict[k],datasetdict[k],10,i,k)
        if len(temp)>0:
            Dict_neighboring_genes_range10[k]=temp
    singleinstance_neighboringdict_combined_range10[i]=Dict_neighboring_genes_range10


emptykeyslist=[]
for i,j in singleinstance_neighboringdict_combined_range10.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del singleinstance_neighboringdict_combined_range10[i]



singleinstance_neighboringdict_combined_range14={}

for i,j in Dict_singlegene_instances.items():
    Dict_neighboring_genes_range14={}
    for k,l in j.items():
        temp=find_neighbor(j[k],uniquedict[k],datasetdict[k],14,i,k)
        if len(temp)>0:
            Dict_neighboring_genes_range14[k]=temp
    singleinstance_neighboringdict_combined_range14[i]=Dict_neighboring_genes_range14


emptykeyslist=[]
for i,j in singleinstance_neighboringdict_combined_range14.items():
    if len(j)==0:
        emptykeyslist.append(i)
for i in emptykeyslist:
    del singleinstance_neighboringdict_combined_range14[i]


### function to get separate dicts of locus tags,protein sequences and gene names inorder to compare and write to a fasta file####
def getrequiredgenes(frame,number_of_genes,drug):
    temp_locus_array=[]
    temp_protein_array=[]
    temp_genename=[]
    for i in frame:
        temp_locus_array.append(list(i['Locus_Tag']))
        temp_protein_array.append(i['ProteinSequence'])
        temp_genename.append((i['GeneName']))
 
    a=temp_locus_array
    a=list(itertools.chain.from_iterable(a))
   
    b=temp_protein_array
    b=list(itertools.chain.from_iterable(b))
   
    c=temp_genename
    c=list(itertools.chain.from_iterable(c))
    
  
    if number_of_genes==10:
        if len(frame)>1:
            return(a[:10]+a[-10:], b[:10]+b[-10:],c[:10]+c[-10:])## for more than one card genes in a genomes
        else:
            return a, b,c
        #return(a[:5]+a[-5:], b[:5]+b[-5:],c[:5]+c[-5:])## for more than one card genes in a genomes
    if number_of_genes==14:
        if len(frame)>1:
            return(a[:14]+a[-14:], b[:14]+b[-14:],c[:14]+c[-14:])## for more than one card genes in a genomes
        else:
            return a,b,c



Drug_multiple_instance_range_10_locus_dict10={}
Drug_multiple_instance_range_10_protein_dict10={}
Drug_multiple_instance_range_10_genename_dict10={}

for i,j in multipleinstance_neighboringdict_combined_range_10.items():
    print(i)
    locustags_dict_10={}
    protein_dict_10={}
    genename_dict_10={}
    for k,l in j.items():
        locustags_dict_10[k],protein_dict_10[k],genename_dict_10[k]=getrequiredgenes(j[k],10,i)
    #print(len(locustags_dict_10),len(protein_dict_10),len(genename_dict_10))  
    Drug_multiple_instance_range_10_locus_dict10[i]=locustags_dict_10
    Drug_multiple_instance_range_10_protein_dict10[i]=protein_dict_10
    Drug_multiple_instance_range_10_genename_dict10[i]=genename_dict_10



Drug_multiple_instance_range_14_locus_dict14={}
Drug_multiple_instance_range_14_protein_dict14={}
Drug_multiple_instance_range_14_genename_dict14={}

for i,j in multipleinstance_neighboringdict_combined_range_14.items():
    locustags_dict_14={}
    protein_dict_14={}
    genename_dict_14={}
    for k,l in j.items():
        locustags_dict_14[k],protein_dict_14[k],genename_dict_14[k]=getrequiredgenes(j[k],14,i)
    #print(len(locustags_dict_10),len(protein_dict_10),len(genename_dict_10))  
    Drug_multiple_instance_range_14_locus_dict14[i]=locustags_dict_14
    Drug_multiple_instance_range_14_protein_dict14[i]=protein_dict_14
    Drug_multiple_instance_range_14_genename_dict14[i]=genename_dict_14


Drug_singleinstance_range_10_locus_dict10={}
Drug_singleinstance_range_10_protein_dict10={}
Drug_singleinstance_range_10_genename_dict10={}

for i,j in singleinstance_neighboringdict_combined_range10.items():
    locustags_dict_10={}
    protein_dict_10={}
    genename_dict_10={}
    for k,l in j.items():
        locustags_dict_10[k],protein_dict_10[k],genename_dict_10[k]=getrequiredgenes(j[k],10,i)
    #print(len(locustags_dict_10),len(protein_dict_10),len(genename_dict_10))  
    Drug_singleinstance_range_10_locus_dict10[i]=locustags_dict_10
    Drug_singleinstance_range_10_protein_dict10[i]=protein_dict_10
    Drug_singleinstance_range_10_genename_dict10[i]=genename_dict_10


# In[176]:


Drug_singleinstance_range_14_locus_dict14={}
Drug_singleinstance_range_14_protein_dict14={}
Drug_singleinstance_range_14_genename_dict14={}

for i,j in singleinstance_neighboringdict_combined_range14.items():
    locustags_dict_14={}
    protein_dict_14={}
    genename_dict_14={}
    for k,l in j.items():
        locustags_dict_14[k],protein_dict_14[k],genename_dict_14[k]=getrequiredgenes(j[k],14,i)
    #print(len(locustags_dict_10),len(protein_dict_10),len(genename_dict_10))  
    Drug_singleinstance_range_14_locus_dict14[i]=locustags_dict_14
    Drug_singleinstance_range_14_protein_dict14[i]=protein_dict_14
    Drug_singleinstance_range_14_genename_dict14[i]=genename_dict_14



import os.path

save_path = 'Multiple_instance/'


for i,j in multipleinstance_neighboringdict_combined_range_14.items():
    #print(i)
    #print(main_protein_dict16[i])
    #name=i+".fasta"
    name = os.path.join(save_path, i+".fasta") 
    #print(name)
    filef=open(name,"w")
    for k in j.keys():
        for a,b in zip(Drug_multiple_instance_range_14_locus_dict14[i][k],Drug_multiple_instance_range_14_protein_dict14[i][k]):
            #print(b)
            filef.writelines(">" + '{}'.format(a))
            filef.write("\n")
            filef.writelines('{}'.format(b))
            filef.write("\n")
    filef.close()


# In[164]:
save_path = 'Single_instance/'

for i,j in singleinstance_neighboringdict_combined_range14.items():
    #print(i)
    #print(main_protein_dict16[i])
    #name=i+".fasta"
    name = os.path.join(save_path, i+".fasta") 
    #print(name)
    filef=open(name,"w")
    for k in j.keys():
        for a,b in zip(Drug_singleinstance_range_14_locus_dict14[i][k],Drug_singleinstance_range_14_protein_dict14[i][k]):
            #print(b)
            filef.writelines(">" + '{}'.format(a))
            filef.write("\n")
            filef.writelines('{}'.format(b))
            filef.write("\n")
    filef.close()

