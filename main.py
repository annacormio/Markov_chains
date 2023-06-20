import itertools
import math
import random
import pandas as pd


#given a sequence file extract the subsequences starting and ending at positions indexed at i in the relative start and end lists
def seq_extraction (start,end,file):
    #create list of CpG strings
    l=[] #initialize list od CpG seq
    for i in range(len(start)): #for each iterated index
        l.append(file[start[i]:(end[i]+1)]) #append in l the indexed CpG islands substrings
    return l



#given the CpG island sequences extracted list and the file it computes randomly the start positions and extract seq starting from them of the same length of the CpG islands
def non_CpG_seq_extraction(cpg,file):
    start=[] #list of starting points
    end=[] #list of end points
    for e in cpg: #iterate through each CpG sequence in the list
        s_i=random.randint(0,len(file)-1) #random value for start index position
        start.append(s_i) #append the position to the start indeces list
        end.append(s_i+len(e)) #end position is the start + the length of the current CpG island
    return seq_extraction(start,end,file) #extract sequences from the file given start and end list indexes




#given a list of sequences it returns the model built on those sequences
def model(l):
    nt = ['A', 'C', 'G', 'T'] #nt. list
    df = pd.DataFrame(index=nt, columns=nt)  # intialize dimer codintional frequency matrix
    for p in itertools.product(nt, repeat=2): #p = all possible dimers (heterodimers and homodimers) built with the nt.
        dim = p[0] + p[1] #p is a tuple so i converti it into a string by adding the 2 string cmoponents
        count_dim = 0 #initialize count of the dimer
        count_gen_dim = 0 #initialize count of the first nt. of the dimer
        if p[0]==p[1]: #if we have an homodimer
            c=p[0]
            for s in l: #for each sequence in l
                for i in range(len(s)-1): #for each character in the sequence
                    if s[i]==c and s[i+1]==c: #if the following character is equal to the current one I add one to the count
                        count_dim+=1
                count_gen_dim += s.count(c)  # count all dimers starting with same character as the dimer
            freq = count_dim / count_gen_dim  # calculate frequency of p dimer over all possible dimers starting with same nt.
            df.loc[[p[0]], [p[1]]] = freq  # insert in the df

        else: #dimer is not an homodimer
            for s in l:  # for each sequence in the CpG sequences list
                count_dim += s.count(dim)  # count p specific dimer in the seq
                count_gen_dim +=s.count(p[0])   # count all dimers starting with same character as the dimer
            freq = count_dim / count_gen_dim  # calculate frequency of p dimer over all possible dimers starting with same nt.
            df.loc[[p[0]], [p[1]]] = freq  # insert in the df
    return df



#given the inside model in_model and outside model out_model returns the S score for the query sequence q
def is_CpG(in_model, out_model, q):
    inside=math.log(0.25,10) #log0.25 is the initial probability of the first nt. of the query
    outside=math.log(0.25,10)
    for n in range(len(q)-1): #iterate through each nt. of the query
        dim=q[n:n+2] #extract the dimer starting with n
        inside=inside+math.log(in_model.loc[dim[1],dim[0]],10) #add the log of the conditional proability of that dimer to the inside score
        outside=outside+math.log(out_model.loc[dim[1], dim[0]],10)#add the log of the conditional proability of that dimer to the outside score
    S=inside-outside
    #print('log of probabilities:',S)
    return S


#returns the window w offsets of the genome g when the log ratio computed is >S calculated with the inside and outised model passed in the function
def is_CpG_windowed(g,w,in_model,out_model):
    offsets=[] #initialize list that store CpG island windows
    for o in range(len(g)-w): #iterate through each position in the genome until the last possible window selection
        r=is_CpG(in_model,out_model,g[o:o+w]) #get the S score for that window
        if r>0: #if S>0 then the window is a CpG island
             offsets.append((r,o+1)) #its result and its offset position in the genome are added to the list
    return offsets




if __name__ =="__main__":

    #OPEN CHR22 SEQUENCE FILE
    I = open("chr22.fa","r")  # open file
    f = I.read().replace('\n','')  # read file content concatenating all lines as one string
    I.close() #close file
    file = f[6:].upper()  # slice first 6 characters of the string (header) and convert all characters in uppercase


    #RETIEVE START AND END POSITIONS OF CPG ISLANDS SEQUENCES
    index = pd.read_csv("index.txt", sep='\t')  # convert index txt file (each column separated by a tab \t) in dataframe
    index_chr22 = index[index["chr"] == "chr22"]  # dataframe of chromosome 22 information on CpG islands
    # print(index_chr22)
    start = list(index_chr22["start"])  # store column of start positions of chr22 in a list
    end = list(index_chr22["end"])  # store column of end positions of chr22 in a list
    # print('start position list \n', start)
    # print('end position list \n',end)

    #BUILD CPG AND NON_CPG MODELS
    seq_CpG=seq_extraction(start,end,file) #extract CpG islands sequences from file
    CpG_model=model(seq_CpG) #build conditional frequency matrix for all 16 dimers (model)
    print('\n CpG_model:\n',CpG_model)

    without_N_file = file.replace('N', '') #remove N from the sequence
    seq_non_CpG=non_CpG_seq_extraction(seq_CpG,without_N_file) #extract random sequences from the file of the same length of the CpG island seuqences
    non_CpG_model=model(seq_non_CpG) #build model for non CpG island sequences
    print('\n non_CpG_model:\n',non_CpG_model)

    #QUERY A SEQUENCE WITH THE MODEL CREATED
    query="".join(random.choice('ATCG') for i in range(100))
    #query='GCGCGGCAAG'  #test--> it is a CpG island
    infer=is_CpG(CpG_model,non_CpG_model,query)
    print(f'Is "{query}" a CpG island sequence?\n',infer>0)

    # BONUS: QUERY SUSEQUENT WINDOWS OF THE GENOME WITH THE MODEL CREATED
    #genome="".join(random.choice('ATCG') for i in range(1000))
    s=random.randint(0,len(file)-1)
    genome=without_N_file [s:s+1000]
    #print(genome)
    tot_len=0
    for s in seq_CpG:
        tot_len+= len(s)
    avg=round(tot_len/len(seq_CpG))
    #print(avg)
    print(is_CpG_windowed(genome,avg,CpG_model,non_CpG_model))














