import sys
input_file = sys.argv[1] 

### Collect positive strand ###
header = ''
pos_seq = ''
with open(input_file) as fp:
    for line in fp:
        if line[0] == '>': 
            header = line[1:].replace('\n','')
        else: 
            line = line.replace('\n','')
            pos_seq = pos_seq + line
            
### Generate negative strand ###
neg_seq = ''
def Replication(nucleotide):
    switch = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    return switch.get(nucleotide)
for pos_nucl in pos_seq: 
    neg_seq = Replication(pos_nucl) + neg_seq #Note: reverse complement!!!


### Find the ORFs in Pos_seq ###
three_mers_pos1 = []
three_mers_pos2 = []
three_mers_pos3 = []
for i in range(0, len(pos_seq)): 
    if len(pos_seq[i:i+3]) == 3: 
        if i%3 == 0: 
            three_mers_pos1.append(pos_seq[i:i+3])
        elif i%3 == 1: 
            three_mers_pos2.append(pos_seq[i:i+3])
        elif i%3 == 2: 
            three_mers_pos3.append(pos_seq[i:i+3])

### Find the ORFs in Neg_seq ###
three_mers_neg1 = []
three_mers_neg2 = []
three_mers_neg3 = []
for i in range(0, len(neg_seq)): 
    if len(neg_seq[i:i+3]) == 3: 
        if i%3 == 0: 
            three_mers_neg1.append(neg_seq[i:i+3])
        elif i%3 == 1: 
            three_mers_neg2.append(neg_seq[i:i+3])
        elif i%3 == 2: 
            three_mers_neg3.append(neg_seq[i:i+3])

### Define a function for extracting genes >=198 ###
def Extract(three_mers):
    #record indices of start codon
    start_ind = []
    ind = 0
    for mer in three_mers: 
        if mer == 'ATG':
            start_ind.append(ind)
        ind = ind +1
    #record indices of end codon
    end_ind = []
    ind = 0
    for mer in three_mers: 
        if mer == 'TAA' or mer == 'TAG' or mer == 'TGA': #Important or format that I always mess up
            end_ind.append(ind)
        ind = ind +1
    #now extract genes based on "start_ind" and "end_ind"
    ext_genes = []
    for i in start_ind: 
        for j in end_ind: 
            if i<j: 
                gene_l = three_mers[i:j+1]
                if len(gene_l) >= 66: 
                    gene_str = ''
                    for k in gene_l: 
                        gene_str = gene_str + k
                    ext_genes.append(gene_str)
                break
    return ext_genes

print("This program generates nested ORFs.")
print("+3:" + str(len(Extract(three_mers_pos3))) + "; +2:" + str(len(Extract(three_mers_pos2))) + "; +1:" + str(len(Extract(three_mers_pos1))) + "; -1:" + str(len(Extract(three_mers_neg1))) + "; -2:" + str(len(Extract(three_mers_neg2))) + "; -3:" + str(len(Extract(three_mers_neg3))))

### Write in genes_output.txt with all possible genes of the 6 ORFs ###
    
p_res = open("genes_output.txt", "w")

ORF_genes = {
    '+3':Extract(three_mers_pos3),
    '+2':Extract(three_mers_pos2),
    '+1':Extract(three_mers_pos1),
    '-1':Extract(three_mers_neg1),
    '-2':Extract(three_mers_neg2),
    '-3':Extract(three_mers_neg3),
}

for each in ORF_genes: 
    gene_list = ORF_genes[each]
    count = 1
    for gene in gene_list: 
        p_res.write(">Gene" + str(count) + "; " + each + " reading frame\n")
        p_res.write(gene + "\n")
        count = count +1

p_res.close()


