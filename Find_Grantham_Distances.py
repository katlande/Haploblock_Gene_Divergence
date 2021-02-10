#Python2.7
#This script finds the Grantham distance between two aligned peptide fastas. 
#This script does not allow for Gaps.

#Function that reads in a peptide FASTA
def extract_sequences(source):
    source_fastaFile = open(source, 'r')
    source_sequences = []
    for line in source_fastaFile:
line = line.strip()
if ">" in line:
    line = "END"
    source_sequences.append(line)
elif line == "\n" or line == "\r\n" or line == "\r" or line == "":#don't include empty lines
    next
else:
    source_sequences.append(line)
    #print line
    source_fastaFile.close() #create a list called "sequences" with sequence data

    source_sequences.append("END")
    #print source_sequences

    #Merge individual sequences into one line:
    pos = 1
    source_sequences2 = []
    while pos < (len(source_sequences)):
#print "pos", pos
if source_sequences[pos] == "END":
    pos = pos + 1
    #print "skipped index", pos
else:
    if source_sequences[pos+1] == "END":
#print "single index", pos
q = source_sequences[pos]
source_sequences2.append(q)
pos = pos + 1
#print pos
    else:
minicounter = pos+1
temp_rounds = 1
while not (source_sequences[minicounter] == "END"):
    if temp_rounds == 1:
single_line_fa = source_sequences[pos] + source_sequences[minicounter]
#print single_line_fa
    else:
single_line_fa = single_line_fa + source_sequences[minicounter]
#print single_line_fa
    #print "merged index", pos, minicounter
    temp_rounds = 1 + temp_rounds
    minicounter = minicounter + 1
source_sequences2.append(single_line_fa)
pos = minicounter + 1
    return(source_sequences2)

#Function that reads in tab-separated text files:
def read_tsv(input):
    tab = open(input, 'r')
    tab_list = []
    for line in tab:
line=line.strip()
line=line.split("\t")
tab_list.append(line)
    tab.close()
    return(tab_list)

#read in the Grantham Matrix (matrix of grantham distances across all)
grantham_table = read_tsv("~/Grantham_Matrix.txt")
#print grantham_table

#Read in two peptide fastas, one for each haplotype:
Hap2 = extract_sequences("~/gene_Hap2.prot.fa")
Hap0 = extract_sequences("~/gene_Hap0.prot.fa")

#Reformat peptide strings:
Hap2 = str(Hap2)
Hap2 = Hap2.replace("'", "")
Hap2 = Hap2.replace("]", "")
Hap2 = Hap2.replace("[", "")

Hap0 = str(Hap0)
Hap0 = Hap0.replace("'", "")
Hap0 = Hap0.replace("]", "")
Hap0 = Hap0.replace("[", "")

#counters for missense sites, amino acid grantham value, and stop codons
stop_counter_out = 0
grantham_value_out  = 0
missense_counter_out  = 0
for aa in range(0,len((Hap0))):
    #If fastas are the same, move to next amino acid
    if not (Hap0)[aa] == (Hap2)[aa]:
        #If one fasta has a novel stop codon, add one to the stop codon counter
        if (Hap0)[aa] == "*" or (Hap2)[aa] == "*":
        stop_counter_out  = stop_counter_out  + 1
        #If the fastas are different but neither contains a stop codon, extract the grantham distance between the two amino acids and add it to the peptide total, then add 1 to the missense counter.
        else:
            column = (str(grantham_dict[(Hap0)[aa]])).strip("]")
            column = int(column.strip("["))
            row = (str(grantham_dict[(Hap2)[aa]])).strip("]")
            row = int(row.strip("["))
            grantham_value_out  = grantham_value_out + int(grantham_table[column+1][row+1])
            missense_counter_out  = missense_counter_out  + 1

#Output the result into a list called "temp_gene_data"
temp_gene_data = [gene, grantham_value, stop_counter, missense_counter]

#This script can be easily looped for many genes; output each "temp_gene_data" as a new row in a DF




