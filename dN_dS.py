#Python/2.7
#Find dN/dS between two gene variants using aligned, no-gap DNA and peptide FASTAs


def assign_sequences(fasta_names, fasta_seq, file_name, gene, inv):
    k = 0
    search = str(inv) + "_" + str(file_name) + "_" + str(gene)
    while k < len(fasta_names):
        if not fasta_names[k] == search:
            chrom = "nothing"
            next
            k = k +1
        else:
            #print "match"
            chrom = k
            break
    if chrom == "nothing":
        print "no matching sequence name"
    else:
        return(fasta_seq[chrom])

#Extract a sequence from a fasta by inputing a path to the file:
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

#Find the number of possible missense and sense mutations from any given DNA codon:
def find_n_and_s(codon, aa, genetic_code):
    
    #print "main codon:", codon
    codon =codon.replace("a", "A")
    codon = codon.replace("c", "C")
    codon = codon.replace("t", "T")
    codon = codon.replace("g", "G")
    
    b1 = codon[0]
    b2 = codon[1]
    b3 = codon[2]
    n = 0
    s = 0

    for site in range(0, 3):
        current_base = codon[site]
        
        if codon[0] == "N" or codon[0] == "n" or codon[1] == "N" or codon[1] == "n" or codon[2] == "N" or codon[2] == "n":
            print codon, "contains N value"
            break
        
        #prevent self-mutations:
        if current_base == 'A' or current_base == 'a':
            subs = ['G', 'C', 'T']
        elif current_base == 'T' or current_base == 't':
            subs = ['G', 'C', 'A']
        elif current_base == 'C' or current_base == 'c':
            subs = ['G', 'A', 'T']
        elif current_base == 'G' or current_base == 'g':
            subs = ['A', 'C', 'T']
        else:
            print "non-IUPAC base:", current_base, "in", codon
            break

        temp_codon = [b1,b2,b3]
        for sub in range(0, 3):
            #replace each codon pairwise
            temp_codon[site] = subs[sub]
            #look up temp_codon in gc
            temp_codon2 = temp_codon[0] + temp_codon[1] + temp_codon[2]
            #print "codon =", temp_codon2, "aa =", genetic_code[temp_codon2]
            temp_aa = genetic_code[temp_codon2]
            temp_aa = str(temp_aa)
            
            #print "subbing", subs[sub], "at site", site, "in", codon, "to make", temp_codon2, "\noriginal aa was", aa, "new aa is", temp_aa, "\n"
            
            if temp_aa == aa:
                s = s + 1
            elif not temp_aa == aa:
                n = n + 1
            else:
                print "unknown error"
    
    n_and_s = [n,s]
    return(n_and_s)

#Run the dnds for one gene:
# Reference_Genome = ancestral DNA sequence
# Reference_Genome_p = ancestral PEP sequence
# Fixed_sites = Haplotype 1 DNA sequence
# Fixed_sites_p = Haplotype 1 PEP sequence
# Fixed_sites_pair = Haplotype 2 DNA sequence
# Fixed_sites_pair_p = Haplotype 2 PEP sequence
def get_dnds(Reference_Genome, Reference_Genome_p, Fixed_sites, Fixed_sites_p, Fixed_sites_pair, Fixed_sites_pair_p):
    function_output = []
    
    #find n and s of gene:
    s = 0
    n = 0
    j = 0
    while(j < len(Reference_Genome)/3):
        start = (0 + j*3)
        end = (3 + j*3)
        n_and_s = find_n_and_s(Reference_Genome[start:end], Reference_Genome_p[j], genetic_code)
        n_out = int(n_and_s[0])#potential n sites for one codon
        s_out = int(n_and_s[1])#potential 2 sites for one codon
        n = n + n_out
        s = s + s_out
        j = j + 1
        
    #Find DNA differences in fixed sample:
    j = 0
    mut = 0#fixed mutant count
    pairwise_mut = 0#mutant overlap between haplotypes
    missense = 0#change proteins
    pairwise_mut_p = 0
    while(j < len(Reference_Genome)/3):
        start = (0 + j*3)
        end = (3 + j*3)
        
        #if the ancestral and input are equal, next
        if Fixed_sites[start:end] == Reference_Genome[start:end]:
            #print "no change"
            j = j + 1
            #skip: same as reference
        else:
            #if the input is different from the ancestor, but the same as the alternate, pairwise mutant tally + 1
            if Fixed_sites[start:end] == Fixed_sites_pair[start:end]:
                pairwise_mut = pairwise_mut + 1#haplotypes are the same codon
                #print "pairwise"
                #if the ancestral genome and input genome have different aa, pairwise missense tally + 1
                if Reference_Genome_p[j] == Fixed_sites_p[j]:
                    pairwise_mut_p  = pairwise_mut_p + 1
                    #print "misense"
                    #haplotypes are the same aa
                j = j + 1
            #if the haplotypes do NOT have the same codon:
            else:
                #Add 1 to the unique mutation counter
                mut = mut + 1 #haplotypes have different codons
            #check aa seq:
                #if different from ancestral aa:
                if Fixed_sites_p[j] == Reference_Genome_p[j]:
                    j = j + 1
                    #skip: synonymous
                else: 
                    #add 1 to the missense counter
                    missense = missense + 1
                    #haplotypes have different AAs
                    #print "misense"
                    j = j + 1
    
    #total codons in gene:
    total_codons = len(Reference_Genome)/3
    #total Dn and Ds
    Dn = missense#Hap specific misense counter
    Ds = mut - missense#Hap specific ALL mutations - missense mutations
    pairwise_Dn = pairwise_mut_p#Hap non-specific misense counter
    pairwise_Ds = pairwise_mut - pairwise_mut_p#Hap non-specific mutations - hap non-specific misense mutations
    
    function_output.append(str(total_codons))#0
    function_output.append(str(Dn))#1
    function_output.append(str(Ds))#2
    function_output.append(str(pairwise_Dn))#3
    function_output.append(str(pairwise_Ds))#4
    function_output.append(str(n))#5
    function_output.append(str(s))#6
    
    return(function_output)
# For each gene it returns a list:
# [codon number, missense SNPs between haplotypes, sense SNPs between haplotypes, missense SNPs across haplotypes, sense SNPs across haplotypes, number of possible novel missense SNPs in the gene, number of possible novel sense SNPs in the gene]




