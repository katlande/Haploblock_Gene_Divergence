# Haploblock_Gene_Divergence
Scripts and files produced to analyze SNP divergence between haplotype-specific gene variants of Helianthus haploblocks


### Files

**All_Gene_Variants_Dn.Ds.txt:** All haplotype-specific gene variants, listed with sense and nonsense mutation counts, Ka/Ks for SNPs with AF>0.9, TAIR homolog descriptions, and Ka/Ks p-values.

**All_Gene_Variants_Grantham_Distances.txt:** Grantham distances between haplotype-specific predicted protein products for all haploblocks, listed with missense and stop codon counts, and TAIR homolog descriptions.

**All_GO_Terms_With_Grantham_Sums.txt:** Grantham distances of individual GO terms for each haploblock.

**All_Gene_Variants_MK_Test_Results.txt:** Mcdonald-Kreitman test results for for all hapltype-specific gene variants with at least 4 polymorphic and 4 fixed SNPs.

### Scripts & Supporting Files

**Find_Grantham_Distances.py:** For finding the grantham distance between two aligned, no-gap peptide fastas.
*Find_Grantham_Distances.py requires Grantham_Matrix.txt to run*

**dN_dS.py:** Functions used to find dN/dS between haploblocks, using an ancestral lineage as the background. Excludes fixed differences that are common across haplotypes, as these likely occured before lineage divergence. Requires 6 inputs per gene: DNA & peptide FASTAs for ancestor and both haplotypes. Fastas must be aligned and contain no gaps.

