#**************************************************#
#   MultiMSOAR 2.0: An Accurate Tool to Identify   #
#      Ortholog Groups among Multiple Genomes      #
#                                                  #
#**************************************************# 

-------------------
Usage
-------------------
 
 The software MultiMSOAR 2.0 is an executable file
 included in the current directory. It requires the 
 following files as input:

 Input:
 -----------
		speciesTree	 -   binary species tree in Newick format
    geneFamily   -   gene families for all input genomes
    Si_Sj        -   ortholog pairs and their similarity 
										 scores between genomes Si and Sj

 Run MultiMSOAR2.0
 -------------------
    Move all files Si_Sj to the same directory as the
    the executable file MultiMSOAR2.0, then execute the
    following command:

    ./MultiMSOAR2.0 <#species> <speciesTree> <geneFamily> 
                    <-o GeneInfo> <-o OrthoGroup>

    The <#species> is an integer, representing the number
    of species in comparison. Suppose there are n species, 
    then we need n*(n-1)/2 files (Si_Sj, 0<=i<j<n), which
    contains the orthology information between any pair of
    the n genomes.

 Output:
 -------------
		GeneInfo    -   the file contains information about 
                    gene births, losses and duplications
                    for every input genomes
    OrthoGroup  -   the file contains all ortholog groups 
                    identified by MultiMSOAR 2.0


-------------------
File Format
-------------------

 Si_Sj (0<=i<j<n)
 -----------------
      Si_gene			Sj_gene		SimilarityScore
 Ex:  S0_gene1    S1_gene2     105


 speciesTree
 -----------------
		We use the Newick format to represent a species tree
 Ex: ((S0,S1),(S2,(S3,S4)));


 geneFamily
 -----------------
		Each line contains the gene names from all species 
    within a gene family, separated by '\t'
 Ex: S0_gene1 S1_gene1 S2_gene1
     S1_gene2 S2_gene2


 GeneInfo
 -----------------
		Report the new born gene and duplicated genes in 
    current genomes, and the number of lost genes in
    every genome in the species tree


 OrthoGroup
 -----------------
    Each line reports genes within the same ortholog
    group


--------------------
Sample Input
--------------------
		This directory also contains sample input files:
		
		Si_Sj(0<=i<j<5)
		sampleSpeciesTree
		sampleGeneFamily

		With these sample files, you can run MultiMSOAR2.0
		as follows:

		./MultiMSOAR2.0	5 sampleSpeciesTree sampleGeneFamily
                      sampleGeneInfo sampleOrthoGroup
		
		You will get the results saved in the files 
		sampleGeneInfo and sampleOrthoGroup 

