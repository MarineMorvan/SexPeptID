# SexPeptID

This open-sourced script enable sex identification from dental enamel based on untargeted paleoproteomics data (DDA).

# Typescript tutorial

1. Download, install and open MaxQuant (https://maxquant.org/)
2. In the first tab "Raw data":
   * "Load" your MS files
3. In the second tab "Group specific parameters":
   * In "Digestion" sub-tab, select the enzyme that you used during your samples preparation or select "unspecific" if you did not any enzyme
   * In "Modifications" sub-tab, select at least "oxidation (MP)" and "deamidation (NQ)" for "variable modifications", "fixed modifications" can be removed
4. In the third tab "Global parameters":
   * In "Sequences" sub-tab, add your FASTA file including at least AMELX and AMELY proteins, "Min. peptiden lenght" 6 (Default is 7),"Min. peptide lenght for unspecific search" 6 (Default is 8)
   * In "Indentification" sub-tab, "PSM FDR" 1 (Default is 0.01), "Protein FDR" 1 (Default is 0.01), "Min. score for unmodified peptides" 40 (Default is 0), "Dependent peptides" check box if desired, "De-novo sequencing" check box if desired
   * In "Protein quantification" sub-tab, add your "Modification used in protein quantification" or remove it if you do not use the quantification
5. On the bottom, adjust according to your computer capacity (Default is 4), and "Start"


6. Download and install both R and Rstudio (https://posit.co/download/rstudio-desktop/)
7. Open Rstudio
8. Copy and Paste the scritp from GitHub to Rstudio, and press "enter"
9. Import your MaxQuant output files: peptides.txt and msms.txt
10. In Peptide list tab, uncheck undesired peptides
11. In Graph and Table tab, add a label if desired
12. Click on "Download the Graph" and "Download the Table"
