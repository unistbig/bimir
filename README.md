# bimir
Biclustering analysis of transcriptome big data identifies condition-specific miRNA targets

# Introduction
‘bimir’ is an R package developed to generate biclusters of cell condition-specific microRNA targets from large mRNA transcriptome fold change (FC) matrix (a.k.a., FC table; it contains log2 FC values of 20,639 human genes under 5,158 experimental conditions). Besides the primary purpose of this package, it can be applied to create any constant biclusters allowing small number of noises using the progressive bicluster extension (PBE) algorithm. This package is also downloadable  from the BIMIR web site (http://btool.org/bimir_dir/).

# Installation
##### install.packages('devtools') (skip this step if it is already installed.)
##### library('devtools')
##### install_github('unistbig/bimir')
##### library('bimir')

#### It will take some minutes to install the package because of the large fold change table (~147MB). Using my PC, it took ~16 minutes. If there are problems with installing this package, please leave a message on this page or send me an email (yoonsora@unist.ac.kr)

# Quick start with an example
Let’s generate biclusters of ‘hsa-miR-1-3p’ where its target genes are up-regulated under common experimental conditions. To do this, type following lines.
##### FC = load_FCtable() # To load FC table
##### MP = getMIRprofile(miRNA = 'hsa-miR-1-3p', FCtable = FC, FCcutoff = log2(1.3)) # To generate MIR profile.
##### PBE_MERGE(MIR_profile = MP, mir = 'hsa-miR-1-3p', biclust.path = './', FCcutoff = log2(1.3))

- Then check the directory assigned to ‘biclust.path’ and see the biclusters (biclust_up_*.txt) and corresponding experimental condition (Experimental_condition_up_*.txt) and gene list (Targetlist_up_*.txt).
- Down-bicluster: To generate biclusters of targets commonly down-regulated under multiple conditions, just modify the ‘FCcutoff’ parameter to be negative value (e.g., -log2(1.3)).
- The available miRNAs are listed in the ‘getmiRNAlist()’ function.
