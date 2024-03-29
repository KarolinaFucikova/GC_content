# GC_content
Evolution of GC content in Green algae

The R project contains the data on a selection of Chlorophyte green algae: the table (csv) contains information on the min and max temperatures at the sites of the algae's origins, their GC content in chloroplast genomes and in nuclear ribosomal genes. The tree file contains a Bayesian consensus tree based on chloroplast protein-coding genes, onto which the GC content and temperature data can be mapped.
The R script allows for drawing the tree, mapping the temperature and GC data on it (after dropping tips with no data for that variable), and also graphs correlations among selected variables using ggplot.
Further, the code also calculates phylogenetically independent contrasts, allowing for a phylogeny-corrected regression between temperature and GC content. 
Several plotting options and color schemes are explored and can be recreated using the script. Example plots are saved in the plots folder. 
