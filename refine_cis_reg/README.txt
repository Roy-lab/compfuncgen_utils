
1. To compile, do 
g++ *.C -o refineCisReg

2. To refine the targets based on evolutionary conservation, here is an example for the cichlid species using one motif.
The code works for multiple motifs as well. Just make sure they are in the files listed in the motifcollection file. The logic of the program is
to show the instances of motifs at the leaf nodes only if the motif is conserved at the ancestor. 

Usage:
./refineCisReg speciesorder orthogroups spectiestree motifcollection scorefile outputprefix
speciesorder: this is just a list of species and species the order of the species in the ogid file (the next argument)

orthogroups: this is a file which specifies the orthology relationships across species. Each line has the format of the OGID with a duplication level. The subsequent column is a comma-separated list of genes that are believed to be orthologs of each other.

speciestree: specis the species tree one line per branch of the tree. The format of the line is
<child> left/right <parent>
left/right is arbitrary as this is a binary tree, just use it consistently

motifcollection: this is a file with a collection of motif files per species. The file has two columns, one for the motif file and the second for the name of the species. The motif file is the output of fimo.

scorefile: This file specifies the cost of gain, loss or maintaining the motif assignment. The format is like a 2X2 matrix. See example in the usage below.

outputprefix: The directory which will store the refined motifs.

Example usage:
./refineCisReg cichlidexample/cichlids_specorder.txt cichlidexample/cichlids_ogids.txt cichlidexample/cichlids_tree.txt cichlidexample/motifcollection2.txt score_gainloss.txt cichlidexample/refinedmotifs

