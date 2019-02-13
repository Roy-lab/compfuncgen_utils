This code is the simplest implementation of gain and loss rates of binding sites. There is only one pair of parameters. The same code samples data and learns. We are currently hard-coding this.

To compile the code type make based on the included Makefile information. The gsl directory contains an installation of the GNU Scientific Library, and forms an important include to compile this code, as indicated in the Makefile. 

Usage:
./estimateRate data/blen_hog1spec.txt data/inputmapping.txt data/species_prob_osr6spec.txt samples_1_3.txt 2 3
