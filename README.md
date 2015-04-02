# PhyloMeasures
Fast Computations of Phylogenetic Biodiversity Measures

PhyloMeasures is an open source software package that provides fast computations of phylogenetic biodiversity measures. The package is released under the GNU Public License version 3, and is both available as an R package and as a C++ library. The current GitHub repository stores the code, programs, and documentation of the C++ library. 

PhyloMeasures provides functions for computing the value and the statistical moments of several phylogenetic biodiversity measures. The measures which are supported in the current version of the package are: the Phylogenetic Diversity (PD), the Mean Pairwise Distance (MPD), the Mean Nearest Taxon Distance (MNTD), the Core Ancestor Cost (CAC), the Common Branch Length (CBL), the Community Distance (CD), the Community Distance Nearest Taxon (CDNT), the Phylogenetic Sorensen's Similarity (PhyloSor), and the Unique Fraction (UniFrac).

The main advantages of PhyloMeasures are: 

    1) The package provides exact computation of the statistical moments and the standardised values for several phylogenetic biodiversity measures.

    2) All functions provided in PhyloMeasures are designed to be very efficient in practice, and they can process in a few seconds even trees that consist of a few hundred thousand tips.  The package functions are developed based on standard principles in Algorithms Design, hence their performance scales ideally as the size of the input trees increases.

