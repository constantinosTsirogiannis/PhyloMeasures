////////////////////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2015,  Constantinos Tsirogiannis and Brody Sandel.
//
//    Email: analekta@gmail.com and brody.sandel@bios.au.dk
//
//    This file is part of PhyloMeasures.
//
//    PhyloMeasures is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PhyloMeasures is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with PhyloMeasures.  If not, see <http://www.gnu.org/licenses/>
////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "Phylogenetic_measures_kernel.h"

typedef Phylogenetic_measures_kernel<>     Kernel;
typedef Kernel::Number_type                Number_type;
typedef Kernel::Phylogenetic_diversity     Phylogenetic_diversity;
typedef Phylogenetic_diversity::Tree_type  Tree_type;
   
int main(void)
{
  Tree_type tree;
  
  tree.construct_from_file("tree.tre");

  Phylogenetic_diversity pd(tree);

  int sample_size = 5;
  Number_type mean = pd.compute_expectation(sample_size);

  std::cout << std::endl << " The mean of the PD for sample size 5 is: " << mean 
            << std::endl << std::endl;

  return 0;  
}
