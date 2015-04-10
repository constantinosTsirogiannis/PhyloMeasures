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
typedef Kernel::Mean_pairwise_distance     Mean_pairwise_distance;
typedef Mean_pairwise_distance::Tree_type  Tree_type;
   
int main(void)
{
  Tree_type tree;
  
  tree.construct_from_file("tree.tre");

  Mean_pairwise_distance mpd(tree);

  std::vector<Number_type> results;
  char matrix_filename[] = "example_matrix.csv";
  mpd.csv_matrix_query_basic(matrix_filename, std::back_inserter(results));

  std::cout << std::endl << " The MPD values are the following: " << std::endl;
  std::cout << "----------------------------------" << std::endl;

  for( size_t i=0; i<results.size(); i++)
    std::cout << results[i] << std::endl;

  std::cout << std::endl;

  return 0;  
}
