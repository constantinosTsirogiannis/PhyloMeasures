//////////////////////////////////////////////////////////////////////////////////
//    Copyright (C) 2015,  Constantinos Tsirogiannis.  Email: analekta@gmail.com
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
//////////////////////////////////////////////////////////////////////////////////

#include<fstream>
#include<iomanip>
#include<Phylogenetic_measures_kernel.h>

typedef PhylogeneticMeasures::Numeric_traits_double      Numeric_traits;
typedef Phylogenetic_measures_kernel<Numeric_traits>     Kernel;
typedef Kernel::Number_type                              Number_type;
typedef Kernel::Mean_nearest_taxon_distance              Mean_nearest_taxon_distance;
typedef Mean_nearest_taxon_distance::Tree_type           Tree_type;

int main(int argc, char *argv[])
{
  char *filename1, *filename2;
  std::vector<double> res;
  Tree_type tree;

  if(argc < 5 || argc > 6)
  {
    std::cout << std::endl << " Wrong number of arguments, should be either four or five ... aborting. " 
              << std::endl << std::endl;
    return -1;
  } 

  filename1 = NULL;

  bool standardised = false;

  for(int i=1; i<argc; i++)
    if(std::string(argv[i]) == std::string("-tree"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -tree flag should be followed by a path that points to" 
                  << std::endl << " a text file representing a tree in Newick format ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      filename1 = argv[i+1];
      i++;
    }
    else if(std::string(argv[i]) == std::string("-csv"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -csv flag should be followed by one or two paths " 
                  << std::endl << " that point to .csv matrix files ... aborting." << std::endl << std::endl;
        return -1;
      }

      filename2 = argv[i+1];
      i++;

    } // else if(std::string(argv[i]) == std::string("-csv"))
    else if(std::string(argv[i]) == std::string("-standardised"))
      standardised = true;
    else
    {
      std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
      return -1;
    }

  if(filename1 == NULL)
  {
    std::cout << std::endl << " No input tree file specified ... aborting. " << std::endl << std::endl;
    return -1;
  }

  if(filename2 == NULL)
  {
    std::cout << std::endl << " No input .csv matrix specified ... aborting. " << std::endl<< std::endl;
    return -1;
  }

  std::cout << std::setprecision(5)<< std::fixed;

  std::cout << std::endl;

  // Construct tree 
  tree.construct_from_file(filename1);
  Mean_nearest_taxon_distance mntd(tree);
  
  if(standardised == true)
  {
    if(!tree.is_ultrametric())
    {
      std::cout << " Unfortunately, the current release of the package does not" << std::endl
                << " support the computation of any statistical moments for "<< std::endl
                << " the Mean Nearest Taxon Distance when the input tree " << std::endl
                << " is not ultrametric ... aborting. " << std::endl<< std::endl;
      return -1;
    }

    mntd.csv_matrix_query_standardised(filename2, std::back_inserter(res));
  }
  else
    mntd.csv_matrix_query_basic(filename2, std::back_inserter(res));

  for( int i=0; i<res.size(); i++ )
    std::cout << res[i] << std::endl;    

  std::cout << std::endl;  

  return 0;
  
} // main(int argc, char *argv[])

