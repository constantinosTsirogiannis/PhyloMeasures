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
typedef Kernel::Core_ancestor_cost                       Core_ancestor_cost;
typedef Core_ancestor_cost::Tree_type                    Tree_type;

int main(int argc, char *argv[])
{
  char *filename1, *filename2;
  Number_type chi;
  std::vector<Number_type> res;
  Tree_type tree;
  bool received_chi = false;

  if(argc < 7 || argc > 8)
  {
    std::cout << std::endl << " Wrong number of arguments, should be either seven or eight ... aborting. " 
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
    else if(std::string(argv[i]) == std::string("-chi"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -chi flag should be followed by a real value " 
                  << std::endl << " that belongs to the interval (0.5,1] ... aborting." << std::endl << std::endl;
        return -1;
      }

      chi = Number_type(atof(argv[i+1]));
      i++;

      received_chi = true;

    } // else if(std::string(argv[i]) == std::string("-chi"))
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

  if(received_chi == false)
  {
    std::cout << std::endl << " No input chi value specified ... aborting. " << std::endl<< std::endl;
    return -1;
  }

  if(chi <= Number_type(0.5) || chi > Number_type(1.0) )
  {
    std::cout << std::endl << " Error: Invalid value of parameter chi." 
              << std::endl << " The value of chi must belong to the interval (0.5,1] ... aborting." 
              << std::endl << std::endl;
    return -1;
  } 

  std::cout << std::endl;

  // Construct tree
  tree.construct_from_file(filename1);
  Core_ancestor_cost cac(tree,chi);
  
  if(standardised == true)
    cac.csv_matrix_query_standardised(filename2, std::back_inserter(res));
  else
    cac.csv_matrix_query_basic(filename2, std::back_inserter(res));

  std::cout << std::setprecision(5) << std::fixed;

  for( int i=0; i<res.size(); i++ )
    std::cout << res[i] << std::endl;    

  std::cout << std::endl;  

  return 0;
  
} // main(int argc, char *argv[])

