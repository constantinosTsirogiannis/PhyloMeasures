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
typedef Kernel::Unique_fraction                          Unique_fraction;
typedef Unique_fraction::Tree_type                       Tree_type;

int main(int argc, char *argv[])
{
  char *filename1, *filename2, *filename3, *filename4;
  std::pair<int,int> rows_cols;
  std::vector<double> res;
  Tree_type tree;

  if(argc < 5 || argc > 8)
  {
    std::cout << std::endl << " Wrong number of arguments, should be from four to seven ... aborting. " 
              << std::endl << std::endl;
    return -1;
  } 

  filename1 = NULL;

  bool specific_pairs = false,
       single_matrix = false;

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

      if(i+1 < argc && std::string(argv[i+1])[0] != '-' )
      {
        filename3 = argv[i+1];
        i++;
      }
      else
        single_matrix = true;

    } // else if(std::string(argv[i]) == std::string("-csv"))
    else if(std::string(argv[i]) == std::string("-specific_pairs"))
    {
      specific_pairs = true;

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -specific_pairs flag should be followed by a path that points" 
                  << std::endl << " to a file that describes the query pairs ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      filename4 = argv[i+1];
      i++;

    } // else if(std::string(argv[i]) == std::string("-specific_pairs"))
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
  Unique_fraction unifrac(tree);
  
  if(single_matrix == true && specific_pairs == true)
  {
    int rows = unifrac.csv_matrix_query_specific_pairs_basic(filename2, filename4, std::back_inserter(res));
    rows_cols = std::make_pair(rows,rows);     
  }
  else if(specific_pairs == true)
    rows_cols = unifrac.csv_matrix_query_specific_pairs_basic(filename2, filename3, 
                                                              filename4, std::back_inserter(res));
  else if(single_matrix == true)
  {
    int rows = unifrac.csv_matrix_query_basic(filename2, std::back_inserter(res));
        rows_cols = std::make_pair(rows,rows);
  }
  else 
    rows_cols = unifrac.csv_matrix_query_basic(filename2, filename3, std::back_inserter(res));


  if(specific_pairs == true)
  {
    for( int i=0; i<res.size(); i++ )
        std::cout << res[i] << "  "; 

    std::cout << std::endl;
  }   
  else
    for( int i=0; i<rows_cols.first; i++ )
    {
      for( int j=0; j<rows_cols.second; j++ )
        std::cout << res[i*rows_cols.second +j] << "  "; 

      std::cout << std::endl; 
    }

  std::cout << std::endl;  

  return 0;
  
} // main(int argc, char *argv[])

