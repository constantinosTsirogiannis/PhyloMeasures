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
  char *tree_filename=NULL, *sample_sizes_filename = NULL;
  Tree_type tree;
  int sample_size=0, k;
  Number_type chi;
  bool set_sample_size = false, 
       set_sample_sizes_file = false;
  bool received_chi = false,
       received_k = false;

  if( argc != 9 )
  {
    std::cout << std::endl << " Invalid number of arguments, should be eight ... aborting. " 
              << std::endl << std::endl;
    return -1;
  } 

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

      tree_filename = argv[i+1];
      i++;

    } // if(std::string(argv[i]) == std::string("-tree"))
    else if(std::string(argv[i]) == std::string("-sample_sizes_file"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_sizes_file flag should be followed by a path that points to" 
                  << std::endl << " a text file storing the samples sizes for the which the moments "
                  << std::endl << " of CAC shall be evaluated ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      set_sample_sizes_file = true;
      sample_sizes_filename = argv[i+1];
      i++;

    } // if(std::string(argv[i]) == std::string("-sample_sizes_file"))
    else if(std::string(argv[i]) == std::string("-sample_size"))
    {
      if( argc != 9 )
      {
        std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
        return -1;
      }

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_size flag should be followed by a positive integer that indicates the size" 
                  << std::endl << " of the samples among which the moment of the CAC is evaluated ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      std::string str_size = argv[i+1];

      for(int j=0; j<str_size.size(); j++)
        if( isdigit(str_size[j]) == false )
        {
          std::cout << std::endl << " -sample_size flag should be followed by a positive integer that indicates the size" 
                    << std::endl << " of the samples among which the moment of the CAC is evaluated ... aborting." 
                    << std::endl << std::endl;
          return -1;
        }

      sample_size = atoi(argv[i+1]);
      set_sample_size = true;
      i++;

    } // else if(std::string(argv[i]) == std::string("-sample_size"))
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
    else if(std::string(argv[i]) == std::string("-k"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -k flag should be followed by a positive integer ... aborting." 
        << std::endl << std::endl;
        return -1;
      }

      std::string str_k = argv[i+1]; 

      for(int j=0; j<str_k.size(); j++)
        if( isdigit(str_k[j]) == false )
        {
          std::cout << std::endl << " -k flag should be followed by a positive integer that indicates which" 
                    << std::endl << " moments of the CAC are to be evaluated ... aborting." 
                    << std::endl << std::endl;
          return -1;
        }

      k = Number_type(atoi(argv[i+1]));
      i++;

      received_k = true;

    } // else if(std::string(argv[i]) == std::string("-k"))
    else
    {
      std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
      return -1;
    }


  if(tree_filename == NULL)
  {
    std::cout << std::endl << " No input tree file specified ... aborting. " << std::endl << std::endl;
    return -1;
  }

  if( (set_sample_size == true && sample_size < 0) || 
      (set_sample_sizes_file == true && set_sample_size == true) ||
      (set_sample_sizes_file == false && set_sample_size == false) )
  {
    std::cout << std::endl << " Invalid values of sample size parameters." 
              << std::endl << std::endl;
    return -1;
  } 

  if(chi <= Number_type(0.5) || chi > Number_type(1.0) )
  {
    std::cout << std::endl 
              << " Invalid value of parameter chi. The value of chi must belong to the interval (0.5,1] ." 
              << std::endl << std::endl;
    return -1;
  } 

  if(k < 1)
  {
    std::cout << std::endl << " Invalid value of parameter k. Parameter k must be an integer greater than zero." 
              << std::endl << std::endl;
    return -1;
  }  
 
  std::cout << std::endl;
  std::cout << std::setprecision(5) << std::fixed;

  // Construct tree
  tree.construct_from_file(tree_filename);
  Core_ancestor_cost cac(tree,chi);


  if(set_sample_size == true)
  {
    std::vector<Number_type> moments;
  
    cac.compute_first_k_centralised_moments(k,sample_size,std::back_inserter(moments));

    for( int i=0; i<k; i++ )
      std::cout << moments[i] << " " ; 

    std::cout << std::endl;
  }
  else if(set_sample_sizes_file == true)
  {
    std::vector<int> sample_sizes;

    cac.read_sample_sizes_from_file(sample_sizes_filename, std::back_inserter(sample_sizes));

    for(int i=0; i<sample_sizes.size(); i++)
    {
      std::vector<Number_type> moments;
  
      cac.compute_first_k_centralised_moments(k,sample_sizes[i],std::back_inserter(moments));

      for( int j=0; j<k; j++ )
        std::cout << moments[j] << " " ; 

      std::cout << std::endl;
    }

  }

  std::cout << std::endl; 

  return 0;
  
} // main(...)
