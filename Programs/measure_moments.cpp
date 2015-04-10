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

typedef PhylogeneticMeasures::Numeric_traits_double                Numeric_traits;
typedef Phylogenetic_measures_kernel<Numeric_traits>               Kernel;
typedef Kernel::Number_type                                        Number_type;

typedef Kernel::Phylogenetic_diversity            Phylogenetic_diversity;
typedef Kernel::Mean_pairwise_distance            Mean_pairwise_distance;
typedef Kernel::Mean_nearest_taxon_distance       Mean_nearest_taxon_distance;
typedef Kernel::Core_ancestor_cost                Core_ancestor_cost;
typedef Kernel::Community_distance                Community_distance;
typedef Kernel::Common_branch_length              Common_branch_length;

typedef Kernel::Unimodal_tree                     Unimodal_tree;
typedef Kernel::Bimodal_tree                      Bimodal_tree;
typedef Kernel::Mean_nearest_taxon_distance_tree  Mean_nearest_taxon_distance_tree;


int main(int argc, char *argv[])
{
  char *tree_filename=NULL, *sample_sizes_filename, 
       *measure=NULL, *moment=NULL;
  Unimodal_tree unimode_tree;
  Bimodal_tree bimode_tree;
  Mean_nearest_taxon_distance_tree  mntd_tree;
  int sample_size, sample_size_a, sample_size_b;
  bool set_sample_size = false,
       set_sample_size_a = false,
       set_sample_size_b = false,
       set_sample_sizes_file = false,
       set_moment_type = false;

  if( argc != 8 && argc != 10 )
  {
    std::cout << std::endl << " Invalid number of arguments, should be either seven or nine ... aborting. " 
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
      if( argc != 8 )
      {
        std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
        return -1;
      }

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_sizes_file flag should be followed by a path that points to" 
                  << std::endl << " a text file storing the samples sizes for the which the moment "
                  << std::endl << " of a measure shall be evaluated ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      set_sample_sizes_file = true;
      sample_sizes_filename = argv[i+1];
      i++;

    } // if(std::string(argv[i]) == std::string("-sample_sizes_file")) 
    else if(std::string(argv[i]) == std::string("-sample_size"))
    {
      if( argc != 8 )
      {
        std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
        return -1;
      }

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_size flag should be followed by a positive integer that indicates the size" 
                  << std::endl << " of the samples among which the moment of a measure is evaluated ... aborting." 
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
    else if(std::string(argv[i]) == std::string("-sample_size_a"))
    {
      if( argc != 10 )
      {
        std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
        return -1;
      }

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_size_a flag should be followed by a positive integer that indicates the " 
                  << std::endl << " size of the A samples among which the moment of a measure is evaluated ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      std::string str_size = argv[i+1];

      for(int j=0; j<str_size.size(); j++)
        if( isdigit(str_size[j]) == false )
        {
          std::cout << std::endl << " -sample_size_a flag should be followed by a positive integer that indicates the " 
                    << std::endl << " size of the A samples for which the moment of a measure is evaluated ... aborting."  
                    << std::endl << std::endl;
          return -1;
        }

      sample_size_a = atoi(argv[i+1]);
      set_sample_size_a = true;
      i++;

    } // else if(std::string(argv[i]) == std::string("-sample_size_a"))
    else if(std::string(argv[i]) == std::string("-sample_size_b"))
    {
      if( argc != 10)
      {
        std::cout << std::endl << " Wrong input arguments ... aborting. " << std::endl << std::endl;
        return -1;
      }

      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -sample_size_b flag should be followed by a positive integer that indicates the size" 
                  << std::endl << " of the B samples among which the moment of a measure is evaluated ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      std::string str_size = argv[i+1];

      for(int j=0; j<str_size.size(); j++)
        if( isdigit(str_size[j]) == false )
        {
          std::cout << std::endl << " -sample_size_b flag should be followed by a positive integer that indicates the " 
                    << std::endl << " size of the B samples for which the moment of a measure is evaluated ... aborting."  
                    << std::endl << std::endl;
          return -1;
        }

      sample_size_b = atoi(argv[i+1]);
      set_sample_size_b = true;
      i++;

    } // else if(std::string(argv[i]) == std::string("-sample_size_b"))
    else if(std::string(argv[i]) == std::string("-measure"))
    {
      if(i+1 >= argc || std::string(argv[i+1])[0] == '-')
      {
        std::cout << std::endl << " -measure flag should be followed by a string that indicates" 
                  << std::endl << " which phylogenetic measure shall be computed ... aborting." 
                  << std::endl << std::endl;
        return -1;
      }

      measure = argv[i+1];
      i++;

    } // else if(std::string(argv[i]) == std::string("-measure"))
    else if( std::string(argv[i]) == std::string("-expectation") ||
             std::string(argv[i]) == std::string("-deviation") )
    {
      if(set_moment_type == true)
      {
        std::cout << std::endl << " Can only compute only one moment type at a time ... aborting. " 
                  << std::endl << std::endl;
        return -1;
      }

      moment = argv[i];
  
      set_moment_type = true;
    }
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

  if(measure == NULL)
  {
    std::cout << std::endl << " No measure specified ... aborting. " << std::endl << std::endl;
    return -1;
  }

  if(moment == NULL)
  {
    std::cout << std::endl << " No statistical moment specified ... aborting. " << std::endl << std::endl;
    return -1;
  }

  if( std::string(measure) != std::string("PD") &&
      std::string(measure) != std::string("MPD") && std::string(measure) != std::string("MNTD") && 
      std::string(measure) != std::string("CD") && std::string(measure) != std::string("CBL") )
  {
    std::cout << std::endl << " Invalid measure specified." << std::endl;
    std::cout << " (Should be either \"-MPD\", or \"-PD\"," 
              << " or \"-MNTD\", or \"-CD\", or \"-CBL\" ) " << std::endl << std::endl;
    return -1;
  }

  if( (std::string(measure) == std::string("PD") && 
           ( (set_sample_sizes_file == false && set_sample_size == false) || 
             set_sample_size_a == true || set_sample_size_b == true)) ||
      (std::string(measure) == std::string("MPD")&& 
           ( (set_sample_sizes_file == false && set_sample_size == false) ||
              set_sample_size_a == true || set_sample_size_b == true)) ||
      (std::string(measure) == std::string("MNTD")&&
           ( (set_sample_sizes_file == false && set_sample_size == false) || 
             set_sample_size_a == true || set_sample_size_b == true)) ||
      (std::string(measure) == std::string("CBL") && set_sample_size == true ) ||
      (std::string(measure) == std::string("CD") && set_sample_size == true ) ||
      (set_sample_size == true && sample_size < 0 ) ||
      (set_sample_size_a == true && sample_size_a < 0) ||
      (set_sample_size_b == true && sample_size_b < 0 )  ||
      (set_sample_sizes_file == true && (set_sample_size == true || 
                                         set_sample_size_a == true || 
                                         set_sample_size_b == true ) ) ||
      (set_sample_sizes_file == false && set_sample_size == false && 
       set_sample_size_a == false && set_sample_size_b == false ) ||
      (set_sample_size_a == true && set_sample_size_b == false) || 
      (set_sample_size_a == false && set_sample_size_b == true) )
  {
    std::cout << std::endl << " Invalid sample size ranges specified ... aborting. " << std::endl << std::endl;
    return -1;
  }

  if( std::string(moment) != std::string("-expectation") && std::string(moment) != std::string("-deviation") )
  {
    std::cout << std::endl << " Invalid moment specified." << std::endl;
    std::cout << std::endl << " (Should be either \"-expectation\" or \"-deviation\" ) " << std::endl << std::endl;
    return -1;
  }

  std::cout << std::setprecision(5) << std::fixed;

  if( std::string(measure) == "PD"  || std::string(measure) == "MPD" )
  {
    std::cout << std::endl;
  
    // Construct tree
    unimode_tree.construct_from_file(tree_filename);
  
    Mean_pairwise_distance  mpd(unimode_tree);
    Phylogenetic_diversity  pd(unimode_tree);

    if(set_sample_size)
    {
      if( std::string(measure) == "MPD" && std::string(moment) == "-expectation" )
        std::cout << Number_type(mpd.compute_expectation(sample_size)) << std::endl;
      else if( std::string(measure) == "PD" && std::string(moment) == "-expectation" )
        std::cout << Number_type(pd.compute_expectation(sample_size)) << std::endl;
      else if( std::string(measure) == "MPD" && std::string(moment) == "-deviation" )
        std::cout << Number_type(mpd.compute_deviation(sample_size)) << std::endl;
      else if( std::string(measure) == "PD" && std::string(moment) == "-deviation" )
        std::cout << Number_type(pd.compute_deviation(sample_size)) << std::endl;
    }
    else if(set_sample_sizes_file)
    {
      std::vector<int> sample_sizes;
 
      pd.read_sample_sizes_from_file(sample_sizes_filename, std::back_inserter(sample_sizes));

      for( int i=0; i<sample_sizes.size(); i++ )
      {
        if( std::string(measure) == "MPD" && std::string(moment) == "-expectation" )
          std::cout << Number_type(mpd.compute_expectation(sample_sizes[i])) << std::endl;
        else if( std::string(measure) == "PD" && std::string(moment) == "-expectation" )
          std::cout << Number_type(pd.compute_expectation(sample_sizes[i])) << std::endl;
        else if( std::string(measure) == "MPD" && std::string(moment) == "-deviation" )
          std::cout << Number_type(mpd.compute_deviation(sample_sizes[i])) << std::endl;
        else if( std::string(measure) == "PD" && std::string(moment) == "-deviation" )
          std::cout << Number_type(pd.compute_deviation(sample_sizes[i])) << std::endl;
      }
    }


  } // if( std::string(measure) == "PD" || std::string(measure) == "MPD" )


  if( std::string(measure) == "MNTD" )
  {  
    std::cout << std::endl;

    // Construct tree
    mntd_tree.construct_from_file(tree_filename);

    if(!mntd_tree.is_ultrametric())
    {
      std::cout << std::endl << " Unfortunately, the current release of the package does not" << std::endl
                             << " support the computation of any statistical moments for "<< std::endl
                             << " the Mean Nearest Taxon Distance when the input tree " << std::endl
                             << " is not ultrametric ... aborting. " << std::endl<< std::endl;
      return -1;
    }

    Mean_nearest_taxon_distance mntd(mntd_tree);

    if(set_sample_size)
    {
      if( std::string(moment) == "-expectation" )
        std::cout << Number_type(mntd.compute_expectation(sample_size)) << std::endl;
      else if( std::string(moment) == "-deviation" )
        std::cout << Number_type(mntd.compute_deviation(sample_size)) << std::endl;
    }
    else if(set_sample_sizes_file)
    {
      std::vector<int> sample_sizes;
      mntd.read_sample_sizes_from_file(sample_sizes_filename,std::back_inserter(sample_sizes));

      for(int i=0; i<sample_sizes.size() ; i++)
        if( std::string(moment) == "-expectation" )
          std::cout << Number_type(mntd.compute_expectation(sample_sizes[i])) << std::endl;
        else if( std::string(moment) == "-deviation" )
          std::cout << Number_type(mntd.compute_deviation(sample_sizes[i])) << std::endl;
    }


  } // if( std::string(measure) == "MNTD" )

  if( std::string(measure) == "CD" || std::string(measure) == "CBL" )
  {
    std::cout << std::endl;

    bimode_tree.construct_from_file(tree_filename);
    Community_distance cd(bimode_tree);
    Common_branch_length cbl(bimode_tree);

    if(set_sample_size_a == true)
    {
        if( std::string(measure) == "CD" && std::string(moment) == "-expectation" )
          std::cout << cd.compute_expectation(sample_size_a, sample_size_b) << std::endl;
        else if( std::string(measure) == "CBL" && std::string(moment) == "-expectation" )
          std::cout << cbl.compute_expectation(sample_size_a, sample_size_b) << std::endl;
        else if( std::string(measure) == "CD" && std::string(moment) == "-deviation" )
          std::cout << cd.compute_deviation(sample_size_a, sample_size_b) << std::endl;
        else if( std::string(measure) == "CBL" && std::string(moment) == "-deviation" )
          std::cout << cbl.compute_deviation(sample_size_a, sample_size_b) << std::endl;
    }
    else if(set_sample_sizes_file == true)
    {
      std::vector< std::pair<int,int> > sample_sizes;
      cd.read_sample_size_pairs_from_file(sample_sizes_filename,std::back_inserter(sample_sizes));

      for( int i=0; i<sample_sizes.size(); i++ )
        if( std::string(measure) == "CD" && std::string(moment) == "-expectation" )
          std::cout << cd.compute_expectation(sample_sizes[i].first, sample_sizes[i].second) << std::endl;
        else if( std::string(measure) == "CBL" && std::string(moment) == "-expectation" )
          std::cout << cbl.compute_expectation(sample_sizes[i].first, sample_sizes[i].second) << std::endl;
        else if( std::string(measure) == "CD" && std::string(moment) == "-deviation" )
          std::cout << cd.compute_deviation(sample_sizes[i].first, sample_sizes[i].second) << std::endl;
        else if( std::string(measure) == "CBL" && std::string(moment) == "-deviation" )
          std::cout << cbl.compute_deviation(sample_sizes[i].first, sample_sizes[i].second) << std::endl;
    }

  } // if( std::string(measure) == "CD" || std::string(measure) == "CBL" )

  std::cout << std::endl;  

  return 0;
  
} // main(...)
