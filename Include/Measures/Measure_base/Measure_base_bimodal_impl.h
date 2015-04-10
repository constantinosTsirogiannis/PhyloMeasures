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

#ifndef MEASURE_BASE_BIMODAL_IMPL_H
#define MEASURE_BASE_BIMODAL_IMPL_H

#include<vector>
#include<cmath>
#include<cctype>

  template < class KernelType >
  std::pair<int,int> PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _get_pair_from_string(std::string &str)
  {
    std::pair<int,int> pr;

    size_t res = str.find_first_of('/');

    if(res == std::string::npos || res == str.size()-1 )
    {
      std::string exception_msg(" There was a mistake in the syntax of the sample sizes file.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string str_a = str.substr(0,res),
                str_b = str.substr(res+1);

    if(str_a.size() == 0 || str_b.size() == 0)
    {
      std::string exception_msg(" There was a mistake in the syntax of the sample sizes file.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    for(int i=0; i<str_a.size(); i++)
      if(isdigit(str_a[i]) == false)
      {
        std::string exception_msg;
        exception_msg += " There is an error in the syntax of the sample sizes file.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    for(int i=0; i<str_b.size(); i++)
      if(isdigit(str_b[i]) == false)
      {
        std::string exception_msg;
        exception_msg += " There is an error in the syntax of the sample sizes file.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    pr.first = atoi(str_a.c_str());
    pr.second = atoi(str_b.c_str());
 
    return pr;
  }


  template < class KernelType >
  template < class TreeType, class OutputIterator >
  void PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _read_sample_size_pairs_from_file(char *filename, TreeType &tree, OutputIterator ot)
  {
    std::ifstream in(filename);
    std::vector<std::pair<int,int> > sample_sizes;

    // Reading first file with queries
    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem while opening the file with the sample sizes.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string line;
    std::getline(in,line);
 
    int prev_index=-1, current_index =0;

    while( current_index< line.size()-1 )
    {
      do
      {
        current_index++;
      }
      while(line[current_index]!=',' && current_index < line.size() );

      if(current_index -prev_index < 4)
      {
        std::string exception_msg(" There is a mistake in the syntax of the sample sizes file.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      std::string substring = line.substr(prev_index+1,current_index -prev_index-1);

      std::pair<int,int> sizes = _get_pair_from_string(substring);

      if( sizes.first < 0 || sizes.first > tree.number_of_leaves() ||
          sizes.second < 0 || sizes.second > tree.number_of_leaves()    )
      {
        std::string exception_msg(" One of the sample sizes in the file has an invalid value.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      sample_sizes.push_back(sizes);

      prev_index = current_index;

    } // while( current_index< line.size()-1 )

    for(int i=0; i<sample_sizes.size(); i++)
      *ot++ = sample_sizes[i];

  } // _read_sample_size_pairs_from_file(...)


  template < class KernelType >
  template< class TreeType, class OutputIterator1, class OutputIterator2 >
  void
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _extract_samples_from_matrix( TreeType &tree, 
                                const std::vector<std::string> &names,
                                const std::vector<std::vector<bool> > &matrix,
                                OutputIterator1 ot1, OutputIterator2 ot2 )
  {
    typedef TreeType                              Tree_type;
    typedef typename Tree_type::Leaves_iterator   Leaves_iterator;

    std::vector<int>   column_to_node_vec;
  
    if(names.size() < tree.number_of_leaves())
    {
      std::string warning(" Warning: one of the input matrices has fewer columns than the number of species in the tree.");
      Exception_functor().issue_warning(warning);
    }

    std::vector<bool> checked_names;
    checked_names.assign(tree.number_of_nodes(),false);

    for( int i=0; i<names.size(); i++ )
    {
      Leaves_iterator lv_it = tree.find_leaf(names[i]);

      if( lv_it == tree.leaves_end() )
      {
        std::string exception_msg;
        exception_msg += " One of the species names in input the matrix was not found in the tree (";
        exception_msg += names[i];
        exception_msg += ") \n";  
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

      if(checked_names[(*lv_it).second] == true)
      {
        std::string exception_msg;
        exception_msg += " Two or more columns of the input matrix share the same species name (";
        exception_msg += (*lv_it).first;
        exception_msg += ") \n";   
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 
      else
        checked_names[(*lv_it).second] = true;

      column_to_node_vec.push_back((*lv_it).second);
    }



    // Read the matrix, extracting a query per line.

    for(int i=0; i<matrix.size(); i++)
    {
      std::vector<int> query_nodes;
      int min = tree.number_of_nodes(), max = -1;

      for( int j=0; j<matrix[i].size(); j++)
        if( matrix[i][j] == true )
        {
          query_nodes.push_back( column_to_node_vec[j] );

          if( query_nodes.back() < min )
            min = query_nodes.back();

          if( query_nodes.back() > max )
            max = query_nodes.back();
        }

      *ot1++ = query_nodes;
      *ot2++ = std::make_pair(min,max);

    } // for(int i=0; i<matrix.size(); i++)

  } // _extract_samples_from_matrix( ... )


  // TODO: This is the second function that reads a csv matrix,
  // (given also the one that is used by unimodal measures).
  // It should be better to unify these functions in some way.

  template < class KernelType >
  template< class TreeType, class OutputIterator1, class OutputIterator2 >
  void
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _extract_samples_from_file( TreeType &tree, char *filename,
                              OutputIterator1 ot1, OutputIterator2 ot2 )
  {
    typedef TreeType                              Tree_type;
    typedef typename Tree_type::Leaves_iterator   Leaves_iterator;
  
    std::ifstream in(filename);

    // Reading first file with queries
    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg;
      exception_msg += " There was a problem with opening the file with the matrix.\n"; 
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    // Read the first row, the one that contains the species names

    std::vector<std::string> names;
    std::vector<int>   column_to_node_vec;
    std::string line;
    char a;

    std::getline(in,line);

    int c=0;

    while(c < line.size() )
    {
      std::string str;
      a = line[c];

      while(a != ',' && a != ' ' && a != '\r' && c<line.size() )
      {
        str.push_back(a);
        c++;
        a = line[c];
      }

      if( str.size() > 0 )
        names.push_back(str);

      c++;
    }


    if(names.size() < tree.number_of_leaves())
    {
      std::string warning(" Warning: one of the input matrices has fewer columns than the number of species in the tree.");
      Exception_functor().issue_warning(warning);
    }

    std::vector<bool> checked_names;
    checked_names.assign(tree.number_of_nodes(),false);

    for( int i=0; i<names.size(); i++ )
    {
      Leaves_iterator lv_it = tree.find_leaf(names[i]);

      if( lv_it == tree.leaves_end() )
      {
        std::string exception_msg;
        exception_msg += " One of the species names in the input matrix was not found in the tree (";
        exception_msg += names[i];
        exception_msg += ") \n";
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

      if(checked_names[(*lv_it).second] == true)
      {
        std::string exception_msg;
        exception_msg += " Two or more columns of the input matrix share the same species name (";
        exception_msg += (*lv_it).first;
        exception_msg += ") \n";   
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 
      else
        checked_names[(*lv_it).second] = true;

      column_to_node_vec.push_back((*lv_it).second);
    }

    // Read the rest of the matrix, extracting a query per line.

    while( in.good() )
    {
      line.clear();

      int count=0, start=0;
      std::vector<int> query_nodes;
      int min = tree.number_of_nodes(), max = -1;

      std::getline(in,line);

      // Exclude the first word of the line if first character is not zero or one.

      if( line[start] != '0' && line[start] != '1' )
        while( start < line.size() && line[start] != ',' )
          start++;

      for( int i=start; i<line.size(); i++)
      {
        a = line[i];

        if( a == '1' )
        {
          if(count >=names.size())
          {
            std::string exception_msg;
            exception_msg += " One of the matrix files has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }  

          query_nodes.push_back( column_to_node_vec[count] );

          if( query_nodes.back() < min )
            min = query_nodes.back();

          if( query_nodes.back() > max )
            max = query_nodes.back();

          count++;
        }
        else if ( a == '0' )
        {
          if(count >=names.size())
          {
            std::string exception_msg;
            exception_msg += " One of the matrix files has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }

          count++;
        }
        else if( a != ',' && a !=' ' && a !='\r' && a !='\n')
        {
          std::string exception_msg;
          exception_msg += " One of the matrix files has wrong syntax.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      } // for( int i=start; i<line.size(); i++)

      if( in.good() )
      {
        *ot1++ = query_nodes;
        *ot2++ = std::make_pair(min,max);
      }

    } // while( in.good() )

    in.close();

  } // _extract_samples_from_file( ... )


  template < class KernelType >
  void PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _extract_query_pairs( std::string &str,std::pair<int,int> &pair_a, std::pair<int,int> &pair_b)
  {
    size_t res = str.find_first_of('/');

    if(res == std::string::npos || res == str.size()-1 )
    {
      std::string exception_msg;
      exception_msg += " There is a mistake in the syntax of the sample pairs file.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string str_a = str.substr(0,res),
                str_b = str.substr(res+1);

    if(str_a.size() == 0 || str_b.size() == 0)
    {
      std::string exception_msg;
      exception_msg += " There is a mistake in the syntax of the sample pairs file.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    size_t res_a = str_a.find_first_of('-'),
           res_b = str_b.find_first_of('-');

    if(res_a == std::string::npos)
    {
      for(int i=0; i<str_a.size(); i++)
        if(isdigit(str_a[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      pair_a.first = atoi(str_a.c_str());
      pair_a.second = pair_a.first;
    } 
    else
    {
      std::string str_a_1 = str_a.substr(0,res_a),
                  str_a_2 = str_a.substr(res_a+1);

      for(int i=0; i<str_a_1.size(); i++)
        if(isdigit(str_a_1[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }


      for(int i=0; i<str_a_2.size(); i++)
        if(isdigit(str_a_2[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      pair_a.first = atoi(str_a_1.c_str());
      pair_a.second = atoi(str_a_2.c_str());
    }     

    if(res_b == std::string::npos)
    {
      for(int i=0; i<str_b.size(); i++)
        if(isdigit(str_b[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      pair_b.first = atoi(str_b.c_str());
      pair_b.second = pair_b.first;
    } 
    else
    {
      std::string str_b_1 = str_b.substr(0,res_b),
                  str_b_2 = str_b.substr(res_b+1);

      for(int i=0; i<str_b_1.size(); i++)
        if(isdigit(str_b_1[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }


      for(int i=0; i<str_b_2.size(); i++)
        if(isdigit(str_b_2[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample pairs file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      pair_b.first = atoi(str_b_1.c_str());
      pair_b.second = atoi(str_b_2.c_str());
    }   

    pair_a.first = pair_a.first-1;
    pair_a.second = pair_a.second-1;

    pair_b.first = pair_b.first-1;
    pair_b.second = pair_b.second-1;

  } // _extract_query_pairs(...)


  // The following function reads a file that indicated which pairs of
  // samples between two csv matrices should be used by as input
  // for computing a bimodal measure. The input file should be a txt file
  // that contains a string of the following format:

  // 3-6/7-10,5/2-4,7/8,65-71/12

  // The string in the above example indicates that we want to compute 
  // the value of a two-sample measure between all pairs of samples 
  // such that the first element of the pair is any of the samples
  // from line 3 up to line 6 in the first csv matrix and the second
  // element of the pair is any of the samples that appear from 
  // line 7 to line 10 in the second matrix, then all pairs where the first
  // element is the sample in line 5 of the first matrix and the second
  // element is any of the samples of lines from 2 to 4 in the second matrix etc.
  // Conceptually, the first line of each matrix is the line with index 0, the
  // second is the one with index 1, and so on and so forth.   

  template < class KernelType >
  template< class OutputIterator1, class OutputIterator2 >
  void PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _extract_queries_from_file( char *filename, OutputIterator1 ot_a, OutputIterator2 ot_b )
  {
     std::ifstream in(filename);

    // Reading first file with queries
    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg;
      exception_msg += " There was a problem with opening the file with the sample pairs.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string line;
    std::getline(in,line);
 
    int prev_index=-1, current_index =0;

    // Check if string contains unwanted characters

    for(int i=0; i<line.size(); i++)
      if( isdigit(line[i]) == false && line[i] != ' ' && line[i] != '-' && 
          line[i] != '/' && line[i] != ',' && line[i] != '\r' && line[i] != '\n' )
      {
        std::string exception_msg;
        exception_msg += " The sample pairs file contains unexpected characters.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }
         
    while( current_index< line.size()-1 )
    {
      do
      {
        current_index++;
      }
      while(line[current_index]!=',' && current_index < line.size() );

      if(current_index -prev_index-1 < 3)
      {
        std::string exception_msg;
        exception_msg += " There is a mistake in the syntax of the sample pairs file.\n";     
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      std::string substring = line.substr(prev_index+1,current_index -prev_index-1);

      std::pair<int,int> pair_a, pair_b;
   
      _extract_query_pairs(substring, pair_a, pair_b);

      *ot_a++ = pair_a;
      *ot_b++ = pair_b;  

      prev_index = current_index;

    } // while( current_index< line.size()-1 )

  } // _extract_queries_from_file(...)

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _matrix_query_internal_bimodal( TreeType &tree,
                                  std::vector< std::vector<int> >  &samples_a,
                                  std::vector< std::pair<int,int> > &min_max_a,
                                  std::vector< std::vector<int> >  &samples_b,
                                  std::vector< std::pair<int,int> > &min_max_b,
                                  bool is_double_matrix, 
                                  Measure &msr, bool standardised, OutputIterator ot)
  {
    if(is_double_matrix)
    {
      if(standardised == false)
      {
        for( int i=0; i<samples_a.size(); i++ )
          for( int j=0; j<samples_b.size(); j++ )
            *ot++ = msr(samples_a[i].begin(), samples_a[i].end(),
                        samples_b[j].begin(), samples_b[j].end(),
                        min_max_a[i].first, min_max_a[i].second,
                        min_max_b[j].first, min_max_b[j].second  );
      }
      else
      {
        for( int i=0; i<samples_a.size(); i++ )
          for( int j=0; j<samples_b.size(); j++ )
          {
            Number_type res  = msr(samples_a[i].begin(), samples_a[i].end(),
                               samples_b[j].begin(), samples_b[j].end(),
                               min_max_a[i].first, min_max_a[i].second,
                               min_max_b[j].first, min_max_b[j].second  ),
                        mean = msr.compute_expectation(samples_a[i].size(), samples_b[j].size() ),
                        deviation = msr.compute_deviation(samples_a[i].size(), samples_b[j].size() );

            if(deviation != Number_type(0.0))
              *ot++ = (res-mean)/deviation;
            else
              *ot++ = res-mean;

          } // for( int j=0; j<samples_b.size(); j++ )

      } // else of if(standardised == false)

      return std::make_pair(samples_a.size(), samples_b.size());
    } // if(is_double_matrix)
    else
    {
      if(standardised == false)
      {
        for( int i=0; i<samples_a.size(); i++ )
          for( int j=0; j<samples_a.size(); j++ )
            *ot++ = msr(samples_a[i].begin(), samples_a[i].end(),
                        samples_a[j].begin(), samples_a[j].end(),
                        min_max_a[i].first, min_max_a[i].second,
                        min_max_a[j].first, min_max_a[j].second  );
      }
      else
      {
        for( int i=0; i<samples_a.size(); i++ )
          for( int j=0; j<samples_a.size(); j++ )
          {
            Number_type res = msr(samples_a[i].begin(), samples_a[i].end(),
                                  samples_a[j].begin(), samples_a[j].end(),
                                  min_max_a[i].first, min_max_a[i].second,
                                  min_max_a[j].first, min_max_a[j].second  ),
                        mean = msr.compute_expectation(samples_a[i].size(), samples_a[j].size() ),
                        deviation = msr.compute_deviation(samples_a[i].size(), samples_a[j].size() );

            if(deviation != Number_type(0.0))
              *ot++ = (res-mean)/deviation;
            else
              *ot++ = res-mean;

          } // for( int j=0; j<samples_a.size(); j++ )
      }
 
      return std::make_pair(samples_a.size(), samples_a.size());

    }	// else of if(is_double_matrix)

  } // _matrix_query_internal_bimodal(...)

  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _matrix_query_internal_bimodal_specific_pairs( TreeType &tree,  
                                                 std::vector< std::vector<int> >  &samples_a,
                                                 std::vector< std::pair<int,int> > &min_max_a,
                                                 std::vector< std::vector<int> >  &samples_b,
                                                 std::vector< std::pair<int,int> > &min_max_b,
                                                 std::vector< std::pair<int,int> > &query_intervals_a,
                                                 std::vector< std::pair<int,int> > &query_intervals_b,
                                                 bool is_double_matrix, 
                                                 Measure &msr, bool standardised, OutputIterator ot )
  {
    for(int i=0; i<query_intervals_a.size(); i++)
      if( query_intervals_a[i].first < 0 || query_intervals_a[i].second < 0 ||
          query_intervals_a[i].first >= samples_a.size() || query_intervals_a[i].second >= samples_a.size() )
      {
        std::string exception_msg;
        exception_msg += " An input query pair is out of range.\n";
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    if(is_double_matrix)
    {
      for(int i=0; i<query_intervals_b.size(); i++)
        if( query_intervals_b[i].first < 0 || query_intervals_b[i].second < 0 ||
            query_intervals_b[i].first >= samples_b.size() || query_intervals_b[i].second >= samples_b.size() )
        {
          std::string exception_msg;
          exception_msg += " An input query pair is out of range.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }
    }
    else
      for(int i=0; i<query_intervals_b.size(); i++)
        if( query_intervals_b[i].first < 0 || query_intervals_b[i].second < 0 ||
            query_intervals_b[i].first >= samples_a.size() || query_intervals_b[i].second >= samples_a.size() )
        {
          std::string exception_msg;
          exception_msg += " An input query pair is out of range.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

    if(is_double_matrix)
    {
      if(standardised == false)
      {
        for( int i=0; i<query_intervals_a.size(); i++ )
          for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
            for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
              *ot++ = msr(samples_a[j].begin(), samples_a[j].end(),
                          samples_b[k].begin(), samples_b[k].end(),
                          min_max_a[j].first, min_max_a[j].second,
                          min_max_b[k].first, min_max_b[k].second  );
      }
      else
      {
        for( int i=0; i<query_intervals_a.size(); i++ )
          for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
            for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
            {
              Number_type res = msr(samples_a[j].begin(), samples_a[j].end(),
                                    samples_b[k].begin(), samples_b[k].end(),
                                    min_max_a[j].first, min_max_a[j].second,
                                    min_max_b[k].first, min_max_b[k].second  ),
                          mean = msr.compute_expectation(samples_a[j].size(), samples_b[k].size() ),
                          deviation = msr.compute_deviation(samples_a[j].size(), samples_b[k].size() );

              if(deviation != Number_type(0.0))
                *ot++ = (res-mean)/deviation;
              else
                *ot++ = res-mean;
                                
            } // for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )


      } // else of if(standardised == false)

      return std::make_pair(samples_a.size(), samples_b.size());
    } // if(is_double_matrix)
    else
    {
      if(standardised == false)
      {
        for( int i=0; i<query_intervals_a.size(); i++ )
          for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
            for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
              *ot++ = msr(samples_a[j].begin(), samples_a[j].end(),
                          samples_a[k].begin(), samples_a[k].end(),
                          min_max_a[j].first, min_max_a[j].second,
                          min_max_a[k].first, min_max_a[k].second  );
      }
      else
      {
        for( int i=0; i<query_intervals_a.size(); i++ )
          for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
            for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
            {
              Number_type res = msr(samples_a[j].begin(), samples_a[j].end(),
                                    samples_a[k].begin(), samples_a[k].end(),
                                    min_max_a[j].first, min_max_a[j].second,
                                    min_max_a[k].first, min_max_a[k].second  ),
                          mean = msr.compute_expectation(samples_a[j].size(), samples_a[k].size() ),
                          deviation = msr.compute_deviation(samples_a[j].size(), samples_a[k].size() );

              if(deviation != Number_type(0.0))
                *ot++ = (res-mean)/deviation;
              else
                *ot++ = res-mean;
            }

      } // else of if(standardised == false)

      return std::make_pair(samples_a.size(), samples_a.size());

    } // else of if(is_double_matrix)

  } // _matrix_query_internal_bimodal_specific_pairs(...)

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _matrix_query_bimodal( TreeType &tree, const std::vector<std::string> &names_a,
                         const std::vector<std::vector<bool> > &matrix_a,
                         const std::vector<std::string> &names_b,
                         const std::vector<std::vector<bool> > &matrix_b, 
                         Measure &msr, bool standardised, OutputIterator ot)
  {  
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (&matrix_b != &matrix_a);

    _extract_samples_from_matrix(tree, names_a, matrix_a, std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      _extract_samples_from_matrix( tree, names_b, matrix_b, 
                                    std::back_inserter(samples_b), std::back_inserter(min_max_b) );

    return _matrix_query_internal_bimodal( tree, samples_a, min_max_a,
                                           samples_b, min_max_b, is_double_matrix, 
                                           msr, standardised, ot);

  } // _matrix_query_bimodal(...)

  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _matrix_query_bimodal_specific_pairs( TreeType &tree,  
                                        const std::vector<std::string> &names_a,
                                        const std::vector<std::vector<bool> > &matrix_a,
                                        const std::vector<std::string> &names_b,
                                        const std::vector<std::vector<bool> > &matrix_b, 
                                        const std::vector<std::pair<int, int> > &queries,
                                        Measure &msr, bool standardised, OutputIterator ot )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (&matrix_a != &matrix_b);

    _extract_samples_from_matrix(tree, names_a, matrix_a, std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      _extract_samples_from_matrix(tree, names_b, matrix_b, 
                                   std::back_inserter(samples_b), 
                                   std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    for(int i=0; i<queries.size(); i++)
    {
      query_intervals_a.push_back(std::make_pair(queries[i].first,queries[i].first));
      query_intervals_b.push_back(std::make_pair(queries[i].second,queries[i].second));
    }  


    return _matrix_query_internal_bimodal_specific_pairs(tree, samples_a, min_max_a,
                                                         samples_b, min_max_b, query_intervals_a,
                                                         query_intervals_b, is_double_matrix, 
                                                         msr, standardised, ot );

  } // _matrix_query_bimodal_specific_pairs(...)


  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _csv_matrix_query_bimodal( TreeType &tree, char *filename1, 
                             char *filename2, Measure &msr, bool standardised, OutputIterator ot)
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(filename1) != std::string(filename2));

    _extract_samples_from_file(tree, filename1, std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      _extract_samples_from_file(tree, filename2, std::back_inserter(samples_b), std::back_inserter(min_max_b) );

    return _matrix_query_internal_bimodal( tree, samples_a, min_max_a,
                                           samples_b, min_max_b, is_double_matrix, 
                                           msr, standardised, ot);

  } // _csv_matrix_query_bimodal( TreeType &tree, char *filename1, 
    //                            char *filename2, Measure &msr, OutputIterator ot)


  template < class KernelType >
  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _csv_matrix_query_bimodal_specific_pairs( TreeType &tree,  char *matrix_filename_a, 
                                            char *matrix_filename_b, char *queries_filename,
                                            Measure &msr, bool standardised, OutputIterator ot )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(matrix_filename_a) != std::string(matrix_filename_b));

    _extract_samples_from_file(tree, matrix_filename_a, std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      _extract_samples_from_file(tree, matrix_filename_b, std::back_inserter(samples_b), std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    _extract_queries_from_file( queries_filename, 
                                std::back_inserter(query_intervals_a), 
                                std::back_inserter(query_intervals_b) );

    if(query_intervals_a.size() != query_intervals_b.size())
    {
      std::string exception_msg;
      exception_msg += " The number of query intervals for the first csv matrix is\n";
      exception_msg += " different than the number of intervals specified for the second matrix.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    return _matrix_query_internal_bimodal_specific_pairs(tree, samples_a, min_max_a,
                                                         samples_b, min_max_b, query_intervals_a,
                                                         query_intervals_b, is_double_matrix, 
                                                         msr, standardised, ot);

  } // csv_matrix_query_bimodal_specific_pairs( TreeType &tree, char *filename1, 
    //                                          char *filename2, Measure &msr, OutputIterator ot)
  

  // Input:  Two ranges of iterators that indicate two lists of leaf species names (in std::string format).
  // Output: The value of the current measure for these two sets of species.
  
   template < class KernelType >  
   template < class TreeType, class RangeIterator, class Measure >
   typename KernelType::Number_type 
   PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
   _list_query( TreeType &tree, 
                RangeIterator rbegin_a, RangeIterator rend_a,
                RangeIterator rbegin_b, RangeIterator rend_b,
                Measure& msr )
  {
    typedef typename TreeType::Leaves_iterator  Leaves_iterator;
  
    RangeIterator it;
    std::string str;
    std::vector<int> leaf_indices_a, leaf_indices_b;
    int min_index_a = tree.number_of_nodes(), max_index_a = -1, 
	    min_index_b = tree.number_of_nodes(), max_index_b = -1;

    // Find the indices of all the leaf nodes that correspond to the query species of sample A.
    // Find also the minimum and maximum of those indices (here represented as min_index_a & max_index_a).
    for( it = rbegin_a; it != rend_a; it++ )
    {
      str = *it;	  
      Leaves_iterator lv_it = tree.find_leaf(str);

      if(  lv_it != tree.leaves_end() )
      {
         leaf_indices_a.push_back(lv_it->second);

         if( leaf_indices_a.back() < min_index_a )
           min_index_a = leaf_indices_a.back();

         if( leaf_indices_a.back() > max_index_a )
           max_index_a = leaf_indices_a.back();
      }

    } // for( it = rbegin_a; it != rend_a; it++ )

    // Find the indices of all the leaf nodes that correspond to the query species of sample B.
    // Find also the minimum and maximum of those indices (here represented as min_index_b & max_index_b).
    for( it = rbegin_b; it != rend_b; it++ )
    {
      str = *it;
      Leaves_iterator lv_it = tree.find_leaf(str);

      if( lv_it != tree.leaves_end() )
      {
         leaf_indices_b.push_back(lv_it->second);

         if( leaf_indices_b.back() < min_index_b )
           min_index_b = leaf_indices_b.back();

         if( leaf_indices_b.back() > max_index_b )
           max_index_b = leaf_indices_b.back();
      }

    } // for( it = rbegin_b; it != rend_b; it++ )

   // if( leaf_indices_a.size() + leaf_indices_b.size()  < 2 ||
   //     leaf_indices_a.size() == 0 || leaf_indices_b.size() == 0 )
   //   return Number_type(0.0);

    return msr( leaf_indices_a.begin(), leaf_indices_a.end(),
                leaf_indices_b.begin(), leaf_indices_b.end(),
                min_index_a, max_index_a, min_index_b, max_index_b );
				   
  } //    list_query( TreeType &tree,
    //                RangeIterator rbegin_a, RangeIterator rend_a,
    //                RangeIterator rbegin_b, RangeIterator rend_b )

	
  // Input: A tree and two txt files that each stores a list of species names, each constituting a subset 
  // of the leaf species in the tree and which appear in random order.
  // Output: The value of the current measure for these two sets of species.
  
  template < class KernelType > 
  template < class TreeType, class Measure>  
  typename KernelType::Number_type 
  PhylogeneticMeasures::Measure_base_bimodal<KernelType>::
  _list_query(TreeType &tree, char* filename1, char* filename2, Measure& msr)
  {
    std::vector<std::string> vec1, vec2;

    // Process first file

    std::ifstream in1(filename1);

    if( !( in1.is_open() && in1.good() ) )
    {
      std::string exception_msg;
      exception_msg += " There was a problem with opening the first file with the species names.\n";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string name;
    char a;

    in1 >> a;

    while(in1.good())
    {
      if( a == ',')
      {
        vec1.push_back(name);
        name.clear();
      }
      else if( a != '\n' && a!='\0' )
        name.push_back(a);

      in1 >> a;
    }

    vec1.push_back(name);

    // Process second file

    name.clear();

    std::ifstream in2(filename2);

    if( !( in2.is_open() && in2.good() ) )
    {
      std::string exception_msg;
      exception_msg += " There was a problem with opening the second file with the species names.\n";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    in2 >> a;

    while(in2.good())
    {
      if( a == ',')
      {
        vec2.push_back(name);
        name.clear();
      }
      else if( a != '\n' && a!='\0' )
        name.push_back(a);

      in2 >> a;
    }

    vec2.push_back(name); 

    in1.close();
    in2.close();

    return _list_query(tree, vec1.begin(), vec1.end(), vec2.begin(), vec2.end(), msr);

  }	// list_query(TreeType &tree, char* filename1, char* filename2)
  	
#endif //MEASURE_BASE_BIMODAL_IMPL_H
