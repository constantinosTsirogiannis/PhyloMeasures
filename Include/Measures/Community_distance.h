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

#ifndef COMMUNITY_DISTANCE_H
#define COMMUNITY_DISTANCE_H

#include<vector>

namespace PhylogeneticMeasures {


template< class KernelType >
struct Community_distance: public KernelType::Measure_base_bimodal,
                           public KernelType::template Mean_pairwise_distance_base<typename KernelType::Bimodal_tree>
{
  typedef KernelType                             Kernel;
  typedef typename Kernel::Measure_base_bimodal  Base;
  typedef Community_distance<Kernel>             Self;
  typedef typename Kernel::Number_type           Number_type;
  typedef typename Kernel::Numeric_traits        Numeric_traits;
  typedef typename Numeric_traits::Square_root   Square_root;
  typedef typename Kernel::Bimodal_tree          Tree_type;
  typedef typename Tree_type::Node_type          Node_type;

  typedef typename Kernel::Exception_type        Exception_type;
  typedef typename Kernel::Exception_functor     Exception_functor; 

 public:

  Community_distance(Tree_type &tree):_expectation(-1.0)
  { p_tree = &tree; }

  template< class RangeIterator >
  Number_type operator()( RangeIterator rbegin_a, RangeIterator rend_a,
                          RangeIterator rbegin_b, RangeIterator rend_b,
                          int min_index_a, int max_index_a,
                          int min_index_b, int max_index_b );

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_standardised( const std::vector<std::string> &names_a,
                             const std::vector<std::vector<bool> > &matrix_a,
                             const std::vector<std::string> &names_b,
                             const std::vector<std::vector<bool> > &matrix_b, 
                             OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names_a, matrix_a, names_b, matrix_b, *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_basic( const std::vector<std::string> &names_a,
                      const std::vector<std::vector<bool> > &matrix_a,
                      const std::vector<std::string> &names_b,
                      const std::vector<std::vector<bool> > &matrix_b, 
                      OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names_a, matrix_a, names_b, matrix_b, *this, false, ot); }

  template< class OutputIterator >
  int matrix_query_standardised( const std::vector<std::string> &names,
                                 const std::vector<std::vector<bool> > &matrix, OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names, matrix, *this, true, ot).first; }

  template< class OutputIterator >
  int matrix_query_basic( const std::vector<std::string> &names,
                          const std::vector<std::vector<bool> > &matrix, OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names, matrix, *this, false, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_specific_pairs_standardised( const std::vector<std::string> &names_a,
                                            const std::vector<std::vector<bool> > &matrix_a,
                                            const std::vector<std::string> &names_b,
                                            const std::vector<std::vector<bool> > &matrix_b, 
                                            const std::vector<std::pair<int, int> > &queries, 
                                            OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names_a, matrix_a, 
                                                      names_b, matrix_b, queries, *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_specific_pairs_basic( const std::vector<std::string> &names_a,
                                     const std::vector<std::vector<bool> > &matrix_a,
                                     const std::vector<std::string> &names_b,
                                     const std::vector<std::vector<bool> > &matrix_b, 
                                     const std::vector<std::pair<int, int> > &queries, 
                                     OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names_a, matrix_a, 
                                                      names_b, matrix_b, queries, *this, false, ot); }

  template< class OutputIterator >
  int matrix_query_specific_pairs_standardised( const std::vector<std::string> &names,
                                                const std::vector<std::vector<bool> > &matrix,
                                                const std::vector<std::pair<int, int> > &queries, 
                                                OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names, matrix, queries, *this, true, ot).first; }

  template< class OutputIterator >
  int matrix_query_specific_pairs_basic( const std::vector<std::string> &names,
                                         const std::vector<std::vector<bool> > &matrix,
                                         const std::vector<std::pair<int, int> > &queries, 
                                         OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names, matrix, queries, *this, false, ot).first; }
  
  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_standardised( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename_a, matrix_filename_b, *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_basic( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename_a, matrix_filename_b, *this, false, ot); }

  template< class OutputIterator >
  int csv_matrix_query_standardised( char *matrix_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename, *this, true, ot).first; }

  template< class OutputIterator >
  int csv_matrix_query_basic( char *matrix_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename, *this, false, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_specific_pairs_standardised( char *matrix_filename_a, char *matrix_filename_b, 
                                                char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs(*p_tree, matrix_filename_a, 
                                                          matrix_filename_b, queries_filename,
                                                          *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_specific_pairs_basic( char *matrix_filename_a, char *matrix_filename_b, 
                                         char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs(*p_tree, matrix_filename_a, 
                                                          matrix_filename_b, queries_filename, 
                                                          *this, false, ot); }

  template< class OutputIterator >
  int csv_matrix_query_specific_pairs_standardised( char *matrix_filename, char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs(*p_tree, matrix_filename, 
                                                          queries_filename, *this, true, ot).first; }

  template< class OutputIterator >
  int csv_matrix_query_specific_pairs_basic( char *matrix_filename, char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs(*p_tree, matrix_filename, 
                                                          queries_filename, *this, false, ot).first; }
  
  
  // Input:  Two ranges of iterators that indicate two lists of species names (in std::string format).
  // Output: The value of the current measure for these two sets of species.
  
  template < class RangeIterator >
  Number_type list_query( RangeIterator rbegin_a, RangeIterator rend_a,
                          RangeIterator rbegin_b, RangeIterator rend_b )
  { return this->_list_query(*p_tree, rbegin_a, rend_a, rbegin_b, rend_b, *this); }
  
  
  // Input: A tree and two txt files that each stores a list of species names, each constituting a subset 
  // of the species in the tree and which appear in random order.
  // Output: The value of the current measure for these two sets of species.
  
  Number_type list_query(char* filename1, char* filename2)
  { return this->_list_query(*p_tree, filename1, filename2, *this); }


  // Function that reads a list of pairs of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class OutputIterator >
  void read_sample_size_pairs_from_file(char *filename, OutputIterator ot)
  { this->_read_sample_size_pairs_from_file(filename, *p_tree, ot); }


  Number_type compute_expectation( int sample_size_a, int sample_size_b )
  {
    if( sample_size_a < 0 || sample_size_b < 0 ||
        sample_size_a > p_tree->number_of_leaves() || sample_size_b > p_tree->number_of_leaves()  )
    {
      std::string exception_msg;
      exception_msg += " Request to compute expectation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size_a == 0 || sample_size_b == 0 )
      return Number_type(0.0);

    if(_expectation != Number_type(-1.0) )
      return _expectation;

   int s(p_tree->number_of_leaves());

   _expectation = Number_type(2.0)*(this->total_path_costs(*p_tree))/(Number_type(s)*Number_type(s));

    return _expectation;
  } // csv_matrix_distance_query( char *filename, OutputIterator ot)

  Number_type compute_variance( int sample_size_a, int sample_size_b );

  Number_type compute_deviation( int sample_size_a, int sample_size_b )
  { 
    if( sample_size_a < 0 || sample_size_b < 0 ||
        sample_size_a > p_tree->number_of_leaves() || sample_size_b > p_tree->number_of_leaves()  )
    {
      std::string exception_msg;
      exception_msg += " Request to compute deviation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    Number_type variance = compute_variance(sample_size_a, sample_size_b);   

    if( variance < Number_type(0.0) ) 
      return Number_type(0.0);

    return Square_root()(variance); 
  }

  private:

   Tree_type *p_tree; // Stores a pointer to a phylogenetic tree object.
   Number_type _expectation;

}; // struct Community_distance

} // PhylogeneneticMeasures

#include "Community_distance_impl.h"

#endif // COMMUNITY_DISTANCE_H
