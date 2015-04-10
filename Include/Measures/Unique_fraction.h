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

#ifndef UNIQUE_FRACTION_H
#define UNIQUE_FRACTION_H

#include<vector>

namespace PhylogeneticMeasures {

template< class KernelType >
struct Unique_fraction: public KernelType::Measure_base_bimodal
{
  typedef KernelType                               Kernel;
  typedef typename Kernel::Measure_base_bimodal    Base;
  typedef Unique_fraction<Kernel>                  Self;
  typedef typename Kernel::Number_type             Number_type;
  typedef typename Kernel::Numeric_traits          Numeric_traits;
  typedef typename Numeric_traits::Square_root     Square_root;
  typedef typename Kernel::Bimodal_tree            Bimodal_tree_type;
  typedef Bimodal_tree_type                        Tree_type;
  typedef typename Bimodal_tree_type::Node_type    Node_type;
  typedef typename Kernel::Edge_relation_type      Edge_relation_type;

  typedef typename Kernel::Phylogenetic_diversity  Phylogenetic_diversity;
  typedef typename Kernel::Common_branch_length    Common_branch_length;

  typedef typename Kernel::Exception_type          Exception_type;
  typedef typename Kernel::Exception_functor       Exception_functor; 


 public:

  Unique_fraction(Tree_type &tree)
  { p_tree = &tree; }


  template< class RangeIterator >
  Number_type operator()( RangeIterator rbegin_a, RangeIterator rend_a,
                          RangeIterator rbegin_b, RangeIterator rend_b,
                          int min_index_a, int max_index_a,
                          int min_index_b, int max_index_b );


  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_basic( const std::vector<std::string> &names_a,
                      const std::vector<std::vector<bool> > &matrix_a,
                      const std::vector<std::string> &names_b,
                      const std::vector<std::vector<bool> > &matrix_b, 
                      OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names_a, matrix_a, names_b, matrix_b, *this, false, ot); }

  template< class OutputIterator >
  int matrix_query_basic( const std::vector<std::string> &names,
                          const std::vector<std::vector<bool> > &matrix, OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names, matrix, *this, false, ot).first; }


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
  int matrix_query_specific_pairs_basic( const std::vector<std::string> &names,
                                         const std::vector<std::vector<bool> > &matrix,
                                         const std::vector<std::pair<int, int> > &queries, 
                                         OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names, matrix, queries, *this, false, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_basic( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename_a, matrix_filename_b, *this, false, ot); }


  template< class OutputIterator >
  int csv_matrix_query_basic( char *matrix_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename, *this, false, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_specific_pairs_basic( char *matrix_filename_a, char *matrix_filename_b, 
                                         char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs(*p_tree, matrix_filename_a, 
                                                          matrix_filename_b, queries_filename, 
                                                          *this, false, ot); }


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
  
  Number_type list_query(char* filename_a, char* filename_b)
  { return this->_list_query(*p_tree, filename_a, filename_b, *this); }


  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  Number_type compute_expectation( int sample_size_a, int sample_size_b )
  {  
    std::string exception_msg;
    exception_msg += " The computation of the expectation of UniFrac is not available.\n";     
    Exception_type excp;
    excp.get_error_message(exception_msg);
    Exception_functor excf;
    excf(excp);
    return Number_type(-1.0); 
  }

  Number_type compute_variance( int sample_size_a, int sample_size_b,
                                Number_type expect = Number_type(-1.0) )
  {
    std::string exception_msg;
    exception_msg += " The computation of the variance of UniFrac is not available.\n";     
    Exception_type excp;
    excp.get_error_message(exception_msg);
    Exception_functor excf;
    excf(excp);
    return Number_type(-1.0); 
  } 

  Number_type compute_deviation( int sample_size_a, int sample_size_b,
                                 Number_type expect = Number_type(-1.0) )
  {
    std::string exception_msg;
    exception_msg += " The computation of the deviation of UniFrac is not available.\n";     
    Exception_type excp;
    excp.get_error_message(exception_msg);
    Exception_functor excf;
    excf(excp); 
    return Number_type(-1.0); 
  }


 private:

  Bimodal_tree_type   *p_tree; // Stores a pointer to a bimodal phylogenetic tree object.

}; // struct Unique_fraction

} // PhylogeneneticMeasures

#include "Unique_fraction_impl.h"

#endif // UNIQUE_FRACTION_H
