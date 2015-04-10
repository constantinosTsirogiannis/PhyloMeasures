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

#ifndef COMMON_BRANCH_LENGTH_H
#define COMMON_BRANCH_LENGTH_H

#include<vector>


namespace PhylogeneticMeasures {

template< class KernelType >
struct Common_branch_length: public KernelType::Measure_base_bimodal
{
  typedef KernelType                             Kernel;
  typedef typename Kernel::Measure_base_bimodal  Base;
  typedef Common_branch_length<Kernel>           Self;
  typedef typename Kernel::Number_type           Number_type;
  typedef typename Kernel::Numeric_traits        Numeric_traits;
  typedef typename Numeric_traits::Square_root   Square_root;
  typedef typename Kernel::Bimodal_tree          Tree_type;
  typedef typename Tree_type::Node_type          Node_type;
  typedef typename Kernel::Edge_relation_type    Edge_relation_type;

  typedef typename Kernel::Exception_type        Exception_type;
  typedef typename Kernel::Exception_functor     Exception_functor; 

  Common_branch_length(Tree_type &tree)
  { p_tree = &tree;}

 private:

  template< class OutputIterator >
  void _compute_subtree_sums( int &index, Number_type &subtree_distances,
                              Number_type &h_products_a, Number_type &h_products_b,
                              Number_type &double_products, OutputIterator ot,
                              Number_type &sum_subtree, Number_type &sum_subtract );

  void _compute_subtree_sums(Number_type &sum_subtree, Number_type &sum_subtract)
  {
    Node_type root = p_tree->root();

    for( int i = 0; i < root.children.size(); i++ )
    {
      std::vector< std::pair<Number_type, int> > subtree_leaves;
      Number_type subtree_distances(0.0), h_products_a(0.0),
                  h_products_b(0.0), double_products(0.0);

      _compute_subtree_sums( root.children[i], subtree_distances, h_products_a,
                             h_products_b, double_products,
                             std::back_inserter(subtree_leaves), sum_subtree, sum_subtract );

      subtree_leaves.clear();
    }

  } // void _compute_subtree_sums( ... )


  public:

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
  { return this->_matrix_query_bimodal_specific_pairs( *p_tree, names_a, matrix_a, 
                                                       names_b, matrix_b, queries, *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_specific_pairs_basic( const std::vector<std::string> &names_a,
                                     const std::vector<std::vector<bool> > &matrix_a,
                                     const std::vector<std::string> &names_b,
                                     const std::vector<std::vector<bool> > &matrix_b, 
                                     const std::vector<std::pair<int, int> > &queries, 
                                     OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs( *p_tree, names_a, matrix_a, 
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
  { return this->_csv_matrix_query_bimodal_specific_pairs( *p_tree, matrix_filename_a, 
                                                           matrix_filename_b, queries_filename,
                                                           *this, true, ot); }

  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_specific_pairs_basic( char *matrix_filename_a, char *matrix_filename_b, 
                                         char *queries_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal_specific_pairs( *p_tree, matrix_filename_a, 
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
  
  Number_type list_query(char* filename_a, char* filename_b)
  { return this->_list_query(*p_tree, filename_a, filename_b, *this); }

  // Function that reads a list of pairs of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class OutputIterator >
  void read_sample_size_pairs_from_file(char *filename, OutputIterator ot)
  { this->_read_sample_size_pairs_from_file(filename, *p_tree, ot); }


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  void compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves,
                                                 std::vector<Number_type> &hypergeom, bool is_a );



  // Simulates f(x) = \binom{x}{a}/\binom{s}{a}, x \in [ _sample_size_a , _number_of_leaves ]

  Number_type hypergeom_a( int x )
  {
     if( x < _sample_size_a || x > _number_of_leaves )
       return Number_type(0.0);

     if( x == _number_of_leaves )
       return Number_type(1.0);

     return _hypergeom_a[x-_sample_size_a];
  }


  // Simulates f(x) = \binom{x}{b}/\binom{s}{b}, x \in [ _sample_size_b , _number_of_leaves ]

  Number_type hypergeom_b( int x )
  {
     if( x < _sample_size_b || x > _number_of_leaves )
       return Number_type(0.0);

     if( x == _number_of_leaves )
       return Number_type(1.0);

     return _hypergeom_b[x-_sample_size_b];
  }


  Number_type compute_expectation( int sample_size_a, int sample_size_b );

  Number_type two_edge_pr( int se, int sl , Edge_relation_type er, bool is_a );

  Number_type compute_variance( int sample_size_a, int sample_size_b);
  
  Number_type compute_variance_slow( int sample_size_a, int sample_size_b );
								
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

  Number_type compute_deviation_slow( int sample_size_a, int sample_size_b)
  { 
    Number_type variance = compute_variance_slow(sample_size_a, sample_size_b);   

    if( variance < Number_type(0.0) ) 
      return Number_type(0.0);

    return Square_root()(variance); 
  }

  private:

  Tree_type                *p_tree; // Stores a pointer to a phylogenetic tree object.
  std::vector<Number_type> _hypergeom_a,// Stores f(x) = \binom{x}{a}/\binom{s}{a}, 
                                                         // x \in [ _sample_size_a, _number_of_leaves ]

                           _hypergeom_b;// Stores f(x) = \binom{x}{b}/\binom{s}{b}, 
                                        // x \in [ _sample_size_b, _number_of_leaves ]

  int _sample_size_a, _sample_size_b;
  int _number_of_leaves;

}; // struct Commmon_branch_length

} // PhylogeneneticMeasures

#include "Common_branch_length_impl.h"

#endif // COMMON_BRANCH_LENGTH_H
