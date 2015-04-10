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

#ifndef MEASURE_BASE_BIMODAL_H
#define MEASURE_BASE_BIMODAL_H

#include<vector>
#include<cmath>

namespace PhylogeneticMeasures {

template< class KernelType >
struct Measure_base_bimodal
{
  typedef KernelType                          Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Exception_type     Exception_type;
  typedef typename Kernel::Exception_functor  Exception_functor; 

  Measure_base_bimodal(){}

protected:

  std::pair<int,int> _get_pair_from_string(std::string &str);

  // Function that reads a list of pairs of integers from a file,
  // assumed to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given tree. 
  template < class TreeType, class OutputIterator >
  void _read_sample_size_pairs_from_file(char *filename, TreeType &tree, OutputIterator ot);

  template< class TreeType, class OutputIterator1, class OutputIterator2 >
  void _extract_samples_from_matrix( TreeType &tree, 
                                     const std::vector<std::string> &names,
                                     const std::vector<std::vector<bool> > &matrix,
                                     OutputIterator1 ot1, OutputIterator2 ot2 );

  template< class TreeType, class OutputIterator1, class OutputIterator2 >
  void 
  _extract_samples_from_file( TreeType &tree, char *filename, 
                              OutputIterator1 ot1, OutputIterator2 ot2 );

  void _extract_query_pairs( std::string &str,std::pair<int,int> &pair_a, std::pair<int,int> &pair_b);

  template< class OutputIterator1, class OutputIterator2 >
  void _extract_queries_from_file( char *filename, OutputIterator1 ot_a, OutputIterator2 ot_b );

  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_internal_bimodal( TreeType &tree,
                                  std::vector< std::vector<int> >  &samples_a,
                                  std::vector< std::pair<int,int> > &min_max_a,
                                  std::vector< std::vector<int> >  &samples_b,
                                  std::vector< std::pair<int,int> > &min_max_b,
                                  bool is_double_matrix, 
                                  Measure &msr, bool standardised, OutputIterator ot);


  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_internal_bimodal_specific_pairs( TreeType &tree,  
                                                 std::vector< std::vector<int> >  &samples_a,
                                                 std::vector< std::pair<int,int> > &min_max_a,
                                                 std::vector< std::vector<int> >  &samples_b,
                                                 std::vector< std::pair<int,int> > &min_max_b,
                                                 std::vector< std::pair<int,int> > &query_intervals_a,
                                                 std::vector< std::pair<int,int> > &query_intervals_b,
                                                 bool is_double_matrix, 
                                                 Measure &msr, bool standardised, OutputIterator ot );


  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_bimodal( TreeType &tree, const std::vector<std::string> &names_a,
                         const std::vector<std::vector<bool> > &matrix_a,
                         const std::vector<std::string> &names_b,
                         const std::vector<std::vector<bool> > &matrix_b, 
                         Measure &msr, bool standardised, OutputIterator ot);


  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_bimodal( TreeType &tree, const std::vector<std::string> &names,
                         const std::vector<std::vector<bool> > &matrix,
                         Measure &msr, bool standardised, OutputIterator ot)
  { return _matrix_query_bimodal(tree,names,matrix,names,matrix,msr,standardised,ot);}

  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_bimodal_specific_pairs( TreeType &tree,  
                                        const std::vector<std::string> &names_a,
                                        const std::vector<std::vector<bool> > &matrix_a,
                                        const std::vector<std::string> &names_b,
                                        const std::vector<std::vector<bool> > &matrix_b, 
                                        const std::vector<std::pair<int, int> > &queries,
                                        Measure &msr, bool standardised, OutputIterator ot );


  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _matrix_query_bimodal_specific_pairs( TreeType &tree,  
                                        const std::vector<std::string> &names,
                                        const std::vector<std::vector<bool> > &matrix,
                                        const std::vector<std::pair<int, int> > &queries,
                                        Measure &msr, bool standardised, OutputIterator ot )
  {  return _matrix_query_bimodal_specific_pairs( tree, names, matrix, names, matrix, 
                                                  queries, msr, standardised, ot); }

  template< class TreeType, class Measure, class OutputIterator >
  std::pair<int, int>
  _csv_matrix_query_bimodal( TreeType &tree, char *filename1, 
                            char *filename2, Measure &msr, bool standardised, OutputIterator ot);

  template< class TreeType, class Measure, class OutputIterator >
  std::pair<int,int>
  _csv_matrix_query_bimodal( TreeType &tree,  char *filename, Measure &msr, bool standardised, OutputIterator ot)
  { return _csv_matrix_query_bimodal(tree,filename,filename,msr, standardised, ot); }

  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _csv_matrix_query_bimodal_specific_pairs( TreeType &tree,  char *matrix_filename_a, 
                                            char *matrix_filename_b, char *queries_filename,
                                            Measure &msr, bool standardised, OutputIterator ot );

  template< class TreeType, class Measure, class OutputIterator>
  std::pair<int, int>
  _csv_matrix_query_bimodal_specific_pairs( TreeType &tree,  char *matrix_filename, 
                                            char *queries_filename,
                                            Measure &msr, bool standardised, OutputIterator ot )
  { return _csv_matrix_query_bimodal_specific_pairs( tree, matrix_filename, matrix_filename, 
                                                     queries_filename, msr, standardised, ot); }
  
  // Input:  Two ranges of iterators that indicate two lists of species names (in std::string format).
  // Output: The value of the current measure for these two sets of species.
  
  template < class TreeType, class RangeIterator, class Measure >
  Number_type _list_query( TreeType &tree,  
                           RangeIterator rbegin_a, RangeIterator rend_a,
                           RangeIterator rbegin_b, RangeIterator rend_b,
                           Measure &msr );
  
  
  // Input: A tree and two txt files that each stores a list of species names, each constituting a subset 
  // of the species in the tree and which appear in random order.
  // Output: The value of the current measure for these two sets of species.
  
  template< class TreeType, class Measure >
  Number_type _list_query(TreeType &tree, char* filename1, char* filename2, Measure &msr);
  
}; // Measure_base_bimodal

} // namespace PhylogeneticMeasures

#include "Measure_base_bimodal_impl.h"

#endif // MEASURE_BASE_BIMODAL_H

