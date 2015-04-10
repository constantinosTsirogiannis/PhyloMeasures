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

#ifndef COMMUNITY_DISTANCE_NEAREST_TAXON_H
#define COMMUNITY_DISTANCE_NEAREST_TAXON_H

#include<vector>

namespace PhylogeneticMeasures {

template< class KernelType >
struct Community_distance_nearest_taxon: public KernelType::Measure_base_bimodal
{
  typedef KernelType                                              Kernel;
  typedef typename Kernel::Measure_base_bimodal                   Base;
  typedef Community_distance_nearest_taxon<Kernel>                Self;
  typedef typename Kernel::Number_type                            Number_type;
  typedef typename Kernel::Numeric_traits                         Numeric_traits;
  typedef typename Kernel::Community_distance_nearest_taxon_tree  Tree_type;
  typedef typename Tree_type::Node_type                           Node_type;
  typedef typename Kernel::Mean_pairwise_distance                 Mean_pairwise_distance;

  typedef typename Kernel::Exception_type                         Exception_type;
  typedef typename Kernel::Exception_functor                      Exception_functor; 

  
  Community_distance_nearest_taxon(Tree_type &tree)
  { p_tree = &tree;}

  private:

  Number_type _compute_subtree_min_values_a( Tree_type &tree, int current_index );
  void _propagate_min_values_a( Tree_type &tree, int current_index );
  void _compute_rest_tree_min_values_a( Tree_type &tree, int current_index );

  Number_type _compute_subtree_min_values_b( Tree_type &tree, int current_index );
  void _propagate_min_values_b( Tree_type &tree, int current_index );
  void _compute_rest_tree_min_values_b( Tree_type &tree, int current_index );

  void _clear_auxiliary_data( Tree_type &tree, int index );

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_query_internal_directed( std::vector<std::vector<int> > &samples_a,
                                   std::vector<std::pair<int,int> > &min_max_a, 
                                   std::vector<std::vector<int> > &samples_b,
                                   std::vector<std::pair<int,int> > &min_max_b, 
                                   bool is_double_matrix,
                                   OutputIterator ot_a_to_b, OutputIterator ot_b_to_a );

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_query_internal_directed_specific_pairs( std::vector<std::vector<int> > &samples_a,
                                                  std::vector<std::pair<int,int> > &min_max_a, 
                                                  std::vector<std::vector<int> > &samples_b,
                                                  std::vector<std::pair<int,int> > &min_max_b,
                                                  std::vector<std::pair<int,int> > &query_intervals_a,
                                                  std::vector<std::pair<int,int> > &query_intervals_b, 
                                                  bool is_double_matrix,
                                                  OutputIterator ot_a_to_b, OutputIterator ot_b_to_a );

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_distance_query_directed( const std::vector<std::string> &names_a, 
                                   const std::vector<std::vector<bool> > &matrix_a,
                                   const std::vector<std::string> &names_b, 
                                   const std::vector<std::vector<bool> > &matrix_b, 
                                   OutputIterator ot_a_to_b, OutputIterator ot_b_to_a );

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_distance_query_directed_specific_pairs( const std::vector<std::string> &names_a, 
                                                  const std::vector<std::vector<bool> > &matrix_a,
                                                  const std::vector<std::string> &names_b, 
                                                  const std::vector<std::vector<bool> > &matrix_b, 
                                                  const std::vector<std::pair<int, int> > &queries, 
                                                  OutputIterator ot_a_to_b, 
                                                  OutputIterator ot_b_to_a);

  template< class OutputIterator >
  std::pair<int, int>
  _csv_matrix_distance_query_directed( char *matrix_filename_a, char *matrix_filename_b,  
                                       OutputIterator ot_a_to_b, OutputIterator ot_b_to_a);

  template< class OutputIterator >
  std::pair<int, int>
  _csv_matrix_distance_query_directed_specific_pairs( char *matrix_filename_a, char *matrix_filename_b,
                                                      char *queries_filename, OutputIterator ot_a_to_b, 
                                                      OutputIterator ot_b_to_a);


  template< class OutputIterator >
  std::pair<int, int>
  _matrix_query_internal_averaged( std::vector< std::vector<int> >  &samples_a, 
                                   std::vector< std::pair<int,int> > &min_max_a,
                                   std::vector< std::vector<int> >  &samples_b, 
                                   std::vector< std::pair<int,int> > &min_max_b,
                                   bool is_double_matrix,
                                   OutputIterator ot );

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_query_internal_averaged_specific_pairs( std::vector< std::vector<int> >   &samples_a, 
                                                  std::vector< std::pair<int,int> > &min_max_a,
                                                  std::vector< std::vector<int> >   &samples_b, 
                                                  std::vector< std::pair<int,int> > &min_max_b,
                                                  std::vector< std::pair<int,int> > &query_intervals_a,
                                                  std::vector< std::pair<int,int> > &query_intervals_b,
                                                  bool is_double_matrix,
                                                  OutputIterator ot);

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_distance_query_averaged( const std::vector<std::string> &names_a, 
                                   const std::vector<std::vector<bool> > &matrix_a,
                                   const std::vector<std::string> &names_b, 
                                   const std::vector<std::vector<bool> > &matrix_b, 
                                   OutputIterator ot);

  template< class OutputIterator >
  std::pair<int, int>
  _matrix_distance_query_averaged_specific_pairs( const std::vector<std::string> &names_a, 
                                                  const std::vector<std::vector<bool> > &matrix_a,
                                                  const std::vector<std::string> &names_b, 
                                                  const std::vector<std::vector<bool> > &matrix_b, 
                                                  const std::vector<std::pair<int, int> > &queries, 
                                                  OutputIterator ot);

  template< class OutputIterator >
  std::pair<int, int>
  _csv_matrix_distance_query_averaged( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot );


  template< class OutputIterator >
  std::pair<int, int>
  _csv_matrix_distance_query_averaged_specific_pairs( char *matrix_filename_a, char *matrix_filename_b, 
                                                      char *queries_filename, OutputIterator ot );


  public:


  // The following function returns two distances: the sum of nearest distances between any
  // species in sample A and all the species in sample B, and the sum of distances between any
  // species in sample B and all the species in sample A.

  template< class RangeIterator >
  std::pair<Number_type, Number_type> directed_distances( RangeIterator rbegin_a, RangeIterator rend_a,
                                                          RangeIterator rbegin_b, RangeIterator rend_b,
                                                          int min_index_a, int max_index_a,
                                                          int min_index_b, int max_index_b );

  template< class RangeIterator >
  Number_type operator()( RangeIterator rbegin_a, RangeIterator rend_a,
                          RangeIterator rbegin_b, RangeIterator rend_b,
                          int min_index_a, int max_index_a,
                          int min_index_b, int max_index_b )
  {
    Number_type n_a(rend_a-rbegin_a), n_b(rend_b-rbegin_b);

    if(n_a == Number_type(0.0) || n_b == Number_type(0.0))
      return Number_type(0.0);

    std::pair<Number_type,Number_type> dpr = directed_distances( rbegin_a, rend_a, rbegin_b, rend_b,
                                                                 min_index_a,  max_index_a,
                                                                 min_index_b,  max_index_b );

    Number_type distance_a = dpr.first/n_a,
                distance_b = dpr.second/n_b;

    if(distance_a > distance_b)
      return distance_a;

    return distance_b;
	
  } // operator()(...)


  template< class RangeIterator >
  Number_type averaged( RangeIterator rbegin_a, RangeIterator rend_a,
                        RangeIterator rbegin_b, RangeIterator rend_b,
                        int min_index_a, int max_index_a,
                        int min_index_b, int max_index_b )
  {
    Number_type n_a(rend_a-rbegin_a), n_b(rend_b-rbegin_b);

    if(n_a == Number_type(0.0) || n_b == Number_type(0.0))
      return Number_type(0.0);


    std::pair<Number_type,Number_type> dpr = directed_distances( rbegin_a, rend_a, rbegin_b, rend_b,
                                                                 min_index_a,  max_index_a,
                                                                 min_index_b,  max_index_b );
    return (dpr.first + dpr.second)/(n_a+n_b);
  }

  template< class RangeIterator >
  Number_type distance_a_to_b( RangeIterator rbegin_a, RangeIterator rend_a,
                               RangeIterator rbegin_b, RangeIterator rend_b,
                               int min_index_a, int max_index_a,
                               int min_index_b, int max_index_b )
  {
    Number_type n_a(rend_a-rbegin_a), n_b(rend_b-rbegin_b);

    if(n_a == Number_type(0.0) || n_b == Number_type(0.0))
      return Number_type(0.0);

    std::pair<Number_type,Number_type> dpr = directed_distances( rbegin_a, rend_a, rbegin_b, rend_b,
                                                                 min_index_a,  max_index_a,
                                                                 min_index_b,  max_index_b );

    return dpr.first/n_a;
  }

  template< class RangeIterator >
  Number_type distance_b_to_a( RangeIterator rbegin_a, RangeIterator rend_a,
                               RangeIterator rbegin_b, RangeIterator rend_b,
                               int min_index_a, int max_index_a,
                               int min_index_b, int max_index_b )
  {
    Number_type n_a(rend_a-rbegin_a), n_b(rend_b-rbegin_b);

    if(n_a == Number_type(0.0) || n_b == Number_type(0.0))
      return Number_type(0.0);

    std::pair<Number_type,Number_type> dpr = directed_distances( rbegin_a, rend_a, rbegin_b, rend_b,
                                                                 min_index_a,  max_index_a,
                                                                 min_index_b,  max_index_b );

    return dpr.second/n_b;
  }

  template< class RangeIterator >
  std::pair<Number_type, Number_type>
  slow_directed_distances( RangeIterator rbegin_a, RangeIterator rend_a,
                           RangeIterator rbegin_b, RangeIterator rend_b,
                           int min_index_a, int max_index_a,
                           int min_index_b, int max_index_b );


  template< class RangeIterator >
  Number_type slow_operator( RangeIterator rbegin_a, RangeIterator rend_a,
                             RangeIterator rbegin_b, RangeIterator rend_b,
                             int min_index_a, int max_index_a,
                             int min_index_b, int max_index_b )
  {
    std::pair<Number_type, Number_type> res = slow_directed_distances( rbegin_a, rend_a,
                                                                       rbegin_b, rend_b,
                                                                       min_index_a, max_index_a,
                                                                       min_index_b, max_index_b );
    int count_a=0,count_b=0;

    for(RangeIterator r=rbegin_a; r!=rend_a; r++)
      count_a++;

    for(RangeIterator r=rbegin_b; r!=rend_b; r++)
      count_b++;;

    Number_type dist_a = res.first/Number_type(count_a),
                dist_b = res.second/Number_type(count_b);
 
    return std::max(dist_a, dist_b);
  } 

  template< class RangeIterator >
  Number_type slow_averaged( RangeIterator rbegin_a, RangeIterator rend_a,
                             RangeIterator rbegin_b, RangeIterator rend_b,
                             int min_index_a, int max_index_a,
                             int min_index_b, int max_index_b )
  {
    std::pair<Number_type, Number_type> res = slow_directed_distances( rbegin_a, rend_a,
                                                                       rbegin_b, rend_b,
                                                                       min_index_a, max_index_a,
                                                                       min_index_b, max_index_b );

    int count_a=0,count_b=0;

    for(RangeIterator r=rbegin_a; r!=rend_a; r++)
      count_a++;

    for(RangeIterator r=rbegin_b; r!=rend_b; r++)
      count_b++;

    return Number_type(res.first+res.second)/Number_type(count_a+count_b);
  } 

  template< class RangeIterator >
  Number_type slow_distance_a_to_b( RangeIterator rbegin_a, RangeIterator rend_a,
                                    RangeIterator rbegin_b, RangeIterator rend_b,
                                    int min_index_a, int max_index_a,
                                    int min_index_b, int max_index_b )
  {
    std::pair<Number_type, Number_type> res = slow_directed_distances( rbegin_a, rend_a,
                                                                       rbegin_b, rend_b,
                                                                       min_index_a, max_index_a,
                                                                       min_index_b, max_index_b );

    int count_a=0;

    for(RangeIterator r=rbegin_a; r!=rend_a; r++)
      count_a++;

    return res.first/Number_type(count_a);
  } 

  template< class RangeIterator >
  Number_type slow_distance_b_to_a( RangeIterator rbegin_a, RangeIterator rend_a,
                                    RangeIterator rbegin_b, RangeIterator rend_b,
                                    int min_index_a, int max_index_a,
                                    int min_index_b, int max_index_b )
  {
    std::pair<Number_type, Number_type> res = slow_directed_distances( rbegin_a, rend_a,
                                                                       rbegin_b, rend_b,
                                                                       min_index_a, max_index_a,
                                                                       min_index_b, max_index_b );

    int count_b=0;

    for(RangeIterator r=rbegin_b; r!=rend_b; r++)
      count_b++;

    return res.second/Number_type(count_b);
  }


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

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_basic(const std::vector<std::string> &names_a, 
                     const std::vector<std::vector<bool> > &matrix_a, 
                     const std::vector<std::string> &names_b, 
                     const std::vector<std::vector<bool> > &matrix_b, OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names_a, matrix_a, names_b, matrix_b, *this, false, ot); }

  template< class OutputIterator >
  int matrix_query_basic( const std::vector<std::string> &names, 
                          const std::vector<std::vector<bool> > &matrix, OutputIterator ot)
  { return this->_matrix_query_bimodal(*p_tree, names, matrix, *this,  false, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_specific_pairs_basic( const std::vector<std::string> &names_a, 
                                     const std::vector<std::vector<bool> > &matrix_a, 
                                     const std::vector<std::string> &names_b, 
                                     const std::vector<std::vector<bool> > &matrix_b, 
                                     const std::vector<std::pair<int,int> > &queries, 
                                     OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names_a, matrix_a,
                                                      names_b, matrix_b, queries, *this, false, ot); }

  template< class OutputIterator >
  int matrix_query_specific_pairs_basic( const std::vector<std::string> &names, 
                                         const std::vector<std::vector<bool> > &matrix,  
                                         const std::vector<std::pair<int,int> > &queries, OutputIterator ot)
  { return this->_matrix_query_bimodal_specific_pairs(*p_tree, names, matrix, queries, *this, false, ot).first; }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_basic( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename_a, matrix_filename_b, *this, false, ot); }

  template< class OutputIterator >
  int csv_matrix_query_basic( char *matrix_filename, OutputIterator ot)
  { return this->_csv_matrix_query_bimodal(*p_tree, matrix_filename, *this,  false, ot).first; }


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

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_directed_basic(const std::vector<std::string> &names_a, 
                              const std::vector<std::vector<bool> > &matrix_a, 
                              const std::vector<std::string> &names_b, 
                              const std::vector<std::vector<bool> > &matrix_b, 
                              OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  { return this->_matrix_distance_query_directed( names_a, matrix_a, names_b, matrix_b, ot_a_to_b, ot_b_to_a); }

  template< class OutputIterator >
  int matrix_query_directed_basic(const std::vector<std::string> &names, 
                                  const std::vector<std::vector<bool> > &matrix,
                                  OutputIterator ot)  
  {
    std::vector<Number_type> foo; 
    return this->_matrix_distance_query_directed(names, matrix, names, matrix, ot, std::back_inserter(foo)).first; 
  }


  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_directed_specific_pairs_basic( const std::vector<std::string> &names_a, 
                                              const std::vector<std::vector<bool> > &matrix_a, 
                                              const std::vector<std::string> &names_b, 
                                              const std::vector<std::vector<bool> > &matrix_b, 
                                              const std::vector<std::pair<int,int> > &queries, 
                                              OutputIterator ot_a_to_b,
                                              OutputIterator ot_b_to_a )
  { return this->_matrix_distance_query_directed_specific_pairs(names_a, matrix_a,
                                                                names_b, matrix_b, queries, ot_a_to_b, ot_b_to_a); }

  template< class OutputIterator >
  int matrix_query_directed_specific_pairs_basic( const std::vector<std::string> &names, 
                                                  const std::vector<std::vector<bool> > &matrix,  
                                                  const std::vector<std::pair<int,int> > &queries,
                                                  OutputIterator ot_a_to_b,
                                                  OutputIterator ot_b_to_a )
  { return this->_matrix_distance_query_directed_specific_pairs(names, matrix, names, matrix, 
                                                                queries, ot_a_to_b, ot_b_to_a).first; }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////


  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_directed_basic( char *matrix_filename_a, char *matrix_filename_b, 
                                   OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  { return this->_csv_matrix_distance_query_directed( matrix_filename_a, matrix_filename_b, ot_a_to_b, ot_b_to_a);}


  template< class OutputIterator >
  int csv_matrix_query_directed_basic( char *matrix_filename, OutputIterator ot )
  {
    std::vector<Number_type> foo;  
    return this->_csv_matrix_distance_query_directed( matrix_filename, matrix_filename, ot, std::back_inserter(foo)).first;
  }


  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_directed_specific_pairs_basic( char *matrix_filename_a, 
                                                  char *matrix_filename_b, char *queries_filename,
                                                  OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  { return this->_csv_matrix_distance_query_directed_specific_pairs(matrix_filename_a, matrix_filename_b, 
                                                                    queries_filename, ot_a_to_b, ot_b_to_a);}

  template< class OutputIterator >
  int csv_matrix_query_directed_specific_pairs_basic( char *matrix_filename, 
                                                      char *queries_filename, OutputIterator ot_a_to_b, 
                                                      OutputIterator ot_b_to_a )
  { return this->_csv_matrix_distance_query_directed_specific_pairs(matrix_filename, matrix_filename, 
                                                                    queries_filename, ot_a_to_b, ot_b_to_a).first;}

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_averaged_basic( const std::vector<std::string> &names_a, 
                               const std::vector<std::vector<bool> > &matrix_a, 
                               const std::vector<std::string> &names_b, 
                               const std::vector<std::vector<bool> > &matrix_b, 
                               OutputIterator ot )
  { return this->_matrix_distance_query_averaged(names_a, matrix_a, names_b, matrix_b, ot); }

  template< class OutputIterator >
  int matrix_query_averaged_basic( const std::vector<std::string> &names, 
                                   const std::vector<std::vector<bool> > &matrix,
                                   OutputIterator ot )
  { return this->_matrix_distance_query_averaged(names, matrix, names, matrix, ot).first; }


  template< class OutputIterator >
  std::pair<int, int>
  matrix_query_averaged_specific_pairs_basic( const std::vector<std::string> &names_a, 
                                              const std::vector<std::vector<bool> > &matrix_a, 
                                              const std::vector<std::string> &names_b, 
                                              const std::vector<std::vector<bool> > &matrix_b, 
                                              const std::vector<std::pair<int,int> > &queries, 
                                              OutputIterator ot )
  { return this->_matrix_distance_query_averaged_specific_pairs( names_a, matrix_a, names_b, matrix_b, queries, ot); }

  template< class OutputIterator >
  int matrix_query_averaged_specific_pairs_basic( const std::vector<std::string> &names, 
                                                  const std::vector<std::vector<bool> > &matrix,  
                                                  const std::vector<std::pair<int,int> > &queries,
                                                  OutputIterator ot )
  { return this->_matrix_distance_query_averaged_specific_pairs(names, matrix, names, matrix, queries, ot).first; }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  template< class OutputIterator >
  std::pair<int, int>
  csv_matrix_query_averaged_basic( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot )
  { return this->_csv_matrix_distance_query_averaged(matrix_filename_a, matrix_filename_b, ot); }

  template< class OutputIterator >
  int csv_matrix_query_averaged_basic( char *matrix_filename, OutputIterator ot )
  { return this->_csv_matrix_distance_query_averaged(matrix_filename, matrix_filename, ot).first; }

  template< class OutputIterator>
  std::pair<int,int>
  csv_matrix_query_averaged_specific_pairs_basic( char *matrix_filename_a, char *matrix_filename_b,
                                                  char *queries_filename, OutputIterator ot )
  { return this->_csv_matrix_distance_query_averaged_specific_pairs(matrix_filename_a,
                                                                    matrix_filename_b,
                                                                    queries_filename,ot); }

  template< class OutputIterator>
  int csv_matrix_query_averaged_specific_pairs_basic( char *matrix_filename,
                                                      char *queries_filename, OutputIterator ot )
  { return this->_csv_matrix_distance_query_averaged_specific_pairs(matrix_filename, matrix_filename, 
                                                                    queries_filename, ot).first; }

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////

  Number_type compute_expectation( int sample_size_a, int sample_size_b )
  {  
    std::string exception_msg;
    exception_msg += " The computation of the expectation of CDNT is not available.\n";     
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
    exception_msg += " The computation of the variance of CDNT is not available.\n";     
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
    exception_msg += " The computation of the deviation of the CDNT is not available.\n";     
    Exception_type excp;
    excp.get_error_message(exception_msg);
    Exception_functor excf;
    excf(excp); 
    return Number_type(-1.0); 
  }

  private:

  Tree_type    *p_tree; // Stores a pointer to a phylogenetic tree object.
  int          _sample_size_a, _sample_size_b;
  int          _number_of_leaves;

}; // struct Community_distance_nearest_taxon

} // namespace PhylogeneneticMeasures

#include "Community_distance_nearest_taxon_impl.h"

#endif // COMMUNITY_DISTANCE_NEAREST_TAXON_H
