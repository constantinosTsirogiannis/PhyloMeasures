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

#ifndef MEAN_PAIRWISE_DISTANCE_IMPL_H
#define MEAN_PAIRWISE_DISTANCE_IMPL_H

#include<vector>

  template < class KernelType>
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index )
  {
    if(rend-rbegin < 2)
      return Number_type(0.0);

    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index).mark = true;

    Number_type total_distance(0.0);
    int count_sample_nodes= int(rend-rbegin);

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);
    p_tree->assign_marked_subtree_leaves(intersection_index);


    // Trace all edges in the marked subtree. Any edge that connects
    // a parent node v to child node u in this subtree appears as many
    // times in the final solution as k * (n-k) where k = u.marked_subtree_leaves
    // and n is the total number of nodes in the marked subtree.

    for( int i=1; i < p_tree->number_of_marked_nodes(); i++)
    {
      Node_type v = p_tree->node(p_tree->marked_node(i));

      total_distance += Number_type(v.distance)*
                        Number_type(v.marked_subtree_leaves)*
                        Number_type(count_sample_nodes - v.marked_subtree_leaves);
    }

    // Unmark all marked nodes and their 'marked_subtree_leaves' fields.

    p_tree->unmark_Steiner_tree_of_sample(rbegin, rend);

    return total_distance*Number_type(2.0)/
           ( Number_type( count_sample_nodes)*Number_type(count_sample_nodes-1) );

  } // Mean_pairwise_distance<KernelType>::operator()(...)

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_expectation( int sample_size )
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute expectation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size <= 1)
      return Number_type(0.0);

    if(_expectation != Number_type(-1.0) )
      return _expectation;

   int s(p_tree->number_of_leaves());

   _expectation = Number_type(2.0)*(this->total_path_costs(*p_tree))/(Number_type(s)*Number_type(s-1));

    return _expectation;
	
  } // compute_expectation( int sample_size )

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance<KernelType>::
  compute_variance( int sample_size )
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute variance with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size <= 1 || sample_size == p_tree->number_of_leaves() )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    int s(p_tree->number_of_leaves());

    Number_type c1 = Number_type(Number_type(4*(sample_size-2))*Number_type(sample_size-3))/
                                 (Number_type(s)*Number_type(sample_size)*
                                  Number_type(sample_size-1)*Number_type(s-1)*
                                  Number_type(s-2)*Number_type(s-3)),
                c2 = Number_type(4*(sample_size-2))/
                                  (Number_type(s)*Number_type(sample_size)*Number_type(sample_size-1)*Number_type(s-1)*
                                   Number_type(s-2)),
                c3 = Number_type(4)/ (Number_type(s)*Number_type(sample_size)*Number_type(sample_size-1)*
                                      Number_type(s-1));

    if( this->sum_all_edges_costs() == Number_type(-1.0) )
      this->compute_all_costs_values(*p_tree);

    return ( c1 * (this->total_path_costs(*p_tree))*(this->total_path_costs(*p_tree)) +
             (c2-c1)*(this->sum_all_leaf_costs()) +
             (c1-(Number_type(2.0)*c2)+c3)*(this->sum_all_edges_costs()) -
             compute_expectation(sample_size)*compute_expectation(sample_size) );
			 
  } // compute_variance( ... )

#endif //MEAN_PAIRWISE_DISTANCE_IMPL_H
