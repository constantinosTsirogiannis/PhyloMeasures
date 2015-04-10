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

#ifndef COMMUNITY_DISTANCE_IMPL_H
#define COMMUNITY_DISTANCE_IMPL_H

  template < class KernelType>
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Community_distance<KernelType>::
  operator()( RangeIterator rbegin_a, RangeIterator rend_a,
              RangeIterator rbegin_b, RangeIterator rend_b,
              int min_index_a, int max_index_a,
              int min_index_b, int max_index_b )
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    if(rbegin_a == rend_a || rbegin_b == rend_b)
      return Number_type(0.0);

    int min_index = std::min( min_index_a, min_index_b ),
        max_index = std::max( max_index_a, max_index_b );

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index).mark = true;
    p_tree->node(intersection_index).mark_b = true;

    Number_type total_distance(0.0);

    int count_sample_nodes_a = int(rend_a-rbegin_a), 
        count_sample_nodes_b = int(rend_b-rbegin_b);

    // Mark all the nodes which fall on a path that connects leaves of the same sample. 

    p_tree->mark_Steiner_trees_of_samples(rbegin_a, rend_a, rbegin_b, rend_b);
    p_tree->assign_marked_subtree_leaves_both_sets(intersection_index);


    // Trace all edges in the marked subtree. Any edge that connects
    // a parent node v to child node u in this subtree appears as many
    // times in the final solution as k * (n-k) where k = u.marked_subtree_leaves
    // and n is the total number of nodes in the marked subtree.

    for( int i=1; i < p_tree->number_of_marked_nodes(); i++)
    {
      Node_type v = p_tree->node(p_tree->marked_node(i));
      total_distance += Number_type(v.distance)*Number_type(v.marked_subtree_leaves)*
                        Number_type(count_sample_nodes_b - v.marked_subtree_leaves_b);
    }

    for( int i=1; i < p_tree->number_of_marked_nodes_b(); i++)
    {
      Node_type v = p_tree->node(p_tree->marked_node_b(i));
      total_distance += Number_type(v.distance)*Number_type(v.marked_subtree_leaves_b)*
                        Number_type(count_sample_nodes_a - v.marked_subtree_leaves);
    }

    p_tree->unmark_Steiner_trees_of_samples(rbegin_a, rend_a, rbegin_b, rend_b);

    return total_distance/( Number_type(count_sample_nodes_a)*Number_type(count_sample_nodes_b) );

  } // Community_distance<KernelType>::operator()(...)

  

  template < class KernelType>
  typename KernelType::Number_type PhylogeneticMeasures::Community_distance<KernelType>::
  compute_variance( int sample_size_a, int sample_size_b )
  {
    if( sample_size_a < 0 || sample_size_a > p_tree->number_of_leaves() ||
        sample_size_b < 0 || sample_size_b > p_tree->number_of_leaves()   )
    {
      std::string exception_msg;
      exception_msg += " Request to compute variance with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size_a < 1 || sample_size_b < 1 )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    int s(p_tree->number_of_leaves());

    Number_type k1 = Number_type((Number_type(4*(sample_size_a-1)*Number_type(sample_size_b-1)))/
                                 (Number_type(sample_size_a)*
                                  Number_type(sample_size_b)*
                                  Number_type(s)*Number_type(s)*
                                  Number_type(s-1)*Number_type(s-1)) ),

                k2 = ( (Number_type(2*sample_size_a-2)*Number_type(sample_size_b-1))/
                                  (Number_type(sample_size_a)*
                                   Number_type(sample_size_b)*
                                   Number_type(s)*Number_type(s)*
                                   Number_type(s-1)*Number_type(s-1)) ) +
                     (  Number_type(sample_size_a+sample_size_b-2)/
                        (Number_type(sample_size_a)*
                                   Number_type(sample_size_b)*
                                   Number_type(s)*Number_type(s)*
                                   Number_type(s-1)) ) ,

                k3 = (( Number_type(2*sample_size_a-2)*Number_type(sample_size_b-1))/
                                  (Number_type(sample_size_a)*
                                   Number_type(sample_size_b)*
                                   Number_type(s)*Number_type(s)*
                                   Number_type(s-1)*Number_type(s-1)) ) +
                     (  Number_type(2)/
                        (Number_type(sample_size_a)*
                         Number_type(sample_size_b)*
                         Number_type(s)*Number_type(s)));

    if( this->sum_all_edges_costs() == Number_type(-1.0) )
      this->compute_all_costs_values(*p_tree);

    return (  k1 * this->total_path_costs(*p_tree)*this->total_path_costs(*p_tree) +
              (k2-k1)*this->sum_all_leaf_costs() +
              (k1-(Number_type(2.0)*k2)+k3)*this->sum_all_edges_costs() -
              (compute_expectation(sample_size_a,sample_size_b)*
              compute_expectation(sample_size_a,sample_size_b)) );

  } // Number_type compute_variance( ... )

#endif // COMMUNITY_DISTANCE_IMPL_H
