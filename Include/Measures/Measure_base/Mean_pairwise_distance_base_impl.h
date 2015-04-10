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

#ifndef MEAN_PAIRWISE_DISTANCE_BASE_IMPL_H
#define MEAN_PAIRWISE_DISTANCE_BASE_IMPL_H

  template < class KernelType, class TreeType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance_base<KernelType,TreeType>::
  _compute_single_edge_path_costs( Tree_type &tree, int node_index, 
                                   Number_type sum_anc1, Number_type sum_anc2,
                                   Number_type& all_weights )
  {
    Node_type u = tree.node(node_index);
    Number_type sum_off(0.0);

    for( int i=0; i<u.number_of_children(); i++ )
    {
      Node_type v = tree.node(u.children[i]);

      Number_type tmp_sum_anc1 = sum_anc1 + Number_type(v.distance)*Number_type(tree.number_of_leaves()
                                 - v.all_subtree_leaves),
                  tmp_sum_anc2 = sum_anc2 + Number_type(v.distance)*Number_type(v.all_subtree_leaves);

      sum_off += _compute_single_edge_path_costs( tree, u.children[i], tmp_sum_anc1, tmp_sum_anc2, all_weights );
    }

    Number_type so = Number_type(tree.number_of_leaves() - u.all_subtree_leaves)*sum_off,
                sa = Number_type(u.all_subtree_leaves)*sum_anc1,
                si = Number_type(u.all_subtree_leaves)*(all_weights-sum_anc2-sum_off);

    _edge_path_costs[node_index] = so+sa+si;

    return sum_off+(Number_type(u.distance)*Number_type(u.all_subtree_leaves));

  } // _compute_single_edge_path_costs(...)


  // Computes the sum of the costs of all possible simple
  // paths between pairs of leaves in the tree. 
  template < class KernelType, class TreeType>
  typename KernelType::Number_type PhylogeneticMeasures::Mean_pairwise_distance_base<KernelType,TreeType>::
  total_path_costs(Tree_type &tree)
  {
    if( _total_path_costs!= Number_type(-1.0) )
      return _total_path_costs;

    Number_type _total_path_costs = Number_type(0.0);

    tree.assign_all_subtree_leaves(tree.root_index());

    int s(tree.number_of_leaves());

    for( int i=0; i<tree.number_of_nodes()-1; i++ ) // We exclude the edge of the root node.
      _total_path_costs += Number_type(tree.node(i).distance)*Number_type(tree.node(i).all_subtree_leaves)
                                                              *Number_type(s-tree.node(i).all_subtree_leaves);

    return _total_path_costs;
	
  } // total_path_costs()

  // Helper function that computes all the sums-of-path-costs
  // values that are stored for future use.

  template <class KernelType, class TreeType>
  void PhylogeneticMeasures::Mean_pairwise_distance_base<KernelType,TreeType>::
  compute_all_costs_values(Tree_type &tree)
  {
    _sum_all_edges_costs = Number_type(0.0);
    _sum_all_leaf_costs   = Number_type(0.0);

    if(_edge_path_costs.size() == 0 )
      compute_all_edge_path_costs(tree);

    for( int i=0; i<tree.number_of_nodes()-1; i++)
    {
      Node_type u = tree.node(i);

      _sum_all_edges_costs += Number_type(u.distance)*_edge_path_costs[i];

      if( u.number_of_children() == 0 )
        _sum_all_leaf_costs += _edge_path_costs[i]*_edge_path_costs[i];
    }

  } // compute_all_costs_values()


  template <class KernelType, class TreeType>
  void PhylogeneticMeasures::Mean_pairwise_distance_base<KernelType,TreeType>::
  compute_all_edge_path_costs(Tree_type &tree)
  {
    _edge_path_costs.assign(tree.number_of_nodes()-1, Number_type(0.0) );

    Number_type all_weights(0.0);

    for( int i=0; i<tree.number_of_nodes()-1; i++)
      all_weights += Number_type(tree.node(i).distance)*Number_type(tree.node(i).all_subtree_leaves);

    Node_type root = tree.root();

    for( int i=0; i<root.number_of_children(); i++)
    {
      Node_type child = tree.node(root.children[i]);

      _compute_single_edge_path_costs( tree, root.children[i], Number_type(child.distance)*
                                       Number_type(tree.number_of_leaves() - child.all_subtree_leaves) ,
                                       Number_type(child.distance)*Number_type(child.all_subtree_leaves) , all_weights );
    }

  } // compute_all_edge_path_costs()

#endif //MEAN_PAIRWISE_DISTANCE_BASE_IMPL_H
