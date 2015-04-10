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

#ifndef COMMUNITY_DISTANCE_NEAREST_TAXON_IMPL_H
#define COMMUNITY_DISTANCE_NEAREST_TAXON_IMPL_H

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf in its subtree that belongs to sample A.
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _compute_subtree_min_values_a( Tree_type &tree, int current_index )
  {

    for( int i=0; i<tree.node(current_index).marked_children.size(); i++ )
    {
      Number_type cmin = _compute_subtree_min_values_a(tree, tree.node(current_index).marked_children[i]);

      if( tree.node(current_index).first_min_a == Number_type(-1.0) || 
          cmin < tree.node(current_index).first_min_a)
      {
        tree.node(current_index).second_min_a = tree.node(current_index).first_min_a;
        tree.node(current_index).first_min_a = cmin;
      }
      else if( tree.node(current_index).second_min_a == Number_type(-1.0) || 
               cmin < tree.node(current_index).second_min_a )   
        tree.node(current_index).second_min_a = cmin;

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

    if( tree.node(current_index).marked_children.size() == 0 )
    {
      tree.node(current_index).first_min_a = Number_type(0.0);
      tree.node(current_index).second_min_a = Number_type(0.0);
    }
  

    return tree.node(current_index).first_min_a + Number_type(tree.node(current_index).distance);

  } // _compute_subtree_min_values_a( ... )


  // Helper function for adjusting the .rest_tree_min_a values in the
  // the Steiner tree of sample B nodes.
  template< class KernelType >
  void PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _propagate_min_values_a( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children_b.size(); i++ )
    {
      Node_type ch_node = tree.node(current_node.marked_children_b[i]);

      tree.node(current_node.marked_children_b[i]).rest_tree_min_a =
        current_node.rest_tree_min_a + Number_type(ch_node.distance);

      _propagate_min_values_a(tree,current_node.marked_children_b[i]);

    } // for( int i=0; i<current_node.marked_children_b.size(); i++ )

  } // void _propagate_min_values_a( ... )


  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf which belongs to sample A
  // and which lies outside its subtree.
  // This function produces a meaningful result only if it is ran after a call
  // of the function _compute_subtree_min_values().
  template< class KernelType >
  void PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _compute_rest_tree_min_values_a( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Node_type child = tree.node(current_node.marked_children[i]);

      if(current_node.rest_tree_min_a == Number_type(-1.0) &&  current_node.second_min_a == Number_type(-1.0))
      {
        if(tree.node(current_node.marked_children[i]).number_of_children() != 0)
          tree.node(current_node.marked_children[i]).rest_tree_min_a = Number_type(-1.0);
        else
          tree.node(current_node.marked_children[i]).rest_tree_min_a = Number_type(0.0);
      }
      else
      {
        Number_type cmin = child.first_min_a + Number_type(child.distance),
                    comp_min;

        if( current_node.first_min_a == cmin )
          comp_min = current_node.second_min_a;
        else
          comp_min = current_node.first_min_a;

        if( (current_node.rest_tree_min_a < comp_min && current_node.rest_tree_min_a != Number_type(-1.0))
            || comp_min == Number_type(-1.0) )
          tree.node(current_node.marked_children[i]).rest_tree_min_a = 
                current_node.rest_tree_min_a + Number_type(child.distance);
        else
          tree.node(current_node.marked_children[i]).rest_tree_min_a = comp_min + Number_type(child.distance);
      }

      _compute_rest_tree_min_values_a(tree, current_node.marked_children[i]);

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

    // Extract the children nodes that are part of the Steiner tree of sample B
    // and which do not appear in the Steiner tree of sample A.
    // For each such child, update the .rest_tree_min_a distances in its subtree.

    for( int i=0; i<current_node.marked_children_b.size(); i++ )
    {
      Node_type ch_node = tree.node(current_node.marked_children_b[i]);

      if(ch_node.first_min_a==Number_type(-1.0))
      {
        if( (current_node.rest_tree_min_a < current_node.first_min_a
            && current_node.rest_tree_min_a != Number_type(-1.0))
            || current_node.first_min_a == Number_type(-1.0) )
          tree.node(current_node.marked_children_b[i]).rest_tree_min_a =
                   current_node.rest_tree_min_a + Number_type(ch_node.distance);
        else
          tree.node(current_node.marked_children_b[i]).rest_tree_min_a =
                   current_node.first_min_a + Number_type(ch_node.distance);

        _propagate_min_values_a(tree,current_node.marked_children_b[i]);
      }

    } // for( int i=0; i<current_node.marked_children_b.size(); i++ )

  } // _compute_rest_tree_min_values_a( ... )

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf in its subtree that belongs to sample B.
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _compute_subtree_min_values_b( Tree_type &tree, int current_index )
  {
    for( int i=0; i<tree.node(current_index).marked_children_b.size(); i++ )
    {
      Number_type cmin = _compute_subtree_min_values_b(tree, tree.node(current_index).marked_children_b[i]);

      if( tree.node(current_index).first_min_b == Number_type(-1.0) || 
          cmin < tree.node(current_index).first_min_b)
      {
        tree.node(current_index).second_min_b = tree.node(current_index).first_min_b;
        tree.node(current_index).first_min_b = cmin;
      }
      else if( tree.node(current_index).second_min_b == Number_type(-1.0) || 
               cmin < tree.node(current_index).second_min_b )
        tree.node(current_index).second_min_b = cmin;

    } // for( int i=0; i<current_node.marked_children_b.size(); i++ )

    if( tree.node(current_index).marked_children_b.size() == 0 )
    {
      tree.node(current_index).first_min_b = Number_type(0.0);
      tree.node(current_index).second_min_b = Number_type(0.0);
    }

    return tree.node(current_index).first_min_b + Number_type(tree.node(current_index).distance);

  } // _compute_subtree_min_values_b( ... )

  // Helper function for adjusting the .rest_tree_min_b values in the
  // the Steiner tree of sample A nodes.
  template< class KernelType >
  void PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _propagate_min_values_b( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Node_type ch_node = tree.node(current_node.marked_children[i]);

      tree.node(current_node.marked_children[i]).rest_tree_min_b =
        current_node.rest_tree_min_b + Number_type(ch_node.distance);

      _propagate_min_values_b(tree,current_node.marked_children[i]);

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

  } // void _propagate_min_values_b( ... )

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf which belongs to sample B
  // and which lies outside its subtree.
  // This function produces a meaningful result only if it is ran after a call
  // of the function _compute_subtree_min_values().
  template< class KernelType >
  void PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _compute_rest_tree_min_values_b( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children_b.size(); i++ )
    {
      Node_type child = tree.node(current_node.marked_children_b[i]);

  
      if(current_node.rest_tree_min_b == Number_type(-1.0) &&  current_node.second_min_b == Number_type(-1.0))
      {
        if(tree.node(current_node.marked_children_b[i]).number_of_children() != 0)
          tree.node(current_node.marked_children_b[i]).rest_tree_min_b = Number_type(-1.0);
        else
          tree.node(current_node.marked_children_b[i]).rest_tree_min_b = Number_type(0.0);
      }
      else
      {
        Number_type cmin = child.first_min_b + Number_type(child.distance),
                    comp_min;

        if( current_node.first_min_b == cmin )
          comp_min = current_node.second_min_b;
        else
          comp_min = current_node.first_min_b;

        if( (current_node.rest_tree_min_b < comp_min && current_node.rest_tree_min_b != Number_type(-1.0))
            || comp_min == Number_type(-1.0) )
          tree.node(current_node.marked_children_b[i]).rest_tree_min_b = current_node.rest_tree_min_b + 
                                                                         Number_type(child.distance);
        else
          tree.node(current_node.marked_children_b[i]).rest_tree_min_b = comp_min + Number_type(child.distance);
       }

      _compute_rest_tree_min_values_b(tree, current_node.marked_children_b[i]);

    } // for( int i=0; i<current_node.marked_children.size(); i++ )


    // Extract the children nodes that are part of the Steiner tree of sample A
    // and which do not appear in the Steiner tree of sample B.
    // For each such child, update the .rest_tree_min_b distances in its subtree.

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Node_type ch_node = tree.node(current_node.marked_children[i]);

      if(ch_node.first_min_b==Number_type(-1.0))
      {
        if( (current_node.rest_tree_min_b < current_node.first_min_b
            && current_node.rest_tree_min_b != Number_type(-1.0))
            || current_node.first_min_b == Number_type(-1.0) )
          tree.node(current_node.marked_children[i]).rest_tree_min_b =
                   current_node.rest_tree_min_b + Number_type(ch_node.distance);
        else
          tree.node(current_node.marked_children[i]).rest_tree_min_b =
                   current_node.first_min_b + Number_type(ch_node.distance);

        _propagate_min_values_b(tree,current_node.marked_children[i]);
      }

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

  } // _compute_rest_tree_min_values_b( ... )


  // The next function clears the auxiliary data that are stored in
  // the tree nodes.

  template< class KernelType >
  void PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
       _clear_auxiliary_data( Tree_type &tree, int index )
  {
      tree.node(index).first_min_a = Number_type(-1.0);
      tree.node(index).second_min_a = Number_type(-1.0);
      tree.node(index).rest_tree_min_a = Number_type(-1.0);
      tree.node(index).first_min_b = Number_type(-1.0);
      tree.node(index).second_min_b = Number_type(-1.0);
      tree.node(index).rest_tree_min_b = Number_type(-1.0);

      tree.node(index).mark = false;
      tree.node(index).mark_b = false;
      tree.node(index).marked_subtree_leaves = 0;

      for( int i=0; i<tree.node(index).marked_children.size(); i++ )
        _clear_auxiliary_data(tree,tree.node(index).marked_children[i]);

      for( int i=0; i<tree.node(index).marked_children_b.size(); i++ )
        _clear_auxiliary_data(tree,tree.node(index).marked_children_b[i]);

      tree.node(index).marked_children.clear();
      tree.node(index).marked_children_b.clear();
	  
  } //_clear_auxiliary_data( Tree_type &tree, int index )

  template< class KernelType >
  template < class RangeIterator >
  std::pair< typename KernelType::Number_type, typename KernelType::Number_type >
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  directed_distances( RangeIterator rbegin_a, RangeIterator rend_a,
                      RangeIterator rbegin_b, RangeIterator rend_b,
                      int min_index_a, int max_index_a,
                      int min_index_b, int max_index_b)
  {
    if(p_tree->number_of_nodes()<2)
      return std::make_pair(Number_type(0.0), Number_type(0.0));

    if(rbegin_a == rend_a || rbegin_b == rend_b)
      return std::make_pair(Number_type(0.0),Number_type(0.0));

    int min_index = std::min( min_index_a, min_index_b ),
        max_index = std::max( max_index_a, max_index_b );

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two extremal paths coincide (intersection_index
    // designates a leaf node) then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return std::make_pair(Number_type(0.0), Number_type(0.0));

    // Next we traverse the Steiner tree of all leaf nodes, and compute which are the nodes
    // that stand on the Steiner tree of the leaves of sample A (marked_children).
    // In fact, the code below will mark also the nodes that stand on the
    // path from the intersection_index to the root of the Steiner tree of the
    // sample A leaves.

    p_tree->node(intersection_index).mark = true;
    p_tree->node(intersection_index).mark_b = true;

    p_tree->mark_Steiner_trees_of_samples(rbegin_a, rend_a, rbegin_b, rend_b);

    _compute_subtree_min_values_a( *p_tree, intersection_index );
    _compute_rest_tree_min_values_a( *p_tree, intersection_index );
    _compute_subtree_min_values_b( *p_tree, intersection_index );
    _compute_rest_tree_min_values_b( *p_tree, intersection_index );

    // Nodes that appear in both samples have nearest taxon distance = 0

    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
      p_tree->node(*rit).rest_tree_min_a = Number_type(0.0);

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
      p_tree->node(*rit).rest_tree_min_b = Number_type(0.0);

    Number_type total_dist_a(0.0), total_dist_b(0.0);

    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
      total_dist_a = total_dist_a + p_tree->node((*rit)).rest_tree_min_b;

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
      total_dist_b = total_dist_b + p_tree->node((*rit)).rest_tree_min_a;


    // Clear auxiliary data.

    _clear_auxiliary_data(*p_tree,intersection_index);

    return std::make_pair(total_dist_a, total_dist_b);

  } //directed_distances(...)


  template< class KernelType >
  template < class RangeIterator >
  std::pair< typename KernelType::Number_type, typename KernelType::Number_type >
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  slow_directed_distances( RangeIterator rbegin_a, RangeIterator rend_a,
                           RangeIterator rbegin_b, RangeIterator rend_b,
                           int min_index_a, int max_index_a,
                           int min_index_b, int max_index_b)
  {
    std::vector<Number_type> min_costs_a, min_costs_b;

    int count_a = 0;

    for(RangeIterator r1 = rbegin_a; r1 != rend_a; r1++)
    {
      min_costs_a.push_back(Number_type(-1.0));

      int count_b=0; 

      for(RangeIterator r2 = rbegin_b; r2 != rend_b; r2++)
      {
        std::vector<int> node_pair;

        node_pair.push_back(*r1);
        node_pair.push_back(*r2);

        if( r1 == rbegin_a )
          min_costs_b.push_back(Number_type(-1.0));

        int steiner_root = p_tree->compute_intersection_node_index(*r1,*r2);

        Number_type dist(0.0);

        if( p_tree->node(steiner_root).number_of_children() != 0)
        {
          p_tree->node(steiner_root).mark = true;

          dist = p_tree->mark_Steiner_tree_of_sample(node_pair.begin(), node_pair.end());
          p_tree->unmark_Steiner_tree_of_sample(node_pair.begin(), node_pair.end());
        }

        if(min_costs_a[count_a] > dist || min_costs_a[count_a] == Number_type(-1.0) )
          min_costs_a[count_a] = dist;

        if(min_costs_b[count_b] > dist || min_costs_b[count_b] == Number_type(-1.0) )
          min_costs_b[count_b] = dist;

        count_b++; 
      }  

      count_a++;

    } // for(RangeIterator r1 = rbegin_a; r1 != rend_a; r1++)

    Number_type total_dist_a = Number_type(0.0),
                total_dist_b = Number_type(0.0);

    int i=0; 

    for(RangeIterator r1 = rbegin_a; r1 != rend_a; r1++)
    {
      if(min_costs_a[i] != Number_type(-1.0))
        total_dist_a += min_costs_a[i];

      i++;
    }

    i=0;

    if( rbegin_a != rend_a)
      for(RangeIterator r1 = rbegin_b; r1 != rend_b; r1++)
      { 
 
        if(min_costs_b[i] != Number_type(-1.0))
          total_dist_b += min_costs_b[i];

        i++;
      }

    return std::make_pair(total_dist_a, total_dist_b);

  } // slow_directed_distances(...)

  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_query_internal_directed( std::vector<std::vector<int> > &samples_a,
                                   std::vector<std::pair<int,int> > &min_max_a, 
                                   std::vector<std::vector<int> > &samples_b,
                                   std::vector<std::pair<int,int> > &min_max_b, 
                                   bool is_double_matrix,
                                   OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  {
    std::vector< std::vector<Number_type> > tmp_res_a_to_b;

    if(!is_double_matrix)
    {
      tmp_res_a_to_b.assign(samples_a.size(),std::vector<Number_type>());

      for( int i=0; i<samples_a.size(); i++ )
        tmp_res_a_to_b[i].assign(samples_a.size(),Number_type(0.0));
    }

    if(is_double_matrix)
    {
      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<samples_b.size(); j++ )
        {
          std::pair<Number_type, Number_type > dpr = directed_distances(samples_a[i].begin(), samples_a[i].end(),
                                                                        samples_b[j].begin(), samples_b[j].end(),
                                                                        min_max_a[i].first, min_max_a[i].second,
                                                                        min_max_b[j].first, min_max_b[j].second  );


          if(samples_a[i].size() != 0)
            *ot_a_to_b++ = dpr.first/Number_type(samples_a[i].size());
          else
            *ot_a_to_b++ = Number_type(0.0);

          if(samples_b[j].size() != 0)
            *ot_b_to_a++ = dpr.second/Number_type(samples_b[j].size());
          else
            *ot_b_to_a++ = Number_type(0.0);

        }

      return std::make_pair(samples_a.size(), samples_b.size());

    } // if(is_double_matrix)
    else
    {
      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<=i; j++ )
        {
          std::pair<Number_type, Number_type > dpr = directed_distances(samples_a[i].begin(), samples_a[i].end(),
                                                                        samples_a[j].begin(), samples_a[j].end(),
                                                                        min_max_a[i].first, min_max_a[i].second,
                                                                        min_max_a[j].first, min_max_a[j].second  );
          
          if(samples_a[i].size() != 0)
            tmp_res_a_to_b[i][j] = dpr.first/Number_type(samples_a[i].size());
          else
            tmp_res_a_to_b[i][j] = Number_type(0.0);

          if(samples_a[j].size() != 0)
            tmp_res_a_to_b[j][i] = dpr.second/Number_type(samples_a[j].size());
          else
            tmp_res_a_to_b[j][i] = Number_type(0.0);

        } // for( int j=0; j<=i; j++ )


      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<samples_a.size(); j++ )
          *ot_a_to_b++ = tmp_res_a_to_b[i][j];

      return std::make_pair(samples_a.size(), samples_a.size());

    } // else of if(is_double_matrix)

  } // _matrix_query_internal_directed(...) 


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_query_internal_directed_specific_pairs( std::vector<std::vector<int> > &samples_a,
                                                  std::vector<std::pair<int,int> > &min_max_a, 
                                                  std::vector<std::vector<int> > &samples_b,
                                                  std::vector<std::pair<int,int> > &min_max_b,
                                                  std::vector<std::pair<int,int> > &query_intervals_a,
                                                  std::vector<std::pair<int,int> > &query_intervals_b, 
                                                  bool is_double_matrix,
                                                  OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
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
      for( int i=0; i<query_intervals_a.size(); i++ )
        for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
          for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
          {
            std::pair<Number_type, Number_type> res = directed_distances( samples_a[j].begin(), samples_a[j].end(),
                                                                          samples_b[k].begin(), samples_b[k].end(),
                                                                          min_max_a[j].first, min_max_a[j].second,
                                                                          min_max_b[k].first, min_max_b[k].second  );

             if(samples_a[j].size() != 0)
               *ot_a_to_b++ = res.first/Number_type(samples_a[j].size());
             else
               *ot_a_to_b++ = Number_type(0.0);

             if(samples_b[k].size() != 0)
               *ot_b_to_a++ = res.second/Number_type(samples_b[k].size());
             else
               *ot_b_to_a++ = Number_type(0.0);

          } // for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )

      return std::make_pair(samples_a.size(), samples_b.size());

    } // if(is_double_matrix)
    else
    {
      for( int i=0; i<query_intervals_a.size(); i++ )
        for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
          for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
          {
            std::pair<Number_type, Number_type> res = directed_distances( samples_a[j].begin(), samples_a[j].end(),
                                                                          samples_a[k].begin(), samples_a[k].end(),
                                                                          min_max_a[j].first, min_max_a[j].second,
                                                                          min_max_a[k].first, min_max_a[k].second  );

             if(samples_a[j].size() != 0)
               *ot_a_to_b++ = res.first/Number_type(samples_a[j].size());
             else
               *ot_a_to_b++ = Number_type(0.0);

             if(samples_a[k].size() != 0)
               *ot_b_to_a++ = res.second/Number_type(samples_a[k].size());
             else
               *ot_b_to_a++ = Number_type(0.0);

          }  // for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )

      return std::make_pair(samples_a.size(), samples_a.size());

    } // else of if(is_double_matrix)

  } // _matrix_query_internal_directed_specific_pairs(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_distance_query_directed( const std::vector<std::string> &names_a, 
                                   const std::vector<std::vector<bool> > &matrix_a,
                                   const std::vector<std::string> &names_b, 
                                   const std::vector<std::vector<bool> > &matrix_b, 
                                   OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    this->_extract_samples_from_matrix(*p_tree, names_a, matrix_a, 
                                       std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(&matrix_a != &matrix_b)
      this->_extract_samples_from_matrix(*p_tree, names_b, matrix_b,
                                         std::back_inserter(samples_b), std::back_inserter(min_max_b) );

    bool is_double_matrix = (&matrix_a != &matrix_b);

    return _matrix_query_internal_directed( samples_a, min_max_a, samples_b, min_max_b,
                                            is_double_matrix, ot_a_to_b, ot_b_to_a );
  }

  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_distance_query_directed_specific_pairs( const std::vector<std::string> &names_a, 
                                                  const std::vector<std::vector<bool> > &matrix_a,
                                                  const std::vector<std::string> &names_b, 
                                                  const std::vector<std::vector<bool> > &matrix_b, 
                                                  const std::vector<std::pair<int, int> > &queries, 
                                                  OutputIterator ot_a_to_b, 
                                                  OutputIterator ot_b_to_a)
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (&matrix_a != &matrix_b);

    this->_extract_samples_from_matrix(*p_tree, names_a, matrix_a, 
                                       std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_matrix(*p_tree, names_b, matrix_b, 
                                         std::back_inserter(samples_b), std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    for(int i=0; i<queries.size(); i++)
    {
      query_intervals_a.push_back(std::make_pair(queries[i].first, queries[i].first));
      query_intervals_b.push_back(std::make_pair(queries[i].second, queries[i].second));
    } 

    return _matrix_query_internal_directed_specific_pairs( samples_a, min_max_a, samples_b,
                                                           min_max_b, query_intervals_a, 
                                                           query_intervals_b, is_double_matrix,
                                                           ot_a_to_b, ot_b_to_a);
  }


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _csv_matrix_distance_query_directed( char *matrix_filename_a, char *matrix_filename_b, 
                                       OutputIterator ot_a_to_b, OutputIterator ot_b_to_a )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(matrix_filename_a) != std::string(matrix_filename_b));

    this->_extract_samples_from_file(*p_tree, matrix_filename_a, 
                                     std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_file(*p_tree, matrix_filename_b, 
                                       std::back_inserter(samples_b), std::back_inserter(min_max_b) );

   return _matrix_query_internal_directed( samples_a, min_max_a, samples_b, min_max_b,
                                           is_double_matrix, ot_a_to_b, ot_b_to_a );

  } // _csv_matrix_distance_query_directed(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _csv_matrix_distance_query_directed_specific_pairs( char *matrix_filename_a, char *matrix_filename_b,
                                                      char *queries_filename, OutputIterator ot_a_to_b, 
                                                      OutputIterator ot_b_to_a)
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(matrix_filename_a) != std::string(matrix_filename_b));

    this->_extract_samples_from_file(*p_tree, matrix_filename_a, 
                                     std::back_inserter(samples_a), 
                                     std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_file(*p_tree, matrix_filename_b, 
                                       std::back_inserter(samples_b), 
                                       std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    this->_extract_queries_from_file( queries_filename, 
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

    return _matrix_query_internal_directed_specific_pairs(samples_a, min_max_a, samples_b, min_max_b, 
                                                          query_intervals_a, query_intervals_b, 
                                                          is_double_matrix, ot_a_to_b, ot_b_to_a);  

  } // _csv_matrix_distance_query_directed_specific_pairs(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_query_internal_averaged( std::vector< std::vector<int> >  &samples_a, 
                                   std::vector< std::pair<int,int> > &min_max_a,
                                   std::vector< std::vector<int> >  &samples_b, 
                                   std::vector< std::pair<int,int> > &min_max_b,
                                   bool is_double_matrix,
                                   OutputIterator ot )
  {
    std::vector< std::vector<Number_type> > tmp_res;

    tmp_res.assign(samples_a.size(),std::vector<Number_type>());

    if(is_double_matrix)
      for( int i=0; i<samples_a.size(); i++ )
        tmp_res[i].assign(samples_b.size(),Number_type(0.0));
    else
      for( int i=0; i<samples_a.size(); i++ )
        tmp_res[i].assign(samples_a.size(),Number_type(0.0));

    if(is_double_matrix)
      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<samples_b.size(); j++ )
        {
          Number_type dpr = this->averaged(samples_a[i].begin(), samples_a[i].end(),
                                           samples_b[j].begin(), samples_b[j].end(),
                                           min_max_a[i].first, min_max_a[i].second,
                                           min_max_b[j].first, min_max_b[j].second  );
          tmp_res[i][j] = dpr;
        }
    else
      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<=i; j++ )
        {
          Number_type dpr = averaged(samples_a[i].begin(), samples_a[i].end(),
                                      samples_a[j].begin(), samples_a[j].end(),
                                      min_max_a[i].first, min_max_a[i].second,
                                      min_max_a[j].first, min_max_a[j].second  );

          tmp_res[i][j] = tmp_res[j][i] = dpr;
        }

    if(is_double_matrix)
    {
      for( int i=0; i<samples_a.size(); i++ )
        for( int j=0; j<samples_b.size(); j++ )
          *ot++ = tmp_res[i][j];

      return std::make_pair(samples_a.size(), samples_b.size());
    }

    for( int i=0; i<samples_a.size(); i++ )
      for( int j=0; j<samples_a.size(); j++ )
        *ot++ = tmp_res[i][j];

    return std::make_pair(samples_a.size(), samples_a.size());

  } // _matrix_query_internal_averaged(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_query_internal_averaged_specific_pairs( std::vector< std::vector<int> >   &samples_a, 
                                                  std::vector< std::pair<int,int> > &min_max_a,
                                                  std::vector< std::vector<int> >   &samples_b, 
                                                  std::vector< std::pair<int,int> > &min_max_b,
                                                  std::vector< std::pair<int,int> > &query_intervals_a,
                                                  std::vector< std::pair<int,int> > &query_intervals_b,
                                                  bool is_double_matrix,
                                                  OutputIterator ot)
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
      for( int i=0; i<query_intervals_a.size(); i++ )
        for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
          for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
            *ot++ = this->averaged( samples_a[j].begin(), samples_a[j].end(),
                                    samples_b[k].begin(), samples_b[k].end(),
                                    min_max_a[j].first, min_max_a[j].second,
                                    min_max_b[k].first, min_max_b[k].second  );

      return std::make_pair(samples_a.size(), samples_b.size());
    } // if(is_double_matrix)
    else
    {
      for( int i=0; i<query_intervals_a.size(); i++ )
        for( int j=query_intervals_a[i].first; j<=query_intervals_a[i].second; j++ )
          for( int k=query_intervals_b[i].first; k<=query_intervals_b[i].second; k++ )
            *ot++ = this->averaged( samples_a[j].begin(), samples_a[j].end(),
                                    samples_a[k].begin(), samples_a[k].end(),
                                    min_max_a[j].first, min_max_a[j].second,
                                    min_max_a[k].first, min_max_a[k].second  );

      return std::make_pair(samples_a.size(), samples_a.size());

    } // else of if(is_double_matrix)

  } // _matrix_query_internal_averaged_specific_pairs(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_distance_query_averaged( const std::vector<std::string> &names_a, 
                                   const std::vector<std::vector<bool> > &matrix_a,
                                   const std::vector<std::string> &names_b, 
                                   const std::vector<std::vector<bool> > &matrix_b, 
                                   OutputIterator ot)
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (&matrix_a != &matrix_b);

    this->_extract_samples_from_matrix(*p_tree, names_a, matrix_a,
                                       std::back_inserter(samples_a), 
                                       std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_matrix(*p_tree, names_b, matrix_b, 
                                         std::back_inserter(samples_b), 
                                         std::back_inserter(min_max_b) );

    return _matrix_query_internal_averaged( samples_a, min_max_a, samples_b, 
                                            min_max_b, is_double_matrix, ot );
  }


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _matrix_distance_query_averaged_specific_pairs( const std::vector<std::string> &names_a, 
                                                  const std::vector<std::vector<bool> > &matrix_a,
                                                  const std::vector<std::string> &names_b, 
                                                  const std::vector<std::vector<bool> > &matrix_b, 
                                                  const std::vector<std::pair<int, int> > &queries, 
                                                  OutputIterator ot)
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (&matrix_a != &matrix_b);

    this->_extract_samples_from_matrix(*p_tree, names_a, matrix_a, std::back_inserter(samples_a), 
                                       std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_matrix(*p_tree,  names_b, matrix_b, std::back_inserter(samples_b), 
                                         std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    for(int i=0; i<queries.size(); i++)
    {
      query_intervals_a.push_back(std::make_pair(queries[i].first, queries[i].first));
      query_intervals_b.push_back(std::make_pair(queries[i].second, queries[i].second));
    } 

    return   _matrix_query_internal_averaged_specific_pairs( samples_a, min_max_a, 
                                                             samples_b, min_max_b, 
                                                             query_intervals_a, query_intervals_b,
                                                             is_double_matrix, ot);

  } // _matrix_distance_query_averaged_specific_pairs(...)


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _csv_matrix_distance_query_averaged( char *matrix_filename_a, char *matrix_filename_b, OutputIterator ot )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(matrix_filename_a) != std::string(matrix_filename_b));

    this->_extract_samples_from_file(*p_tree, matrix_filename_a,
                                     std::back_inserter(samples_a), std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_file(*p_tree, matrix_filename_b, 
                                       std::back_inserter(samples_b), 
                                       std::back_inserter(min_max_b) );


    return _matrix_query_internal_averaged( samples_a, min_max_a, samples_b, 
                                            min_max_b, is_double_matrix, ot );
  }


  template< class KernelType >
  template< class OutputIterator >
  std::pair<int, int>
  PhylogeneticMeasures::Community_distance_nearest_taxon<KernelType>::
  _csv_matrix_distance_query_averaged_specific_pairs( char *matrix_filename_a, char *matrix_filename_b, 
                                                      char *queries_filename, OutputIterator ot )
  {
    std::vector< std::vector<int> >  samples_a, samples_b;
    std::vector< std::pair<int,int> > min_max_a, min_max_b;

    bool is_double_matrix = (std::string(matrix_filename_a) != std::string(matrix_filename_b));

    this->_extract_samples_from_file(*p_tree, matrix_filename_a, std::back_inserter(samples_a), 
                                     std::back_inserter(min_max_a) );

    if(is_double_matrix)
      this->_extract_samples_from_file(*p_tree, matrix_filename_b, std::back_inserter(samples_b), 
                                       std::back_inserter(min_max_b) );

    std::vector< std::pair<int,int> > query_intervals_a, query_intervals_b;

    this->_extract_queries_from_file( queries_filename, 
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

    return   _matrix_query_internal_averaged_specific_pairs( samples_a, min_max_a, 
                                                             samples_b, min_max_b, 
                                                             query_intervals_a, query_intervals_b,
                                                             is_double_matrix, ot);

  } // _csv_matrix_distance_query_averaged_specific_pairs(...)

#endif //COMMUNITY_DISTANCE_NEAREST_TAXON_IMPL_H
