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

#ifndef MEAN_NEAREST_TAXON_DISTANCE_IMPL_H
#define MEAN_NEAREST_TAXON_DISTANCE_IMPL_H

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf in its subtree.
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_subtree_min_values( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Number_type cmin = _compute_subtree_min_values(tree, current_node.marked_children[i]);

      if( current_node.first_min == Number_type(-1.0) || cmin < current_node.first_min)
      {
        tree.node(current_index).second_min = tree.node(current_index).first_min;
        tree.node(current_index).first_min = cmin;
      }
      else if(current_node.second_min == Number_type(-1.0) || cmin < current_node.second_min )
        tree.node(current_index).second_min = cmin;

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

    if( current_node.marked_children.size() == 0 )
    {
      tree.node(current_index).first_min = Number_type(0.0);
      tree.node(current_index).second_min = Number_type(0.0);
    }

    return tree.node(current_index).first_min + Number_type(tree.node(current_index).distance);

  } // _compute_subtree_min_values( ... )

  // The following function computes the two smallest costs of paths
  // that connect node with index `current_index' to a leaf outside its subtree.
  // This function produces a meaningful result only if it is ran after a call
  // of the function _compute_subtree_min_values().
  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_rest_tree_min_values( Tree_type &tree, int current_index )
  {
    Node_type current_node = tree.node(current_index);

    int first_min_index = -1, second_min_index = -1;
    Number_type current_first_min(-1.0), current_second_min(-1.0);

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      Node_type child = tree.node(current_node.marked_children[i]);
      int child_index = current_node.marked_children[i];      

      if(first_min_index == -1 || child.first_min + Number_type(child.distance) < current_first_min)
      {
        second_min_index = first_min_index;
        current_second_min = current_first_min; 
        first_min_index = child_index;
        current_first_min = child.first_min + Number_type(child.distance);
      }
      else
      if(second_min_index == -1 || child.first_min + Number_type(child.distance) < current_second_min)
      {
        second_min_index = child_index;
        current_second_min = child.first_min + Number_type(child.distance);
      }

    } // for( int i=0; i<current_node.marked_children.size(); i++ )   

    for( int i=0; i<current_node.marked_children.size(); i++ )
    {
      int child_index = current_node.marked_children[i];
      Node_type child = tree.node(current_node.marked_children[i]);

      Number_type cmin = std::min(current_first_min, current_node.rest_tree_min);

      if(current_node.rest_tree_min == Number_type(-1.0))
        cmin = current_first_min;

      if( child_index == first_min_index)
      { 
        if(current_second_min == Number_type(-1.0) || 
           (current_node.rest_tree_min < current_second_min && current_node.rest_tree_min != Number_type(-1.0)) )
          tree.node(current_node.marked_children[i]).rest_tree_min = current_node.rest_tree_min 
                                                                     + Number_type(child.distance);   
        else
          tree.node(current_node.marked_children[i]).rest_tree_min = current_second_min + Number_type(child.distance);   
      }
      else
        tree.node(current_node.marked_children[i]).rest_tree_min = cmin + Number_type(child.distance);

      _compute_rest_tree_min_values(tree, current_node.marked_children[i]);

    } // for( int i=0; i<current_node.marked_children.size(); i++ )

  } // _compute_rest_tree_min_values( ... )


  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index )
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    if(int(rend-rbegin) < 2)
      return Number_type(0.0);

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    // If the two paths coincide (intersection_index designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index).mark = true;

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);

    _compute_subtree_min_values( *p_tree, intersection_index );
    _compute_rest_tree_min_values( *p_tree, intersection_index );

    Number_type total_dist(0.0);

    // Unmark all marked nodes.

    for( RangeIterator rit = rbegin; rit != rend; rit++ )
    {
      total_dist = total_dist + p_tree->node((*rit)).rest_tree_min;

      p_tree->node(*rit).mark = false;
      p_tree->node(*rit).marked_children.clear();
      p_tree->node(*rit).first_min = Number_type(-1.0);
      p_tree->node(*rit).second_min = Number_type(-1.0);
      p_tree->node(*rit).rest_tree_min = Number_type(-1.0);
      Node_type v = p_tree->node(*rit);

      while( (!p_tree->is_root(v)) && p_tree->node(v.parent).mark == true )
      {
        p_tree->node(v.parent).mark = false;
        p_tree->node(v.parent).marked_children.clear();
        p_tree->node(v.parent).first_min = Number_type(-1.0);
        p_tree->node(v.parent).second_min = Number_type(-1.0);
        p_tree->node(v.parent).rest_tree_min = Number_type(-1.0);
        v = p_tree->node(v.parent);
      }
    }

    return total_dist/Number_type(rend - rbegin);

  } // Mean_nearest_taxon_distance<KernelType>::operator()(...)


  template< class KernelType >
  template< class OutputIterator >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  _compute_subtree_sums( int index, Number_type& sum_of_products, OutputIterator ot,
                         Number_type &sum_subtree, Number_type &sum_subtract )
  {
    Node_type v = p_tree->node(index);

    Number_type mhyperg = this->hypergeom_minus_one(_number_of_leaves-v.all_subtree_leaves);

    std::vector<int> subtree_leaves;

    for( int i=0; i<v.children.size(); i++ )
    {
      Number_type sum_of_pr(0.0);
      std::vector< std::pair<Number_type,int> > c_subtree_leaves;

      _compute_subtree_sums( v.children[i], sum_of_pr,
                             std::back_inserter(c_subtree_leaves),
                             sum_subtree, sum_subtract );

      sum_subtree  += Number_type(v.distance)*sum_of_pr*mhyperg;

      for( int j=0; j< c_subtree_leaves.size(); j++ )
      {
        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        Number_type(v.all_subtree_leaves)*c_subtree_leaves[j].second*
                        hypergeom_minus_two(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second);

        *ot++ = c_subtree_leaves[j];
      }

      sum_of_products += sum_of_pr;
    }

    sum_subtree  += Number_type(v.distance)*Number_type(v.distance)*Number_type(v.all_subtree_leaves)*mhyperg;

    sum_subtract += Number_type(v.distance)*Number_type(v.distance)*
                    Number_type(v.all_subtree_leaves)*Number_type(v.all_subtree_leaves)*
                    two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT);

    sum_of_products += Number_type(v.distance)*Number_type(v.all_subtree_leaves);
    *ot++ = std::make_pair(Number_type(v.distance), v.all_subtree_leaves);

  } // _compute_subtree_sums( int index, ... )


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves    = number_of_leaves;

    if( !_hypergeom.empty() )
      _hypergeom.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i= _number_of_leaves-1; i>=_sample_size; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x-Number_type(_sample_size)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom.push_back( tempgeom[i] );
	  
  } // compute_all_hypergeometric_probabilities(...)


  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD expectation on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size <= 1)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());
    compute_all_hypergeometric_probabilities( sample_size, p_tree->number_of_leaves());

    Number_type expectation(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ ) // We exclude the edge of the root node.
        expectation += Number_type(p_tree->node(i).distance)*Number_type(p_tree->node(i).all_subtree_leaves)*
                       hypergeom_minus_one(p_tree->number_of_leaves()-p_tree->node(i).all_subtree_leaves);

    return expectation*Number_type(2.0)/Number_type(sample_size);
	
  } // compute_expectation( int sample_size )

  
  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_variance( int sample_size, Number_type expect )
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD variance on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(sample_size <= 1)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    Number_type exp;

    if( expect != Number_type(-1.0) )
      exp = expect;
    else
      exp = compute_expectation( sample_size );

    Number_type sum_subtree(0.0), sum_subtract(0.0), sum_third_case(0.0),
                sum_self(0.0), sum_self_third_case(0.0), sum_same_class_third_case(0.0), total_sum(0.0);

    _compute_subtree_sums(sum_subtree, sum_subtract);

    std::map< int, Number_type > edge_classes;

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type v = p_tree->node(i);

      if( edge_classes.find(v.all_subtree_leaves) == edge_classes.end() )
        edge_classes[v.all_subtree_leaves] = Number_type(v.distance)*Number_type(v.all_subtree_leaves);
      else
        edge_classes[v.all_subtree_leaves] += Number_type(v.distance)*Number_type(v.all_subtree_leaves);

      sum_self += Number_type(v.distance)*Number_type(v.distance)*Number_type(v.all_subtree_leaves)*
                 two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR);

      sum_self_third_case += Number_type(v.distance)*Number_type(v.distance)*
                             Number_type(v.all_subtree_leaves)*Number_type(v.all_subtree_leaves)*
                             two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT);
    }

    typename std::map< int, Number_type >::iterator cit1, cit2;

    for(cit1 = edge_classes.begin(); cit1 != edge_classes.end(); cit1++ )
    {
      sum_same_class_third_case += cit1->second*cit1->second*
                                   two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT);



      for(cit2 = edge_classes.begin(); cit2 != cit1; cit2++ )
        sum_third_case += cit1->second*cit2->second*
                          two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT);
    }

    sum_same_class_third_case = ((sum_same_class_third_case - sum_self_third_case)/Number_type(2.0))
                                  + sum_self_third_case;

    sum_third_case += sum_same_class_third_case;

    total_sum = (Number_type(2.0)*(sum_third_case - sum_subtract + sum_subtree)) - sum_self;
    total_sum = Number_type(4.0)*total_sum/(Number_type(_sample_size)*Number_type(_sample_size));

    return total_sum-(exp*exp);

  } //compute_variance( ... )

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Mean_nearest_taxon_distance<KernelType>::
  compute_variance_slow( int sample_size, Number_type expect )
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

    if(!p_tree->is_ultrametric())
    {
      std::string exception_msg;
      exception_msg += " Request to compute MNTD variance on a non-ultrametric tree.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size <= 1 )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(*p_tree, p_tree->root_index());

    Number_type exp;

    if( expect != Number_type(-1.0) )
      exp = expect;
    else
      exp = compute_expectation( sample_size );

    if(p_tree->subtree_edges_size() == 0)
      p_tree->compute_subtree_edges(p_tree->root_index());

    Number_type sum(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type u = p_tree->node(i);

      for( int j=0; j<p_tree->number_of_nodes()-1; j++ )
      {
        Node_type v = p_tree->node(j);

        if( j < i && j >= i - p_tree->subtree_edges(i) )
          sum += Number_type(u.distance)*Number_type(v.distance)*
                 Number_type(v.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::OFFSPRING);
        else if( i <= j && i >= j - p_tree->subtree_edges(j) )
          sum += Number_type(u.distance)*Number_type(v.distance)*
                 Number_type(u.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::ANCESTOR);
        else
          sum += Number_type(u.distance)*Number_type(v.distance)*Number_type(u.all_subtree_leaves)*
                 Number_type(v.all_subtree_leaves)*two_edge_pr(u.all_subtree_leaves, 
                                                               v.all_subtree_leaves, Kernel::INDEPENDENT);

      } // for( int j=0; j<p_tree->number_of_nodes()-1; j++ )

    } // for( int i=0; i<p_tree->number_of_nodes()-1; i++ )


    return (Number_type(4.0)*sum/((sample_size)*(sample_size)))-(exp*exp);

  } // compute_variance_slow( ... )


#endif //MEAN_NEAREST_TAXON_DISTANCE_IMPL_H
