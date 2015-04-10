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

#ifndef PHYLOGENETIC_DIVERSITY_IMPL_H
#define PHYLOGENETIC_DIVERSITY_IMPL_H

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
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

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the sample and the node indicated by intersection_index.
    // For each edge that we traverse, add the corresponding weight to the solution.

    Number_type total_dist = p_tree->mark_Steiner_tree_of_sample(rbegin, rend);


    // Unmark all marked nodes.

    p_tree->unmark_Steiner_tree_of_sample(rbegin, rend);

    return total_dist;
  } // Phylogenetic_diversity<KernelType>::operator()(...)


  template< class KernelType >
  template< class OutputIterator >
  void PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
  _compute_subtree_sums( int index, Number_type &subtree_distances,
                         Number_type &h_products, OutputIterator ot,
                         Number_type &sum_subtree, Number_type &sum_subtract )
  {
    Node_type v = p_tree->node(index);

    std::vector<int> subtree_leaves;

    for( int i=0; i<v.children.size(); i++ )
    {
      Number_type sub_dist(0.0), h_prod(0.0);
      std::vector< std::pair<Number_type,int> > c_subtree_leaves;

      _compute_subtree_sums( v.children[i], sub_dist, h_prod,
                             std::back_inserter(c_subtree_leaves),
                             sum_subtree, sum_subtract );

      sum_subtree  += (Number_type(v.distance)*Number_type(sub_dist))
                      -(Number_type(v.distance)*Number_type(sub_dist)*hypergeom(v.all_subtree_leaves))
                      -(Number_type(v.distance)*h_prod);

      sum_subtract += (Number_type(v.distance)*Number_type(sub_dist))
                      -(Number_type(v.distance)*Number_type(sub_dist)*hypergeom(_number_of_leaves - v.all_subtree_leaves))
                      -(Number_type(v.distance)*h_prod);

      for( int j=0; j< c_subtree_leaves.size(); j++ )
      {
        sum_subtree  += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom(v.all_subtree_leaves-c_subtree_leaves[j].second);

        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second);

        *ot++ = c_subtree_leaves[j];
      }

      subtree_distances += sub_dist;
      h_products += h_prod;
    }

    sum_subtree  += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0) -
                    two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR) );

    sum_subtract += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0) -
                    two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT) );

    subtree_distances += Number_type(v.distance);
    h_products += Number_type(v.distance)*hypergeom(_number_of_leaves-v.all_subtree_leaves);
    *ot++ = std::make_pair(Number_type(v.distance), v.all_subtree_leaves);

  } // _compute_subtree_sums( int index, ... )



  template< class KernelType >
  void PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
  _compute_subtree_sums(Number_type &sum_subtree, Number_type &sum_subtract)
  {
    Node_type root = p_tree->root();

    for( int i = 0; i < root.children.size(); i++ )
    {
      std::vector< std::pair<Number_type, int> > subtree_leaves;
      Number_type subtree_distances(0.0), h_products(0.0);

      _compute_subtree_sums( root.children[i], subtree_distances, h_products,
                             std::back_inserter(subtree_leaves), sum_subtree, sum_subtract );

      subtree_leaves.clear();
    }

  } // _compute_subtree_sums( ... )


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
  compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves = number_of_leaves;

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
  typename KernelType::Number_type PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
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

    if( sample_size <= 1)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());
    compute_all_hypergeometric_probabilities( sample_size, p_tree->number_of_leaves());

    Number_type expectation(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ ) // We exclude the edge of the root node.
      expectation += Number_type(p_tree->node(i).distance)*(Number_type(1.0) 
                     - hypergeom(p_tree->node(i).all_subtree_leaves)
                     - hypergeom(p_tree->number_of_leaves() - p_tree->node(i).all_subtree_leaves) );

    return expectation;
	
  } //compute_expectation( int sample_size )

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
  compute_variance( int sample_size)
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

    if(sample_size <= 1)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    Number_type exp = compute_expectation( sample_size );

    Number_type sum_subtree(0.0), sum_subtract(0.0), sum_third_case(0.0),
                sum_self(0.0), sum_self_third_case(0.0), sum_same_class_third_case(0.0), total_sum(0.0);

    _compute_subtree_sums(sum_subtree, sum_subtract);

    std::map< int, Number_type > edge_classes;

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type v = p_tree->node(i);

      if( edge_classes.find(v.all_subtree_leaves) == edge_classes.end() )
        edge_classes[v.all_subtree_leaves] = Number_type(v.distance);
      else
        edge_classes[v.all_subtree_leaves] += Number_type(v.distance);

      sum_self += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0) -
                 two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR) );

      sum_self_third_case += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0) -
                             two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT) );
    }

    typename std::map< int, Number_type >::iterator cit1, cit2;

    for(cit1 = edge_classes.begin(); cit1 != edge_classes.end(); cit1++ )
    {
      sum_same_class_third_case += cit1->second*cit1->second*( Number_type(1.0) -
                                   two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT) );

      for(cit2 = edge_classes.begin(); cit2 != cit1; cit2++ )
        sum_third_case += cit1->second*cit2->second*( Number_type(1.0) -
                          two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT) );
    }

    sum_same_class_third_case = ((sum_same_class_third_case - sum_self_third_case)/Number_type(2.0))
                                  + sum_self_third_case;

    sum_third_case += sum_same_class_third_case;

    total_sum = (Number_type(2.0)*(sum_third_case - sum_subtract + sum_subtree)) - sum_self;

    return total_sum-(exp*exp);

  } // compute_variance( ... )
  

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Phylogenetic_diversity<KernelType>::
  compute_variance_slow( int sample_size )
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

    if(sample_size <= 1)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(*p_tree, p_tree->root_index());

    Number_type exp = compute_expectation( sample_size );

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
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0) -
                                         two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, Kernel::OFFSPRING) );
        else if( i <= j && i >= j - p_tree->subtree_edges(j) )
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0) -
                                         two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR) );
        else
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0) -
                                         two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT) );
      }

    } // for( int i=0; i<p_tree->number_of_nodes()-1; i++ )

    return sum-(exp*exp);

  } // compute_variance_slow( ... )

#endif //PHYLOGENETIC_DIVERSITY_IMPL_H
