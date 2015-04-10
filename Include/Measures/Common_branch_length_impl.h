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

#ifndef COMMON_BRANCH_LENGTH_IMPL_H
#define COMMON_BRANCH_LENGTH_IMPL_H

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Common_branch_length<KernelType>::
  operator()( RangeIterator rbegin_a, RangeIterator rend_a,
              RangeIterator rbegin_b, RangeIterator rend_b,
              int min_index_a, int max_index_a,
              int min_index_b, int max_index_b )
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    if(rbegin_a == rend_a || rbegin_b == rend_b)
      return Number_type(0.0);

    int intersection_index_a = p_tree->compute_intersection_node_index(min_index_a, max_index_a),
        intersection_index_b = p_tree->compute_intersection_node_index(min_index_b, max_index_b);

    // Set A: If the extremal paths coincide (intersection_index_a designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index_a).children.size() == 0 )
      return Number_type(0.0);

    // Set B: If the extremal paths coincide (intersection_index_b designates a leaf node)
    // then return zero distance.
    if( p_tree->node(intersection_index_b).children.size() == 0 )
      return Number_type(0.0);

    p_tree->node(intersection_index_a).mark = true;
    p_tree->node(intersection_index_b).mark_b = true;


    // Mark all the nodes which fall on a path that connects a 
    // leaf node of set A and the node indicated by intersection_index.

    Number_type total_dist(0.0);

    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
    {
      int index = *rit;
      Node_type v = p_tree->node(index);

      p_tree->node(index).mark = true;

      while( (!p_tree->is_root(v)) )
      {
        p_tree->node(v.parent).marked_children.push_back(index);

        if(p_tree->node(v.parent).mark == false)
        {
          p_tree->node(v.parent).mark = true;
          v = p_tree->node(v.parent);
        }
        else
          break;

      } // while( (!p_tree->is_root(v)) )

    } // for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )

    // Mark all the nodes which fall on a path that connects a 
    // leaf node of set B and the node indicated by intersection_index.
    // For each edge which is marked for both sets, add the associated
    // weight to the solution.

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
    {
      int index = *rit;
      Node_type v = p_tree->node(index);

      p_tree->node(index).mark_b = true;

      if( index != p_tree->root_index() && p_tree->node(index).mark == true )
        total_dist += Number_type(v.distance);

      while( (!p_tree->is_root(v)) )
      {
        p_tree->node(v.parent).marked_children_b.push_back(index);

        if(p_tree->node(v.parent).mark_b == false)
        {
          int ind_new = v.parent;
          p_tree->node(v.parent).mark_b = true;
          v = p_tree->node(ind_new);

          if(v.mark == true && ind_new != intersection_index_a  )
            total_dist += Number_type(v.distance);
        }
        else
          break;

      } // while( (!p_tree->is_root(v)) )

    } // for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )


    // Unmark all marked nodes.
    p_tree->unmark_Steiner_trees_of_samples(rbegin_a, rend_a, rbegin_b, rend_b);

    return total_dist;
  } // Common_branch_length<KernelType>::operator()(...)

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Common_branch_length<KernelType>::
  compute_variance( int sample_size_a, int sample_size_b)
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

    if( sample_size_a < 2 || sample_size_b < 2 )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    Number_type exp = this->compute_expectation( sample_size_a, sample_size_b );

    Number_type sum_subtree(0.0), sum_subtract(0.0), sum_third_case(0.0),
                sum_self(0.0), sum_self_third_case(0.0), sum_same_class_third_case(0.0), total_sum(0.0);

    this->_compute_subtree_sums(sum_subtree, sum_subtract);


    std::map< int, Number_type > edge_classes;


    for( int i=0; i<p_tree->number_of_nodes()-1; i++ )
    {
      Node_type v = p_tree->node(i);


      if( edge_classes.find(v.all_subtree_leaves) == edge_classes.end() )
        edge_classes[v.all_subtree_leaves] = Number_type(v.distance);
      else
        edge_classes[v.all_subtree_leaves] += Number_type(v.distance);


      sum_self += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, true)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, false)
                    + ( two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, true)*
                        two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR,false)));


      sum_self_third_case += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, true)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, false)
                    + ( two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, true)*
                        two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT,false)));

    }

    typename std::map< int, Number_type >::iterator cit1, cit2;


    for(cit1 = edge_classes.begin(); cit1 != edge_classes.end(); cit1++ )
    {

      sum_same_class_third_case += Number_type(cit1->second)*Number_type(cit1->second)*( Number_type(1.0)
                    - two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT, true)
                    - two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT, false)
                    + ( two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT, true)*
                        two_edge_pr(cit1->first, cit1->first, Kernel::INDEPENDENT,false)));

      for(cit2 = edge_classes.begin(); cit2 != cit1; cit2++ )
        sum_third_case += Number_type(cit1->second)*Number_type(cit2->second)*( Number_type(1.0)
                    - two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT, true)
                    - two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT, false)
                    + ( two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT, true)*
                        two_edge_pr(cit1->first, cit2->first, Kernel::INDEPENDENT,false)));

    }


    sum_same_class_third_case = ((sum_same_class_third_case - sum_self_third_case)/Number_type(2.0))
                                  + sum_self_third_case;

    sum_third_case += sum_same_class_third_case;

    total_sum = (Number_type(2.0)*(sum_third_case - sum_subtract + sum_subtree)) - sum_self;

    return total_sum-(exp*exp);

  } // compute_variance( ... )


  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Common_branch_length<KernelType>::
  compute_variance_slow( int sample_size_a, int sample_size_b)
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

    if( sample_size_a < 2 || sample_size_b < 2 )
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    Number_type exp = this->compute_expectation( sample_size_a, sample_size_b );

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
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0)
                                         - two_edge_pr(u.all_subtree_leaves, 
                                                       v.all_subtree_leaves, Kernel::OFFSPRING, true)
                                         - two_edge_pr(u.all_subtree_leaves, 
                                                       v.all_subtree_leaves, Kernel::OFFSPRING, false)
                                         + ( two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                             Kernel::OFFSPRING, true)
                                        *two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                                     Kernel::OFFSPRING, false)));
        else if( i <= j && i >= j - p_tree->subtree_edges(j) )
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0)
                                         - two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, true)
                                         - two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, false)
                                         + ( two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                             Kernel::ANCESTOR, true)*
                                           two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                                       Kernel::ANCESTOR,false)));
        else
          sum += Number_type(u.distance)*Number_type(v.distance)*( Number_type(1.0)
                                         - two_edge_pr(u.all_subtree_leaves, 
                                                       v.all_subtree_leaves, Kernel::INDEPENDENT,true)
                                         - two_edge_pr(u.all_subtree_leaves, 
                                                       v.all_subtree_leaves, Kernel::INDEPENDENT,false)
                                         + (two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                            Kernel::INDEPENDENT,true)
                                       *two_edge_pr(u.all_subtree_leaves, v.all_subtree_leaves, 
                                                    Kernel::INDEPENDENT,false)));
      } // for( int j=0; ... )
    } // for( int i=0; ... )

    return sum-(exp*exp);

  } // compute_variance_slow( ... )

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Common_branch_length<KernelType>::
  two_edge_pr( int se, int sl , Edge_relation_type er, bool is_a )
  {

    if( is_a )
      switch(er)
      {
        case Kernel::OFFSPRING:    return hypergeom_a(se)+hypergeom_a(_number_of_leaves-sl)-hypergeom_a(se-sl);
        case Kernel::ANCESTOR:     return hypergeom_a(sl)+hypergeom_a(_number_of_leaves-se)-hypergeom_a(sl-se);
        case Kernel::INDEPENDENT:  return hypergeom_a(_number_of_leaves-se)+hypergeom_a(_number_of_leaves-sl)
                                          -hypergeom_a(_number_of_leaves-se-sl);
      }
    else
      switch(er)
      {
        case Kernel::OFFSPRING:    return hypergeom_b(se)+hypergeom_b(_number_of_leaves-sl)-hypergeom_b(se-sl);
        case Kernel::ANCESTOR:     return hypergeom_b(sl)+hypergeom_b(_number_of_leaves-se)-hypergeom_b(sl-se);
        case Kernel::INDEPENDENT:  return hypergeom_b(_number_of_leaves-se)+hypergeom_b(_number_of_leaves-sl)
                                          -hypergeom_b(_number_of_leaves-se-sl);
      }

    return Number_type(-1.0);
	
  } // two_edge_pr( int se, int sl , Edge_relation_type er, bool is_a )

  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Common_branch_length<KernelType>::
  compute_expectation( int sample_size_a, int sample_size_b )
  {
    if( sample_size_a < 0 || sample_size_a > p_tree->number_of_leaves() ||
        sample_size_b < 0 || sample_size_b > p_tree->number_of_leaves()   )
    {
      std::string exception_msg;
      exception_msg += " Request to compute expectation with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if( sample_size_a < 2 || sample_size_b < 2)
      return Number_type(0.0);

    p_tree->assign_all_subtree_leaves(p_tree->root_index());
    compute_all_hypergeometric_probabilities( sample_size_a, p_tree->number_of_leaves(), _hypergeom_a,true);
    compute_all_hypergeometric_probabilities( sample_size_b, p_tree->number_of_leaves(), _hypergeom_b,false);

    Number_type expectation(0.0);

    for( int i=0; i<p_tree->number_of_nodes()-1; i++ ) // We exclude the edge of the root node.
    {
      Number_type a1 = hypergeom_a(p_tree->node(i).all_subtree_leaves),
                  a2 = hypergeom_a(p_tree->number_of_leaves() - p_tree->node(i).all_subtree_leaves),
                  b1 = hypergeom_b(p_tree->node(i).all_subtree_leaves),
                  b2 = hypergeom_b(p_tree->number_of_leaves() - p_tree->node(i).all_subtree_leaves);

      expectation += Number_type(p_tree->node(i).distance)*(Number_type(1.0)-a1-a2-b1-b2+(a1+a2)*(b1+b2));
    }

    return expectation;
	
  } // compute_expectation( int sample_size_a, int sample_size_b )


  template< class KernelType >
  void PhylogeneticMeasures::Common_branch_length<KernelType>::
  compute_all_hypergeometric_probabilities( int sample_size, int number_of_leaves,
                                            std::vector<Number_type> &hypergeom, bool is_a )
  {
    _number_of_leaves = number_of_leaves;

    if(is_a)
      _sample_size_a = sample_size;
    else
      _sample_size_b = sample_size;

    if( !hypergeom.empty() )
      hypergeom.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i= _number_of_leaves-1; i>=sample_size; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x-Number_type(sample_size)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      hypergeom.push_back( tempgeom[i] );
 
  } // compute_all_hypergeometric_probabilities(...)
  

  template< class KernelType >
  template< class OutputIterator >
  void PhylogeneticMeasures::Common_branch_length<KernelType>::
  _compute_subtree_sums( int &index, Number_type &subtree_distances,
                         Number_type &h_products_a, Number_type &h_products_b,
                         Number_type &double_products, OutputIterator ot,
                         Number_type &sum_subtree, Number_type &sum_subtract )
  {
    Node_type v = p_tree->node(index);

    std::vector<int> subtree_leaves;

    Number_type hpgm_a = hypergeom_a(v.all_subtree_leaves), 
                hpgm_b = hypergeom_b(v.all_subtree_leaves),
                hpgm_a_minus = hypergeom_a(_number_of_leaves - v.all_subtree_leaves),
                hpgm_b_minus = hypergeom_b(_number_of_leaves - v.all_subtree_leaves);

    for( int i=0; i<v.children.size(); i++ )
    {
      Number_type sub_dist(0.0), h_prod_a(0.0), h_prod_b(0.0), double_prod(0.0);
      std::vector< std::pair<Number_type,int> > c_subtree_leaves;

      _compute_subtree_sums( v.children[i], sub_dist, h_prod_a, h_prod_b,
                             double_prod, std::back_inserter(c_subtree_leaves),
                             sum_subtree, sum_subtract );


      sum_subtree  += (Number_type(v.distance)*sub_dist)-
                      (Number_type(v.distance)*sub_dist*hpgm_a)
                      -(Number_type(v.distance)*h_prod_a);

      sum_subtree  += -(Number_type(v.distance)*sub_dist*hpgm_b)
                      -(Number_type(v.distance)*h_prod_b);

      sum_subtree  += (Number_type(v.distance)*sub_dist*hpgm_a*hpgm_b) +
                      (Number_type(v.distance)*hpgm_a*h_prod_b)+
                      (Number_type(v.distance)*h_prod_a*hpgm_b) +
                      (Number_type(v.distance)*double_prod);


      sum_subtract += (Number_type(v.distance)*sub_dist)-
                      (Number_type(v.distance)*sub_dist*hpgm_a_minus)
                      -(Number_type(v.distance)*h_prod_a);

      sum_subtract += -(Number_type(v.distance)*sub_dist*hpgm_b_minus)
                      -(Number_type(v.distance)*h_prod_b);

      sum_subtract += (Number_type(v.distance)*sub_dist*hpgm_a_minus*hpgm_b_minus) +
                      (Number_type(v.distance)*hpgm_a_minus*h_prod_b)+
                      (Number_type(v.distance)*h_prod_a*hpgm_b_minus) +
                      (Number_type(v.distance)*double_prod);

      for( int j=0; j< c_subtree_leaves.size(); j++ )
      {

        sum_subtree  += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_a(v.all_subtree_leaves-c_subtree_leaves[j].second);

        sum_subtree  += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_b(v.all_subtree_leaves-c_subtree_leaves[j].second);

        sum_subtree  += Number_type(-1.0)*Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_b(v.all_subtree_leaves-c_subtree_leaves[j].second)*
                        (hpgm_a+hypergeom_a(_number_of_leaves-c_subtree_leaves[j].second));

        sum_subtree  += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_a(v.all_subtree_leaves-c_subtree_leaves[j].second)*
                        (hypergeom_b(v.all_subtree_leaves-c_subtree_leaves[j].second)
                         -hpgm_b-hypergeom_b(_number_of_leaves-c_subtree_leaves[j].second));

        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_a(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second);

        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_b(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second);

        sum_subtract += Number_type(-1.0)*Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_b(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second)*
                        (hpgm_a_minus+hypergeom_a(_number_of_leaves-c_subtree_leaves[j].second));

        sum_subtract += Number_type(v.distance)*c_subtree_leaves[j].first*
                        hypergeom_a(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second)*
                        (hypergeom_b(_number_of_leaves-v.all_subtree_leaves-c_subtree_leaves[j].second)
                         -hpgm_b_minus-hypergeom_b(_number_of_leaves-c_subtree_leaves[j].second));

        *ot++ = c_subtree_leaves[j];
      }

      subtree_distances += sub_dist;
      h_products_a += h_prod_a;
      h_products_b += h_prod_b;
      double_products += double_prod;
    }


    sum_subtree  += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, true)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, false)
                    + ( two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR, true)*
                        two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::ANCESTOR,false)));

    sum_subtract += Number_type(v.distance)*Number_type(v.distance)*( Number_type(1.0)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, true)
                    - two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, false)
                    + ( two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT, true)*
                        two_edge_pr(v.all_subtree_leaves, v.all_subtree_leaves, Kernel::INDEPENDENT,false)));

    subtree_distances += Number_type(v.distance);
    h_products_a += Number_type(v.distance)*hpgm_a_minus;
    h_products_b += Number_type(v.distance)*hpgm_b_minus;
    double_products += Number_type(v.distance)*hpgm_a_minus*hpgm_b_minus;

    *ot++ = std::make_pair(Number_type(v.distance), v.all_subtree_leaves);

  } // _compute_subtree_sums( int index, ... )




#endif // COMMON_BRANCH_LENGTH_IMPL_H
