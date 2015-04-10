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

#ifndef CORE_ANCESTOR_COST_IMPL_H
#define CORE_ANCESTOR_COST_IMPL_H

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  operator()( RangeIterator rbegin, RangeIterator rend,
              int min_index, int max_index)
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(0.0);

    int number_of_input_leaves = rend-rbegin;
    Number_type temp_rchi = Ceiling()(this->chi()*Number_type(number_of_input_leaves)); 
    int rchi = int(To_double()(temp_rchi));

    if(rchi == 0)
	  return Number_type(0.0);

    // Mark all the nodes which fall on a path that connects
    // a leaf node in the input sample and the root node.

    p_tree->mark_Steiner_tree_of_sample(rbegin, rend);
    p_tree->assign_marked_subtree_leaves(p_tree->root_index() );

    // Starting from the root, traverse the tree to find
    // the core ancestor of the input leaf sample.

    Number_type total_dist = Number_type(0.0);
    Node_type v = p_tree->root();
    bool found = false;

    int last_heir = p_tree->root_index();

    do
    {
      if(v.number_of_marked_children() > 0)
      {
        int heir=-1;

        for(int i=0; i<v.number_of_marked_children(); i++)
          if( p_tree->node(v.marked_children[i]).marked_subtree_leaves >= rchi )
          {
            heir = v.marked_children[i];
	    break;
          }

	if(heir == -1)
	  found = true;
	else
	{
	  total_dist += Number_type(p_tree->node(heir).distance);
	  v = p_tree->node(heir);
	  last_heir = heir;
	}

      } // if(v.number_of_children() > 0)

    }while(found == false && v.number_of_children() > 0);

    p_tree->unmark_Steiner_tree_of_sample(rbegin, rend);

    return total_dist;

  } // Core_ancestor_cost<KernelType>::operator()(...)


  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  slow_operator( RangeIterator rbegin, RangeIterator rend,
                 int min_index, int max_index)
  {
    int number_of_input_leaves = rend-rbegin;
    Number_type temp_rchi = Ceiling()(this->chi()*Number_type(number_of_input_leaves));
	int rchi = int(To_double(temp_rchi));

    if(rchi == 0)
	  return Number_type(0.0);

    for(int i=0; i<p_tree->number_of_nodes(); i++)
    {
      p_tree->node(i).marked_subtree_leaves = 0;
      p_tree->node(i).mark= false;
    }

    for(RangeIterator rit=rbegin; rit != rend; rit++)
    {
      p_tree->node(*rit).mark=true;

      Node_type u = p_tree->node(*rit);

      while(u.parent >=0)
      {
        p_tree->node(u.parent).mark=true;
        u = p_tree->node(u.parent);
      }
    }

    p_tree->assign_marked_subtree_leaves(p_tree->root_index());

    int largest_index=-1;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      if(p_tree->node(i).marked_subtree_leaves >= rchi && i>largest_index  )
      {
        int new_index= -1;

        for(int j=0; j<p_tree->node(i).number_of_children(); j++)
        {
          int child = p_tree->node(i).children[j];

          if(p_tree->node(child).marked_subtree_leaves >= rchi)
            new_index = child;
        }

        if(new_index==-1)
          largest_index=i;

      } // if(p_tree->node(i).marked_subtree_leaves >= ... )

    Node_type v= p_tree->node(largest_index);
    Number_type dist(0.0);

    while(v.parent >=0)
    {
      dist = dist + Number_type(v.distance);
      v = p_tree->node(v.parent);
    }

    for(int i=0; i<p_tree->number_of_nodes(); i++)
    {
      p_tree->node(i).marked_subtree_leaves = 0;
      p_tree->node(i).mark= false;
    }

    return dist;

  } //slow_operator(...)


  // Computes all together the probability values f(x) = \binom{x}{r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_hypergeometric_probabilities_a( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves = number_of_leaves;

    if( !_hypergeom_a.empty() )
      _hypergeom_a.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i= _number_of_leaves-1; i>=_sample_size; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x-Number_type(_sample_size)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_a.push_back( tempgeom[i] );

  } //compute_all_hypergeometric_probabilities_a(...)


  // Computes all together the probability values f(x) = \binom{x}{s-r}/\binom{s}{r}

  template< class KernelType >
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_hypergeometric_probabilities_b( int sample_size, int number_of_leaves)
  {
    _sample_size = sample_size;
    _number_of_leaves = number_of_leaves;

    if( !_hypergeom_b.empty() )
      _hypergeom_b.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(1.0));

    for( int i=_number_of_leaves-1; i>=_number_of_leaves-_sample_size; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()*
                          (Number_type(_sample_size)+x-Number_type(_number_of_leaves))/x  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_b.push_back( tempgeom[i] );

  } //compute_all_hypergeometric_probabilities_b(...)


  template< class KernelType >
  typename KernelType::Number_type PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_node_probability(int number_of_subtree_leaves, int sample_size)
  {
    int number_of_all_leaves = p_tree->number_of_leaves();
    
    Number_type temp_rchi = Ceiling()(_chi*Number_type(sample_size));
    int leaf_fraction = int(To_double()(temp_rchi));

    if(leaf_fraction < 1)
      return Number_type(0.0);

    if(std::max(1,sample_size+number_of_subtree_leaves-number_of_all_leaves) > 
                  std::min(number_of_subtree_leaves,sample_size))
      return Number_type(0.0);

    if(leaf_fraction > number_of_subtree_leaves)
      return Number_type(0.0);

    if(number_of_subtree_leaves == number_of_all_leaves)
      return Number_type(1.0);

    Number_type basic_hypergeom, extra_hypergeom(0.0),
		probability_product(1.0),
		cumulative_probability(0.0);

    int i;

    if( number_of_all_leaves-sample_size-number_of_subtree_leaves >= 0 )
    {
      basic_hypergeom = hypergeom_a(number_of_all_leaves - number_of_subtree_leaves);
      i=1;
    }
    else
    {
      basic_hypergeom = hypergeom_b(number_of_subtree_leaves);
      i=sample_size+number_of_subtree_leaves-number_of_all_leaves+1;
    }

    if(i-1>=leaf_fraction)
      extra_hypergeom = Number_type(1.0);

    for( /*We settled i already*/ ; i<=std::min(number_of_subtree_leaves,sample_size); i++ )
    {
      Number_type f1(sample_size-i+1), f2(number_of_subtree_leaves-i+1),
	          f3(std::max(number_of_all_leaves-number_of_subtree_leaves-sample_size+i,1)), 
                  f4(i), numerator, denominator;

      numerator = f1*f2;
      denominator = f3*f4;

      probability_product = probability_product*(numerator/denominator);

      if(i>=leaf_fraction)
        cumulative_probability+=probability_product;
    }

    return basic_hypergeom*(cumulative_probability+extra_hypergeom);

  } // compute_node_probability(...)


  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_all_root_path_costs(OutputIterator ot)
  {
    std::vector<Number_type> path_costs;
    path_costs.assign( p_tree->number_of_nodes() , Number_type(0.0) );

    std::queue< std::pair<int,Number_type> > node_cost_pairs;
	node_cost_pairs.push( std::make_pair(p_tree->root_index(), Number_type(0.0)) );

    while(!node_cost_pairs.empty())
    {
      int index = node_cost_pairs.front().first;
	  Number_type cost = node_cost_pairs.front().second;

      node_cost_pairs.pop();
      path_costs[index] = cost;

      Node_type v = p_tree->node(index);

      for(int i=0; i<v.number_of_children(); i++)
      {
	int child_index = v.children[i];
        Number_type new_cost = cost + Number_type(p_tree->node(child_index).distance);
        node_cost_pairs.push(std::make_pair(child_index,new_cost));
      }

    } // while(!node_indices.empty())

    for(int i=0; i<path_costs.size(); i++)
      *ot++ = path_costs[i];

  } // compute_all_root_path_costs(OutputIterator ot)


  template< class KernelType >
  template<class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_first_k_raw_moments(int k, int sample_size, OutputIterator ot)
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute moments with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(p_tree->number_of_leaves()<=1 || sample_size == 0)
    {
      for(int i=0; i<k; i++)
        *ot++ = Number_type(0.0);

      return;
    }

    std::vector<Number_type> root_path_costs, subtree_size_probabilities,
	                     node_probabilities, costs_kth_exponent;

    p_tree->assign_all_subtree_leaves(p_tree->root_index());

    this->compute_all_hypergeometric_probabilities_a(sample_size, p_tree->number_of_leaves());
    this->compute_all_hypergeometric_probabilities_b(sample_size, p_tree->number_of_leaves());
    this->compute_all_root_path_costs(std::back_inserter(root_path_costs));
    subtree_size_probabilities.assign(p_tree->number_of_leaves(), Number_type(0.0));
    node_probabilities.assign(p_tree->number_of_nodes(), Number_type(0.0));
    costs_kth_exponent.assign(p_tree->number_of_nodes(), Number_type(1.0));

    std::set<int> unique_subtree_sizes;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      unique_subtree_sizes.insert(p_tree->node(i).all_subtree_leaves);

    typename std::set<int>::iterator sit;

    for( sit=unique_subtree_sizes.begin(); sit!=unique_subtree_sizes.end(); sit++)
      subtree_size_probabilities[(*sit)-1] = this->compute_node_probability(*sit,sample_size);

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      node_probabilities[i] = subtree_size_probabilities[p_tree->node(i).all_subtree_leaves-1];

    std::vector<Number_type> final_probabilities = node_probabilities;

    for(int i=0; i<p_tree->number_of_nodes(); i++)
      for( int j=0; j<p_tree->node(i).number_of_children(); j++ )
        final_probabilities[i] = final_probabilities[i] - node_probabilities[p_tree->node(i).children[j]];

    for(int i=0; i<k; i++)
    {
      Number_type ith_raw_moment(0.0);

      for(int j=0; j < p_tree->number_of_nodes(); j++)
      {
        costs_kth_exponent[j] = costs_kth_exponent[j]*root_path_costs[j];
        ith_raw_moment += costs_kth_exponent[j]*final_probabilities[j];
      }

      *ot++ = ith_raw_moment;

    } // for(int i=0; i<k; i++)

  } //compute_first_k_raw_moments(int k, OutputIterator ot)


  template< class KernelType >
  template <class OutputIterator>
  void PhylogeneticMeasures::Core_ancestor_cost<KernelType>::
  compute_first_k_centralised_moments( int k, int sample_size, OutputIterator ot)
  {
    if(sample_size < 0 || sample_size > p_tree->number_of_leaves())
    {
      std::string exception_msg;
      exception_msg += " Request to compute moments with sample size which is out of range.\n";     
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(p_tree->number_of_leaves()<=1 || sample_size == 0)
    {
      for(int i=0; i<k; i++)
        *ot++ = Number_type(0.0);

      return;
    }

    if( k <= 0 )
    {
      std::string msg(" Invalid value of parameter k. Value of k must be > 0 .\n");
      Exception_type excp;
      excp.get_error_message(msg);
      Exception_functor excf;
      excf(excp);
    }

    if(k==1)
    {
      std::vector<Number_type>  expec;
      compute_first_k_raw_moments(1, sample_size, ot);
      return;
    }

    // Check which fixes a case where a floating point error occurs.
    bool ultrametric = p_tree->is_ultrametric();
    
    if(k==2)
    {
      std::vector<Number_type>  raw_moments;
      this->compute_first_k_raw_moments(2, sample_size, std::back_inserter(raw_moments));

      *ot++ = raw_moments[0];

      if(sample_size != 1 || ultrametric == false)
        *ot++ = raw_moments[1]-(raw_moments[0]*raw_moments[0]);
      else    
        *ot++ = Number_type(0.0);         

      return;
    }

    std::vector<Number_type>  raw_moments;

    this->compute_first_k_raw_moments(k, sample_size, std::back_inserter(raw_moments));

    Number_type variance = raw_moments[1]-(raw_moments[0]*raw_moments[0]), variance_powers(1.0),
	        minus_mean(Number_type(-1.0)*raw_moments[0]), result(0.0);

    if(sample_size == 1 && ultrametric == true)
      variance = Number_type(0.0);

    *ot++ = raw_moments[0];
    *ot++ = variance;

    std::vector<Number_type> minus_mean_powers;
    std::vector<Number_type> binom_coeff;

    minus_mean_powers.push_back(Number_type(1.0));

    for(int i=0; i<k; i++)
      minus_mean_powers.push_back(minus_mean_powers.back()*minus_mean);

    variance_powers= variance*variance;

    for(int i=3; i<=k; i++)
    {
      std::vector<Number_type> coefficients;
      this->compute_k_binomial_coefficients(i,std::back_inserter(coefficients));
        
      result=Number_type(0.0);

      for(int j=0; j<i; j++)
        result = result + (minus_mean_powers[j]*(coefficients[j]*raw_moments[i-j-1]));

      result = result + (minus_mean_powers[i]*coefficients[i]);

      variance_powers = variance_powers*variance;
 
      if(variance != Number_type(0.0))
        *ot++ = result;
      else
        *ot++ = Number_type(0.0);
    }

  } // _compute_first_k_centralised_moments( int k, OutputIterator ot)


#endif //CORE_ANCESTOR_COST_IMPL_H
