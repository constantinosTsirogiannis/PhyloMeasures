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

#ifndef UNIQUE_FRACTION_IMPL_H
#define UNIQUE_FRACTION_IMPL_H

  template< class KernelType >
  template < class RangeIterator >
  typename KernelType::Number_type PhylogeneticMeasures::Unique_fraction<KernelType>::
  operator()( RangeIterator rbegin_a, RangeIterator rend_a,
              RangeIterator rbegin_b, RangeIterator rend_b,
              int min_index_a, int max_index_a,
              int min_index_b, int max_index_b )
  {
    if(p_tree->number_of_nodes()<2)
      return Number_type(1.0);

    if(rend_a - rbegin_a == 0 || rend_b - rbegin_b == 0)
      return Number_type(1.0);

    int min_index = std::min(min_index_a, min_index_b),
        max_index = std::max(max_index_a, max_index_b);     

    Number_type union_dist(0.0), intersection_dist(0.0);

    int intersection_index = p_tree->compute_intersection_node_index(min_index, max_index);

    if(p_tree->node(intersection_index).number_of_children() == 0)
      return Number_type(1.0);

    p_tree->node(intersection_index).mark = true;
    p_tree->node(intersection_index).mark_b = true;


    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
    {
      int index = *rit;

      Node_type v = p_tree->node(index);

      p_tree->node(index).mark = true;

      while( (!p_tree->is_root(v)) )
      {
        union_dist = union_dist + Number_type(p_tree->node(index).distance);
        p_tree->node(v.parent).marked_children.push_back(index);

        if( p_tree->node(v.parent).mark == false )
        {
          p_tree->node(v.parent).mark = true;
          index = v.parent;
          v = p_tree->node(v.parent);
        }
        else
          break;

      } // while( (!this->is_root(v)) )

    } // for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )


    // Mark the nodes of the Steiner tree of sample B.

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
    {
      int index = *rit;

      Node_type v = p_tree->node(index);

      p_tree->node(index).mark_b = true;

      while( (!p_tree->is_root(v)) )
      {
        p_tree->node(v.parent).marked_children_b.push_back(index);

        if(p_tree->node(index).mark == false)
          union_dist = union_dist + Number_type(p_tree->node(index).distance);
        else
          intersection_dist = intersection_dist + Number_type(p_tree->node(index).distance);

        if( p_tree->node(v.parent).mark_b == false )
        {
          p_tree->node(v.parent).mark_b = true;
          index = v.parent;
          v = p_tree->node(v.parent);
        }
        else
          break;

      } // while( (!this->is_root(v)) )

    } // for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )

    p_tree->unmark_Steiner_trees_of_samples(rbegin_a, rend_a, rbegin_b, rend_b);

    return (union_dist - intersection_dist)/union_dist;

  } // Unique_fraction<KernelType>::operator()(...)

#endif // UNIQUE_FRACTION_IMPL_H
