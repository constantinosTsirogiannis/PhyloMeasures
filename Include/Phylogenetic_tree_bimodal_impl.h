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

#ifndef PHYLOGENETIC_TREE_BIMODAL_IMPL_H
#define PHYLOGENETIC_TREE_BIMODAL_IMPL_H

  // Consider the minimum cost Steiner tree that connects all leaves of a given sample. 
  // The following function computes for each node v of this tree the number of 
  // sample leaves situated at the subtree of v. This number is assigned to the field
  // 'marked_subtree_leaves' of v.

  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_bimodal<KernelType, NodeType>::
  _assign_marked_subtree_leaves_b( int index )
  {  
    this->node(index).marked_subtree_leaves_b=0;

    if( this->node(index).mark_b == false )
      return 0;

    _marked_nodes_b.push_back(index);

    // Leaf node
    if( this->node(index).number_of_children() == 0 )
      this->node(index).marked_subtree_leaves_b=1;
    else
    {
      Node_type v = this->node(index);

      // Interior node
      for( int i = 0; i < v.number_of_marked_children_b(); i++ )
        this->node(index).marked_subtree_leaves_b += this->_assign_marked_subtree_leaves_b(v.marked_children_b[i]);
    }

    return this->node(index).marked_subtree_leaves_b;

  }  // _assign_marked_subtree_leaves_b( int index )
  
  
  
  // Given a tree and two samples of its leaf nodes that are represented
  // by the iterator ranges [rbegin_a,rend_a) and [rbegin_b,rend_b), the 
  // following function marks all the tree nodes that fall on the two Steiner trees 
  // that connect the leaves of each sample. For this function to work appropriately, 
  // we should mark the roots of the Steiner trees in advance.

  template < class KernelType, class NodeType >
  template< class RangeIterator > 
  void PhylogeneticMeasures::Phylogenetic_tree_bimodal<KernelType, NodeType>::
  mark_Steiner_trees_of_samples( RangeIterator rbegin_a, RangeIterator rend_a, 
                                 RangeIterator rbegin_b, RangeIterator rend_b )
  {  
    // Mark the nodes of the Steiner tree for sample A.

    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
    {
      int index = *rit;

      Node_type v = this->node(index);

      this->node(index).mark = true;

      while( (!this->is_root(v)) )
      {
        this->node(v.parent).marked_children.push_back(index);

        if( this->node(v.parent).mark == false )
        {
          this->node(v.parent).mark = true;
          index = v.parent;
          v = this->node(v.parent);
        }
        else
          break;

      } // while( (!this->is_root(v)) )

    } // for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )


    // Mark the nodes of the Steiner tree of sample B.

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
    {
      int index = *rit;

      Node_type v = this->node(index);

      this->node(index).mark_b = true;

      while( (!this->is_root(v)) )
      {
        this->node(v.parent).marked_children_b.push_back(index);

        if( this->node(v.parent).mark_b == false )
        {
          this->node(v.parent).mark_b = true;
          index = v.parent;
          v = this->node(v.parent);
        }
        else
          break;

      } // while( (!this->is_root(v)) )

    } // for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )

  } // mark_Steiner_trees_of_samples( ... )


  // Given a tree and two samples of its leaf nodes that are represented
  // by the iterator ranges [rbegin_a,rend_a) and [rbegin_b,rend_b), the 
  // following function unmarks the marked nodes that fall on paths between
  // leaves of the same sample. However, if there exist any isolated marked nodes
  // these are not going to get unmarked by this function.

  template < class KernelType, class NodeType >
  template< class RangeIterator > 
  void PhylogeneticMeasures::Phylogenetic_tree_bimodal<KernelType, NodeType>::
  unmark_Steiner_trees_of_samples( RangeIterator rbegin_a, RangeIterator rend_a,
                                   RangeIterator rbegin_b, RangeIterator rend_b )
  {  
    // Unmark all marked nodes of the first sample and their 'marked_subtree_leaves' fields.

    for( RangeIterator rit = rbegin_a; rit != rend_a; rit++ )
    {
      this->node(*rit).mark = false;
      this->node(*rit).marked_subtree_leaves = 0;
      Node_type v = this->node(*rit);

      while( (!this->is_root(v)) && this->node(v.parent).mark == true )
      {
        this->node(v.parent).mark = false;
        this->node(v.parent).marked_subtree_leaves = 0;
        this->node(v.parent).marked_children.clear();
        v = this->node(v.parent);
      }
    }

    // Unmark all marked nodes of the second sample and their 'marked_subtree_leaves' fields.

    for( RangeIterator rit = rbegin_b; rit != rend_b; rit++ )
    {
      this->node(*rit).mark_b = false;
      this->node(*rit).marked_subtree_leaves_b = 0;
      Node_type v = this->node(*rit);

      while( (!this->is_root(v)) && this->node(v.parent).mark_b == true )
      {
        this->node(v.parent).mark_b = false;
        this->node(v.parent).marked_subtree_leaves_b = 0;
        this->node(v.parent).marked_children_b.clear();
        v = this->node(v.parent);
      }
    }

    this->clear_marked_nodes();
    this->clear_marked_nodes_b();

  } // unmark_Steiner_trees_of_samples(...)

#endif // PHYLOGENETIC_TREE_BIMODAL_IMPL_H

