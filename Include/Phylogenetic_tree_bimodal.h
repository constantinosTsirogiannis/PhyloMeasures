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

#ifndef PHYLOGENETIC_TREE_BIMODAL_H
#define PHYLOGENETIC_TREE_BIMODAL_H

namespace PhylogeneticMeasures 
{

// Definition of a Phylogenetic tree structure
template < class KernelType, class NodeType = typename KernelType::Bimodal_node >
class Phylogenetic_tree_bimodal: public PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>
{
 public:

  typedef KernelType                                                                     Kernel;
  typedef NodeType                                                                       Node_type;
  typedef typename Kernel::Numeric_traits                                                Numeric_traits;
  typedef typename Kernel::Number_type                                                   Number_type;
  typedef typename PhylogeneticMeasures::Phylogenetic_tree_base<Kernel, Node_type>       Base;
  typedef Phylogenetic_tree_bimodal<Kernel, Node_type>                                   Self;
  typedef Self                                                                           Tree_type;

  Phylogenetic_tree_bimodal():Base(){}

 private:
  
  // Consider the minimum cost Steiner tree that connects all leaves of a given sample. 
  // The following function computes for each node v of this tree the number of 
  // sample leaves situated at the subtree of v. This number is assigned to the field
  // 'marked_subtree_leaves_b' of v.

  int _assign_marked_subtree_leaves_b( int index );    
  
  
 public:
 
  std::pair<int,int> assign_marked_subtree_leaves_both_sets( int index )
  { return std::make_pair( this->assign_marked_subtree_leaves(index), this->_assign_marked_subtree_leaves_b(index) ); }
 
  // Given a tree and two samples of its leaf nodes that are represented
  // by the iterator ranges [rbegin_a,rend_a) and [rbegin_b,rend_b), the 
  // following function marks all the tree nodes that fall on the two Steiner trees 
  // that connect the leaves of each sample. For this function to work appropriately, 
  // we should mark the roots of the Steiner trees in advance.

  template< class RangeIterator > 
  void mark_Steiner_trees_of_samples( RangeIterator rbegin_a, RangeIterator rend_a,
                                      RangeIterator rbegin_b, RangeIterator rend_b );

 
  // Given a tree and two samples of its leaf nodes that are represented
  // by the iterator ranges [rbegin_a,rend_a) and [rbegin_b,rend_b), the 
  // following function unmarks the marked nodes that fall on paths between
  // leaves of the same sample. However, if there exist any isolated marked nodes
  // these are not going to get unmarked by this function.

  template< class RangeIterator > 
  void unmark_Steiner_trees_of_samples( RangeIterator rbegin_a, RangeIterator rend_a,
                                        RangeIterator rbegin_b, RangeIterator rend_b );

  int number_of_marked_nodes_b()
  { return _marked_nodes_b.size(); }

  int marked_node_b(int i)
  { return _marked_nodes_b[i]; }

  void clear_marked_nodes_b()
  { _marked_nodes_b.clear(); }
  
 private:

   std::vector<int> _marked_nodes_b; // Stores the indices of the nodes that
                                     // constitute the Steiner tree that connects the leaves of the sample. 

}; // class Phylogenetic_tree_bimodal

} // namespace PhylogeneticMeasures

#include "Phylogenetic_tree_bimodal_impl.h"

#endif // PHYLOGENETIC_TREE_BIMODAL_H

