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

#ifndef TREE_NODE_UNIMODAL_H
#define TREE_NODE_UNIMODAL_H

#include<vector>
#include<string>

namespace PhylogeneticMeasures 
{
  // Definitions of a phylogenetic tree node used for measures
  // that involve a single sample of the tree species.

  template <class KernelType>
  struct Tree_node_unimodal
  {
    typedef KernelType   Kernel;

    std::string taxon; // Species/taxonomy name.
    double distance; // Distance from father node.
    std::vector<int> children; // Indices to child nodes.
    std::vector<int> marked_children; // Indices to child nodes that have been marked.
    int parent; // Index to parent node.
    bool mark; // Auxiliary flag.
    int all_subtree_leaves; // Number of all leaf nodes.
                        // in the subtree of this node.
    int marked_subtree_leaves; // Number of leaves in the subtree
                           // of this node which are also elements
                           // of the query sample.

    Tree_node_unimodal():distance(-1.0),parent(-1),mark(false),marked_subtree_leaves(0)
    {}

    int number_of_children()
    { return children.size(); }

    int number_of_marked_children()
    { return marked_children.size(); }

  }; // struct Tree_node_unimodal


  // Tree node which allows for
  // extra data stored with the node.

  template< class AUXILIARY_TYPE >
  struct Tree_node_unimodal_augmented: 
  public Tree_node_unimodal<typename AUXILIARY_TYPE::Kernel>, public AUXILIARY_TYPE
  {
    typedef AUXILIARY_TYPE                   Auxiliary_type;
    typedef typename Auxiliary_type::Kernel  Kernel;
    typedef Tree_node_unimodal<Kernel>       Base_type;

    Tree_node_unimodal_augmented():Base_type(), Auxiliary_type(){}
  };

} // namespace PhylogeneticMeasures 

#endif //PHYLOGENETIC_TREE_NODE_UNIMODAL_H
