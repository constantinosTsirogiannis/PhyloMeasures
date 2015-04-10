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

#ifndef PHYLOGENETIC_TREE_BASE_H
#define PHYLOGENETIC_TREE_BASE_H

#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<set>
#include<map>
#include<cmath>

namespace PhylogeneticMeasures 
{

// Definition of a Phylogenetic tree structure
template < class KernelType, class NodeType = typename KernelType::Unimodal_node >
class Phylogenetic_tree_base
{
  public:

  typedef KernelType                                       Kernel;
  typedef NodeType                                         Node_type;
  typedef typename Kernel::Numeric_traits                  Numeric_traits;
  typedef typename Kernel::Number_type                     Number_type;
  typedef typename Numeric_traits::Absolute_value          Absolute_value;
  typedef typename Numeric_traits::Is_exact                Is_exact;
  typedef Phylogenetic_tree_base<Node_type, Kernel>        Self; 

  typedef typename Kernel::Exception_type                  Exception_type;
  typedef typename Kernel::Exception_functor               Exception_functor;  

private:

 // Minor predicate used for efficient sorting of alpharithmetic strings.
 struct _Is_placed_before
 {
   _Is_placed_before(){};

   bool operator()( const std::string &s1, const std::string &s2 ) const
   {
     if( s1.size() < s2.size() )
       return true;

     if( s1.size() > s2.size() )
       return false;

     int comp = s1.compare(s2);

     if(comp < 0)
       return true ;

     return false;
  }
 }; // Is_placed_before

 typedef std::map< std::string, int, _Is_placed_before >  Leaves_container_type;

public:

 typedef typename Leaves_container_type::iterator  Leaves_iterator;

private:

  int _assign_subtree_leaves( int index, bool marked );

  int _post_order_traversal( int index, std::vector<int> &swap_vector, int next_index_value);

  int _recursive_parse_tree_string( std::string &tree_string, int index );

  int _extract_taxon_and_distance( std::string &tree_string, int index, Node_type &v);

  template <class OUTSTRM >  
  void _print_tree(OUTSTRM &outs, int index );

  
  void _print_tree(Node_type &root);
  Number_type _check_if_ultrametric(int index);

  bool _check_if_ultrametric();

  int _compute_subtree_edges( int index )
  {
    if(this->root_index() == index)
      _subtree_edges.assign(this->number_of_nodes(),0);

    Node_type u = this->node(index);

    for(int i=0; i<u.number_of_children(); i++ )
      _subtree_edges[index] += _compute_subtree_edges( u.children[i] );

    return _subtree_edges[index]+1;
  }

public:

  Phylogenetic_tree_base():_root_index(-1),_is_ultrametric(-1),
                           _assigned_all_subtree_leaves(false),
                           _is_ultrametric_predicate_precision(0.001){}

  Node_type& node( int index)
  { return _container[index]; }

  int number_of_nodes()
  { return _container.size(); }

  int number_of_leaves()
  { return _leaves.size(); }

  void set_is_ultrametric_predicate_precision(Number_type &is_ultrametric_predicate_precision)
  { 
    _is_ultrametric_predicate_precision = is_ultrametric_predicate_precision;
    _is_ultrametric=-1;
  }

  bool is_ultrametric()
  {
    if(_is_ultrametric == -1)
      _check_if_ultrametric();

    if(_is_ultrametric == 1)
      return true;
    
    return false;
  }

  // Returns the index of the root node, duh.
  int root_index()
  { return _root_index; }

  // Returns the root node of the tree.
  Node_type root()
  {
    if( _root_index < 0)
      return Node_type();

    return _container[_root_index];
  }

  // Checks if the given node has no parent node.
  bool is_root( Node_type& node)
  { return node.parent == -1; }

  Leaves_iterator leaves_begin()
  { return _leaves.begin(); }

  Leaves_iterator leaves_end()
  { return _leaves.end(); }

  Leaves_iterator find_leaf(std::string &name)
  { return _leaves.find(name); }

  Leaves_iterator find_leaf(const std::string &name)
  { return _leaves.find(name); }

  int subtree_edges( int i )
  { 
    if(_subtree_edges.size()==0 && _container.size() > 0)    
      _compute_subtree_edges(_root_index);

    return _subtree_edges[i]; 
  }

  // This function computes the Sackin's index of a phylogenetic tree;
  // this index is equal to the average number of leaf nodes that
  // appear in the subtree of a node in the tree.
  Number_type sackins_index();
  
  // The next function computes for every edge e how many leaves
  // fall in the subtree of e, and returns the number of
  // distinct values for this measure among all edges in the tree.
  // This is the index which we call "DSSI" in the following paper:
  //
  // C. Tsirogiannis, B. Sandel and A. Kalvisa.
  // New Algorithms for Computing Phylogenetic Biodiversity.
  // In Proc. Workshop of Algorithms in Biology (WABI), 2014. To appear.

  int distinct_subspecies_sizes_index();
  
  // Computes the sum of the costs of all possible simple
  // paths between pairs of leaves in the tree. 
  Number_type total_path_costs();
  
  ///////////////////////////////////////////////////////////
  // Functions that are used for marking nodes in the tree //
  ///////////////////////////////////////////////////////////
  
  int assign_marked_subtree_leaves( int index )
  { return _assign_subtree_leaves(index,true); }

  int assign_all_subtree_leaves(int index )
  { 
    if( _assigned_all_subtree_leaves )
      return this->node(index).all_subtree_leaves;

    _assigned_all_subtree_leaves = true;

    return _assign_subtree_leaves(index,false);
  }

  int number_of_marked_nodes()
  { return _marked_nodes.size(); }

  int marked_node(int i)
  { return _marked_nodes[i]; }

  void clear_marked_nodes()
  { _marked_nodes.clear(); }

  bool& assigned_all_subtree_leaves()
  { return _assigned_all_subtree_leaves;}
  
 
  // Starting from the nodes indicated by index1 and index2, trace the paths
  // from those nodes up to the root of the tree. Find the deepest node at
  // which the two paths intersect (intersection_index) and mark this node.
  int compute_intersection_node_index( int index1, int index2 );

  template< class RangeIterator >
  Number_type mark_Steiner_tree_of_sample(RangeIterator rbegin, RangeIterator rend );

  template< class RangeIterator >
  void unmark_Steiner_tree_of_sample(RangeIterator rbegin, RangeIterator rend );

  template< class OutputIterator >
  Number_type compute_shortest_path_cost_to_subtree_leaf(int index, OutputIterator ot);
  
  template< class OutputIterator >
  Number_type compute_longest_path_cost_to_subtree_leaf(int index, OutputIterator ot);

  void compute_shortest_path_cost_to_outside_subtree_leaf( int index, 
                                                           std::vector<Number_type> &inside_subtree_costs,
                                                           std::vector<Number_type> &outside_subtree_costs);

  int compute_number_of_distance_violations();
  
  ///////////////////////////////////////////////////////
  // Functions that are used for handling strings that // 
  // represent trees in the Newick format.             //
  ///////////////////////////////////////////////////////
  
  // Constructs a new tree from the input string 'tree_str'.
  // Precondition: 'tree_str' is written in Newick format.
  void construct_from_string( std::string &tree_str );
  
  // Constructs a new tree from the input text file 'filename'.
  // Precondition: 'filename' stores a text string in Newick format.
  void construct_from_file(const char *filename );
  
  // Constructs a new tree from a given set of arrays which .
  // This is a practical way to build a tree from the data of 
  // the standard tree object that is used in the R package ape.
  void construct_from_edge_data(std::vector<int> &parents, std::vector<int> &children, 
                                std::vector<Number_type> &edge_weights, 
                                std::vector<std::string> &species_names);

  template <class OUTSTRM >  
  void print_tree(OUTSTRM &outs);

  void print_tree()
  { print_tree(std::cout); }

  void clear()
  {    
    _container.clear();
    _leaves.clear();
    _root_index=-1;
    _is_ultrametric = -1; 
    _assigned_all_subtree_leaves = false;
  }

private:

  std::vector<Node_type> _container;  // Vector that stores all the tree nodes in post-order
                                      // (except the root which is the last element of the vector).

  Leaves_container_type   _leaves; // Structure that maps taxon names of leaves
                                     // to indices of nodes in the tree.
                                     // TODO: Check if standard library string
                                     // comparison predicate is faster.

  int  _root_index; // Index to root node;
  int  _is_ultrametric; // Indicates if the tree is ultrametric: 
                        // that means all simple paths from the root to leaf nodes 
						// in the tree have the same cost.

  std::vector<int> _subtree_edges; // this vector stores for every tree edge
                                   // the number of edges that appear in its
                                   // subtree. This allows for computing in
                                   // constant time if an edge e appears in the
                                   // subtree of an edge l given their post-order indices.

   bool _assigned_all_subtree_leaves; // Used to determine if the number of subtree_leaves
                                  // has already been computed for every tree node.
								  
   std::vector<int> _marked_nodes; // Stores the indices of the nodes that
                                   // constitute the Steiner tree of a sample of leaves. 

   Number_type _is_ultrametric_predicate_precision; // Stores the precision error that
                                                    // we allow for deciding whether two
                                                    // paths from the root to leaf nodes have
                                                    // the same length.
								   
}; // class Phylogenetic_tree_base

} // namespace PhylogeneticMeasures

#include "Phylogenetic_tree_base_impl.h"

#endif // PHYLOGENETIC_TREE_H

