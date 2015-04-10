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

#ifndef PHYLOGENETIC_TREE_BASE_IMPL_H
#define PHYLOGENETIC_TREE_BASE_IMPL_H

#include<string>
#include<stack>
#include<vector>
#include<set>
#include<cstdlib>
#include<sstream>
#include<fstream>
  
  template < class KernelType, class NodeType >
  typename PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::Number_type
  PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::_check_if_ultrametric(int index)
  {
    Node_type v = node(index);

    if( v.children.size() == 0 )
      return Number_type(0.0);

    Number_type res(-1.0);

    for(int i=0; i<v.children.size(); i++)
    {
      Number_type leaf_path_cost = _check_if_ultrametric(v.children[i]);

      if(leaf_path_cost < Number_type(0.0))
        return Number_type(-1.0);

      leaf_path_cost = leaf_path_cost + Number_type(node(v.children[i]).distance);

      if(res == Number_type(-1.0))
        res = leaf_path_cost;
      else if( Absolute_value()(res - leaf_path_cost) > _is_ultrametric_predicate_precision ) 
        return Number_type(-1.0); // Numeric error when checking equality

    } // for(int i=0; i<v.children.size(); i+

    return res;

  } // Number_type _check_if_ultrametric(int index)


  template < class KernelType, class NodeType >
  bool PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::_check_if_ultrametric()
  {
    if( _is_ultrametric == 1 || _is_ultrametric == 0 )
      return bool(_is_ultrametric);

    Number_type res = _check_if_ultrametric(root_index());

    if( res >= Number_type(0.0) )
    {
      _is_ultrametric = 1;
      return true; 
    }
    
    _is_ultrametric = 0;
    return false;

  }

  // This function computes the Sackin's index of a phylogenetic tree;
  // this index is equal to the average number of leaf nodes that
  // appear in the subtree of each node in the tree.
  template < class KernelType, class NodeType >
  typename PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::Number_type 
  PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::sackins_index()
  {
    this->assign_all_subtree_leaves(this->root_index());

    Number_type sum_subtree_leaves(0.0);

    for( int i=0 ; i < this->number_of_nodes(); i++ )
      sum_subtree_leaves += Number_type(this->node(i).all_subtree_leaves);

    sum_subtree_leaves = sum_subtree_leaves - this->number_of_leaves();

    return sum_subtree_leaves;
  }
  

  // This function computes for every edge e how many leaves
  // fall in the subtree of e, and returns the number of
  // distinct values for this measure among all edges in the tree.
  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::distinct_subspecies_sizes_index()
  {
    this->assign_all_subtree_leaves(this->root_index());

    std::set<int> distinct_values;

    for( int i=0 ; i < this->number_of_nodes(); i++ )
      distinct_values.insert(this->node(i).all_subtree_leaves);

    return distinct_values.size();
  }
  
  // Consider the min cost Steiner tree that spans a given sample of leaf nodes. 
  // The following function computes for each node v of this subtree the number of sample leaf nodes
  // situated at the subtree of v. This number is assigned to the field
  // 'marked_subtree_leaves' of v. The field '.mark' of v is assumed to have already
  // been set to 'true' if the subtree of v contains at least one sample leaf.

  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::
  _assign_subtree_leaves( int index, bool marked )
  {  
    if( marked == true )
      this->node(index).marked_subtree_leaves=0;
    else
      this->node(index).all_subtree_leaves=0;

    if( marked == true && this->node(index).mark == false )
      return 0;

    if( marked == true )
      _marked_nodes.push_back(index);

    // Leaf node
    if( this->node(index).number_of_children() == 0 )
    {
      if( marked == true )
        this->node(index).marked_subtree_leaves=1;
      else       
        this->node(index).all_subtree_leaves=1;
    }
    else
    {
      Node_type v = this->node(index);

      // Interior node
      
      if(marked==true)
        for( int i = 0; i < v.number_of_marked_children(); i++ )
          this->node(index).marked_subtree_leaves += _assign_subtree_leaves(v.marked_children[i],true);
      else
        for( int i = 0; i < v.number_of_children(); i++ )
          this->node(index).all_subtree_leaves += _assign_subtree_leaves(v.children[i],false);
   

    } // else of if( tree.node(index).number_of_children() == 0 )

    if( marked == true )
      return this->node(index).marked_subtree_leaves;
    else
      return this->node(index).all_subtree_leaves;
	  
  } // _assign_subtree_leaves(...)

  // Given a tree and a sample of its leaf nodes that are represented
  // by the iterator range [rbegin,rend), the following function 
  // marks all the tree nodes that fall on the Steiner tree 
  // that connects the sample leaves. The function returns
  // the cost (i.e. the sum of the weights of the edges) of this
  // Steiner tree. For this function to work appropriately, 
  // we should mark the root of the Steiner tree in advance.

  template < class KernelType, class NodeType >
  template< class RangeIterator >
  typename PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::Number_type 
  PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::mark_Steiner_tree_of_sample
  ( RangeIterator rbegin, RangeIterator rend )
  {
    Number_type total_dist(0.0);

    for( RangeIterator rit = rbegin; rit != rend; rit++ )
    {
      int index = *rit;
      Node_type v = this->node(index);

      this->node(index).mark = true;

      if( index != this->root_index() ) // Maybe redundant
        total_dist += Number_type(v.distance);

      while( (!this->is_root(v)) )
      {
        this->node(v.parent).marked_children.push_back(index);

        if(this->node(v.parent).mark == false)
        {
          this->node(v.parent).mark = true;
          index = v.parent;
          v = this->node(v.parent);
          total_dist += Number_type(v.distance);
        }
        else
          break;

      } // while( (!tree.is_root(v)) )

    } // for( RangeIterator rit = rbegin; rit != rend; rit++ )
  
    return total_dist;

  } // mark_Steiner_tree_of_sample(...)


  // Given a tree and a sample of its leaf nodes that are represented
  // by the iterator range [rbegin,rend), the following function
  // unmarks the marked nodes that fall on paths between sample
  // leaves. However, if there exist any isolated marked nodes
  // these are not going to get unmarked by this function.

  template < class KernelType, class NodeType >
  template< class RangeIterator >
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType,NodeType>::unmark_Steiner_tree_of_sample
  (RangeIterator rbegin, RangeIterator rend )
   {
    for( RangeIterator rit = rbegin; rit != rend; rit++ )
    {
      this->node(*rit).mark = false;
      this->node(*rit).marked_subtree_leaves = 0;
      Node_type v = this->node(*rit);
     

      while( (!this->is_root(v)) && this->node(v.parent).mark == true )
      {
        this->node(v.parent).mark = false;
        this->node(v.parent).marked_children.clear();
        this->node(v.parent).marked_subtree_leaves = 0;
        v = this->node(v.parent);
      }
    }

    clear_marked_nodes();

  } // unmark_Steiner_tree_of_nodes(...)
  
  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  compute_intersection_node_index( int index1, int index2 )
  {
    if(index1 == index2)
      return index1;

    int current_index=index1;

    while( current_index != root_index() )
    {
      if( current_index >= index2 && index2 >= current_index - this->subtree_edges(current_index) )
        return current_index;
      else 
        current_index = node(current_index).parent;
    }

    return current_index;

  } // compute_intersection_node_index(...)

    
  template < class KernelType, class NodeType >
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::construct_from_file( const char *filename )
  {
    std::ifstream in(filename);
    std::string tree_str;

    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file that describes the tree.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    // Read tree string

    char a;

    in >> a;

    while(in.good())
    {
      tree_str.push_back(a);
      in >> a;
    }

    construct_from_string( tree_str );
	
  } // construct_from_file( const char *filename ) 


  template < class KernelType, class NodeType >
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  construct_from_string( std::string &tree_string)
  { 
    // Clear the tree from previously stored information
    this->clear();

    _recursive_parse_tree_string(tree_string, 0);
 
    _root_index = _container.size()-1;

    // Create the data structure that maps each leaf node name
    // to the index of the respective leaf node.
    for( int i=0; i<_container.size(); i++ )
      if( _container[i].number_of_children() == 0 )
        _leaves[_container[i].taxon] = i;
  }


  template < class KernelType, class NodeType >
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  construct_from_edge_data(std::vector<int> &parents, std::vector<int> &children, 
                           std::vector<Number_type> &edge_weights, 
                           std::vector<std::string> &species_names)
  {
    // Clear the tree from previously stored information
    this->clear();

    int number_of_edges = children.size();
    int number_of_leaves = species_names.size();
    int number_of_nodes = number_of_edges+1;

    if(number_of_nodes<1)
      return;

    _container.assign(number_of_nodes, Node_type());

    // Populate the node container with 
    // the current set of input data.

    for(int i=0; i<number_of_edges; i++)
    {
      int parent_index = parents[i]-1,
          child_index  = children[i]-1;

      _container[parent_index].children.push_back(child_index);
      _container[child_index].parent = parent_index;

      _container[child_index].distance = edge_weights[i];

      if(child_index<number_of_leaves)
        _container[child_index].taxon = species_names[child_index];
      else
      {
        std::stringstream sstrm;
        sstrm << child_index;
        _container[child_index].taxon += '\'';
        _container[child_index].taxon += sstrm.str();
        _container[child_index].taxon += '\'';
      }

    } // for(int i=0; i<number_of_edges; i++)


    // Next: positions of nodes into the container will get
    // swapped so that they match a post-order traversal of the tree.

    std::vector<int> swap_vector;

    swap_vector.assign(_container.size(),-1);

    // Find root index

    int index = 0;

    while(_container[index].parent != -1)
      index = _container[index].parent;

    _root_index = index;
    std::stringstream sstrm;
    sstrm << _root_index;
    _container[_root_index].taxon += '\'';
    _container[_root_index].taxon += sstrm.str();
    _container[_root_index].taxon += '\'';

    _post_order_traversal(_root_index,swap_vector,0);

    bool swapped_zero = false;

    for(int i=0; i<number_of_nodes; i++)
      if(swap_vector[i] > 0 || (swap_vector[i] == 0 && swapped_zero == false) )
      {
        int temp_index = i;
        Node_type swapped_node = _container[temp_index], temp_node;

        while( swap_vector[temp_index] > 0 || 
               ( swap_vector[temp_index] == 0 && swapped_zero==false )  )
        {
          Node_type temp_node = _container[swap_vector[temp_index]];
          _container[swap_vector[temp_index]] = swapped_node;
          swapped_node = temp_node; 

          int new_index = swap_vector[temp_index]; 
          swap_vector[temp_index] = -swap_vector[temp_index];
          
          if(swap_vector[temp_index] == 0)
            swapped_zero = true;

          temp_index = new_index;
        }

        _container[-swap_vector[temp_index]] = swapped_node;

      } // if(swap_vector[i] != -1)


    for(int i=0; i<_container.size(); i++)
    {
      if(_container[i].parent !=-1)
        _container[i].parent = -swap_vector[_container[i].parent];

      for(int j=0; j<_container[i].number_of_children(); j++)
        _container[i].children[j] = -swap_vector[_container[i].children[j]];
    }


    for(int i=0; i<number_of_nodes; i++)
      if(_container[i].number_of_children() == 0)
        _leaves[_container[i].taxon]=i;

    _root_index = _container.size()-1;
    swap_vector.clear();

    //std::ofstream of("temp_tree.txt");
    //print_tree(of);
  
  }// construct_from_edge_data(...)

  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  _post_order_traversal( int index, std::vector<int> &swap_vector, int next_index_value)
  {
    int updated_value = next_index_value;

    for(int i=0; i<_container[index].number_of_children(); i++)
      updated_value = _post_order_traversal(_container[index].children[i], swap_vector, updated_value);

    swap_vector[index] = updated_value++;

    return updated_value;
  }

  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  _recursive_parse_tree_string( std::string &tree_string, int index )  
  {
    bool end_of_node_found = false;
    Node_type v; 

    if(tree_string[index] == '(')
    {
      while(!end_of_node_found)
      {

        index = _recursive_parse_tree_string(tree_string, ++index );
        v.children.push_back(_container.size()-1);

        do
        {
          index++;
        }while(tree_string[index]==' ' || tree_string[index]=='\n' );       
     
        if(tree_string[index]==')' )
          end_of_node_found = true;
        else if(tree_string[index]!=',') 
        {
          std::string exception_msg(" The input tree string does not follow a valid Newick format.\n");
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      } // while(!end_of_node_found)

      index++;

    } // if(tree_string[index] == '(')
   
    index = _extract_taxon_and_distance(tree_string, index, v);
    _container.push_back(v);

    for(int i=0; i<v.number_of_children(); i++)
      _container[v.children[i]].parent = _container.size()-1;

    return index;

  } // _recursive_parse_tree_string( ... )

  template < class KernelType, class NodeType >
  int PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  _extract_taxon_and_distance( std::string &tree_string, int index, Node_type &v)
  {
    std::string temp, taxon, distance;

    while( index < tree_string.length() && tree_string[index] != ')'  && 
           tree_string[index] != ',' && tree_string[index] != ';' )
    {
      temp.push_back(tree_string[index]);
      index++;
    }
       
    size_t split_index = temp.find(':');

    if(split_index != std::string::npos)
    {
      v.taxon = temp.substr(0,split_index);
      v.distance = double(atof(temp.substr(split_index+1).c_str()));
    }
    else
      v.taxon = temp;

    index--;

    return index;

  } // _extract_taxon_and_distance(...)

  template < class KernelType, class NodeType >
  template <class OUTSTRM >  
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  print_tree(OUTSTRM &outs)
  {
    if(this->number_of_nodes() == 0)
      return; 

    _print_tree(outs, this->root_index());

    outs << ";" << std::endl;
  }


  template < class KernelType, class NodeType >
  template <class OUTSTRM >  
  void PhylogeneticMeasures::Phylogenetic_tree_base<KernelType, NodeType>::
  _print_tree(OUTSTRM &outs, int index )  
  {
    Node_type v = this->node(index);

    if(v.number_of_children()>0)
    {
      outs << "(";

      for(int i=0; i<v.number_of_children(); i++)
      {
        _print_tree(outs,v.children[i]);

        if(i < v.number_of_children()-1)
          outs << ",";
      }

      outs << ")";
    }

    if( v.distance >= double(0.0) )
      outs << v.taxon << ":" << v.distance ;
    else
      outs << v.taxon;

  } // _print_tree(...)

#endif // PHYLOGENETIC_TREE_BASE_IMPL_H
