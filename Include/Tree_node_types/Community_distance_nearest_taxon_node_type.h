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

#ifndef COMMUNITY_DISTANCE_NEAREST_TAXON_NODE_TYPE_H
#define COMMUNITY_DISTANCE_NEAREST_TAXON_NODE_TYPE_H

namespace PhylogeneticMeasures 
{
  template< class KernelType>
  struct Community_distance_nearest_taxon_auxiliary_data
  {
    typedef KernelType                    Kernel;
    typedef typename Kernel::Number_type  Number_type;

    Number_type first_min_a, second_min_a, rest_tree_min_a,
                first_min_b, second_min_b, rest_tree_min_b;

    Community_distance_nearest_taxon_auxiliary_data():first_min_a(Number_type(-1.0)), second_min_a(Number_type(-1.0)), 
	                                                  rest_tree_min_a(Number_type(-1.0)), 
                                                          first_min_b(Number_type(-1.0)), 
	                                                  second_min_b(Number_type(-1.0)), 
                                                          rest_tree_min_b(Number_type(-1.0)){}
  };

  template< typename KernelType>
  struct Community_distance_nearest_taxon_node_type:
  public KernelType::template Bimodal_node_augmented<Community_distance_nearest_taxon_auxiliary_data<KernelType> >
  {};

} // namespace PhylogeneticMeasures 

#endif //COMMUNITY_DISTANCE_NEAREST_TAXON_NODE_TYPE_H
