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

#ifndef NUMERIC_TRAITS_DOUBLE_H
#define NUMERIC_TRAITS_DOUBLE_H

#include<cmath>

namespace PhylogeneticMeasures 
{
  struct Numeric_traits_double
  {
    typedef double Number_type;

    class Is_exact
    {
      public:
 
        bool operator()(void)
        { return false;}
    };

    class To_double
    {
      public:

        double operator()(double x)
        { return x; }
    };

    class Ceiling
    {
      public:

        double operator()(double x)
        { return ceil(x); }
    };

    class Square_root
    {
      public:

        double operator()(double x)
        { return std::sqrt(x); }
    };

    class Absolute_value
    {
      public:
 
        double operator()(double x)
        { return std::abs(x);}
    };

  }; // struct Numeric_traits_double

} //namespace PhylogeneticMeasures

#endif // NUMERIC_TRAITS_DOUBLE_H
