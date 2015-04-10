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

#ifndef EXCEPTION_FUNCTOR_LOUD_H
#define EXCEPTION_FUNCTOR_LOUD_H

#include<string>
#include<cstdlib>
#include<exception>
#include"Exception_type.h"

namespace ExceptionRelatedTypes
{
  class Exception_functor
  {
   public:

    Exception_functor(){}
    
    void operator()(Exception_type excp)
    {
      std::string msg = excp.return_error_message();
      std::cout << " Error:" << msg << std::endl;
      exit(1);
    }

    void issue_warning(std::string str)
    { std::cout << str << std::endl << std::endl;}

  }; // class Exception_functor

} // ExceptionRelatedTypes

#endif // EXCEPTION_FUNCTOR_LOUD_H
