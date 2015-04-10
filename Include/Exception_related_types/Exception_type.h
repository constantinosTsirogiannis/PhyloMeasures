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

#ifndef EXCEPTION_TYPE_H
#define EXCEPTION_TYPE_H

#include<string>
#include<cstdlib>

namespace ExceptionRelatedTypes
{
  class Exception_type
  {    
   public:

    void get_error_message(std::string msg)
    { _msg = msg;}

    std::string return_error_message()
    { return _msg; }

   private:

    std::string _msg;

  }; // class Exception_type


} // ExceptionRelatedTypes

#endif // EXCEPTION_TYPE_H
