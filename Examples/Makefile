#################################################################################
##    Copyright (C) 2015,  Constantinos Tsirogiannis.  Email: analekta@gmail.com
##
##    This file is part of PhyloMeasures.
##
##    PhyloMeasures is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    PhyloMeasures is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhyloMeasures.  If not, see <http://www.gnu.org/licenses/>
#################################################################################

CC = g++
CFLAGS =
INCLUDES = -I../Include/
DEPENDENCIES = ../Include/ ../Include/Tree_node_types/ ../Include/Measures/ ../Include/Measures/Measure_base/ ../Include/Exception_related_types/

all: compute_PD_expectation compute_MPD_matrix_query 

compute_PD_expectation: compute_PD_expectation.o
	$(CC) $(CFLAGS) $(INCLUDES) -o compute_PD_expectation compute_PD_expectation.o
compute_PD_expectation.o: compute_PD_expectation.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c compute_PD_expectation.cpp

compute_MPD_matrix_query: compute_MPD_matrix_query.o
	$(CC) $(CFLAGS) $(INCLUDES) -o compute_MPD_matrix_query compute_MPD_matrix_query.o
compute_MPD_matrix_query.o: compute_MPD_matrix_query.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c compute_MPD_matrix_query.cpp

clean:
	$(RM) *.o *~
