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

all: csv_query_CD csv_query_PD csv_query_CDNT csv_query_CAC csv_query_MNTD csv_query_CBL csv_query_MPD csv_query_PhyloSor csv_query_UniFrac measure_moments CAC_moments 

csv_query_CD: csv_query_CD.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_CD csv_query_CD.o
csv_query_CD.o: csv_query_CD.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_CD.cpp

csv_query_PD: csv_query_PD.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_PD csv_query_PD.o
csv_query_PD.o: csv_query_PD.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_PD.cpp

csv_query_CDNT: csv_query_CDNT.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_CDNT csv_query_CDNT.o
csv_query_CDNT.o: csv_query_CDNT.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_CDNT.cpp

csv_query_CAC: csv_query_CAC.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_CAC csv_query_CAC.o
csv_query_CAC.o: csv_query_CAC.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_CAC.cpp

csv_query_MNTD: csv_query_MNTD.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_MNTD csv_query_MNTD.o
csv_query_MNTD.o: csv_query_MNTD.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_MNTD.cpp

csv_query_CBL: csv_query_CBL.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_CBL csv_query_CBL.o
csv_query_CBL.o: csv_query_CBL.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_CBL.cpp

csv_query_MPD: csv_query_MPD.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_MPD csv_query_MPD.o
csv_query_MPD.o: csv_query_MPD.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_MPD.cpp

csv_query_PhyloSor: csv_query_PhyloSor.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_PhyloSor csv_query_PhyloSor.o
csv_query_PhyloSor.o: csv_query_PhyloSor.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_PhyloSor.cpp

csv_query_UniFrac: csv_query_UniFrac.o
	$(CC) $(CFLAGS) $(INCLUDES) -o csv_query_UniFrac csv_query_UniFrac.o
csv_query_UniFrac.o: csv_query_UniFrac.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c csv_query_UniFrac.cpp

measure_moments: measure_moments.o
	$(CC) $(CFLAGS) $(INCLUDES) -o measure_moments measure_moments.o
measure_moments.o: measure_moments.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c measure_moments.cpp

CAC_moments: CAC_moments.o
	$(CC) $(CFLAGS) $(INCLUDES) -o CAC_moments CAC_moments.o
CAC_moments.o: CAC_moments.cpp $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) -c CAC_moments.cpp

clean:
	$(RM) *.o *~
