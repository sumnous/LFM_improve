#if !defined(LFM_INCLUDED)
#define LFM_INCLUDED



#include <math.h>
#include <iostream>
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <string> 
#include <fstream>
#include <ctime>
#include <iterator>
#include <algorithm>




using namespace std;

#define UNLIKELY -2147483647


#define R2_IM1 2147483563
#define R2_IM2 2147483399
#define R2_AM (1.0/R2_IM1)
#define R2_IMM1 (R2_IM1-1)
#define R2_IA1 40014
#define R2_IA2 40692
#define R2_IQ1 53668
#define R2_IQ2 52774
#define R2_IR1 12211
#define R2_IR2 3791
#define R2_NTAB 32
#define R2_NDIV (1+R2_IMM1/R2_NTAB)
#define R2_EPS 1.2e-7
#define R2_RNMX (1.0-R2_EPS)



#include "readwrite.cpp"

#include "sarray.h"
#include "static_network.h"


	

ofstream arcout ("./Library/Files/archive.dat");
ifstream get_part("./Library/Files/archive.dat");

map <double, int> occurrences;		// this map cointains the fitness of the partition and how many times it has been found
map <double, streampos> archive_pos;




map <double, int>::iterator find_value_in_occurrences (const double & f)  {

	map <double, int>::iterator it_mr;	
	for (it_mr=occurrences.lower_bound(f-2.*R2_EPS); it_mr!=occurrences.end(); it_mr++) {
		if (fabs(it_mr->first-f)<R2_EPS)
			break;
	}
	
	return it_mr;

}



#include "ddi.h"
#include "lfm_network.h"


#endif




