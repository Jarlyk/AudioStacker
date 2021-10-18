#pragma once

#include <math.h>

inline double atodb(double a) {
	if (a <= 0) return -75;
	
	return 20*log10(a);
}

inline double dbtoa(double d) {
	return pow(10, d/20);
}