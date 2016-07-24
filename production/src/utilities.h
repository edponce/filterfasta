#ifndef UTILITIES_H
#define UTILITIES_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <unistd.h>

#define ERROR       -1         // Error code from failed functions

double get_wtime();
int computePartitionOffsets(long long int **, int *, int, long int, char);


#endif
