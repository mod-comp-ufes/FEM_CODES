#ifndef _C_Operations_h_
#define _C_Operations_h_

#include <stdio.h>
#include <string.h>


// copy x to y
void dmemcpy(int n, double *x, double *y);

// set zero to x
void memsetzero(int n, double *x);

#endif