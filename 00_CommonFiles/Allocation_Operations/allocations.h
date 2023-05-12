#ifndef _allocations_h_
#define _allocations_h_

#include "stdio.h"
#include "stdlib.h"


void *mycalloc(char *, int, int);
void myfree(void *ptr);
void list_leaks();

#endif
