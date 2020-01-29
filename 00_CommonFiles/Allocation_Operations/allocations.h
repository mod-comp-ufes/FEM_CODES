# ifndef allocations
# define allocations

#include "stdio.h"
#include "stdlib.h"

void *mycalloc(char *, int, int);
void myfree(void *ptr);
void list_leaks();

#endif
