# ifndef allocations
# define allocations

#include "stdio.h"
#include "stdlib.h"

void *mycalloc(char *, int, int);
void myfree(void *ptr);
void list_leaks();
void free_leaks();
void list_leaks_and_free();

#endif
