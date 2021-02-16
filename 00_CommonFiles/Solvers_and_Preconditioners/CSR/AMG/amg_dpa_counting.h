# ifndef counting_h
# define counting_h

#include <stdlib.h>

/*----------------------------------------------------------------
	This is a specific data structure for Counting Sort algorithm
	for the DPA matching algorithm.

	Every vertex i has a m coefficient. The focus is to order all
	vertexes by m.

	Each vertex i with a m coefficient equals x is appended to
	the x-th list of an array doubly linked lists.
  ----------------------------------------------------------------*/

// doubly linked node
typedef struct Nodec {
    int i, m;
    struct Nodec *next, *prev;
} Nodec;

// doubly linked list
typedef struct {
    int n;
    Nodec *head, *tail;
} listc;

// counting sort structure
typedef struct {
    listc **count;	// list of pointers for counting array elements
    Nodec **nodes;	// counting array
    int *m;		// array of m values for each vertex i
    int min_m;		// minimun index with non empty list
    int n;		// number of vertexes
} counting;

listc * create_lc();

void destroy_lc(listc *l);

/*----------------------------------------------------------------*/

counting * create_c(int n);

void destroy_c(counting * c);

Nodec * insert_c(counting *c, int i, int m) ;

void remove_c(counting *c, int i);

Nodec* move_c(counting *c, int i, int new_m);

# endif
