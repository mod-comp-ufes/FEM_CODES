# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include "amg_fib_heap.h"
# include "../../../Allocation_Operations/allocations.h"
# define nINF -2147483648

fib_heap * make_heap () {
	fib_heap *H = (fib_heap *) mycalloc("fib_heap of make_heap",1,sizeof(fib_heap));
	H->n = 0;
	H->min = NULL;
	return H;
};

node_fh * insert_heap (fib_heap *H, int v, int k) {
	node_fh *x = (node_fh *) mycalloc("x of insert_heap",1,sizeof(node_fh));
	x->v = v;
	x->key = k;
	x->degree = 0;
	x->p = NULL;
	x->child = NULL;
	x->mark = 0;
	if (!H->min) {
		x->left = x;
		x->right = x;
		H->min = x;
	} else {
		x->left = H->min->left;
		x->right = H->min;
		H->min->left->right = x;
		H->min->left = x;
		if (x->key < H->min->key) H->min = x;
	}
	H->n++;
	return x;
};

node_fh * minimum_heap (fib_heap *H) { return H->min; };

void link (fib_heap *H, node_fh *y, node_fh *x) {
	y->left->right = y->right;
	y->right->left = y->left;
	if (!x->child) {
		y->left = y;
		y->right = y;
		x->child = y;
	} else {
		y->left = x->child->left;
		y->right = x->child;
		x->child->left->right = y;
		x->child->left = y;
	}
	y->p = x;
	x->degree++;
	y->mark = 0;
};

void consolidate (fib_heap *H) {
	int i, d, max = floor(log2(H->n))+2;
	node_fh **A, *x, *y, *w;
	A = (node_fh **) mycalloc("A of consolidate",max,sizeof(node_fh *));
	for (i=0; i<max; i++) A[i] = NULL;
	w = H->min;
	do {
		x = w;
		d = x->degree;
		w = w->right;
		while (d<max && A[d]) {
			y = A[d];
			if (x->key > y->key) { y = x; x = A[d]; }
			if (w == y) w = y->right;
			if (H->min == y) H->min = y->right;
			link(H, y, x);
			A[d] = NULL;
			d++;
		}
		A[d] = x;
	} while (w!=H->min);
	H->min = NULL;
	for (i=0; i<max; i++) {
		if (A[i]) {
			if (!H->min) {
				A[i]->left = A[i];
				A[i]->right = A[i];
				H->min = A[i];
			} else {
				A[i]->left = H->min->left;
				A[i]->right = H->min;
				H->min->left->right = A[i];
				H->min->left = A[i];
				if (A[i]->key < H->min->key) H->min = A[i];
			}
		}
	}
	myfree(A);
};

node_fh * extract_min_heap (fib_heap *H) {
	node_fh *z = H->min, *x;
	if (z) {
		if (z->child) {
			x = z->child;
			do { x->p = NULL; x = x->right; } while (x!=z->child);
			H->min->left->right = z->child;
			z->child->left->right = H->min;
			x = z->child->left;
			z->child->left = H->min->left;
			H->min->left = x;
		}
		z->left->right = z->right;
		z->right->left = z->left;
		if (z == z->right) H->min = NULL;
		else { H->min = z->right; consolidate(H); }
		H->n--;
	}
	return z;
};

void cut (fib_heap *H, node_fh *x, node_fh *y) {
	x->left->right = x->right;
	x->right->left = x->left;
	if (y->child == x) y->child = x->right;
	y->degree--; if (!y->degree) y->child = NULL;
	x->left = H->min->left;
	x->right = H->min;
	H->min->left->right = x;
	H->min->left = x;
	x->p = NULL;
	x->mark = 0;
};

void cascading_cut (fib_heap *H, node_fh *y) {
	node_fh *z = y->p;
	if (z) {
		if (!y->mark) y->mark = 1;
		else { cut(H, y, z); cascading_cut(H, z); }
	}
};

void decrease_key_heap (fib_heap *H, node_fh *x, int k) {
	node_fh *y;
	if (k > x->key) printf("Erro! Novo valor de key eh maior que o atual\n");
	x->key = k;
	y = x->p;
	if (y && x->key < y->key) { cut(H, x, y); cascading_cut(H, y); }
	if (x->key < H->min->key) H->min = x;
};

void delete_heap (fib_heap *H, node_fh *x) {
	node_fh *min;
	decrease_key_heap(H, x, nINF);
	min = extract_min_heap(H);
	myfree(min);
};
